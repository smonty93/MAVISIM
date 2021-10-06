import numpy as np
from astropy.io import fits
from tqdm import tqdm

class PSF:
    """
    PSF object class.

    ...
    
    Attributes
    ----------
    fft_data : np.array
        array storing the (minimal) rfft2 data of the given PSF
    xpos : float
        x-position of PSF point source in arcsec
    ypos : float
        y-position of PSF point source in arcsec
    Lambda : float
        wavelength used to capture the PSF

    Methods
    -------
    N/A

    """

    def __init__(self, fits_ext, padto, *, dtype=np.complex128):
        """Init PSF and perform RFFT2"""
        self.fft_data = (np.fft.rfft2(fits_ext.data,s=padto)).astype(dtype)
        self.xpos     = float(fits_ext.header["YPOS"])
        self.ypos     = float(fits_ext.header["XPOS"])
        self.Lambda   = float(fits_ext.header["LAMBDA"])
    
class TileGenerator:
    """
    Object for generating tiles to be sliced into final image.

    ...

    Attributes
    ----------
    source : List
        source list of objects, each with 3 elements:
            [0] -> flux
            [1] -> position [arcsec]
            [2] -> gaussian covariance [arcsec^2]
    psfs : List
        list of PSF objects (see PSF class)
    pixsize : float
        pixel size in arcsec used to build tile
    gauss_width_pix : int
        width of Gaussian square support in pixels
    gauss_width_as : float
        width of Gaussian square support in arcsec (edge to edge)
    gauss_dim : int ndarray
        shape of Gaussian square support
    psf_file : str
        fits filename storing all PSFs + PSF metadata
    psf_pitch : float
        gap between adjacent PFSs on square field grid
    psf_width_pix : int
        width of PSF square support in pixels
    psf_width_as : float
        width of PSF square support in arcsec (edge to edge)
    fourier_tile_dim : int ndarray
        shape of tile before cropping
    static : bool
        True -> use static PSF
        False -> use field variable PSF
    dtype : np.dtype
        dtype of all large arrays

    Methods
    -------

    """
    
    def __init__( self, source, psfs_file, gauss_width_pix, *,
             dtype = np.complex128, which_psf = None):
        """Initialise tile generator object by preparing source list and PSFs

        ...
        
        Parameters
        ----------
        source_list : List
            table of sources as output by Source.py
        psfs_file : str
            fits file containing all PSFs and PSF metadata
        gauss_width_pix: int
            support size in pixels to build the Gaussian star kernel
        dtype : np.dtype, optional
            complex data type to work with in DFT space
        which_psf : int, optional
            default value (None) -> use field varying PSFs
            otherwise -> use fits HDU[which_psf+1] PSF only
        """
        
        self.dtype=dtype

        # Parse source info:
        self.source_flux = source.flux.copy()
        self.source_pos  = source.gauss_pos.copy()
        self.source_cov  = source.gauss_cov.copy()
        
        # Define nominal image geometry:
        # TODO read this from fits Primary HDU metadata
        self.pixsize          = 3.75e-3         
        self.gauss_width_pix  = gauss_width_pix
        self.gauss_width_as   = self.gauss_width_pix*self.pixsize
        self.gauss_dim        = np.array([self.gauss_width_pix,self.gauss_width_pix])

        # Load PSF file and define Fourier Image Geometry:
        self.psf_file          = psfs_file
        psfs_fits              = fits.open(psfs_file, lazy_load_hdus=True)
        # TODO read this from fits Primary HDU metadata
        self.psf_pitch         = 3.0 
        self.psf_width_pix     = psfs_fits[1].data.shape[0]
        self.psf_width_as      = self.psf_width_pix * self.pixsize
        
        self.fourier_tile_dim = self.gauss_dim + np.array(psfs_fits[1].data.shape)
        
        if which_psf is None:
            # Use variable PSFs
            self.static = False
            # Create PSF objects and do their FFTs
            print("Using varible PSFs")
            print("Doing PSF FFTs:")
            self.psfs = []
            for psf in tqdm(psfs_fits[1:],leave=False):
                self.psfs.append(PSF(psf, self.fourier_tile_dim))
            psfs_fits.close()
            print("Done.                          ")
        else:
            # Use static PSFs
            self.static = True
            self.psfs = []
            self.psfs.append(PSF(psfs_fits[1+which_psf], self.fourier_tile_dim))
            psfs_fits.close()
            print(f"Using static PSF: {which_psf:3d} (fits hdu {which_psf+1:3d})" + 
                    f"\nPSF pos: ({self.psfs[0].xpos:0.3f},{self.psfs[0].ypos:0.3f})\"")

        self._init = True
        
    def get_effective_psf_fft(self, s_pos):
        """Takes star information and computes effective PSF.

        From star position, the convex combination PSF is found
        (equivalent to bilinear interpolation since stars are defined
        on square grid). The resulting effective PSF is added to the 
        internal _psf_array to be used in the get_tile pipeline.

        ...
        
        Parameters
        ----------
        star_pos : np.ndarray : position [arcsec]
        
        """
        
        if self.static == True:
            # Use static PSF
            self._psf_array += self.psfs[0].fft_data
        else:
            # Perform convex combination of PSFs based on star position
            for psf in self.psfs:
                gamma_check = \
                    (np.abs(psf.xpos-s_pos[0]) < self.psf_pitch) and \
                    (np.abs(psf.ypos-s_pos[1]) < self.psf_pitch)
                if gamma_check:
                    # PSF is a vertex of convex hull around star.
                    # Compute the PSFs convex combination coefficient:
                    gamma = (
                                (1-np.abs(psf.xpos-s_pos[0])/self.psf_pitch) * \
                                (1-np.abs(psf.ypos-s_pos[1])/self.psf_pitch)
                             ).astype(psf.fft_data.dtype)
                else:
                    continue
                # add result to PSF array
                self._psf_array += gamma*psf.fft_data
    
    def get_tile(self,index):
        """Get the tile corresponding to source[index]


        From the tile_generator object tgen, calling tgen.get_tile(index) will
        generate the tile corresponding to the tgen.source_pos[index] star by 
        interpolating the 4 neighbouring PSFs and convolving this effective
        PSF with a Gaussian kernel defined by tgen.source_cov[index].

        The output of this is a tile which has been trimmed down to the input
        PSF dimensions, as well as the coordinates of the bottom-left-corner
        of the tile so that it may be sliced into the final image properly.
        
        ...

        Parameters
        ----------
        index : int
            index of star in source table to generate tile for

        Returns
        -------
        out : real ndarray
            tile to be sliced into final image
        bottom_left_corner : ndarray
            coordinates of bottom left corner of bottom left pixel in arcsec
        timing : dict
            dictionary with timing for various parts of program
        """
        
        if self._init == True:
            # Initialise internal variables first time around
            
            # new shape is ~ half full tile size due to RFFT2 optimisation
            new_shape = [self.fourier_tile_dim[0],self.fourier_tile_dim[1]//2+1]
            self._psf_array = np.zeros(new_shape,dtype=self.dtype)

            # Coordinates (uu,vv) of DFT space:
            self._nx = self.fourier_tile_dim[0] # square image only
            self._T  = self.dtype(self.gauss_width_as + self.psf_width_as)
            self._Ts = self.dtype(self._T/self._nx)
            uu,vv = np.meshgrid(
                np.linspace(0,1/self._Ts,self._nx+1,dtype=self.dtype)[:-1]-1/(2*self._Ts),
                np.linspace(0,1/self._Ts,self._nx+1,dtype=self.dtype)[:-1]-1/(2*self._Ts)
                )
            uu = np.fft.fftshift(uu)
            vv = np.fft.fftshift(vv)
            self._fft_pos = np.c_[uu[:new_shape[0],:new_shape[1]].flatten(),
                                 vv[:new_shape[0],:new_shape[1]].flatten()].T
            
            # Run einsum calcs once to determine optimal sequence:
            self.optimize_star_kernel()
            
            self._init = False
        
        # Pick star from source table
        star_flux = self.source_flux[index]
        star_pos  = self.source_pos[index]
        star_cov  = self.source_cov[index]
        if self.dtype==np.complex64:
            star_flux = np.float32(star_flux)
            star_pos  = star_pos.astype(np.float32)
            star_cov  = star_cov.astype(np.float32)

        # Clear internal array:
        self._psf_array *= 0.0

        # Compute effective PSF:
        self.get_effective_psf_fft(star_pos)

        # Prepare for FFT Gaussian computation:
        offset = (((((star_pos % self.pixsize)/self.pixsize)+0.5)%1)-1)*self.pixsize # this has to be easier
        _star_pos = (self.psf_width_as+self.gauss_width_as)/2 * np.array([1.0,1.0]) + offset
        bottom_left_corner = star_pos-offset-self.psf_width_as/2 - 0.5*self.pixsize

        # Compute star Gaussian and convolve with PSF:
        self._psf_array *= self.get_star_kernel_fft(star_flux, star_cov, _star_pos)
        
        # Inverse FFT and trimming:
        out = (np.fft.fftshift(
            np.fft.irfft2(
                (self._psf_array)
            ).astype(self.dtype))).real[:self.psf_width_pix,:self.psf_width_pix]

        return out, bottom_left_corner    

    def get_star_kernel_fft(self, flux, cov, mu):
        """Compute star Gaussian based in DFT space.

        Directly computes the FFT of the Gaussian kernel with
        appriate amplitude, width, and offset to suit the tile
        being generated.

        Uses optimised np.einsum so requires running:
            optimize_star_kernel()
        first.

        FFT description of shifted Gaussian given by:
            TBD
        ...

        Inputs
        ------
        flux : float
            flux of star
        cov : np.ndarray
            covariance matrix of star Gaussian
        mu : np.ndarray
            position of star

        Output
        ------
        gaussian_fft : ndarray


        """
        offset = (2*np.pi*1j)*mu
        gaussian_fft = (flux*(self._nx**2/self._T**2)*2*np.pi*np.sqrt(np.linalg.det(cov)))*np.exp(
                (-2*(np.pi)**2)*np.einsum(
                    "ij,ii,ij->j",
                    self._fft_pos,
                    cov,
                    self._fft_pos,
                    optimize=self._esp1
                ) - np.einsum(
                    "ij,i->j",
                    self._fft_pos,
                    offset,
                    optimize=self._esp2
                )).reshape(self._psf_array.shape) 
        # TODO brute force normalise. Was an error of about 0.2% last I checked.
        return gaussian_fft
    
    def optimize_star_kernel(self):
        """Runs star kernel once to optimise np.einsum"""

        cov = self.source_cov[0]
        mu  = self.source_pos[0]
        offset = 2*np.pi*1j*mu.flatten()
        self._esp1 = np.einsum_path(
                    "ij,ii,ij->j",
                    self._fft_pos,
                    cov,
                    self._fft_pos,
                    optimize="optimal"
                )[0]
        self._esp2 = np.einsum_path(
                    "ij,i->j",
                    self._fft_pos,
                    -offset,
                    optimize="optimal"
                )[0]
    

class ImageGenerator:
    """
    Generate image from sliced tiles, one tile per source object.

    Usage:


    ...

    Attributes
    ----------
    pixsize : float
        pixel size in arcsec of image
    fov : float
        FoV of full image
    full_image : real ndarray
        final image at original pixel size (i.e., before rebinning)
    tile_gen : TileGenerator
        tile generator object used to create tiles to slice into final image


    Methods
    -------
    main()
        generates the full image and saves it to attribute self.full_image
    get_rebinned_cropped(rebin_factor, cropped_width_as)
        rebin the image in self.full_image by a rebin factor (rebin_factor), and 
        crop to a specified FoV (cropped_width_as)
    rebin(arr, new_shape)
        rebin an array (arr) into new_shape. new_shape must be achievable
        by an integer rebinning.
    get_source_list()
        getter function for source list as seen by tile generator
    """
    def __init__(self, array_width_pix, pixsize, source, psfs_file, gauss_width_pix,
            which_psf = None):        
        """Contructor for ImageGenerator object

        ...

        Parameters
        ----------
        array_width_pix : int
            width of full image in pixels before rebinning
        pixsize : float
            pixel size in arcsec before rebinning
        source_list : table
            source list table from Source.py format
        psfs_file : str
            filename for fits file containing all PSFs and PSF metadata
        gauss_width_pix : int
            width of Gaussian square support in pixels
        which_psf : int, optional
            default = None -> use variable PSF over field
            otherwise -> use fits HDU[which_psf+1] PSF only
        """

        self.pixsize = pixsize
        self.xx = np.arange(array_width_pix)*pixsize
        self.fov = self.xx[-1]-self.xx[0]
        self.xx -= (self.fov/2+self.pixsize/2)
        self.full_image = np.zeros([self.xx.shape[0],self.xx.shape[0]])
        self.tile_gen = TileGenerator(source, psfs_file, gauss_width_pix, which_psf=which_psf)
        self.nsource = source.flux.shape[0]
        
    def main(self):
        """Loop over all stars and add the tile to the full image."""
        for ni in tqdm(range(self.nsource),leave=False):
            # Generate the tile:
            tile,origin = self.tile_gen.get_tile(ni)
            # Find the location of the tile in the full image:
            xstart = np.abs(self.xx-origin[0]).argmin()
            ystart = np.abs(self.xx-origin[1]).argmin()
            # Slice the tile in:
            self.full_image[ystart:ystart+self.tile_gen.psf_width_pix,
                            xstart:xstart+self.tile_gen.psf_width_pix] += tile
        
    def get_rebinned_cropped(self,rebin_factor,cropped_width_as):
        """Rebin self.full_image after cropping to desired rebin factor."""

        xx_cropped_id = np.abs(self.xx)<=cropped_width_as/2
        cropped_im = self.full_image[xx_cropped_id,:][:,xx_cropped_id]
        rebinned_im = self.rebin(cropped_im,np.array(cropped_im.shape)//rebin_factor)
        return rebinned_im

    def rebin(self, arr, new_shape):
        """Rebin array (arr) into new shape (new_shape)."""
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
                 new_shape[1], arr.shape[1] // new_shape[1])
        return arr.reshape(shape).mean(-1).mean(1)
