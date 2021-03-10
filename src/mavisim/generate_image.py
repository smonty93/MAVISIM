import numpy as np
import pyfftw
import time
from astropy.io import fits
from copy import deepcopy

class PSF:
    """PSF object class"""
    def __init__(self, fits_ext, padto, *, dtype=np.complex128):
        self.fft_data = (np.fft.rfft2(fits_ext.data,s=padto)).astype(dtype)
        self.xpos     = fits_ext.header["YPOS"]
        self.ypos     = fits_ext.header["XPOS"]
        self.Lambda   = fits_ext.header["LAMBDA"]
        print(f"x: {self.xpos:7.3f}\", y: {self.ypos:7.3f}\"",end="\r")
    
class TileGenerator:
    def __init__( self, source_list, psfs_file, gauss_width_pix, *,
             dtype = np.complex128 ):
        
        self.dtype=dtype

        # Parse source info:
        self.source = deepcopy(source_list["Gauss_Info"])
        
        # Define nominal image geometry:
        self.pixsize               = 3.75e-3 # TODO read this from fits Primary HDU metadata
        self.gauss_width_pix  = gauss_width_pix
        self.gauss_width_as   = self.gauss_width_pix*self.pixsize
        self.gauss_dim        = np.array([self.gauss_width_pix,self.gauss_width_pix])

        # Load PSF file and define Fourier Image Geometry:
        self.psf_file          = psfs_file
        psfs_fits              = fits.open(psfs_file, lazy_load_hdus=True)
        self.psf_pitch         = 3.0 # TODO read this from fits Primary HDU metadata
        self.psf_width_pix     = psfs_fits[1].data.shape[0]
        self.psf_width_as      = self.psf_width_pix * self.pixsize
        
        self.fourier_tile_dim = self.gauss_dim + np.array(psfs_fits[1].data.shape)
        
        # Create PSF objects and do their FFTs
        print("Doing PSF FFTs:")
        self.psfs = []
        for psf in psfs_fits[1:]:
            self.psfs.append(PSF(psf, self.fourier_tile_dim))
        psfs_fits.close()
        print("Done.                          ")

        self.init = True
        
    def get_effective_psf_fft(self, star, psf_out):
        # Takes in a single star object (list of [flux, pos, cov])
        # and all PSFs, and returns the effective PSF of that star.
        s_pos = star[1]
        for psf in self.psfs:
            gamma_check = \
                (np.abs(psf.xpos-s_pos[0]) < self.psf_pitch) and \
                (np.abs(psf.ypos-s_pos[1]) < self.psf_pitch)
            if gamma_check:
                gamma = (
                            (1-np.abs(psf.xpos-s_pos[0])/self.psf_pitch) * \
                            (1-np.abs(psf.ypos-s_pos[1])/self.psf_pitch)
                         ).astype(psf.fft_data.dtype)
            else:
                continue
            psf_out += gamma*psf.fft_data
            # print(f"{psf.xpos:0.5f} {psf.ypos:0.5f} {gamma.real:0.4f}")
        return psf_out
    
    def get_tile(self,index):
        t1 = time.time()
        if self.init == True:
            new_shape = [self.fourier_tile_dim[0],self.fourier_tile_dim[1]//2+1]
            self.psf_array = np.zeros(new_shape,dtype=self.dtype)
            self.nx = self.fourier_tile_dim[0] # square image only
            self.T  = self.dtype(self.gauss_width_as + self.psf_width_as)
            self.Ts = self.dtype(self.T/self.nx)
            uu,vv = np.meshgrid(
                np.linspace(0,1/self.Ts,self.nx+1,dtype=self.dtype)[:-1]-1/(2*self.Ts),
                np.linspace(0,1/self.Ts,self.nx+1,dtype=self.dtype)[:-1]-1/(2*self.Ts)
                )
            uu = np.fft.fftshift(uu)
            vv = np.fft.fftshift(vv)
            #self.fft_pos = np.c_[uu.flatten(),vv.flatten()].T
            self.fft_pos = np.c_[uu[:new_shape[0],:new_shape[1]].flatten(),
                                 vv[:new_shape[0],:new_shape[1]].flatten()].T
            self.optimize_star_kernel(self.source[0])
            self.init = False
        t2 = time.time()
        star = self.source[index]
        if self.dtype==np.complex64:
            star = [x.astype(np.float32) for x in star]
        self.psf_array *= 0.0
        t3 = time.time()
        self.get_effective_psf_fft(star, self.psf_array)
        t4 = time.time()
        offset = (((((star[1] % self.pixsize)/self.pixsize)+0.5)%1)-1)*self.pixsize # this has to be easier
        star_pos = (self.psf_width_as+self.gauss_width_as)/2 * np.array([1.0,1.0]) + offset
        bottom_left_corner = star[1]-offset-self.psf_width_as/2 - 0.5*self.pixsize
        t5 = time.time()
        self.psf_array *= star[0] * self.get_star_kernel_fft(star[2], star_pos)
        t6 = time.time()
        out = (np.fft.fftshift(
            np.fft.irfft2(
                (self.psf_array)
            ).astype(self.dtype))).real[:self.psf_width_pix,:self.psf_width_pix]
        t7 = time.time()

        return out, bottom_left_corner, {
            #"init":t2-t1,
            "clear":t3-t2,
            "psf":t4-t3,
            "prep":t5-t4,
            "star":t6-t5,
            "fft":t7-t6,
        }
    
    def get_star_kernel_fft(self, cov, mu):
        offset = 2*np.pi*1j*mu.flatten()
        gaussian_fft = 2*np.pi*np.sqrt(np.linalg.det(cov))*np.exp(
                -2*(np.pi)**2*np.einsum(
                    "ij,ii,ij->j",
                    self.fft_pos,
                    cov,
                    self.fft_pos,
                    optimize=self.esp1
                ) + np.einsum(
                    "ij,i->j",
                    self.fft_pos,
                    -offset,
                    optimize=self.esp2
                )).reshape(self.psf_array.shape)
        gaussian_fft *= self.nx**2/self.T**2
        # TODO brute force normalise. Was an error of about 0.2% last I checked.
        return gaussian_fft
    
    def optimize_star_kernel(self, star):
        cov = star[2]
        mu = star[1]
        offset = 2*np.pi*1j*mu.flatten()
        self.esp1 = np.einsum_path(
                    "ij,ii,ij->j",
                    self.fft_pos,
                    cov,
                    self.fft_pos,
                    optimize="optimal"
                )[0]
        self.esp2 = np.einsum_path(
                    "ij,i->j",
                    self.fft_pos,
                    -offset,
                    optimize="optimal"
                )[0]
    
    def main(self): 
        image = self.get_tile()
        def rebin(arr, new_shape):
            shape = (new_shape[0], arr.shape[0] // new_shape[0],
                     new_shape[1], arr.shape[1] // new_shape[1])
            return arr.reshape(shape).mean(-1).mean(1)
        
        self.image_resampled = rebin(image,np.array(image.shape)//2)
    
class ImageGenerator:
    def __init__(self, array_width_pix, pixsize, source_list, psfs_file, gauss_width_pix):        
        self.pixsize = pixsize
        self.xx = np.arange(array_width_pix)*pixsize
        self.fov = self.xx[-1]-self.xx[0]
        self.xx -= (self.fov/2+self.pixsize/2)
        self.full_image = np.zeros([self.xx.shape[0],self.xx.shape[0]])
        self.tile_gen = TileGenerator(source_list, psfs_file, gauss_width_pix)
        
    def main(self):
        time_profile = None
        print(" | ".join([f"{x:7s}" for x in self.tile_gen.get_tile(0)[2]]) + " ||  ETA")
        for ni in range(len(self.tile_gen.source)):
            tile = self.tile_gen.get_tile(ni)
            xstart = np.abs(self.xx-tile[1][0]).argmin()
            ystart = np.abs(self.xx-tile[1][1]).argmin()
            self.full_image[ystart:ystart+self.tile_gen.psf_width_pix,
                            xstart:xstart+self.tile_gen.psf_width_pix] += tile[0]
            if time_profile is not None:
                time_profile = {x:time_profile[x]+tile[2][x] for x in time_profile}
            else:
                time_profile = tile[2]
            eta = (len(self.tile_gen.source)-ni)*(sum(time_profile.values())/(ni+1))
            print(" | ".join([f"{time_profile[x]/(ni+1):7.5f}" for x in time_profile]) + \
                    f" || {eta:0.2f}          ",end="\r")
        print("\ndone")
        print("Totals:")
        print(" | ".join([f"{x:7s}" for x in time_profile])+" || TOTAL")
        total = sum(time_profile.values())
        print(" | ".join([f"{time_profile[x]:7.5f}" for x in time_profile]) + f" || {total:0.2f}")
        
    def get_rebinned_cropped(self,rebin_factor,cropped_width_as):
        xx_cropped_id = np.abs(self.xx)<=cropped_width_as/2
        xx_cropped = self.xx[xx_cropped_id]
        cropped_im = self.full_image[xx_cropped_id,:][:,xx_cropped_id]
        rebinned_im = self.rebin(cropped_im,np.array(cropped_im.shape)//2)
        return rebinned_im

    def rebin(self, arr, new_shape):
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
                 new_shape[1], arr.shape[1] // new_shape[1])
        return arr.reshape(shape).mean(-1).mean(1)
