import numpy as np
import time
from astropy.io import fits
from copy import deepcopy

class PSF:
    """PSF object class"""
    def __init__(self, fits_ext, padto, *, dtype=np.complex64,
            debug_scale_factor=1):
        self.fft_data = np.fft.fftshift(np.fft.fft2(fits_ext.data,s=padto)).astype(dtype)
        self.xpos     = fits_ext.header["YPOS"]/debug_scale_factor
        self.ypos     = fits_ext.header["XPOS"]/debug_scale_factor
        self.Lambda   = fits_ext.header["LAMBDA"]
        print(f"x: {self.xpos:7.3f}\", y: {self.ypos:7.3f}\"",end="\r")
    
class ImageGenerator:
    def __init__( self, source_list, psfs_file, image_extent, *,
            debug_scale_factor: int = 1, dtype = np.complex128 ):
        
        self.dtype=dtype

        # Parse source info:
        self.source = deepcopy(source_list["Gauss_Info"])
        for ni in range(len(self.source)):
            self.source[ni][1] /= debug_scale_factor
        for ni in range(len(self.source)):
            self.source[ni][0] = dtype(self.source[ni][0])
            self.source[ni][1] = dtype(self.source[ni][1])
            self.source[ni][2] = dtype(self.source[ni][2])
        
        # Define nominal image geometry:
        self.pixsize               = 3.75e-3 # TODO read this from fits Primary HDU metadata
        self.final_image_width_as  = (image_extent[1]-image_extent[0])/debug_scale_factor
        self.final_image_width_pix = int(round(self.final_image_width_as/self.pixsize))
        self.final_image_dim       = np.array([
                self.final_image_width_pix,
                self.final_image_width_pix
            ])

        # Load PSF file and define Fourier Image Geometry:
        self.psf_file          = psfs_file
        psfs_fits              = fits.open(psfs_file, lazy_load_hdus=True)
        self.psf_pitch         = 3.0 # TODO read this from fits Primary HDU metadata
        self.psf_pitch = self.psf_pitch/debug_scale_factor
        self.psf_width_pix     = psfs_fits[1].data.shape[0]
        self.psf_width_as      = self.psf_width_pix * self.pixsize
        self.fourier_image_dim = self.final_image_dim + np.array(psfs_fits[1].data.shape)
        
        # Create PSF objects and do their FFTs
        print("Doing PSF FFTs:")
        self.psfs = []
        for psf in psfs_fits[1:]:
            self.psfs.append(PSF(psf, self.fourier_image_dim, debug_scale_factor=debug_scale_factor))
        psfs_fits.close()
        print("Done.                          ")
        
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
    
    def get_image(self):
        image_fft = np.zeros(self.fourier_image_dim,dtype=self.dtype)
        psf_origin = self.dtype(np.array([1.0,1.0]) * (self.psf_width_as/2))
        storage_array = image_fft.copy()
        time_clearing = 0
        time_star = 0
        time_psf = 0
        time_product = 0
        nx = self.fourier_image_dim[0] # square image only
        T  = self.dtype(self.final_image_width_as + self.psf_width_as)
        Ts = self.dtype(T/nx)
        uu,vv = np.meshgrid(
                np.linspace(0,1/Ts,nx+1,dtype=self.dtype)[:-1]-1/(2*Ts),
                np.linspace(0,1/Ts,nx+1,dtype=self.dtype)[:-1]-1/(2*Ts)
                )
        fft_pos = np.c_[uu.flatten(),vv.flatten()].T
        esp1,esp2 = self.optimize_star_kernel(self.source[0], self.dtype(fft_pos), nx, T)
        print("Building image from sources:")
        for star in self.source:
            print(f"flux: {star[0].real:5.3e}, x: {star[1][0].real:7.3f}\", y: {star[1][1].real:7.3f}\"", end="\r")
            t1 = time.time()
            storage_array *= 0.0
            t2 = time.time()
            x2 = self.get_effective_psf_fft(star, storage_array)
            t3 = time.time()
            star_pos = star[1] - psf_origin
            x1 = star[0] * self.get_star_kernel_fft(star[2], star_pos, fft_pos, nx, T, esp1, esp2)
            t4 = time.time()
            image_fft += x1 * x2
            t5 = time.time()
    
            time_clearing += t2-t1
            time_psf      += t3-t2
            time_star     += t4-t3
            time_product  += t5-t4
        print("Done.                                                 ")
        print(f"time clearing : {time_clearing:f}")
        print(f"time star     : {time_star:f}")
        print(f"time psf      : {time_psf:f}")
        print(f"time product  : {time_product:f}")
        
        return image_fft
    
    def get_star_kernel_fft(self, cov, mu, fft_pos, nx, T, esp1, esp2):
        offset = 2*np.pi*1j*mu.flatten()
        gaussian_fft = 2*np.pi*np.sqrt(np.linalg.det(cov))*np.exp(
                -2*(np.pi)**2*np.einsum(
                    "ij,ii,ij->j",
                    fft_pos,
                    cov,
                    fft_pos,
                    optimize=esp1
                ) + np.einsum(
                    "ij,i->j",
                    fft_pos,
                    -offset,
                    optimize=esp2
                )).reshape([nx,nx])
        gaussian_fft *= nx**2/T**2
        # TODO brute force normalise. Was an error of about 0.2% last I checked.
        return gaussian_fft
    
    def get_star_kernel_fft_new(self, cov, mu, fft_pos, nx, T, esp1, esp2):
        offset = 2*np.pi*1j*mu.flatten()
        gaussian_fft = 2*np.pi*np.sqrt(np.linalg.det(cov))*np.exp(
                np.array([
                    -2*(np.pi)**2 * x.T @ cov @ x - x.T @ offset for x in fft_pos.T
                    ])).reshape([nx,nx])
        gaussian_fft *= nx**2/T**2
        # TODO brute force normalise. Was an error of about 0.2% last I checked.
        return gaussian_fft
    
    def optimize_star_kernel(self, star, fft_pos, nx, T):
        cov = star[2]
        mu = star[1]
        offset = 2*np.pi*1j*mu.flatten()
        esp1 = np.einsum_path(
                    "ij,ii,ij->j",
                    fft_pos,
                    cov,
                    fft_pos,
                    optimize="optimal"
                )
        esp2 = np.einsum_path(
                    "ij,i->j",
                    fft_pos,
                    -offset,
                    optimize="optimal"
                )
        return esp1[0],esp2[0]
    
    def main(self): 
        image = self.get_image()
        image_cropped = (np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(image)).astype(self.dtype))).real[
                self.psf_width_pix//2:-self.psf_width_pix//2,
                self.psf_width_pix//2:-self.psf_width_pix//2
                ]
        
        def rebin(arr, new_shape):
            shape = (new_shape[0], arr.shape[0] // new_shape[0],
                     new_shape[1], arr.shape[1] // new_shape[1])
            return arr.reshape(shape).mean(-1).mean(1)
        
        self.image_resampled = rebin(image_cropped,self.final_image_dim//2)
    
    
