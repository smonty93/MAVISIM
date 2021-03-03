# Standard
import numpy as np
import matplotlib.pyplot as plt; plt.ion()
import sys
import os

# Astropy
from astropy.io import fits, ascii
from astropy.table import Table

from src.mavisim.generate_image import ImageGenerator

# Project-specific/
#sys.path.append('/src')
cwd = os.getcwd()
os.chdir(os.getcwd() + "/src")

#import mavisim.rampup
#import mavisim.addconstantsky
#import mavissimim.addnoise
#import mavissimim.trimimage
import mavisim.input_parameters as input_par
#from rebin import rebin

from mavisim.Source import Source
#from mavissimim.AOGaussGrid import AOGaussGrid
#from mavissimim.SeeingGrid import SeeingGrid
from mavisim.InputCoo import InputCoo

os.chdir(cwd)
# Nbody input plus chosen exposure time
ngc3201 = input_par.nbody_in

exp_time = 3600
source = Source(ngc3201, exp_time, tt_amp_fac=1.0, static_dist=False, stat_amp=1.0).main()
#self, mavis_src, exp_time, static_dist, tt_amp_fac, stat_amp

image_gen = ImageGenerator(source[:100], "src/e2epsfs/e2e_psfs.fits", 
        [-15,15,-15,15], debug_scale_factor=16)

"""
#debug only
source = source[:10]

#debug only
debug_scale_factor = 4

dtype = np.complex128
psf_fits = fits.open("src/e2epsfs/e2e_psfs.fits")
psfs = []
for psf in psf_fits[1:]:
    tmp_metadata = {
            "xpos"    : psf.header["YPOS"]/debug_scale_factor,
            "ypos"    : psf.header["XPOS"]/debug_scale_factor,
            "pixsize" : psf.header["PIXSIZE"],
            "Lambda"  : psf.header["LAMBDA"]
        }
    tmp_data = psf.data
    psfs.append(PSF(tmp_data,tmp_metadata))

#debug only
for ni in range(len(source["Gauss_Info"])):
    source["Gauss_Info"][ni][1] /= debug_scale_factor

for ni in range(len(source["Gauss_Info"])):
    source["Gauss_Info"][ni][0] = dtype(source["Gauss_Info"][ni][0])

pad_for_psf = psfs[0].pixsize*psfs[0].data.shape[0]

psf_pitch = 30/10/debug_scale_factor
target_image_dim = np.array([8000//debug_scale_factor,8000//debug_scale_factor])
image_dim = target_image_dim + np.array(psfs[0].data.shape)
compute_psf_ffts(psfs, image_dim, dtype=dtype)

extent = np.array([psfs[0].xpos,psfs[-1].xpos+pad_for_psf,psfs[0].ypos,psfs[-1].ypos+pad_for_psf])

image = get_image(source["Gauss_Info"], psfs, psf_pitch, extent, image_dim, dtype=dtype)
image_cropped = (np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(image)).astype(dtype))).real[
        psfs[0].data.shape[0]//2:-psfs[0].data.shape[0]//2,
        psfs[0].data.shape[0]//2:-psfs[0].data.shape[0]//2
        ]

def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

image_resampled = rebin(image_cropped,target_image_dim//2)
"""
