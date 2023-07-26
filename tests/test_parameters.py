from astropy.io import ascii
from astropy import units as u
import numpy as np

####
# Input Catalogue info. BYO Catalogue here:
input_cat = ascii.read("tests/test_cat")

####
# Gaussian Point Source Info. E.g., to investigate sensitivity to vibrations
vib_term = 0.48  # 0.48 pixels, FWHM = 8.225 mas extra vibration, accounting for the 10mas min FWHM
cd_term = 0.4  # 0.4 pixels, FWHM = 7.05 pixels to account for charge diffusion, base jitter in the case of a static jitter kernal
gauss_offset = 0.0  # pixels, offset to make sure the rebinning maintains a centred PSF

####
# Static Distortion. E.g., to investigate sensitivity to VLT+MAVIS distortions
plate_scale = 1.29  # arcsec/mm, to convert location in mm to arcseconds
dynamic_amp = 1  # amplification factor to increase the distortion

####
# Sky Background
surf_bright = 21.61   # mag/square arcsecond

####
# MAVIS specific. Unlikely to be modified by the user:

# CCD info
ccd_size = 4000  # Change this to 4096 for FFT speed
ccd_sampling = 0.0075  # Pixel/arcsec, should be the same as the sampled FV PSF prior to convolution

# noise related parameters
psf_wavelength = 550   # nm (monochromatic version)
filt_width = 88        # nm, width of V band
collecting_area = (np.pi * 4.0**2 * (1 - 0.16**2)) * u.m**2

# FoV for trimming input catalogue
MAVIS_fov = 30  # arcseconds
buffer = 10.005  # buffer in arcseconds to add to the MAVIS fov to account for stars outside the science fov

# Detector Characteristics
slow_rdnoise = 3  # electrons, best
fast_rdnoise = 5  # electrons, worst
bit_depth = 16  # standard
sat_point = 0  # could specificy this instead
gain = 1.0  # taken fron Francois ADU/e
noise_mode = "slow"
AOMthruput = 0.63  # for 550 nm --> generalise this!
VLTthruput = 0.75  # approx flat across all wavelengths
QE = 0.89 		   # for 550nm  --> generalise this
