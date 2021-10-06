from astropy.io import ascii
from astropy import units as u
import numpy as np

start_wavelength =  550 # nm, start of wavelength slice (=stop wavelength if monochromatic)
stop_wavelength =  550 # nm, end of wavelength slice
filt_width = 88        # nm, width of V band
collecting_area = (np.pi * 4.0**2 * (1 - 0.16**2)) * u.m**2
MAVIS_fov = 30 # arcseconds
buffer = 10.005 #buffer in arcseconds to add to the MAVIS fov to account for stars outside the science fov

# Input Catalogue info (if relevant)
# TODO: move to data dir
input_cat = ascii.read("data/example/ngc3201_mavisim_input0year_imbh_4000stars2_radecFixedFlux")
nbody_yr = 0

# CCD info
ccd_size =  4000 # Change this to 4096 for FFT speed
ccd_sampling = 0.0075 # Pixel/arcsec, should be the same as the sampled FV PSF prior to convolution
psf_sampling = 0.0075 # Upsampling the PSF to avoid ringing, rebinning happens in code
psf_wavelength = 550

# Gaussian Point Source Info
vib_term = 0.48 #0.48 pixels, FWHM = 8.225 mas extra vibration, accounting for the 10mas min FWHM
cd_term = 0.4 #0.4 pixels, FWHM = 7.05 pixels to account for charge diffusion, base jitter in the case of a static jitter kernal
gauss_offset = -0.5 # pixels, offset to make sure the rebinning maintains a centred PSF

# Static Distortion
static_distort =  ascii.read("data/StaticDistortTransmissive")
plate_scale =  1.29 # arcsec/mm, to convert location in mm to arcseconds
dynamic_amp = 1 # amplification factor to increase the distortion

# Detector Characteristics
slow_rdnoise =  3 #electrons, best
fast_rdnoise =  5 #electrons, worst
bit_depth = 16	  #standard
sat_point = 0 	  #could specificy this instead
gain = 1.0        #taken fron Francois ADU/e
noise_mode = "slow"
AOMthruput = 0.63 # for 550 nm --> generalise this!
VLTthruput = 0.75 # approx flat across all wavelengths
QE = 0.89 		   # for 550nm  --> generalise this

# Sky Background
surf_bright = 21.61   # mag/square arcsecond
# TODO: move this to data dir
sky_bright =  ascii.read("data/sky_emission.dat")
sky_bright["lambda_nm"].unit = u.nm
sky_bright["flux"].unit = u.s**(-1) * u.m**(-2) * u.um**(-1) * u.arcsec**(-2)
airmass =  1.155
