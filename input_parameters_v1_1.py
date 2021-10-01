from astropy.io import fits
from astropy.io import ascii
from astropy import units as u
import numpy as np

# ------------------------------------------------------------------------------------------------------------------------#
# 												Free Parameters									  						  #
# ------------------------------------------------------------------------------------------------------------------------# 
input_file = ascii.read("fornax5_formatted_trim0") # input catalogue to simulate

filter = "V" 		 # specify the closest broadband filter to the monochromatic wavelength being studied (e.g. 550nm -> V)
psf_wavelength = 550 # nm, the wavelength of the PSF database of choice

noise_mode = "slow"	# detector mode (slow or fast readout)

# ------------------------------------------------------------------------------------------------------------------------#
# 									Hardcoded Characteristics (no need to change)										  #
# ------------------------------------------------------------------------------------------------------------------------#

# Telescope and instrument characteristics
collecting_area = (np.pi * 4.0**2 * (1 - 0.16**2)) * u.m**2 # VLT collecting area accounting for the central obscuration
MAVIS_fov = 30 		  # arcseconds
buffer = 10.005 	  #buffer in arcseconds to add to the MAVIS fov to account for stars outside the science fov

# CCD info
ccd_size =  4000 	  # Rounded size to reflect rounded pixel size
ccd_sampling = 0.0075 # Pixel/arcsec, should be the same as the sampled FV PSF prior to convolution
psf_sampling = 0.0075 # Upsampling the PSF to avoid ringing, rebinning happens in code

# Gaussian Point Source Info
vib_term = 0.48 	# FWHM = 7.05 mas extra vibration, accounting for the 10mas min FWHM
cd_term = 0.4 		# FWHM = 7.05 pixels to account for charge diffusion, base jitter in the case of a static jitter kernal
gauss_width = 34 	# pixels, 11 x 11 array to store gaussian with central pixel = xcog, ycog
gauss_wing = 17 	# pixels, size of the wings of the gaussian (centred at xcog, ycog) extending to gauss_width total
gauss_offset = -0.5 # pixels, offset to make sure the rebinning maintains a centred PSF

# Static Distortion
static_distort =  ascii.read("data/StaticDistortTransmissive")      # Static distortion map reflecting the contribution of the MAVIS optics
plate_scale =  1.29                                                 # arcsec/mm, to convert location in mm to arcseconds
dynamic_amp = 1                                                     # NOT IN USE, amplification factor to increase the distortion

# Detector Characteristics
slow_rdnoise =  3 	# electrons, best
fast_rdnoise =  5 	# electrons, worst
bit_depth = 16	  	# standard
sat_point = 0 	  	# could specificy this instead
gain = 1.0        	# taken fron Francois ADU/e
AOMthruput = 0.63 	# for 550 nm --> generalise this!
VLTthruput = 0.75 	# approx flat across all wavelengths
QE = 0.89 		   	# for 550nm  --> generalise this!

# Sky Background
filt_surf_bright = {"U": 22.28, "B": 22.64, "V": 21.61, "R": 20.87, "I": 19.71} # mag/arcsec2, Generalising for future use. Paranal sky brightness from Patat 2003
filt_widths = {"U": 66, "B": 94, "V": 88, "R": 138, "I": 149}					# nm, Filter widths
surf_bright = filt_surf_bright[filter]   										# mag/square arcsecond
filt_width = filt_widths[filter]												# nm, Filter widths
airmass =  1.155																# NOT USED, average airmass estimate for simulated observation, just for info
