from astropy.io import fits
from astropy.io import ascii
from astropy import units as u
import numpy as np

# ------------------------------------------------------------------------------------------------------------------------#
# 												Free Parameters									  						  #
# ------------------------------------------------------------------------------------------------------------------------# 
path_to_mavisim = "/Users/stephanie/Dropbox/ANUResearch/AO/MAVISIM1/" # location of the mavisim directory
path_to_data = path_to_mavisim + "data/"	# location of the data files MAVISIM needs to run (the big PSF, PSF database, the TT jitter map database, the file of weights)

input_file = ascii.read("ngc3201_mavisim_input0year_imbh_4000stars2_radecFixedFlux") # input catalogue to simulate

filter = "V" 		 # specify the closest broadband filter to the monochromatic wavelength being studied (e.g. 550nm -> V)
psf_wavelength = 550 # nm, the wavelength of the PSF database of choice
fv_psf_path = path_to_data + "PSF_Grid_1ArcSecFoV_75masSampling_Jan2020Code_NoTT_550nm/"												 # Database of field variable PSFs (11 x 11 grid)								
static_psf = fits.open(fv_psf_path + "PSF_0_0dir_arcsec_1arcsec_550nm.fits") # Single PSF to use for static case, assumed to be best (central PSF)

tt_residual_map = "TT_jitter_mavis_astrad10_10_10_angles0_120_240_mag15.0_15.0_15.0_121dirs_sr80.0_80.0_80.0_at1650nm.fits" # tip-tilt residual map to use in the case of a field variable TT kernel
tt_kernel = 8.0 																											# generic tip-tilt kernel FWHM in mas if the user wants a static tt kernel not derived from the map

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

# Arcetri PSF + weights to help blend the wings
big_psf =  fits.open(path_to_data + "PSF_0_0dir_arcsec_40arcsec_75mas_550nm.fits") # big PSF to gather info for seeing wings
big_psf_ramp = 3.5 																   # pixels, to help ramp the seeing wings into the AO core
weight_grid = ascii.read(path_to_data + "Seeing_Weights_Aug2020")				   # weights to help with realistically ramping the PSF wings

# TT Jitter
tt_jitter =  fits.open(path_to_data + "TTJitterMaps/" + tt_residual_map)		   # TT residual map used to recover tip-tilt contribution from NGS
tt_theta_conv = "Degrees"
shift_tech = "Matrix"

# Static Distortion
static_distort =  ascii.read(path_to_data + "StaticDistortTransmissive")		   # Static distortion map reflecting the contribution of the MAVIS optics
plate_scale =  1.29 															   # arcsec/mm, to convert location in mm to arcseconds
dynamic_amp = 1 																   # NOT IN USE, amplification factor to increase the distortion

# Field Variable PSF Information
# -------------------------------------------------------------------------------------------------
# PSFs created using 0.8" seeing profile and configuration 4 that includes the most recent updates
# -------------------------------------------------------------------------------------------------
fv_psf_filename = "PSF_%s_%sdir_arcsec_1arcsec_" + str(psf_wavelength) + "nm.fits" 							# Generic name of FV PSF file for
fv_psf_grid = np.array([-15, -12, -9, -6, -3, 0, 3, 6, 9, 12, 15]) 					# Array of positions where the field variable PSFs are sampled (in arcseconds from centre of FoV at 0, 0)
theory_core_rad = (40 * psf_wavelength*1e-9)/(2 * 8.0) 								# Guess at where the theoretical control radius should be, assuming MAVIS has 40 subapertures (d) and ctrl_rad = (d * lambda)/2D - Monochromatic for now!

# Handle the cases of different sampling rates for the big PSF and the small PSFS
psf_core_rad_pix  = int((theory_core_rad  * 206265.0)/psf_sampling)
seeing_core_rad_pix  = int((theory_core_rad  * 206265.0)/ccd_sampling)
ramp_size = 2 															# ramp size in pixels to blend the FV PSF into the seeing wings
                
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
