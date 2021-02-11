from astropy.io import fits
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
nbody_in = ascii.read("../src/mavissimim/ngc3201_mavissimim_input0year_imbh_4000stars2_radecFixedFlux")
nbody_yr = 0

# CCD info
ccd_size =  4000 # Change this to 4096 for FFT speed
ccd_sampling = 0.0075 # Pixel/arcsec, should be the same as the sampled FV PSF prior to convolution
psf_sampling = 0.0075 # Upsampling the PSF to avoid ringing, rebinning happens in code

# Gaussian Point Source Info
vib_term = 0.48 #0.48 pixels, FWHM = 7.05 mas extra vibration, accounting for the 10mas min FWHM
gauss_sigma = 0.62 #0.4 pixels, FWHM = 7.05 pixels to account for charge diffusion, base jitter in the case of a static jitter kernal
gauss_width = 34 # pixels, 11 x 11 array to store gaussian with central pixel = xcog, ycog
gauss_wing = 17 # pixels, size of the wings of the gaussian (centred at xcog, ycog) extending to gauss_width total
gauss_offset = -0.5 # pixels, offset to make sure the rebinning maintains a centred PSF

# Arcetri PSF + weights to help blend the wings
big_psf =  fits.open("../src/mavissimim/PSF_0_0dir_arcsec_40arcsec_75mas_550nm.fits")
big_psf_ramp = 3.5 # pixels_
weight_grid = ascii.read("../src/mavissimim/Seeing_Weights_Aug2020")

# TT Jitter
tt_jitter =  fits.open("../src/mavissimim/TTJitterMaps/TT_jitter_mavis_astrad10_10_10_angles0_120_240_mag15.0_15.0_15.0_121dirs_sr80.0_80.0_80.0_at1650nm.fits")
tt_theta_conv = "Degrees"
shift_tech = "Matrix"

# Static Distortion
static_distort =  ascii.read("../src/mavissimim/StaticDistortTransmissive")
plate_scale =  1.29 # arcsec/mm, to convert location in mm to arcseconds
dynamic_amp = 1 # amplification factor to increase the distortion

# Field Variable PSF Information
# -------------------------------------------------------------------------------------------------
# PSFs created using 0.8" seeing profile and configuration 4 that includes the most recent updates
# -------------------------------------------------------------------------------------------------
fv_psf_path = "../src/mavissimim/PSF_Grid_1ArcSecFoV_75masSampling_Jan2020Code_NoTT/"
static_psf = fits.open("../src/mavissimim/PSF_Grid_1ArcSecFoV_75masSampling_Jan2020Code_NoTT/PSF_0_0dir_arcsec_1arcsec_550nm.fits")
fv_psf_filename = "PSF_%s_%sdir_arcsec_1arcsec_550nm.fits" # Generic name of FV PSF file for
fv_psf_grid = np.array([-15, -12, -9, -6, -3, 0, 3, 6, 9, 12, 15]) # array of positions where the field variable PSFs are sampled (in arcseconds from centre of FoV)
theory_core_rad = (40 * start_wavelength*1e-9)/(2 * 8.0) # Guess at where the theoretical control radius should be, assuming MAVIS has 40 subapertures (d) and ctrl_rad = (d * lambda)/2D - Monochromatic for now!

# Handle the cases of different sampling rates for the big PSF and the small PSFS
psf_core_rad_pix  = int((theory_core_rad  * 206265.0)/psf_sampling)
seeing_core_rad_pix  = int((theory_core_rad  * 206265.0)/ccd_sampling)
ramp_size = 2 # ramp size in pixels to blend the FV PSF into the seeing wings

# Rebinning for PSFs sampled at higher spatial sampling (NECESSARY)
rebin_size = int(1/ccd_sampling)

if rebin_size%2 != 0:
	rebin_size = rebin_size + 1
                
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
sky_bright =  ascii.read("../src/mavissimim/sky_emission.dat")
sky_bright["lambda_nm"].unit = u.nm
sky_bright["flux"].unit = u.s**(-1) * u.m**(-2) * u.um**(-1) * u.arcsec**(-2)
airmass =  1.155

debug_file = "/home/montys/AO/MAVISSimIm/AstrometryPaper/SimI/debug"
