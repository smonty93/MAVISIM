# ----------------------------------------------------------------------------
#
# TITLE - add constant sky
# AUTHOR - Stephanie Monty
# PROJECT - mavisim
# CONTENTS:
#   1. This function determines a constant sky background value by converting from a value in mag/arcsec^2 to a value in photons for each pixel
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Function to add a constant sky background to the images

ONLY MONOCHROMATIC FOR NOW
'''
__author__ = "Stephanie Monty"

### Imports

## Basic
import numpy as np

## Astropy
from astropy import units as u


def add_constant_sky_pixel(input_par, exp_time):
	"""
	Args:
		exp_time = exposure time in seconds to convert from photons/s to photons
		
	Returns:
		sky_value = a global sky value in photons to add to every pixel

	""" 

	square_arcsec_pervoxel = (input_par.ccd_sampling**2) * u.arcsec**2

	# Assuming the surface brightness is passed in as mag/arcsec^2 do the following:
	flux_jy = (3631 * 10**(-1*input_par.surf_bright/2.5))*u.Jy*u.nm**(-1)

	flux_ph_s_nm_cm2 = flux_jy * (1.51e3/input_par.psf_wavelength)

	flux_ph_s_nm_m2 = flux_ph_s_nm_cm2 * 10**4

	flux_ph_s_m2 = flux_ph_s_nm_m2 * input_par.filt_width #in nm

	sky_value = flux_ph_s_m2 * square_arcsec_pervoxel  * exp_time * input_par.collecting_area

	
	return (sky_value.value)
