# ----------------------------------------------------------------------------
#
# TITLE - add constant sky
# AUTHOR - Stephanie Monty
# PROJECT - MAVISSimIm
# CONTENTS:
#   1. Fill this in
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Function to add a constant sky background to the images

ONLY MONOCHROMATIC FOR NOW --> Update this
'''
__author__ = "Stephanie Monty"

### Imports

## Basic
import numpy as np

## Astropy
from astropy import units as u

## Project Specific
import mavissimim.input_parameters as input_par


def add_constant_sky_pixel(exp_time):
	"""
	Add a constant sky to the image

	Result is photons/voxel
	"""
	# Mean wavelength of the study
	mean_wavelength = (input_par.start_wavelength+input_par.stop_wavelength)/2

	square_arcsec_pervoxel = (input_par.ccd_sampling**2) * u.arcsec**2


	# If the surface brightness is passed in as mag/arcsec^2 do the following:
	if input_par.surf_bright != 0:
		flux_jy = (3631 * 10**(-1*input_par.surf_bright/2.5))*u.Jy*u.nm**(-1)

		flux_ph_s_nm_cm2 = flux_jy * (1.51e3/mean_wavelength)

		flux_ph_s_nm_m2 = flux_ph_s_nm_cm2 * 10**4

		flux_ph_s_m2 = flux_ph_s_nm_m2 * input_par.filt_width #in nm

		sky_value = flux_ph_s_m2 * square_arcsec_pervoxel  * exp_time * input_par.collecting_area

	else:
		start_index = np.where(input_par.sky_bright["lambda_nm"] == input_par.start_wavelength)[0][0]
		stop_index = np.where(input_par.sky_bright["lambda_nm"] == input_par.stop_wavelength)[0][0]

		# Slice the sky emission spectrum using our start and stop wavelengths
		
		# Code for monochromatic case
		if start_index == stop_index:
			sky_slice = input_par.sky_bright[start_index]
			lambda_slice = filt_width * 1.e-3 * u.um

		else:
			sky_slice = input_par.sky_bright[start_index:stop_index]
			lambda_slice = (input_par.stop_wavelength - input_par.start_wavelength) * filt_width * 1.e-3 * u.um


		sky_value = np.sum(sky_slice["flux"] * square_arcsec_pervoxel  * exp_time * input_par.collecting_area * lambda_slice)
	
	return (sky_value.value)
