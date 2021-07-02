# ----------------------------------------------------------------------------
#
# TITLE - BigPSF
# AUTHOR - Stephanie Monty
# PROJECT - mavisim
# CONTENTS:
#   1. This class returns a CCD-sized field of weighted seeing wings to create the final image
#
#		a) seeing_field = a CCD-sized array populated with the weighted seeing wings
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Generates the CCD-sized array of seeing wings using the big PSF and the weighted Gaussians

'''
__author__ = "Stephanie Monty"

### Imports

# Standard
import numpy as np
import sys

# Scipy
from scipy import signal

# Project-specific
import mavisim.rampup as rampup

class SeeingGrid:

	"""
	Args:
		input_par = input parameters either hardcoded or altered by the user
		gauss_field = the grid of weighted Gaussians created in AOGaussGrid
		
	Returns:
		seeing_field = a CCD-sized numpy array that contains the weighted seeing wings

	"""

	def __init__(self, input_par, gauss_grid):
		self.input_par = input_par

		self.gauss_grid = gauss_grid


	def main(self):

		big_psf_data = self.input_par.big_psf[0].data

		# Normalize the big PSF
		big_psf_norm = big_psf_data/np.sum(big_psf_data)#(big_psf_data - np.amin(big_psf_data))/(np.amax(big_psf_data) - np.amin(big_psf_data))

		xmin = (self.input_par.seeing_core_rad_pix + 5) - self.input_par.big_psf_ramp/2.0
		xmax = (self.input_par.seeing_core_rad_pix + 5) + self.input_par.big_psf_ramp/2.0

		# Ramp up the seeing limited wings
		big_psf_smooth = big_psf_norm * rampup.ramp_up(xmin, xmax, big_psf_data.shape[0])

		# Create the final grid of seeing wings at the positions of the stellar sources
		seeing_grid = signal.fftconvolve(self.gauss_grid, big_psf_smooth, mode="same")

		return seeing_grid