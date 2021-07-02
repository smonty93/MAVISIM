# ----------------------------------------------------------------------------
#
# TITLE - bilinear interpolation
# AUTHOR - Stephanie Monty
# PROJECT - mavisim
# CONTENTS:
#   1. This class performs a bilinear interpolation of the four nearest PSFs to the source position to create the final field variable PSF for that source.
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Interpolate between the four nearest PSFs to create the FV PSF at each point
'''
__author__ = "Stephanie Monty"

### Imports

# Standard
import numpy as np
import os

# Astropy
from astropy.table import Table
from astropy.table import Column
from astropy.io import fits

# Project-specific
import mavisim.findclosestvalue as find_closest_value
import mavisim.rampdown as rampdown


class BilinearInterpolation:

	"""
	Args:
		x_mu = x position of the source in arcseconds (from the centre of the FoV) convention is -15, 15 in top left corner
		y_mu = y position -"- 
		
	Returns:
		final_psf = the field variable ramped PSF core interpolated to the position of the source

	"""

	def __init__(self, input_par, x_mu, y_mu):
		self.input_par = input_par
		self.x_mu = x_mu
		self.y_mu = y_mu

	def interp_psf(self):

		"""
		Args:
			None
			
		Returns:
			final_psf = the field variable ramped PSF core interpolated to the position of the source

		"""

		# Create the set of possible PSF locations
		fv_psf_set = self.input_par.fv_psf_grid

		# Check first if the star is located exactly where a PSF exists, if so pass that PSF back
		if np.any(np.isin(fv_psf_set, self.x_mu)) and np.any(np.isin(fv_psf_set, self.x_mu)):

			fv_psf = self.find_fv_psf(self.x_mu, self.y_mu)
			final_psf = self.ramp_psf(fv_psf)

		else:
			# Locate where the star would fall within this set to gather the set of four nearest PSFs
			set_x = np.searchsorted(fv_psf_set, self.x_mu)
			set_y = np.searchsorted(fv_psf_set, self.y_mu)
		
			# Deal with the case when the star is at the edge of the image, but within either x or y bounds
			# Note that right and left refer to the position of the star in the fv_psf_set (-15, -12 ... 12, 15)
			if (set_x == 0) and (set_y < len(fv_psf_set)):
				x_left = fv_psf_set[set_x]
				x_right = fv_psf_set[set_x + 1]

				y_left = fv_psf_set[set_y - 1]
				y_right = fv_psf_set[set_y]

			elif (set_y == 0) and (set_x < len(fv_psf_set)):
				x_left = fv_psf_set[set_x - 1]
				x_right = fv_psf_set[set_x]

				y_left = fv_psf_set[set_y]
				y_right = fv_psf_set[set_y + 1]

			elif (set_x == 0) and (set_y == 0):
				x_left = fv_psf_set[set_x]
				x_right = fv_psf_set[set_x + 1]

				y_left = fv_psf_set[set_y]
				y_right = fv_psf_set[set_y + 1]

			elif (set_x == len(fv_psf_set)) and (set_y < len(fv_psf_set)):
				x_left = fv_psf_set[set_x - 2]
				x_right = fv_psf_set[set_x  - 1]

				y_left = fv_psf_set[set_y - 1]
				y_right = fv_psf_set[set_y]

			elif (set_y == len(fv_psf_set)) and (set_x < len(fv_psf_set)):
				x_left = fv_psf_set[set_x - 1]
				x_right = fv_psf_set[set_x]

				y_left = fv_psf_set[set_y - 2]
				y_right = fv_psf_set[set_y - 1]

			elif (set_x == 0) and (set_y == len(fv_psf_set)):
				x_left = fv_psf_set[set_x]
				x_right = fv_psf_set[set_x + 1]

				y_left = fv_psf_set[set_y - 2]
				y_right = fv_psf_set[set_y - 1]

			elif (set_x == len(fv_psf_set)) and (set_y == 0):
				x_left = fv_psf_set[set_x - 2]
				x_right = fv_psf_set[set_x  - 1]

				y_left = fv_psf_set[set_y]
				y_right = fv_psf_set[set_y + 1]


			# Deal with the case when the star is completely beyond both the x and y bounds
			elif (set_x == len(fv_psf_set)) and (set_y == len(fv_psf_set)):
				x_left = fv_psf_set[set_x - 2]
				x_right = fv_psf_set[set_x  - 1]

				y_left = fv_psf_set[set_y - 2]
				y_right = fv_psf_set[set_y - 1]

			else:
				x_left = fv_psf_set[set_x - 1]
				x_right = fv_psf_set[set_x]

				y_left = fv_psf_set[set_y - 1]
				y_right = fv_psf_set[set_y]


			# Sort the corners so that they follow this convention
			#
			# 1----2
			# |    |
			# 3----4
			#
			corners = np.array([[x_left, y_right], [x_right, y_right], [x_left, y_left], [x_right, y_left]])

			psf1 = self.find_fv_psf(corners[0][0], corners[0][1])
			psf2 = self.find_fv_psf(corners[1][0], corners[1][1])
			psf3 = self.find_fv_psf(corners[2][0], corners[2][1])
			psf4 = self.find_fv_psf(corners[3][0], corners[3][1])

			interp_psf = self.bilinear_interpolation(psf1, psf2, psf3, psf4, corners)

			final_psf = self.ramp_psf(interp_psf)

		return (final_psf)

	def bilinear_interpolation(self, psf1, psf2, psf3, psf4, corners):

		"""
		Args:
			psf1, psf2, psf3, psf4 = 2D arrays of the four closest PSFs loaded from fits files
			corners = array of the four corner locations to weight the PSF contributions
			
		Returns:
			PSF_interp = interpolated PSF created from the four weighted PSFs (2D array)

		"""
		# Applying the following formula:
		# PSF_1 = psf1((l-a)/l) + psf2((l-b)/l)
		# PSF_2 = psf3((l-a)/l) + psf4((l-b)/l)
		# PSF_final = PSF_1((l-c)/l) + PSF_2((l-d)/l)
		#
		# Distance conventions are as follows
		#
		# 1-------2
		# |       c
		# | a o b |
		# |       d
		# 3-------4
		#

		# Distance between the corners (side of the box), fixed by the grid resolution
		l = np.absolute(self.input_par.fv_psf_grid[0] - self.input_par.fv_psf_grid[1])

		# Distance between the x position of the star (in array coord) and the left corners
		# VERY IMPORTANT THIS HAS TO FOLLOW CONVENTION OF (0, 0) IN CENTRE OF CCD
		a = np.absolute(corners[0][0] - self.x_mu)

		# Distance between the x position of the star and the right corners
		b = np.absolute(corners[1][0] - self.x_mu)

		PSF_1 = (((l - a)/l) * psf1) + (((l - b)/l) * psf2)
		PSF_2 = (((l - a)/l) * psf3) + (((l - b)/l) * psf4)

		# Distance between the y position of the star and the top corners
		c = np.absolute(corners[1][1] - self.y_mu)

		# Distance between the y position of the star and the bottom corners
		d = np.absolute(corners[3][1] - self.y_mu)

		PSF_interp = (((l - c)/l) * PSF_1) + (((l - d)/l) * PSF_2)


		return PSF_interp


	def find_fv_psf(self, x_closest, y_closest):
		"""
		Args:
			x_closest, y_closest = (x, y) coordinates in arcsec, centre of FoV at (0, 0)
			
		Returns:
			fv_psf_norm = 2D array of the closest PSF to the star loaded from the PSF database and normalised

		"""
		
		# Store the location of the current directory
		cur_dir = os.getcwd()
		
		# Change to the directory where the FV PSFs are stored
		os.chdir(self.input_par.fv_psf_path)
		
		# Open the FV PSF that's closest to the point source
		fv_psf = fits.open(self.input_par.fv_psf_filename % (str(int(x_closest)), str(int(y_closest))))
		fv_psf_truval = fv_psf[0].data
		fv_psf_norm = fv_psf_truval/np.sum(fv_psf_truval)
		fv_psf.close()
		
		# Change back to the current working directory
		os.chdir(cur_dir)
		
		return (fv_psf_norm)

	def ramp_psf(self, fv_psf_norm):
		"""
		Args:
			fv_psf_norm = normalised 2D array representation of the PSF
			
		Returns:
			fv_psf_ramped = normalised PSF with ramped the edges surrounding the AO core to smooth the transition from AO to seeing wings

		"""

		# Find the location of the start and stop of the ramp function
		xmin = (self.input_par.psf_core_rad_pix+5) - self.input_par.ramp_size/2.0
		xmax = (self.input_par.psf_core_rad_pix+5) + self.input_par.ramp_size/2.0

		# Ramp the PSF down to blend with the seeing rings
		fv_psf_ramped = fv_psf_norm * rampdown.ramp_down(xmax, xmin, fv_psf_norm.shape[0])

		return (fv_psf_ramped)