# ----------------------------------------------------------------------------
#
# TITLE - AOGaussGrid
# AUTHOR - Stephanie Monty
# PROJECT - mavisim
# CONTENTS:
#
#	CHANGE AUGUST 20/20
# 	Normalised the product of the convolution of the gauss star with the AO PSF (there was some artificial flux gain during convolution)
#	TT gaussians were contributing an artificial flux gain for every amplification term, the AOGrid is now normalised post convolution and 
#	then multiplied by the flux
#
#   1. This class returns the field of source objects convolved with AO cores (no seeing wings) as well as a field of the raw Gaussians
#
#		a) ao_field = a CCD-sized array populated with the AO cores relevant for each star
#		b) gauss_field = a CCD-sized array populated with the source Gaussians
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Generates two fields, of the AO cores and the source Gaussians for use in building the final image 
Both field variable and static PSFs can be used to build these grids

'''
__author__ = "Stephanie Monty"

### Imports

# Standard
import numpy as np

# Scipy
from scipy import signal

# Project Specific
import mavisim.rampdown as rampdown
from mavisim.BilinearInterpolation import BilinearInterpolation

class AOGaussGrid:

	"""
	Args:
		self.input_par = input parameters either hardcoded or altered by the user
		source = the astropy table that contains the source object generated in Source
		fv_psf = Boolean stating whether we want to use a static or spacially variable PSF
		
	Returns:
		ao_field = a CCD-sized numpy array that contains the product of the convolution of the AO core and the source Gaussian
		gauss_field = A CCD-sized numpy array that contains the raw source gaussians

	"""

	def __init__(self, input_par, source_table, fv_psf=True):
		self.input_par = input_par

		self.source_table = source_table

		# Use a FV PSFs or not (static PSF)
		self.fv_psf = fv_psf

	def main(self):

		(ao_field, gauss_field) = self.make_fields()

		return (ao_field, gauss_field)

	def make_fields(self):
		"""
		Args:
			None
		
		Returns:
			ao_field = a CCD-sized numpy array that contains the product of the convolution of the AO core and the source Gaussian
			gauss_field = A CCD-sized numpy array that contains the raw source gaussians
		"""
	
		# Find the full FoV in pixels, this contains the buffer (which is trimmed later) to contain light from stars outside the FoV
		full_fov_pix = int((self.input_par.MAVIS_fov + self.input_par.buffer)/self.input_par.ccd_sampling)

		# Create the CCD-sized array to store the convolved cores
		ao_field = np.zeros([full_fov_pix, full_fov_pix])

		# Create the CCD-sized array to store the raw Gaussians
		gauss_field = np.zeros([full_fov_pix, full_fov_pix])

		# If FV PSF is not selected, use the static PSF
		if self.fv_psf == False:

			static_psf = self.input_par.static_psf[0].data
			static_norm = static_psf/np.sum(static_psf)
			psf_ramped = self.ramp_psf(static_norm)

		# Loop over all the stars in the source table
		for row in self.source_table:
			
			# If a FV PSF is selected, then update the ramped psf by interpolating the FV PSFs to find the corresponding AO core
			if self.fv_psf == True:
				psf_ramped = BilinearInterpolation(self.input_par, row["X"], row["Y"]).interp_psf()

			# Load the Gaussian for rebinning and the Gaussian for the seeing grid
			gauss = row["Gauss_Src"]

			# Convolve the AO core with the corresponding Gaussian
			ao_conv_core = signal.fftconvolve(psf_ramped, gauss, mode="same")

			# At this point we need to swap x and y to recreate the CCD coordinate system
			x_loc = int(np.around((row["Y"]/self.input_par.ccd_sampling) + ((full_fov_pix)/2.0), 0))
			y_loc = int(np.around((row["X"]/self.input_par.ccd_sampling) + ((full_fov_pix)/2.0), 0))

			# Find the bounds of where the convolved ao core fits into the larger ao field
			(xleft, xright, yleft, yright) = self.find_array_bounds(x_loc, y_loc, ao_conv_core.shape[0])

			# Find the bounds of where the Gaussian fits into the larger Gaussian field
			(xleft2, xright2, yleft2, yright2) = self.find_array_bounds(x_loc, y_loc, self.input_par.gauss_width)

    
    		# Catch the cases where the star is outside the extended FoV
			try: 
				gauss_field[xleft2:xright2, yleft2:yright2] += gauss * row["Wing_Weight"]

				# Store the convolved core and weighted gaussian
				ao_field[xleft:xright, yleft:yright] += ao_conv_core 

			except:
				pass
			

		return (ao_field, gauss_field)

	def find_array_bounds(self, x_loc, y_loc, size):
		"""
		Args:
			x_loc = x position of the centre of the stellar Gaussian
			y_loc = y position of the cetnre of the stellar Gaussian
			size = the size of the AO core or Gaussian
		
		Returns:
			xleft, xright, yleft, yright = positional bounds in the larger array within to store the convolved core or Gaussian
		"""

		xright = int(np.around((x_loc + size/2.0), 0))
		xleft = int(np.around((x_loc - size/2.0), 0))

		yright = int(np.around((y_loc + size/2.0), 0))
		yleft = int(np.around((y_loc - size/2.0), 0))

		return (xleft, xright, yleft, yright)

	def ramp_psf(self, norm_psf):
		"""
		Args:
			norm_psf = normalized psf array
		
		Returns:
			psf_ramped = ramped psf array
		"""

		# Find the location of the start and stop of the ramp function
		xmin = (self.input_par.psf_core_rad_pix+5) - self.input_par.ramp_size/2.0
		xmax = (self.input_par.psf_core_rad_pix+5) + self.input_par.ramp_size/2.0

		# Ramp the PSF down to blend with the seeing rings
		psf_ramped = norm_psf * rampdown.ramp_down(xmax, xmin, norm_psf.shape[0])

		return (psf_ramped)
