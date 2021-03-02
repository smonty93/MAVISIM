# ----------------------------------------------------------------------------
#
# TITLE - Source
# AUTHOR - Stephanie Monty
# PROJECT - MAVISim
# CONTENTS:
#	UPDATE:
#	PSF Must be sampled at atleast 3.75 mas to avoid smoothing from the Gaussian kernel representing each star (FWHM = 2, adds smoothing of 7.5 mas)
#
#	1. NOTE: Very important, this program and all programs assume all projection effects have been taken into account, the cos(Dec) term has been
#			 dealt with
#
#   1. This class returns a source object with the required information to create a simulated image. This information includes:
#
#		a) Star, Flux, X, Dist_X, PM_X, Y, Dist_Y, PM_Y, 
#
#			i)          Star = Number assigned to the star to help with tracking through epochs in PM sims
#			ii)         Flux = Monochromatic flux of the star (photons/s) passed to the program by the input data
#			iii)	  RA/Dec = RA and Dec of input star to aid in velocity dispersion calculations
#			iii)         X/Y = Distance from the centre of the field (in arcseconds)
#			iv)    X_PM/Y_PM = Proper motion of the star (mas)
#			v) X_Dist/Y_Dist = Sub-pixel shifts due to a) the star falling at at a non-integer pixel position, b) static distortion present
#			vi)    Stat_Dist = Static distortion applied [X,Y]
#			vi)    Gauss_Src = Matrix rep. of the MV Gaussian 
#
#names = ('Star', 'Flux', 'RA', 'Dec', 'X', 
#		 'X_PM', 'X_Dist', 'Y', 'Y_PM', 'Y_Dist', 
#		 'Gauss_Info'),
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Generates a source object containing all the point sources used to create the final image. 
A gaussian is created for each point source taking into account sub-pixel shifts due to distortion and non-integer pixel positions.

'''
__author__ = "Stephanie Monty"

### Imports

# Standard
import numpy as np
import sys

# Astropy
from astropy.table import Table

# Project-specific
sys.path.append('../src')
import mavisim.input_parameters as input_par
import mavisim.makestaticdistortmap as make_static_map

class Source:

	"""
	Args:
		mavis_src = the source catalogue input to the program
		exp_time = exposure time in seconds (to multiply with the flux)
		static_dist = Boolean stating whether static distortion should be taken into account (and therefore added) to the image
		stat_amp = Factor to amplify the static distortion (immitate a dynamic term)
		tt_var = Boolean stating whether the TT kernel is spacially variable or is fixed (dealt with already in the creation of the PSFs themselves)
		tt_amp_fac = Factor to amplify the TT kernel
		tt_static = Boolean stating if a static TT kernel should be applied
		tt_kern = Static TT kernel to apply, units are mas!!!

	Returns:
		source = an astropy table containing the above mentioned data

	"""

	def __init__(self, mavis_src, exp_time, tt_amp_fac, static_dist, stat_amp):

		self.mavis_src = mavis_src

		self.exp_time = exp_time

		self.static_dist = static_dist

		self.tt_amp_fac = tt_amp_fac

		# Hard-coded for now
		self.tt_var = False

		if self.tt_var == False:
			self.cov_mat = self.find_covmat()

		# Add a variable TT kernel? If so, create the functions nec. to find the kernel
		# If not, it's assumed to be constant and is included in the PSFs themselves
		if self.tt_var == True:

			(self.semimajor_func, self.semiminor_func, self.theta_func) = make_tt_map.make_ttkernel_map()

		self.stat_amp = stat_amp

		# Track the semi-major and semi-minor axes to get a sense of the scale of the TT distortion
		# Not appending these for now because TT kernel static
		self.track_smajor = []
		self.track_sminor = []

		# Add static field distortion to the image?, If so, create the functions necessary to determine the distortion to add
		if self.static_dist == True:

			(self.dist_x_func_degmm, self.dist_y_func_degmm) = make_static_map.make_static_dist_map()


	def main(self):

		"""
		Args:
			None
			
		Returns:
			source_table = an astropy table containing the output columns (see above)

		"""
		num_stars = self.mavis_src["X"].shape[0]

		# Perform this once for the case of a TT residual 

		# Perform this once storing the information in a table object, thereafter append the table
		star_info = self.mavis_src[0]

		star_row = self.find_star_row(star_info)

		# Hardcode this once... and then never again!
		source_table = Table(data = ([star_row[0]], [star_row[1]], [star_row[2]], [star_row[3]], [star_row[4]],
									[star_row[5]], [star_row[6]], [star_row[7]], [star_row[8]], [star_row[9]],
									[star_row[10]], [star_row[11]]), 
							names=('Star', 'Flux', 'RA', 'Dec', 'X', 
									'X_PM', 'X_Dist', 'Y', 'Y_PM', 'Y_Dist', 
									'Gauss_Info', 'Stat_Dist'),
							meta={"exp_time": self.exp_time})

		source_table["Flux"].unit   = "photons"
		source_table["X"].unit      = "arcseconds"
		source_table["X_PM"].unit   = "pixels/year"
		source_table["X_Dist"].unit = "pixels"
		source_table["Y"].unit      = "arcseconds"
		source_table["Y_PM"].unit   = "pixels/year"
		source_table["Y_Dist"].unit = "pixels"
		source_table["Stat_Dist"].unit = "pixels"


		# Iterate for every star in the input catalogue, create a Gaussian for each point source.
		for star in range(1, num_stars):

			star_info = self.mavis_src[star]

			star_row = self.find_star_row(star_info)

			source_table.add_row(star_row)

		# Print some information on the TT kernel applied (or not)
		try:
			print ("Max semi-major distortion (mas) = %s" % str(np.amax(self.track_smajor)))
			print ("Min semi-major distortion (mas) = %s" % str(np.amin(self.track_smajor)))
			print ("Average semi-major distortion (mas) = %s" % str(np.average(self.track_smajor)))
			print ("	")
			print ("Max semi-minor distortion (mas) = %s" % str(np.amax(self.track_sminor)))
			print ("Min semi-minor distortion (mas) = %s" % str(np.amin(self.track_sminor)))
			print ("Average semi-minor distortion (mas) = %s" % str(np.average(self.track_sminor)))

		except:
			print ("No TT kernel applied")

		# Print information on the static distortion applied

		if self.static_dist == True:
			ave_stat = np.average(np.sqrt(source_table["Stat_Dist"][0]**2 + source_table["Stat_Dist"][1]**2))

			print ("Average X static distortion applied (mas) = %s" % str(ave_stat * (1e3 * input_par.ccd_sampling)))

		else:
			print ("No static distortion applied.")
		
		return source_table


	def find_star_row(self, star_info):

		"""
		Args:
			star_info = row of the astropy table containing all the info required for that star
			
		Returns:
			row_info = the relevant info formatted into a row to return and append the final table with

		"""

		# Save the sub-pixel shift that comes from converting 15.1 degrees to pixels, keep track of this to shift later
		(x_pos_dist, y_pos_dist, stat_dist) = self.find_shift(star_info)
		
		# PM in pixels
		pm_x = (star_info["PM_X"])/input_par.ccd_sampling
		pm_y = (star_info["PM_Y"])/input_par.ccd_sampling
		
		# Multiply the input flux by the exposure time
		total_flux = star_info["Flux"] * self.exp_time

		# Create a column vector for the x,y coordinates of the centre of the Gauss
		# All inputs are in arcseconds, convention is centre of FoV at (0, 0)
		gauss_mux = star_info["X"] + ((input_par.gauss_offset+x_pos_dist+pm_x) * input_par.ccd_sampling)
		gauss_muy = star_info["Y"] + ((input_par.gauss_offset+y_pos_dist+pm_y) * input_par.ccd_sampling)

		gauss_mu = np.r_[gauss_mux, gauss_muy]

		# Create the Gaussian info to pass to Jesse's code, this is a scalar, vector, matrix
		gauss_info = [total_flux, gauss_mu, self.cov_mat]

		# Return the formatted info to append as a row to the final astropy table
		# Swap the x and y for the CCD convention
		row_info = [star_info["Star"], total_flux, star_info["RA"], star_info["Dec"],
					star_info["X"], pm_x, x_pos_dist, star_info["Y"], pm_y, y_pos_dist,
					gauss_info, stat_dist]

		return row_info

	def find_shift(self, star_info):

		"""
		Args:
			star_info = row of the astropy table containing all the info required for that star
			
		Returns:
			x/y_dist = the sub-pixel shift in the x/y-direction for the star

		"""

		# Save the sub-pixel shift that comes from converting e.g. 15.1 degrees to pixels, pass this to the Gaussian later
		x_mu_pix_true = (star_info["X"])/input_par.ccd_sampling
		y_mu_pix_true = (star_info["Y"])/input_par.ccd_sampling

		# Nearest integer pixel
		x_mu_pix = int(np.around(star_info["X"]/input_par.ccd_sampling, 0))
		y_mu_pix = int(np.around(star_info["Y"]/input_par.ccd_sampling, 0))

		# Sub-pixel shift
		delta_xpix = x_mu_pix_true - x_mu_pix
		delta_ypix = y_mu_pix_true - y_mu_pix


		if self.static_dist == True:
			# Retrieve the interpolated distortion at each point and convert from the default output (mm) to pixels

			# Input to functions is in degrees
			x_deg = star_info["X"]/3600.0
			y_deg = star_info["Y"]/3600.0

			# Static distortion for the star given the location, converting from mm to pixels
			x_stat_dist = (((self.dist_x_func_degmm(x_deg, y_deg) * input_par.plate_scale)/input_par.ccd_sampling)[0][0]) * self.stat_amp
			y_stat_dist = (((self.dist_y_func_degmm(x_deg, y_deg) * input_par.plate_scale)/input_par.ccd_sampling)[0][0]) * self.stat_amp

			# Sum the two sub-pixel terms to find the final distortion
			final_x_dist = x_stat_dist + delta_xpix
			final_y_dist = y_stat_dist + delta_ypix

			return (final_x_dist, final_y_dist, np.array([x_stat_dist, y_stat_dist]))

		else:
			# Sub-pixel shift is solely associated with the conversion from arcseconds to pixels
			final_x_dist = delta_xpix
			final_y_dist = delta_ypix

			return (final_x_dist, final_y_dist, np.array([0, 0]))


	def find_covmat(self):
		"""
		At the moment this only needs to be called once, as all the stars have the same covariance matrix.
		This is for the case that the TT residual is already captured in the EtE PSFs
		Args:
			star_info = row of the astropy table containing all the info required for that star
			
		Returns:
			cov_mat = covariance matrix to create the TT kernel

		"""

		# For now the semi-major and semi-minor axes will be identical (assuming only vibration and C.D. terms) 
		# Leave the option to add an additional skewed TT residual kernel in the future
		ttres_semi_major = 0
		ttres_semi_minor = 0

		# Add the fixed vibration and CD terms
		semi_major = np.sqrt(ttres_semi_major**2 + input_par.cd_term**2 + input_par.vib_term**2) * self.tt_amp_fac
		semi_minor = np.sqrt(ttres_semi_minor**2 + input_par.cd_term**2 + input_par.vib_term**2) * self.tt_amp_fac
		
		# Convert to arcseconds
		semi_major_arcsec = semi_major * input_par.ccd_sampling
		semi_minor_arcsec = semi_minor * input_par.ccd_sampling

		# For some reason the array is not oriented the correct way so y and x have to be swapped
		theta = np.pi/2

		# Correct the rotation about the first diagonal by swapping the sine and cosine of the angles
		c, s = np.sin(theta), np.cos(theta)
		R = np.matrix([[c, -s], [s, c]])
		RT = np.matrix([[c, s], [-s, c]])
		S = np.matrix([[semi_major_arcsec**2, 0], [0, semi_minor_arcsec**2]])

		# Generate the covariance matrix following a rotation
		# See https://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
		cov_mat = R @ S @ RT
		
		return (cov_mat)


if __name__ == '__main__':
	main()