# ----------------------------------------------------------------------------
#
# TITLE - Source
# AUTHOR - Stephanie Monty
# PROJECT - MAVISIM
# CONTENTS:
#	V1.1 DIFFERENCES:
# 	1. NOT possible to add NGS tip-tilt jitter term on the fly, the NGS constellation is fixed as the best-case scenario (three stars, spaced evenly throughout the field)
#      and is already built into the e2e PSFs.
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
#			vi)    Gauss_Info = Matrix rep. of the MV Gaussian 
#			vi)    static_dist = Static distortion applied [X,Y]
#
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Generates a source object containing all the point sources used to create the final image. 
The info necessary to build a gaussian is created for each point source taking into account sub-pixel shifts due to distortion and non-integer pixel positions.

'''
__author__ = "Stephanie Monty"

### Imports

# Standard
import numpy as np

# Astropy
from astropy.table import Table

# Project-specific
import mavisim.makestaticdistortmap as make_static_map

class Source:

	"""
	Args:
		input_par = input parameters either hardcoded or altered by the user
		mavis_src = the source catalogue input to the program
		exp_time = exposure time in seconds (to multiply with the flux)
		static_dist = Boolean stating whether static distortion should be taken into account (and therefore added) to the image
		stat_amp = Factor to amplify the static distortion (immitate a dynamic term)
		tt_amp_fac = Factor to amplify the TT kernel

	Returns:
		source = an astropy table containing the above mentioned data

	"""

	def __init__(self, input_par, mavis_src, exp_time, static_dist=False, stat_amp=1.0, tt_amp=1.0):
		self.input_par = input_par

		self.mavis_src = mavis_src

		self.exp_time = exp_time

		self.static_dist = static_dist

		self.stat_amp = stat_amp

		self.tt_amp = tt_amp

		# Add static field distortion to the image?, If so, create the functions necessary to determine the distortion to add
		if self.static_dist == True:

			(self.dist_x_func_degmm, self.dist_y_func_degmm) = make_static_map.make_static_dist_map(self.input_par)

		# Determine the covariance matrix only once (identical for all the stars) because NGS tip-tilt captured in e2e PSFs
		self.cov_mat = self.find_covmat()

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
									'Gauss_Info', 'Static_Dist'),
							meta={"exp_time": self.exp_time})

		source_table["Flux"].unit   = "photons"
		source_table["X"].unit      = "arcseconds"
		source_table["X_PM"].unit   = "pixels/year"
		source_table["X_Dist"].unit = "pixels"
		source_table["Y"].unit      = "arcseconds"
		source_table["Y_PM"].unit   = "pixels/year"
		source_table["Y_Dist"].unit = "pixels"
		source_table["Static_Dist"].unit = "pixels"


		# Iterate for every star in the input catalogue, create a Gaussian for each point source.
		for star in range(1, num_stars):

			star_info = self.mavis_src[star]

			star_row = self.find_star_row(star_info)

			source_table.add_row(star_row)

		# Print information on the static distortion applied
		print (" ")
		print ("Average X static distortion applied (mas) = %s" % str(np.average(np.absolute(source_table["Static_Dist"][:, 0])) * (1e3 * self.input_par.ccd_sampling)))
		
		return source_table

	def find_star_row(self, star_info):

		"""
		Args:
			star_info = row of the astropy table containing all the info required for that star
			
		Returns:
			row_info = the relevant info formatted into a row to return and append the final table with
		"""

		# Save the sub-pixel shift that comes from converting 15.1 degrees to pixels, keep track of this to shift later
		(x_pos_dist, y_pos_dist, static_dist) = self.find_shift(star_info)
		
		# PM in arcseconds
		pm_x = star_info["PM_X"]
		pm_y = star_info["PM_Y"]
		
		# Multiply the input flux by the exposure time
		total_flux = star_info["Flux"] * self.exp_time

		# Create a column vector for the x,y coordinates of the centre of the Gauss
		# All inputs are in arcseconds, convention is centre of FoV at (0, 0)
		gauss_mux = star_info["X"] + ((self.input_par.gauss_offset + x_pos_dist) * self.input_par.ccd_sampling) + pm_x
		gauss_muy = star_info["Y"] + ((self.input_par.gauss_offset + y_pos_dist) * self.input_par.ccd_sampling) + pm_y

		gauss_mu = np.r_[gauss_mux, gauss_muy]

		# Create the Gaussian info to pass to Jesse's code, this is a scalar, vector, matrix
		gauss_info = [total_flux, gauss_mu, self.cov_mat]

		# Return the formatted info to append as a row to the final astropy table
		# Swap the x and y for the CCD convention
		row_info = [star_info["Star"], total_flux, star_info["RA"], star_info["Dec"],
					star_info["X"], pm_x, x_pos_dist, star_info["Y"], pm_y, y_pos_dist,
					gauss_info, static_dist]

		return row_info

	def find_shift(self, star_info):

		"""
		Args:
			star_info = row of the astropy table containing all the info required for that star
			
		Returns:
			x/y_dist = the sub-pixel shift in the x/y-direction for the star

		"""

		# Save the sub-pixel shift that comes from converting e.g. 15.1 degrees to pixels, pass this to the Gaussian later
		x_mu_pix_true = (star_info["X"])/self.input_par.ccd_sampling
		y_mu_pix_true = (star_info["Y"])/self.input_par.ccd_sampling

		# Nearest integer pixel
		x_mu_pix = int(np.around(star_info["X"]/self.input_par.ccd_sampling, 0))
		y_mu_pix = int(np.around(star_info["Y"]/self.input_par.ccd_sampling, 0))

		# Sub-pixel shift
		delta_xpix = x_mu_pix_true - x_mu_pix
		delta_ypix = y_mu_pix_true - y_mu_pix


		if self.static_dist == True:
			# Retrieve the interpolated distortion at each point and convert from the default output (mm) to pixels

			# Input to functions is in degrees
			x_deg = star_info["X"]/3600.0
			y_deg = star_info["Y"]/3600.0

			# Static distortion for the star given the location, converting from mm to pixels
			x_static_dist = (((self.dist_x_func_degmm(x_deg, y_deg) * self.input_par.plate_scale)/self.input_par.ccd_sampling)[0][0]) * self.stat_amp
			y_static_dist = (((self.dist_y_func_degmm(x_deg, y_deg) * self.input_par.plate_scale)/self.input_par.ccd_sampling)[0][0]) * self.stat_amp

			# Sum the two sub-pixel terms to find the final distortion
			final_x_dist = x_static_dist + delta_xpix
			final_y_dist = y_static_dist + delta_ypix

			return (final_x_dist, final_y_dist, np.asarray([x_static_dist, y_static_dist]))

		else:
			# Sub-pixel shift is solely associated with the conversion from arcseconds to pixels
			final_x_dist = delta_xpix
			final_y_dist = delta_ypix

			return (final_x_dist, final_y_dist, np.asarray([0, 0]))


	def find_covmat(self):
		"""
		At the moment this only needs to be called once, as all the stars have the same covariance matrix.
		This is for the case that the TT residual is already captured in the e2e PSFs

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
		semi_major = np.sqrt(ttres_semi_major**2 + self.input_par.cd_term**2 + self.input_par.vib_term**2) * self.tt_amp
		semi_minor = np.sqrt(ttres_semi_minor**2 + self.input_par.cd_term**2 + self.input_par.vib_term**2) * self.tt_amp
		
		# Convert to arcseconds
		semi_major_arcsec = semi_major * self.input_par.ccd_sampling
		semi_minor_arcsec = semi_minor * self.input_par.ccd_sampling

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