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
#   1. This class creates an object with the required information to create a simulated image. This information is accessible via object attributes:
#		         star  =  Number assigned to the star to help with tracking through epochs in PM sims
#		         flux  =  Monochromatic flux of the star (photons/s) passed to the program by the input data
#		       ra/dec  =  RA and Dec of input star to aid in velocity dispersion calculations
#		          x/y  =  Distance from the centre of the field (in arcseconds)
#		    x_pm/y_pm  =  Proper motion of the star (mas)
#		x_dist/y_dist  =  Sub-pixel shifts due to a) the star falling at at a non-integer pixel position, b) static distortion present
#		    gauss_pos  =  Position of the MV Gaussian 
#		    gauss_cov  =  Covariance of the MV Gaussian 
#		  static_dist  =  Static distortion applied [X,Y]
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
from tqdm import tqdm

# Project-specific
from mavisim import make_static_dist_map

class Source:
	"""
	The source object allows access to the parameters of each source (star) required
	to compute the MAVISIM image.

	Args:
		input_par : imported input parameter python file, e.g., `input_parameters.py`.
		exp_time (`float`): exposure time in seconds to simulate.
		static_dist (`bool`, optional): if `True`, add static distortion defined in `input_par`.
		stat_amp (`float`, optional): scaling factor to apply to distortions.
		tt_amp (`float`, optional): scaling factor to apply to tip-tilt blurring kernel width.
	
	Attributes:
		exp_time (`float`): exposure time in seconds to simulate.
		star (`np.ndarray` of `int`): unique ID of each star.
		flux (`np.ndarray` of `float`): flux of each star.
		ra (`np.ndarray` of `float`): RA of each star.
		dec (`np.ndarray` of `float`): Dec of each star.
		x_pos (`np.ndarray` of `float`): X position of each star (in arcsec).
		x_pm (`np.ndarray` of `float`): X proper motion of each star (in pixels).
		x_dist (`np.ndarray` of `float`): sub-pixel X position shift of each star (in pixels).
		y_pos (`np.ndarray` of `float`): Y position of each star (in arcsec).
		y_pm (`np.ndarray` of `float`): Y proper motion of each star (in mas/year).
		y_dist (`np.ndarray` of `float`): sub-pixel Y position shift of each star (in pixels).
		gauss_pos (`np.ndarray` of `float`): X/Y-position of each star (in arcsec).
		gauss_cov (`np.ndarray` of `float`): covariance of Gaussian kernel to simulate tip-tilt blurring (in arcsec^2).
		static_dist	(`np.ndarray` of `float`): static distortion to apply to each source.
	"""

	def __init__(self, input_par, exp_time, static_dist=False, stat_amp=1.0, tt_amp=1.0):
		"""Create a MAVISIM Source data object
		"""
		self._input_par = input_par
		self.exp_time = exp_time
		self._static_dist_flag = static_dist
		self._stat_amp = stat_amp
		self._tt_amp = tt_amp

		# Add static field distortion to the image?, If so, create the functions necessary to determine the distortion to add
		if self._static_dist_flag == True:
			(self._dist_x_func_degmm, self._dist_y_func_degmm) = make_static_dist_map(self._input_par)
		else:
			self._dist_x_func_degmm = self._dist_y_func_degmm = None

		# Determine the covariance matrix only once (identical for all the stars) because NGS tip-tilt captured in e2e PSFs
		self.cov_mat = self._find_covmat()

	def build_source(self):
		""" From the data stored in the object, compute the source data as required
		"""
		dtype = np.float64
		num_stars   = self._input_par.input_cat["X"].shape[0]
		self.star   = np.zeros([num_stars],dtype=int)
		self.flux   = np.zeros([num_stars],dtype=dtype)
		self.ra     = np.zeros([num_stars],dtype=dtype)
		self.dec    = np.zeros([num_stars],dtype=dtype)
		self.x_pos  = np.zeros([num_stars],dtype=dtype)
		self.x_pm   = np.zeros([num_stars],dtype=dtype)
		self.x_dist = np.zeros([num_stars],dtype=dtype)
		self.y_pos  = np.zeros([num_stars],dtype=dtype)
		self.y_pm   = np.zeros([num_stars],dtype=dtype)
		self.y_dist = np.zeros([num_stars],dtype=dtype)
		self.gauss_pos   = np.zeros([num_stars,2],dtype=dtype)
		self.gauss_cov   = np.zeros([num_stars,2,2],dtype=dtype)
		self.static_dist = np.zeros([num_stars,2],dtype=dtype)
		
		for star in tqdm(range(num_stars)):
			star_info = self._input_par.input_cat[star]
			self._compute_row(star,star_info)

	def _compute_row(self, row, star_info):
		"""
		Args:
			star_info = row of the astropy table containing all the info required for that star
			
		Returns:
			row_info = the relevant info formatted into a row to return and append the final table with
		"""

		# Save the sub-pixel shift that comes from converting 15.1 degrees to pixels, keep track of this to shift later
		(x_pos_dist, y_pos_dist, static_dist) = self._find_shift(star_info)
		
		# PM in arcseconds
		pm_x = star_info["PM_X"]
		pm_y = star_info["PM_Y"]
		
		# Multiply the input flux by the exposure time
		total_flux = star_info["Flux"] * self.exp_time

		# Create a column vector for the x,y coordinates of the centre of the Gauss
		# All inputs are in arcseconds, convention is centre of FoV at (0, 0)
		gauss_mux = star_info["X"] + ((self._input_par.gauss_offset + x_pos_dist) * self._input_par.ccd_sampling) + pm_x
		gauss_muy = star_info["Y"] + ((self._input_par.gauss_offset + y_pos_dist) * self._input_par.ccd_sampling) + pm_y
		gauss_mu = np.r_[gauss_mux, gauss_muy]
		
		# Return the formatted info to append as a row to the final astropy table
		# Swap the x and y for the CCD convention

		self.star[row]   = star_info["Star"]
		self.flux[row]   = total_flux
		self.ra[row]     = star_info["RA"]
		self.dec[row]    = star_info["Dec"]
		self.x_pos[row]  = star_info["X"]
		self.x_pm[row]   = pm_x
		self.x_dist[row] = x_pos_dist
		self.y_pos[row]  = star_info["Y"]
		self.y_pm[row]   = pm_y
		self.y_dist[row] = y_pos_dist
		self.gauss_pos[row]   = gauss_mu
		self.gauss_cov[row]   = self.cov_mat
		self.static_dist[row] = static_dist

	def _find_shift(self, star_info):
		"""
		Args:
			star_info = row of the astropy table containing all the info required for that star
			
		Returns:
			x/y_dist = the sub-pixel shift in the x/y-direction for the star

		"""

		# Save the sub-pixel shift that comes from converting e.g. 15.1 degrees to pixels, pass this to the Gaussian later
		x_mu_pix_true = (star_info["X"])/self._input_par.ccd_sampling
		y_mu_pix_true = (star_info["Y"])/self._input_par.ccd_sampling

		# Nearest integer pixel
		x_mu_pix = int(np.around(star_info["X"]/self._input_par.ccd_sampling, 0))
		y_mu_pix = int(np.around(star_info["Y"]/self._input_par.ccd_sampling, 0))

		# Sub-pixel shift
		delta_xpix = x_mu_pix_true - x_mu_pix
		delta_ypix = y_mu_pix_true - y_mu_pix

		if self._static_dist_flag == True:
			# Retrieve the interpolated distortion at each point and convert from the default output (mm) to pixels
			# Input to functions is in degrees
			x_deg = star_info["X"]/3600.0
			y_deg = star_info["Y"]/3600.0

			# Static distortion for the star given the location, converting from mm to pixels
			x_static_dist = (((self._dist_x_func_degmm(x_deg, y_deg) * self._input_par.plate_scale)/self._input_par.ccd_sampling)[0][0]) * self._stat_amp
			y_static_dist = (((self._dist_y_func_degmm(x_deg, y_deg) * self._input_par.plate_scale)/self._input_par.ccd_sampling)[0][0]) * self._stat_amp

			# Sum the two sub-pixel terms to find the final distortion
			final_x_dist = x_static_dist + delta_xpix
			final_y_dist = y_static_dist + delta_ypix

			return (final_x_dist, final_y_dist, np.asarray([x_static_dist, y_static_dist]))

		else:
			# Sub-pixel shift is solely associated with the conversion from arcseconds to pixels
			final_x_dist = delta_xpix
			final_y_dist = delta_ypix

			return (final_x_dist, final_y_dist, np.asarray([0, 0]))

	def decimate(self, nstar):
		"""Decimate the list of objects in the object (e.g., for faster simualtions).

		Args:
			nstar (`int`): number of stars to reduce list down to. The resulting object will
			only contain the first `nstar` stars.
		"""
		if nstar >= self.star.shape[0]:
			raise ValueError(f"Already less than {nstar:d} stars in catalogue ({self.star.shape[0]:d})")
		self.star        =  self.star[:nstar]
		self.flux        =  self.flux[:nstar]
		self.ra          =  self.ra[:nstar]
		self.dec         =  self.dec[:nstar]
		self.x_pos       =  self.x_pos[:nstar]
		self.x_pm        =  self.x_pm[:nstar]
		self.x_dist      =  self.x_dist[:nstar]
		self.y_pos       =  self.y_pos[:nstar]
		self.y_pm        =  self.y_pm[:nstar]
		self.y_dist      =  self.y_dist[:nstar]
		self.gauss_pos   =  self.gauss_pos[:nstar]
		self.gauss_cov   =  self.gauss_cov[:nstar]
		self.static_dist =  self.static_dist[:nstar]

	def _find_covmat(self):
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
		semi_major = np.sqrt(ttres_semi_major**2 + self._input_par.cd_term**2 + self._input_par.vib_term**2) * self._tt_amp
		semi_minor = np.sqrt(ttres_semi_minor**2 + self._input_par.cd_term**2 + self._input_par.vib_term**2) * self._tt_amp
		
		# Convert to arcseconds
		semi_major_arcsec = semi_major * self._input_par.ccd_sampling
		semi_minor_arcsec = semi_minor * self._input_par.ccd_sampling

		# For some reason the array is not oriented the correct way so y and x have to be swapped
		theta = np.pi/2

		# Correct the rotation about the first diagonal by swapping the sine and cosine of the angles
		c, s = np.sin(theta), np.cos(theta)
		R = np.array([[c, -s], [s, c]])
		RT = np.array([[c, s], [-s, c]])
		S = np.array([[semi_major_arcsec**2, 0], [0, semi_minor_arcsec**2]])

		# Generate the covariance matrix following a rotation
		# See https://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
		cov_mat = R @ S @ RT
		
		return (cov_mat)