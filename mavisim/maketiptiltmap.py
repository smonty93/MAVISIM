# ----------------------------------------------------------------------------
#
# TITLE - make_tt_map
# AUTHOR - Stephanie Monty
# PROJECT - mavisim
# CONTENTS:
#   1. 
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Function to create the tip/tilt kernel map, interpolating and creating functions to return the kernel params at every point

NOTE: The convention of the map and the convention that's kept constant throughout is that the x, y positions of the objects are in arcseconds and 
	  relative to the centre of the FoV, (-15, 15) is the top left corner

	Returns: semimajor_func = interpolated function that returns the semi-major axis in arcseconds, given an input cog (in arcseconds)
			 semiminor_func = interpolated function that returns the semi-minor axis is arcseconds, given an input cog (in arcseconds)
			 theta_func = interpolated function that returns the inclination angle of the ellipse relative to the x-axis, given an in input cog
		
'''
__author__ = "Stephanie Monty"

### Imports

## Basic
import numpy as np

## Scipy
from scipy import interpolate

def make_ttkernel_map(input_par):
	"""
	Args:
		input_par = input parameters either hardcoded or altered by the user

	Returns:
		semimajor_func = returns a tip-tilt semi-major axis for a given centre of gravity (location in arcseconds) 
		semiminor_func = returns a tip-tilt semi-minor axis for a given centre of gravity (location in arcseconds) 
		theta_func = returns a inclination angle for a given centre of gravity (location in arcseconds) 

	"""

	tt_jitter_info = input_par.tt_jitter

	# All positions in arcseconds, used as input into the BivariateSpline interpolator (grid-based data)
	x = tt_jitter_info[0].data[0][0]
	y = tt_jitter_info[0].data[1][:, 0]

	# Collect the ellipsoid information required for the interpolation, given the format of the TT kernel maps supplied by Cedric
	semimajor = tt_jitter_info[0].data[2]
	semiminor = tt_jitter_info[0].data[3]

	# Account for the instances where theta is presented in degrees vs radians
	if input_par.tt_theta_conv == "Degrees":
		theta = tt_jitter_info[0].data[4]

	else:
		# Convert it to degrees
		theta = (tt_jitter_info[0].data[4]) * (180.0/np.pi)

	# Correct for the discontinuities in the interpolation
	theta[np.where(theta < 0)] += 180

	#https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html#scipy-interpolate-rectbivariatespline
	semimajor_func = interpolate.RectBivariateSpline(x, y, semimajor)
	semiminor_func = interpolate.RectBivariateSpline(x, y, semiminor)
	theta_func = interpolate.RectBivariateSpline(x, y, theta)

	return (semimajor_func, semiminor_func, theta_func)

