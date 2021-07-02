# ----------------------------------------------------------------------------
#
# TITLE - trim image
# AUTHOR - Stephanie Monty
# PROJECT - mavisim
# CONTENTS:
#   1. This function trims the final image and input coordinate catalogue to the MAVIS FoV, allowing for light from the stars outside the field to be captured
#
#		a) trimmed_image = image array trimmed to the MAVIS FoV
#		b) trimmed_cat = trimmed catalogue 
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Trims the image the requested FoV

'''
__author__ = "Stephanie Monty"

### Imports

# Standard
import numpy as np


def trim_image(input_par, image):

	"""
	Args:
		image = 2D image (with or without added noise)
		
	Returns:
		image_trim = trimmed image 
		
	"""

	# Find the bounds of the new trimmed image (in the FoR of the original image)
	trim_low = int((image.shape[0] - input_par.ccd_size)/2.0)
	trim_high = int(input_par.ccd_size + (image.shape[0] - input_par.ccd_size)/2.0)

	image_trim = image[trim_low:trim_high, trim_low:trim_high]


	return image_trim

