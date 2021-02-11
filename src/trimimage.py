# ----------------------------------------------------------------------------
#
# TITLE - trim image
# AUTHOR - Stephanie Monty
# PROJECT - MAVISSimIm
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
import sys

# Project-specific
sys.path.append('../src')
import mavissimim.input_parameters as input_par


def trim_image(image, input_coo):

	# Find the bounds of the new trimmed image (in the FoR of the original image)
	trim_low = int((image.shape[0] - input_par.ccd_size)/2.0)
	trim_high = int(input_par.ccd_size + (image.shape[0] - input_par.ccd_size)/2.0)

	image_trim = image[trim_low:trim_high, trim_low:trim_high]

	# Trim the input catalogue to only include stars in the MAVIS FoV
	# Set the boundary to include stars just outside of the CCD limit (diagonal)
	r = 3000

	r_coo = np.sqrt((2667 - input_coo["CCD_Mapped_X"])**2 + (2667 - input_coo["CCD_Mapped_Y"])**2)

	trim = np.where(r_coo <= r)[0]
	trimmed_cat = input_coo[trim]

	# Correct for the change from the larger buffered FoV to the FoV of the MAVIS CCD, swap the x and y for the convention of DAOphot
	x = trimmed_cat["CCD_Mapped_Y"]
	y = trimmed_cat["CCD_Mapped_X"]
 
	trimmed_cat["CCD_Mapped_X"] = trimmed_cat["CCD_Mapped_X"] - 667
	trimmed_cat["CCD_Mapped_Y"] = trimmed_cat["CCD_Mapped_Y"] - 667

	return (image_trim, trimmed_cat)

