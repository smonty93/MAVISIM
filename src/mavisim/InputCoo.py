# ----------------------------------------------------------------------------
#
# Fixed May 6/20 disagreement with old code.
#
#
# TITLE - InputCoo
# AUTHOR - Stephanie Monty
# PROJECT - MAVISSimIm
# CONTENTS:
#   1. This class returns an astropy table of the input coordinates of the sources, this can be used for comparisons with the DAOPhot detected positions
#
#		a) Star, Flux, X_True, PM_X, Y_True, PM_Y
#
#			i)	   Star = Number assigned to the star to help with tracking through epochs in PM sims
#			ii)	  Flux = Monochromatic flux of the star (photons/s) passed to the program by the input data
#			iii CCD_Mapped_X/CCD_Mapped_Y	   = True position following sub-pixel shifts due to sub-pixel positions, static distortion and proper motion
#												  Mapped to the FoR of the CCD (swap x and y)
#			iv) CCD_Mapped_PM_X/CCD_Mapped_PM_Y = Proper motion of the star (mas)
#												  Mapped to the FoR of the CCD (swap x and y)
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Generates an input catalogue containing the distorted positions, flux and proper motions of the input source objects. This can be compared to the output
catalogue from DAOPhot FIND

'''
__author__ = "Stephanie Monty"

### Imports

# Standard
import numpy as np
import sys

# Astropy
from astropy.table import Table, Column

# Project-specific
sys.path.append('../src')
import mavisim.input_parameters as input_par

class InputCoo:

	"""
	Args:
		source = the astropy table that contains the source object generated in Source
		
	Returns:
		input_coo = an astropy table containing the TRUE input positions (static distortion, sub-pixel positioning and proper motion taken into account)
					and additional information necessary to compare with the DAOPhot output

	"""

	def __init__(self, source_table, trim_cat=False):

		self.source_table = source_table
		self.trim_cat = trim_cat


	def main(self):

		# Note we need to take two different coordinate conventions into account here, the convention of the CCD (x -> y) and the convention
		# of DAOPhot (starts at position 1, 1)

		# Find the full FoV in pixels, this contains the buffer (which is trimmed later) to contain light from stars outside the FoV
		full_fov_pix = int((input_par.MAVIS_fov + input_par.buffer)/input_par.ccd_sampling)

		# Convert the position to pixels (remove the knowledge of the static distortion)
		x_pos = np.array(np.around(((self.source_table["X"]/input_par.ccd_sampling) + full_fov_pix/2.0), 0), int)
		true_x = x_pos + self.source_table["X_Dist"] - self.source_table["Stat_Dist"][:, 0] + self.source_table["X_PM"]


		y_pos = np.array(np.around(((self.source_table["Y"]/input_par.ccd_sampling) + full_fov_pix/2.0), 0), int)
		true_y = y_pos + self.source_table["Y_Dist"] - self.source_table["Stat_Dist"][:, 1] + self.source_table["Y_PM"]


		# Create the final table, swap the coordinates for the CCD convention
		input_coo = Table(data = ([self.source_table["Star"]]))
		input_coo.add_column(self.source_table["Flux"])
		input_coo.add_column(self.source_table["RA"])
		input_coo.add_column(self.source_table["Dec"])
		input_coo.add_column(Column(true_x), name="CCD_Mapped_X")
		input_coo.add_column(self.source_table["X_PM"], name="CCD_Mapped_PM_X")
		input_coo.add_column(self.source_table["Stat_Dist"][:, 0], name="X Static Dist")
		input_coo.add_column(Column(true_y), name="CCD_Mapped_Y")
		input_coo.add_column(self.source_table["Y_PM"], name="CCD_Mapped_PM_Y")
		input_coo.add_column(self.source_table["Stat_Dist"][:, 1], name="Y Static Dist")

		if self.trim_cat == True:
			# Trim the input catalogue to only contain stars within the MAVIS FoV, a bounding box sides = 30"
			# a,b are the top-left coordinate of the rectangle and (c,d) be its width and height.
			# to judge a point(x0,y0) is in the rectangle, just to check if a < x0 < a+c and b < y0 < b + d

			keep_stars = []

			# Convention is (0, 0) at centre of FoV
			a = -1 * input_par.MAVIS_fov/2.0
			b = input_par.MAVIS_fov/2.0

			width = input_par.MAVIS_fov

			place = 0

			for row in self.source_table:
				if (a <= row["X"]) and (row["X"] <= a + width):
					if (b - width <= row["Y"]) and (row["Y"] <= b):
						keep_stars.append(place)

					else:
						pass
				else:
					pass

				place += 1

			trimmed_coo = input_coo[keep_stars]

			# Correct for the change from the larger buffered FoV to the FoV of the MAVIS CCD
			offset = (full_fov_pix - input_par.ccd_size)/2
			trimmed_coo["CCD_Mapped_X"] = trimmed_coo["CCD_Mapped_X"] - offset
			trimmed_coo["CCD_Mapped_Y"] = trimmed_coo["CCD_Mapped_Y"] - offset

			return (trimmed_coo)

		else:
			return (input_coo)
