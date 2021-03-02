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
#			i)           Star = Number assigned to the star to help with tracking through epochs in PM sims
#			ii)          Flux = Monochromatic flux of the star (photons/s) passed to the program by the input data
#			iii CCD_Mapped_X/CCD_Mapped_Y       = True position following sub-pixel shifts due to sub-pixel positions, static distortion and proper motion
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

	def __init__(self, source_table):

		self.source_table = source_table


	def main(self):

		# Note we need to take two different coordinate conventions into account here, the convention of the CCD (x -> y) and the convention
		# of DAOPhot (starts at position 1, 1)

		# Find the full FoV in pixels, this contains the buffer (which is trimmed later) to contain light from stars outside the FoV
		full_fov_pix = int((input_par.MAVIS_fov + input_par.buffer)/input_par.ccd_sampling)

		# Convert the position to pixels (remove the knowledge of the static distortion)
		x_pos = np.array(np.around(((self.source_table["X"]/input_par.ccd_sampling) + full_fov_pix/2.0), 0), int)
		true_x = x_pos + self.source_table["X_Dist"] - self.source_table["Stat_Dist"][:, 0] + self.source_table["X_PM"] + 1


		y_pos = np.array(np.around(((self.source_table["Y"]/input_par.ccd_sampling) + full_fov_pix/2.0), 0), int)
		true_y = y_pos + self.source_table["Y_Dist"] - self.source_table["Stat_Dist"][:, 1] + self.source_table["Y_PM"] + 1

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

		return (input_coo)
