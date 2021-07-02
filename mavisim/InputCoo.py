# ----------------------------------------------------------------------------
#
# Fixed May 6/20 disagreement with old code.
#
#
# TITLE - InputCoo
# AUTHOR - Stephanie Monty
# PROJECT - mavisim
# CONTENTS:
#   1. This class returns an astropy table of the input coordinates of the sources, this can be used for comparisons with the DAOPhot detected positions
#
#		a) Star, Flux, X_True, PM_X, Y_True, PM_Y
#
#			i)  Star = Number assigned to the star to help with tracking through epochs in PM sims
#			ii) Flux = Monochromatic flux of the star (photons/s) passed to the program by the input data
#			iii CCD_Mapped_X/CCD_Mapped_Y       = True position following sub-pixel shifts due to sub-pixel positions, static distortion and proper motion
#												  Mapped to the FoR of the CCD (swap x and y)
#			iv) CCD_Mapped_PM_X/CCD_Mapped_PM_Y = Proper motion of the star (mas)
#												  Mapped to the FoR of the CCD (swap x and y)
#			v) X Static Dist/Y Static Dist = Distortion applied from the static distortion map in the x- and y-directions
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

# Astropy
from astropy.table import Table, Column


class InputCoo:

	"""
	Args:
		input_par = input parameters either hardcoded or altered by the user
		source = the astropy table that contains the source object generated in Source
		
	Returns:
		trimmed_cat = an astropy table containing the TRUE input positions (static distortion, sub-pixel positioning and proper motion taken into account)
					  and additional information necessary to compare with the DAOPhot output
					  *the catalogue may be trimmed to the 4k x 4k array or not*

	"""

	def __init__(self, input_par, source_table):
		self.input_par = input_par

		self.source_table = source_table


	def main(self):

		# Note we need to take two different coordinate conventions into account here, the convention of the CCD (x -> y) and the convention
		# of DAOPhot (starts at position 1, 1)

		# Find the full FoV in pixels, this contains the buffer (which is trimmed later) to contain light from stars outside the FoV
		full_fov_pix = int((self.input_par.MAVIS_fov + self.input_par.buffer)/self.input_par.ccd_sampling)

		# Convert the position to pixels (remove the knowledge of the static distortion and add the proper motion (if any))
		x_pos = np.array(np.around(((self.source_table["X"]/self.input_par.ccd_sampling) + full_fov_pix/2.0), 0), int)
		true_x = x_pos + self.source_table["X_Dist"] - self.source_table["Static_Dist"][:, 0] + self.source_table["X_PM"] + 1


		y_pos = np.array(np.around(((self.source_table["Y"]/self.input_par.ccd_sampling) + full_fov_pix/2.0), 0), int)
		true_y = y_pos + self.source_table["Y_Dist"] - self.source_table["Static_Dist"][:, 1] + self.source_table["Y_PM"] + 1

		# Create the final table
		input_coo = Table(data = ([self.source_table["Star"]]))
		input_coo.add_column(self.source_table["Flux"])
		input_coo.add_column(self.source_table["RA"])
		input_coo.add_column(self.source_table["Dec"])
		input_coo.add_column(Column(true_x), name="CCD_Mapped_X")
		input_coo.add_column(self.source_table["X_PM"], name="CCD_Mapped_PM_X")
		input_coo.add_column(self.source_table["Static_Dist"][:, 0], name="X Static Dist")
		input_coo.add_column(Column(true_y), name="CCD_Mapped_Y")
		input_coo.add_column(self.source_table["Y_PM"], name="CCD_Mapped_PM_Y")
		input_coo.add_column(self.source_table["Static_Dist"][:, 1], name="Y Static Dist")

		# Roughly trim the input catalogue to only include stars in the MAVIS FoV
		# Set the boundary to include stars just outside of the CCD limit (diagonal)
		image_size = (self.input_par.MAVIS_fov + self.input_par.buffer)/self.input_par.ccd_sampling

		r = np.sqrt(2 * (self.input_par.ccd_size/2.0)**2) + 100

		r_coo = np.sqrt((int(image_size/2.0) - input_coo["CCD_Mapped_X"])**2 + 
						(int(image_size /2.0) - input_coo["CCD_Mapped_Y"])**2)

		trim = np.where(r_coo <= r)[0]
		trimmed_cat = input_coo[trim]

		# Correct for the change from the larger buffered FoV to the FoV of the MAVIS CCD
		trimmed_cat["CCD_Mapped_X"] = trimmed_cat["CCD_Mapped_X"] - (image_size/2 - self.input_par.ccd_size/2)
		trimmed_cat["CCD_Mapped_Y"] = trimmed_cat["CCD_Mapped_Y"] - (image_size/2 - self.input_par.ccd_size/2)

		return (trimmed_cat)
