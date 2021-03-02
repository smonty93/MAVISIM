# ----------------------------------------------------------------------------
#
# TITLE - weight wings
# AUTHOR - Stephanie Monty
# PROJECT - MAVISSimIm
# CONTENTS:
#   1. 
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Function to interpolate the weights for the seeing wings and apply them following convolution with the gaussian grid.

The weights are used to bring the seeing wings up to the height of the AO core and to help with blending.

'''
__author__ = "Stephanie Monty"

### Imports

# Scipy
from scipy import interpolate


def interpolate_weights(weight_table):

	weight_func = interpolate.interp2d(weight_table["X"], weight_table["Y"], weight_table["Seeing_Weight"])

	return weight_func


