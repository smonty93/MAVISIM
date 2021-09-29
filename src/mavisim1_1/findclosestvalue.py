# ----------------------------------------------------------------------------
#
# TITLE - find closest value
# AUTHOR - Stephanie Monty
# PROJECT - MAVISSimIm
# CONTENTS:
#   1. 
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Function to find the closest match between an individual value and an array of values
'''
__author__ = "Stephanie Monty"

### Imports

## Basic
import numpy as np

def find_closest_value(ref, val):
	"""
	find closest value

	Find the closest match between a value and an array of values

	Args:
		ref - array of values to locate match
		val - value to find match for

	Returns:
		index - location of matched value in the reference array
	"""
	index = (np.abs(ref - val)).argmin()
	
	return (index)