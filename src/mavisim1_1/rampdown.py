# ----------------------------------------------------------------------------
#
# TITLE - ramp down
# AUTHOR - Stephanie Monty
# PROJECT - MAVISSimIm
# CONTENTS:
#   1. 
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Function to ramp down the PSF

'''
__author__ = "Stephanie Monty"

## Basics
import numpy as np


def ramp_down(stop, flex_pt, sz):
    """
    A simple ramp down function to blend the AO core with the Moffat wings
    
    Because we're multiplying the functions together to apply the ramp, the ramp down
    then has the following form:
    
    R(r) = - (r - flex_pt) + 1 [flex_pt < r <= 2 * flex_pit]
         = 1                   [r <= flex_pt]
    """
    x = np.arange(sz) - sz/2.0
    xy = np.meshgrid(x,x)
    
    x = np.sqrt(xy[0]**2 + xy[1]**2)
    
    m = -1/(stop - flex_pt)
    
    ramp = np.minimum(1, (m * (x - flex_pt)))
    
    #Remove all the negative values
    ramp[np.where(ramp <= 0)] = 0

    return (ramp)