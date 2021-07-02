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
Function to ramp up the seeing wings
'''

## Basic
import numpy as np

def ramp_up(start, flex_pt, sz):
    """
    A simple ramp up function to blend wings with the AO core
    
    Because we're multiplying the functions together to apply the ramp, the ramp up
    then has the following form:
    
    R(r) = (r - flex_pt) + 1 [r < flex_pt]
         = 1                 [r >= flex_pt]
    """
    x = np.linspace(-1 * (sz/2.0), sz/2.0, sz)
    xy = np.meshgrid(x,x)
    
    x = np.sqrt(xy[0]**2 + xy[1]**2)
    
    m = 1/(flex_pt - start)
    
    ramp = np.maximum(0, (m * (x - start) + 1))
    
    #Replace all values above 1 with 1
    ramp[np.where(ramp >= 1)] = 1
    
    
    return (ramp)
