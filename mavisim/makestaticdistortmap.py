# ----------------------------------------------------------------------------
#
# TITLE - static_distortion_map
# AUTHOR - Stephanie Monty
# PROJECT - mavisim
# CONTENTS:
#   1. 
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
Function to create the static distortion map interpolated functions (to return the distortion in x and y from the location)

dist_x,y = f(position(x), position(y))

	Returns: dist_x_func_degmm = interpolated function that returns the shift in x (in mm, need conversion factor)
			 dist_y_func_degmm = interpolated function that returns the shift in y (in mm, need conversion factor)
		
'''
__author__ = "Stephanie Monty"

### Imports

## Basic
import numpy as np

## Scipy
from scipy import interpolate

# Call this function once to create functions to find the distortion (shift in x & y) at any point

def make_static_dist_map(input_par):
    """
    Args:
        input_par = input parameters either hardcoded or altered by the user

    Returns:
        dist_x_func_degmm = returns shift in x for a given position
        dist_y_func_degmm = returns shift in y for a given position
    """

    # Field_y(deg) Hx Hy  Predicted_x(mm)  Predicted_y(mm)  Real_x(mm)  Real_y(mm)

    field_x = input_par.static_distort["Field_x(deg)"]
    field_y = input_par.static_distort["Field_y(deg)"]

    dist_x = input_par.dynamic_amp * (input_par.static_distort["Predicted_x(mm)"] - input_par.static_distort["Real_x(mm)"])
    dist_y = input_par.dynamic_amp * (input_par.static_distort["Predicted_y(mm)"] - input_par.static_distort["Real_y(mm)"])

    # Create an array of the field positions and the distortion at each pt (the difference)
    dist_all = np.empty([dist_x.shape[0], 4])
    dist_all[:, 0] = field_x
    dist_all[:, 1] = field_y
    dist_all[:, 2] = dist_x
    dist_all[:, 3] = dist_y

    grid_vals = np.unique(field_x)

    # Create grids to save the distortion (difference) in x and y at each point in the grid
    # This is necessary for the interpolation of the distortion (grid-wise interpolation)
    dist_x_grid = np.zeros([len(grid_vals),len(grid_vals)])
    dist_y_grid = np.zeros([len(grid_vals),len(grid_vals)])

    num = 0

    for x in grid_vals:
        sub_array = dist_all[np.where(dist_all[:, 0] == x), :][0]

        for row in np.arange(0, sub_array.shape[0]):

            dist_x_grid[num, row] = sub_array[row, 2]
            dist_y_grid[num, row] = sub_array[row, 3]

        num+=1

    dist_x_func_degmm = interpolate.RectBivariateSpline(grid_vals, grid_vals, dist_x_grid)
    dist_y_func_degmm = interpolate.RectBivariateSpline(grid_vals, grid_vals, dist_y_grid)
    
    return (dist_x_func_degmm, dist_y_func_degmm)