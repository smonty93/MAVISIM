# ----------------------------------------------------------------------------
# 
# TITLE - util
# AUTHOR - Stephanie Monty
# PROJECT - MAVISIM
# 
# ----------------------------------------------------------------------------

### Imports
# Standard
import numpy as np
# Astropy
from astropy.table import Table, Column
from astropy import units as u
# Scipy
from scipy import interpolate

def input_coo(input_par, source):
	""" Take input parameters and parse them as an astropy table for comparison with DAOphot.

	Args:
		input_par: input parameters, e.g., from `input_parameters.py`.
		source (`Source` object): the `Source` object containing all of the sources simulated.
		
	Returns:
		trimmed_cat = an astropy table containing the TRUE input positions (static distortion, 
		sub-pixel positioning and proper motion taken into account) and additional information necessary to compare with the DAOPhot output.
	"""
	# Note we need to take two different coordinate conventions into account here, the convention of the CCD (x -> y) and the convention
	# of DAOPhot (starts at position 1, 1)

	# Find the full FoV in pixels, this contains the buffer (which is trimmed later) to contain light from stars outside the FoV
	full_fov_pix = int(input_par.MAVIS_fov + input_par.buffer)/input_par.ccd_sampling

	# Convert the position to pixels (remove the knowledge of the static distortion and add the proper motion (if any))
	x_pos = np.array(np.around(((source.x_pos/input_par.ccd_sampling) + full_fov_pix/2.0), 0), int)
	true_x = x_pos + source.x_dist - source.static_dist[:, 0] + source.x_pm + 1

	y_pos = np.array(np.around(((source.y_pos/input_par.ccd_sampling) + full_fov_pix/2.0), 0), int)
	true_y = y_pos + source.y_dist - source.static_dist[:, 1] + source.y_pm + 1

	# Create the final table
	input_coo = Table(data = (source.star,),names=("Star",))
	input_coo.add_column(Column(source.flux),name="Flux")
	input_coo.add_column(Column(source.ra),name="RA")
	input_coo.add_column(Column(source.dec),name="Dec")
	input_coo.add_column(Column(true_x), name="CCD_Mapped_X")
	input_coo.add_column(Column(source.x_pm), name="CCD_Mapped_PM_X")
	input_coo.add_column(Column(source.static_dist[:, 0]), name="X Static Dist")
	input_coo.add_column(Column(true_y), name="CCD_Mapped_Y")
	input_coo.add_column(Column(source.y_pm), name="CCD_Mapped_PM_Y")
	input_coo.add_column(Column(source.static_dist[:, 1]), name="Y Static Dist")

	# Roughly trim the input catalogue to only include stars in the MAVIS FoV
	# Set the boundary to include stars just outside of the CCD limit (diagonal)
	image_size = (input_par.MAVIS_fov + input_par.buffer)/input_par.ccd_sampling

	r = np.sqrt(2 * (input_par.ccd_size/2.0)**2) + 100

	r_coo = np.sqrt((int(image_size/2.0) - input_coo["CCD_Mapped_X"])**2 + 
					(int(image_size /2.0) - input_coo["CCD_Mapped_Y"])**2)

	trim = np.where(r_coo <= r)[0]
	trimmed_cat = input_coo[trim]

	# Correct for the change from the larger buffered FoV to the FoV of the MAVIS CCD
	trimmed_cat["CCD_Mapped_X"] = trimmed_cat["CCD_Mapped_X"] - (image_size/2 - input_par.ccd_size/2)
	trimmed_cat["CCD_Mapped_Y"] = trimmed_cat["CCD_Mapped_Y"] - (image_size/2 - input_par.ccd_size/2)

	return trimmed_cat

def add_all_noise(input_par, image, exp_time):
    """ Adds all noise effects onto noiseless image.

    Args:
        input_par: input parameters object.
        image (`np.ndarray`): 2D noise-free image.
        exp_time (`float`): total exposure time in seconds.
        
    Returns:
        image_adu (`np.ndarray`): final image + noise in ADU
    """

    # Create the true image by randomly sampling Poissons with expectation intervals = N* + Nsky
    # Add the readnoise on top by randomly sampling Gaussians with a sigma = readnoise
    
    if input_par.noise_mode == "slow":
        # Read noise in photons
        read_noise = input_par.slow_rdnoise
        
    if input_par.noise_mode == "fast":
        # Read noise in photons
        read_noise = input_par.fast_rdnoise

    # Get the constant sky value in photons/pixel
    sky_value = add_constant_sky_pixel(input_par, exp_time)

    # Remove the chance of negative N photons (only sky in that region), this is due to artifacts from the PSF creation
    image_addsky = image + sky_value
    image_addsky[np.where(image_addsky  < 0)] = sky_value

    # Account for the difference in photons between M1 and the imager
    image_addsky_mirror = image_addsky * input_par.AOMthruput * input_par.VLTthruput

    # Add the shot noise (from Poissonian stats)
    # We are not sky limited on the pixel scale - Poisson statistics
    rng = np.random.default_rng()

    # Draw the image in photons + photon noise
    true_im = rng.poisson(lam=image_addsky_mirror)

    # Convert from photons to electrons
    im_elec = true_im  * input_par.QE

    # Add Read Noise
    im_addnoise = im_elec + np.random.normal(loc=np.zeros([true_im.shape[0], true_im.shape[0]]), scale=read_noise)

    # Convert the final image from electrons to ADU
    # Cap at the saturation point
    if input_par.sat_point  == 0:
        sat_point = np.int(2**input_par.bit_depth - 1) * input_par.gain

    else:
        sat_point = input_par.sat_point

    # Convert from electrons to ADU, round to integers to be physical
    image_adu = (im_addnoise * input_par.gain).astype(np.int)
    image_adu[image_adu > sat_point] = sat_point
	
    return (image_adu)

def add_constant_sky_pixel(input_par, exp_time):
	""" Add constant-valued sky background.
	Args:
        input_par: input parameters object.    
		exp_time (`float`): exposure time in seconds to convert from photons/s to photons.
		
	Returns:
		sky_value (`float`): a global sky value in photons to add to every pixel
	""" 

	square_arcsec_pervoxel = (input_par.ccd_sampling**2) * u.arcsec**2

	# Assuming the surface brightness is passed in as mag/arcsec^2 do the following:
	flux_jy = (3631 * 10**(-1*input_par.surf_bright/2.5))*u.Jy*u.nm**(-1)
	flux_ph_s_nm_cm2 = flux_jy * (1.51e3/input_par.psf_wavelength)
	flux_ph_s_nm_m2 = flux_ph_s_nm_cm2 * 10**4
	flux_ph_s_m2 = flux_ph_s_nm_m2 * input_par.filt_width #in nm
	sky_value = flux_ph_s_m2 * square_arcsec_pervoxel  * exp_time * input_par.collecting_area

	return sky_value.value

def make_static_dist_map(input_par):
    """ Make static distortion map from loaded distortion samples in `input_par`.

    Args:
        input_par: input parameters object.    

    Returns:
        dist_x_func_degmm:  
        dist_y_func_degmm:  
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