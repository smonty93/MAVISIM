# ----------------------------------------------------------------------------
#
# TITLE - add read noise
# AUTHOR - Stephanie Monty
# PROJECT - mavisim
# CONTENTS:
#   1. Add photon noise (Poisson), sky background and read noise (Gaussian) to the image using the detector and sky characteristics included in input_par.
#   2. Handles the throughput of the VLT and Q.E. of the detector in the V band only for now.
#   3. Saturation is also dealt with using the detector characteristics in input_par
#
# ----------------------------------------------------------------------------
### Docstrings and metadata:
"""
Adding all noise to the image and converting from photons to electrons to ADU

"""

__author__ = "Stephanie Monty"

### Imports

## Basic
import numpy as np

## Project Specific
import mavisim.findclosestvalue
from mavisim.addconstantskypixels import add_constant_sky_pixel
from mavisim.trimimage import trim_image



def add_all_noise(input_par, image, exp_time):
    """
    Args:
        image = 2D noise-free image (ao_gauss_grid + seeing_grid)
        exp_time = total exposure time
        
    Returns:
        image_adu = final image + noise in ADU

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

    # Trim the image to 4k x 4k
    final_image = trim_image(input_par, image_adu)
    
    return (final_image)