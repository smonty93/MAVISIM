# ----------------------------------------------------------------------------
#
# TITLE - add read noise
# AUTHOR - Stephanie Monty
# PROJECT - MAVISSimIm
# CONTENTS:
#   1. Fill this in
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
from mavisim.addconstantskypixels import add_constant_sky_pixel
import mavisim.input_parameters as input_par



def add_all_noise(image, exp_time):
    # Create the true image by randomly sampling Poissons with expectation intervals = N* + Nsky
    # Add the readnoise on top by randomly sampling Gaussians with a sigma = readnoise
    
    if input_par.noise_mode == "slow":
        # Read noise in photons
        read_noise = input_par.slow_rdnoise
        
    if input_par.noise_mode == "fast":
        # Read noise in photons
        read_noise = input_par.fast_rdnoise

    # Get the constant sky value in photons/pixel
    sky_value = add_constant_sky_pixel(exp_time)

    print (sky_value)

    # Remove the chance of negative N photons (only sky in that region), this is due to artifacts from the PSF creation
    image_addsky = image + sky_value
    image_addsky[np.where(image_addsky  < 0)] = sky_value

    # Account for the difference in photons between M1 and the imager
    image_addsky_mirror = image_addsky * input_par.AOMthruput * input_par.VLTthruput

    # Add the shot noise (from Poissonian stats)
    # we not sky limited on the pixel scale - Poisson statistics
    rng = np.random.default_rng()

    # Change to Poisson, and add throughput before drawing the Poisson statistics
    true_im = rng.poisson(lam=image_addsky_mirror)

    # Legacy: Gaussian stats
    #true_im = np.random.normal(loc=image_addsky_mirror, scale=np.sqrt(image_addsky))

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

    print (np.average(image_adu[2400:2460, 2100:2160]))
    print (np.std(image_adu[2400:2460, 2100:2160]))

    
    return (image_adu)