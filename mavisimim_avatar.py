# Run the MAVISSimIm On Avatar
#
# Taking the notebook commands and translating them to a scipt to run on Avatar

# Standard
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Astropy
from astropy.io import fits, ascii
from astropy.table import Table

# Project-specific/
sys.path.append('../src')
import mavissimim.rampup
import mavissimim.addconstantsky
import mavissimim.addnoise
import mavissimim.trimimage
import mavissimim.input_parameters as input_par

from mavissimim.Source import Source
from mavissimim.AOGaussGrid import AOGaussGrid
from mavissimim.InputCoo import InputCoo
from mavissimim.SeeingGrid import SeeingGrid

def run_mavissimim():
	exp_time = 30

	n_stack = 0

	# Paths to store the output coordinate files and images
	fits_save_path = ("/home/montys/AO/MAVISSimIm/AstrometryPaper/SimI/") #% str(input_par.nbody_yr))
	coo_save_path = ("/home/montys/AO/MAVISSimIm/AstrometryPaper/SimI/") #% str(input_par.nbody_yr))

	# Name of the final image
	file_name = "StaticPSF_11TT_NoMap_cat2_30G1FFFN.fits" #% (str(input_par.nbody_yr), str(n_stack))#str(input_par.nbody_yr)

	print (n_stack * exp_time)
	print (input_par.fv_psf_path)
	print (input_par.nbody_yr)	

	# Nbody input plus chosen exposure time
	ngc3201 = input_par.nbody_in

	# Make the source
	source = Source(ngc3201, exp_time, static_dist=False, stat_amp=1, tt_var=False, tt_amp_fac=1, tt_static=False, tt_kern=12).main()

	# Create the input catalogue
	input_coo = InputCoo(source).main()

	# Create the ao_field (which will be our final field in the case of the big EtE PSxF)
	(ao_field,gauss_field) = AOGaussGrid(source, fv_psf=False).main()

	# Create the seeing field
	seeing_field = SeeingGrid(gauss_field).main()

	# Make the image
	image = ao_field + seeing_field

	# Code for the option of stacking the images, remake the noise profile each time
	if n_stack != 0:
		# Create a list to store the images to stack (the arrays)
		stack_images = []

		for num in range(0, n_stack):
			# Add shot noise, read noise and convert from electrons to ADU
			image_stack = mavissimim.addnoise.add_all_noise(image, source.meta["exp_time"])

			stack_images.append(image_stack)

		image_allnoise = sum(stack_images)

	else:
		# Add shot noise, read noise and convert from electrons to ADU
		image_allnoise = mavissimim.addnoise.add_all_noise(image, source.meta["exp_time"])

	# Trim the final image
	(image_plusskyrn_trim, trimmed_input_coo) = mavissimim.trimimage.trim_image(image_allnoise, input_coo)

	# Options to save intermediate images (i.e. only the Gaussian field, the AO field or the Sky (wings) field)
	#(gauss_trim, dum) = mavissimim.trimimage.trim_image(gauss_field, input_coo)
	#(ao_trim, dum) = mavissimim.trimimage.trim_image(ao_field, input_coo)
	#(sky_trim, dum) = mavissimim.trimimage.trim_image(seeing_field, input_coo)

	# Save the final image
	final_image = np.array(image_plusskyrn_trim, dtype="float32")
	hdu = fits.PrimaryHDU(final_image)
	hdul = fits.HDUList([hdu])
	hdul.writeto(fits_save_path + file_name + ".fits", overwrite=True)

	 # Save an intermediate image
	#ao_image = np.array(ao_trim, dtype="float32")
	#hdu = fits.PrimaryHDU(ao_image)
	#hdul = fits.HDUList([hdu])
	#hdul.writeto(fits_save_path + file_name + "AO.fits", overwrite=True)

	# Save the coordinates
	ascii.write(trimmed_input_coo, coo_save_path + file_name + '.all', overwrite=True)

if __name__ == "__main__":
	run_mavissimim()





