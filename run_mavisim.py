# Standard
import numpy as np
import matplotlib.pyplot as plt; plt.ion()
import os
import time

t_tot = 0
t1 = time.time()

# Astropy
from astropy.io import fits, ascii
from astropy.table import Table

cwd = os.getcwd()

os.chdir(os.getcwd() + "/src")
import mavisim.addnoise
import mavisim.input_parameters as input_par
from mavisim.generate_image import ImageGenerator
from mavisim.Source import Source
from mavisim.InputCoo import InputCoo

os.chdir(cwd)

ngc3201 = input_par.nbody_in
exp_time = 10
source = Source(ngc3201, exp_time, tt_amp_fac=1.0, static_dist=False, stat_amp=1.0).main()
input_coo = InputCoo(source, trim_cat=True).main()

image_gen = ImageGenerator(array_width_pix=12000, pixsize=3.75e-3, source_list=source,
                            psfs_file="src/e2epsfs/e2e_psfs_100s_lqg.fits", gauss_width_pix=34)
image_gen.main()
image_final = image_gen.get_rebinned_cropped(rebin_factor=2,cropped_width_as=30.0)

extent = 30*np.array([-0.5,0.5,-0.5,0.5])

image_allnoise = mavisim.addnoise.add_all_noise(image_final, source.meta["exp_time"])
image_final_noise = np.array(image_allnoise, dtype="float32")
hdu = fits.PrimaryHDU(image_final_noise)
hdul = fits.HDUList([hdu])
hdul.writeto("output_noisy.fits", overwrite=True)
