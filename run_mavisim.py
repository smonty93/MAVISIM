# Standard
import numpy as np
import matplotlib.pyplot as plt; plt.ion()
import sys
import os

# Astropy
from astropy.io import fits, ascii
from astropy.table import Table

from src.mavisim.generate_image import ImageGenerator

# Project-specific/
#sys.path.append('/src')
cwd = os.getcwd()
os.chdir(os.getcwd() + "/src")

#import mavisim.rampup
#import mavisim.addconstantsky
#import mavissimim.addnoise
#import mavissimim.trimimage
import mavisim.input_parameters as input_par
#from rebin import rebin

from mavisim.Source import Source
#from mavissimim.AOGaussGrid import AOGaussGrid
#from mavissimim.SeeingGrid import SeeingGrid
from mavisim.InputCoo import InputCoo

os.chdir(cwd)
# Nbody input plus chosen exposure time
ngc3201 = input_par.nbody_in

exp_time = 3600
source = Source(ngc3201, exp_time, tt_amp_fac=1.0, static_dist=False, stat_amp=1.0).main()
#self, mavis_src, exp_time, static_dist, tt_amp_fac, stat_amp

image_gen = ImageGenerator(source[:100], "src/e2epsfs/e2e_psfs.fits", 
        [-15,15,-15,15], debug_scale_factor=16)

