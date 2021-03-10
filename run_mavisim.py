# Standard
import numpy as np
import matplotlib.pyplot as plt; plt.ion()
import sys
import os
import time

# Astropy
from astropy.io import fits, ascii
from astropy.table import Table

from src.mavisim.generate_image import TileGenerator

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

array_size = 12000
pixsize = 3.75e-3

# x-position of bottom-left-corner of each pixel:
xx = np.arange(array_size)*pixsize
FoV = xx[-1]-xx[0]
xx -= (FoV/2+pixsize/2)
yy = xx.copy()
full_image = np.zeros([xx.shape[0],yy.shape[0]])

tile_gen = TileGenerator(source, "src/e2epsfs/e2e_psfs.fits", 34)

thresh = 1e-7
t1 = time.time()
t_tot = 0
for ni in range(len(source)):
    tile = tile_gen.get_tile(ni)
    xstart = np.abs(xx-tile[1][0]).argmin()
    ystart = np.abs(yy-tile[1][1]).argmin()
    if np.abs(xx[xstart]-tile[1][0]) > thresh or np.abs(yy[ystart]-tile[1][1]) > thresh:
        print(ni)
        print(xx[xstart])
        print(yy[ystart])
    full_image[ystart:ystart+tile_gen.psf_width_pix,xstart:xstart+tile_gen.psf_width_pix] += tile[0]
    dt = time.time()-t1
    t_tot += dt
    print(f"star: {ni:6d},  took: {dt:7.4f} sec,  eta {(len(source)-ni)*(t_tot/(ni+1)):7.1f} sec",end="\r")
    t1 = time.time()
print("\ndone")
print(f"time: {t2-t1:f}")
print(f"stars: {ni:d}")

cropped_image = full_image[np.abs(yy)<=30/2,:][:,np.abs(xx)<=30/2]
def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)
rebinned_image = rebin(cropped_image,[4000,4000])

