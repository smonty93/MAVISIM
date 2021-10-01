#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 13:34:37 2020

@author: jcranney
"""

import numpy as np
import matplotlib.pyplot as pp

"""
Example run:
M = 5                        # Image height [pixels]
N = 5                        # Image width  [pixels]

W = np.array([M,N])

P0 = (M-1)/2*np.array([1,1])

star = np.array([  0,   0, 1000, 0.5])
                [  y,   x, flux, one_sigma]

I = fillPixelsTT(P0,W,star,cov)
pp.imshow(I)

"""
def fill_pixels(P0, M, N, W,star):
    G = np.meshgrid(range(W[0]),range(W[1]))
    G[0] = G[0].flatten()
    G[1] = G[1].flatten()
    I = np.zeros(W)
    y = G[0]-P0[0]
    x = G[1]-P0[1]
    r = np.sqrt((y-star[0])**2+(x-star[1])**2)
    i = gauss(r,star[2],star[3])
    I += np.reshape(i,[M,N])
    return I
            
    
def gauss(r,flux,width):

    gauss = np.exp(-((r/width)**2)/2)

    # normalise to 0 - 1
    #gauss_scale = (gauss - np.amin(gauss))/(np.amax(gauss) - np.amin(gauss))
    gauss_norm = gauss/np.sum(gauss)

    return (gauss_norm * flux)
