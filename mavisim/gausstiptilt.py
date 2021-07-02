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
cov = np.array([[1.3, 0.9],     # Covariance is defined [[xx, xy],
                [ 0.9, 1.1]])   #                        [yx, yy]]

I = fillPixelsTT(P0,W,star,cov)
pp.imshow(I)

"""
def fill_pixels_tt(P0, M, N, W,star, cov):
    G = np.meshgrid(range(W[0]),range(W[1]))
    G[0] = G[0].flatten()
    G[1] = G[1].flatten()
    I = np.zeros(W)
    y = G[0]-P0[0]-star[0]
    x = G[1]-P0[1]-star[1]

    i = gauss(x,y,cov,star[2])
    I += np.reshape(i,[M,N])

    return I            
    
def gauss(x,y,cov,flux):

    z = np.array([x[:], y[:]]) # (M x N) x 2 array

    # final sigma = np.sqrt(sigma_tip**2 + sigma_tilt**2 + sigma_tt**2)
    gauss = np.exp(-(np.diag(z.T @ np.linalg.solve(cov,z)))/2)
   
    gauss_norm = gauss/np.sum(gauss)

    return (gauss_norm*flux)
