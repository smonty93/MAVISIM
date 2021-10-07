# MAVISIM
Image simulating tool for the next generation ESO instrument MAVIS. 

## Versions
1.0 - Fourier PSF
1.1 - EtE PSF, faster

## Current Capabilities
1. Able to simulate full stellar field (eg. Milky Way Globular Cluster) from an input physical catalogue
2. Optimised for astrometric studies specifically
3. Models three major sources of astrometric error introduced by the AO system:
    1. Tip-tilt residuals originating from uncorrected LO aberrations (dependent on the NGS constellation brightness and geometry) **field variable**
    2. HO aberrations originating from the LGS constellation and characteristics **field variable**
    3. Static field distortion originating from the AO module optics

## Current Limitations
1. Monochromatic images only
2. Currently only compatible with an approximate model of the MAVIS PSF using the Fourier method
3. Only compatible with input **point source catalogues only**

# Features Under Development:
## MAVISIM 1.1: End-to-End PSF
Differences:
1. Uses an end-to-end PSF (1024x1024) for more realistic treatement of PSF and seeing wings
2. Uses analytical expression for Fourier transform of Gaussian to increase the speed
3. Main code loops over every star, calls JC's program that returns ifft(F(EtE PSF) x F(MV Gauss)), placed into large CCD array

## MAVISIM 1.2: Broadband End-to-End PSF
Differences:
1. Broadband images! Starting with the V band

## MAVISIM 2.0
Extended object version, not yet in development.
