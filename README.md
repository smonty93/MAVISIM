# MAVISIM1.0
Image simulating tool for the next generation ESO instrument MAVIS. 

## Current Capabilities
1. Able to simulate full stellar field (eg. Milky Way Globular Cluster) from an input physical catalogue
2. Optimised for astrometric studies specifically
3. Models three major sources of astrometric error introduced by the AO system
  1. Tip-tilt residuals originating from uncorrected LO aberrations (dependent on the NGS constellation brightness and geometry) **field variable**
  2. HO aberrations originating from the LGS constellation and characteristics **field variable**
  3. Static field distortion originating from the AO module optics

## Current Limitations
1. Monochromatic images only
2. Currently only compatible with an approximate model of the MAVIS PSF using the Fourier method
3. Only compatible with input **point source catalogues only**

# MAVISIM 2.0
Currently being developed using pieces of MAVISIM1.0 with several improvements.

## Planned Improvements
1. Broadband images
2. EtE PSF database to replace Fourier PSFs

## Large-Scale Structural Changes
1. PSF database to move to database of wavefronts, taking the wavelength and exposure time information from the user.
Con't Here