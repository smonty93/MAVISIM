# MAVISIM
Image simulating tool for the next generation ESO instrument MAVIS. If you use MAVISIM please cite <a href="https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2192M/abstract" target="_blank"> Monty et al. 2021.</a> 

## Current Version
**1.1 - EtE PSF** - more accurate representation of the PSF, suitable for photometric and astrometric science cases

Documentation can be found <a href="https://mavisim.readthedocs.io/en/latest/" target="_blank"> here</a> 

## Past Versions
1.0 - Fourier PSF - fast, most flexibility for astrometric modeling

Documentation can be found <a href="https://mavisim.readthedocs.io/en/v1.0.3/" target="_blank"> here</a>

## General Features
1. Can simulate full stellar field (eg. Milky Way Globular Cluster) from an input catalogue
2. Models three major sources of astrometric error introduced by the AO system:
    1. Tip-tilt residuals originating from uncorrected LO aberrations (dependent on the NGS constellation brightness and geometry) **field variable**
    2. HO aberrations originating from the LGS constellation and characteristics **field variable**
    3. Static field distortion originating from the AO module optics

## Current Limitations
1. Monochromatic images only
2. Only compatible with input **point source catalogues only**

# Future Versions

## MAVISIM 1.2: Broadband End-to-End PSF
Broadband images! Starting with the V band

## MAVISIM 2.0
Extended object version, not yet in development.
