# MAVISIM

<!--- These are examples. See https://shields.io for others or to customize this set of shields. You might want to include dependencies, project status and licence info here --->

![GitHub contributors](https://img.shields.io/github/contributors/smonty93/mavisim) ![License](https://img.shields.io/github/license/smonty93/mavisim)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mavisim) ![PyPI version](https://img.shields.io/pypi/v/mavisim)
![Tests](https://github.com/smonty93/mavisim/actions/workflows/tests.yml/badge.svg)

Image simulating tool for the next generation ESO instrument MAVIS. If you use MAVISIM please cite <a href="https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2192M/abstract" target="_blank"> Monty et al. 2021.</a>

## Getting started
For full documentation, see <a href="https://mavisim.readthedocs.io/en/latest/" target="_blank">mavisim.readthedocs.io</a>. To get started, install via pip:
```bash
pip install mavisim
```
and use the CLI:
```bash
mavisim --help
```
which should output something like:
```
usage: mavisim [-h] [-o OUTPUT] [-d DISTS] [--static STATIC] [--header HEADER] field psfs

mavis image simulator

positional arguments:
  field                 file containing source field to simulate
  psfs                  specify PSFs fits file on disk to use in simulation

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        specify name for output fits file
  -d DISTS, --dists DISTS
                        specify distortion file on disk to pre-distort sources with
  --static STATIC, -s STATIC
                        use only 1 psf (not field varying), at the given psf fits hdu index
  --header HEADER       dict to be parsed as header for saved fits file

```

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
