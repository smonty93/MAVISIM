#


## PSF
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/generate_image.py/#L5)
```python 
PSF(
   fits_ext, padto, *, dtype = np.complex128
)
```


---
PSF object class.

Helper class to handle PSF data.


**Args**

* **fits_ext** (`astropy.io.fits`) : opened fits file containing PSF data and required
* **padto** (`list` of `int`) : shape of desired Fourier transform array.
* **dtype** (`np.dtype`, optional) : desired complex dtype of Fourier array.
header data (XPOS,YPOS,LAMBDA).


**Attributes**

* **fft_data** (`np.ndarray`) : array storing the (minimal) rfft2 data of the given PSF.
* **xpos** (`float`) : x-position of PSF point source in arcsec.
* **ypos** (`float`) : y-position of PSF point source in arcsec.
* **Lambda** (`float`) : wavelength used to capture the PSF.


----


## TileGenerator
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/generate_image.py/#L31)
```python 
TileGenerator(
   source, psfs_file, gauss_width_pix, *, dtype = np.complex128, which_psf = None
)
```


---
Object for generating tiles to be sliced into final image.


**Args**

* **source** (`Source` object) : Object containing all of the source data, as defined in `Source.py`.
* **psfs_file** (`str`) : path to `.fits` file containing all PSFs and metadata.
* **gauss_width_pix** (`int`) : support size in pixels to build the Gaussian star kernel.
* **dtype** (`np.dtype`, optional) : complex data type to work with in DFT space.
* **which_psf** (`int`, optional) : If specified, use fits HDU[which_psf+1] PSF only.


**Attributes**

* **source** (`Source` object) : collection of sources in `Source`-type object.
* **psfs** (`list` of `PSF` objects) : list of PSF objects (see PSF class).
* **pixsize** (`float`) : pixel size in arcsec used to build tile.
* **gauss_width_pix** (`int`) : width of Gaussian square support in pixels
* **static** (`bool`) : if `True`, use a static PSF, otherwise use field variable PSF.



**Methods:**


### .get_effective_psf_fft
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/generate_image.py/#L108)
```python
.get_effective_psf_fft(
   s_pos
)
```

---
Takes star information and computes effective PSF.

From star position, the convex combination PSF is found
(equivalent to bilinear interpolation since stars are defined
on square grid). The resulting effective PSF is added to the 
internal _psf_array to be used in the get_tile pipeline.


**Args**

* **star_pos** (`np.ndarray`) : position [arcsec].


### .get_tile
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/generate_image.py/#L141)
```python
.get_tile(
   index
)
```

---
Get the tile corresponding to source[index]

From the tile_generator object tgen, calling tgen.get_tile(index) will
generate the tile corresponding to the tgen.source_pos[index] star by 
interpolating the 4 neighbouring PSFs and convolving this effective
PSF with a sub-pixel shifted Dirac-delta function, and if requested, a 
Gaussian kernel defined by tgen.cov_mat .

The output of this is a tile which has been trimmed down to the input
PSF dimensions, as well as the coordinates of the bottom-left-corner
of the tile so that it may be sliced into the final image properly.


**Args**

* **index** (`int`) : index of star in source table to generate tile for.


**Returns**

* **out** (real-valued `np.ndarray`) : tile to be sliced into final image
* **bottom_left_corner** (`np.ndarray`) : coordinates of bottom left corner of bottom left pixel in arcsec


### .get_star_kernel_fft
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/generate_image.py/#L217)
```python
.get_star_kernel_fft(
   flux, mu
)
```

---
Compute star Gaussian based in DFT space.

Directly computes the FFT of the Gaussian kernel with appriate amplitude, 
width, and offset to suit the tile being generated.

Uses optimised np.einsum so requires running `optimize_star_kernel()` first.


**Args**

* **flux** (`float`) : flux of star.
* **mu** (`np.ndarray`) : position of star.


**Returns**

* **gaussian_fft** (`np.ndarray`) : star Gaussian kernel in FFT space.


### .optimize_star_kernel
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/generate_image.py/#L257)
```python
.optimize_star_kernel()
```

---
Runs star kernel once to optimise `np.einsum`


----


## ImageGenerator
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/generate_image.py/#L279)
```python 
ImageGenerator(
   array_width_pix, source, psfs_file, pixsize = 0.00375, gauss_width_pix = 34,
   which_psf = None
)
```


---
Generate image from sliced tiles, one tile per source object.

This is the core object to work with when generating a MAVISIM image from
a `Source` object.


**Args**

* **array_width_pix** (`int`) : width of full image in pixels before rebinning.
* **pixsize** (`float`) : pixel size in arcsec before rebinning.
* **source** (`Source` object) : source list as `Source`-type object
* **psfs_file** (`str`) : filename for fits file containing all PSFs and PSF metadata.
* **gauss_width_pix** (`int`) : width of Gaussian square support in pixels.
* **which_psf** (`int`, optional) : If specified, use fits HDU[which_psf+1] PSF only.


**Attributes**

* **pixsize** (`float`) : pixel size in arcsec of image.
* **fov** (`float`) : FoV of full image.
* **full_image** (real-valued `np.ndarray`) : final image at original pixel size (i.e., before rebinning).
* **tile_gen** (`TileGenerator` object) : tile generator object used to create tiles to slice into final image.



**Methods:**


### .main
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/generate_image.py/#L312)
```python
.main()
```

---
Loop over all stars and add the tile to the full image.


### .get_rebinned_cropped
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/generate_image.py/#L325)
```python
.get_rebinned_cropped(
   rebin_factor, cropped_width_as
)
```

---
Rebin self.full_image after cropping to desired rebin factor.


**Args**

* **rebin_factor** (int) : rebinning factor from high-res image to rebinned image. 
* **cropped_width_as** (float) : desired width of final image in arcsec.
Note that no checking is done on the validity of this, so use with care.

**Returns**

* **rebinned_im** (real-valued `np.ndarray`) : complete image, rebinned and cropped. 

