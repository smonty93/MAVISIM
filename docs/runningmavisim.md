# Running MAVISIM 1.0: The Long Version
Version 1.0 of MAVISIM is designed primarily to explore astrometric science cases as was done in Monty et al. 2021 (link paper). This is reflected in the functionality with many different options for the user to control the astrometric error terms. As such, we will go through the main classes, or objects, in MAVISIM 1.0 to help the user understand how to control each error term and grasp the basics of each component of the code.

## Parameter File
The parameter file primarily contains hardcoded instrument, telescope, detector and sky background characteristics. It's not neccesary to change most of them. However, there are a few parameters at the beginning of the file that allow the user to control the simulation output. These are as follows; 

- `path_to_data`: full path to the location of the mavisim data files (e.g. the PSF database)
- `input_file`: name of the input source catalogue (e.g. nbody file)
- `fv_psf_path`: name of the PSF database (assuming it's stored in the data directory)
- `filter`: the closest broad band filter to the simulation wavelength (recall Version 1.0 is monochromatic), choices are UBVRI
- `psf_wavelength`: PSF database wavelength
- `static_psf`: name of the preferred PSF to use for the static PSF case (no field variability)
- `tt_residual_map`: name of the tip-tilt residual map to use for the spatially variable tip-tilt case
- `tt_kernel`: FWHM in mas of a user-specified tip-tilt kernel

## Input File
MAVISIM 1.0 can **only** parse an input list of point sources at this time (stay tuned for MAVISIM 2.0 to test extended object science cases). The input file is parsed as a Astropy Table and must contain the following columns. Columns with essential information are bolded. We have provided two example catalogues on GitHub in the /example directory to test MAVISIM.

| Column    | Units | Description |
| --------- | ----- | ----------- |
| Star | n/a | Not essential (can set to 0) - helps track the stars across multiple epochs during proper motions studies|
| RA   | degrees | Not essential (can set to 0) - helps with deprojection effects if desired |
| Dec  | degrees | Not essential (can set to 0) - helps with deprojection effects if desired |
| **X**    | **arcseconds** | **X-distance from the centre of the field set as (0", 0")** |
| PM_X | mas/year | Not essential (can set to 0) - can be projected or de-projected proper motion |
| **Y**    | **arcseconds** | **Y-distance from the centre of the field set as (0", 0")** |
| PM_Y |  mas/year | Not essential (can set to 0) - can be projected or de-projected proper motion |
| **Flux** | **photons/s** | **Used to scale the Gaussian representation of each star, multipled by the exposure time for final flux** |

### Relevant Command
```python
#Loading the input catalogue and specifying the exposure time
glob_clust = input_par.input_file
exp_time = 30
```
## Source Class

This where the user can specify the bulk of the astrometric error terms along with the exposure time. Both the tip-tilt residual errors associated with the natural guide star constellation and static distortion introduced by the MAVIS optics are tuned in this object. The remaining high-order PSF spatial variability is dealt with in the AOGaussGrid object. The `Source` class handles all the user specifications and returns a source object that is used in tandem with the MAVIS Fourier PSF to create the final image. The user has the option to change the following keywords:  `static_dist`, `stat_amp`, `tt_var`, `user_tt` and `tt_amp`. The following is a description of each keyword, including the default values.

- **`static_dist`** controls whether static distortion from the MAVIS optics is included (=`True`) in the simulation or not (=`False`). The distortion is included as a shift in the (x, y) position of the centroid of each star. A provisional distortion map from the MAVIS optical design is used to generate a map of the distortion across the field. This is used to recover the specific (x, y) distortion for each star. Note that MAVISIM deals with sub-pixel shifts. **Default value is `False`.**<p>&nbsp;</p>
- **`stat_amp`** controls the magnitude of the static distortion. Change this value to increase the distortion by a scalar in both the x- and y-direction. **Default value is 1.0 (no amplification)**<p>&nbsp;</p>
- **`tt_var`** controls whether the *spatial variability* of the tip-tilt residual error is included (=`True`) or not (=`False`). A map of the spatial dependency of the tip-tilt residual error (expressed as a multivariate Gaussian) is used to extract the tip-tilt for each star. MAVISIM 1.0 comes with three tip-tilt residual maps created using different natural guide star constellations and characteristics. The choice of map is specified in the `input_parameters` file. **Default value is `False`.**<p>&nbsp;</p>
- **`user_tt`** specifies whether the user would like **fix the tip-tilt residual error** as a specific value described by a single Gaussian kernel. This option only works when the tip-tilt residual is *static* (`tt_var=False`). The user specifies the value of the residual in the `input_parameter` file as the `tt_kernel` parameter. The `tt_kernel` is assumed to be the FWHM of the tip-tilt residual kernel in units of milli-arcseconds. By default MAVISIM 1.0 includes a charge diffusion and vibrational term in the creation of the tip-tilt residual error, if this option is used these terms will be neglected. **Default value is `False`.**<p>&nbsp;</p>
- **`tt_amp`** controls the magnitude of the tip-tilt residual error. Change this value to scale the tip-tilt residual error by some fixed amount. If this value is set to 0, only the charge diffusion and vibrational terms will be included in the determination of the tip-tilt residual. **Default value is 1.0 (no amplification)**

### Relevant Command
```python
# Creating a source object with the default built-in static distortion and 
# tip-tilt residual errors
source = Source(input_par, glob_clust, exp_time, static_dist=True, tt_var=True).main()
```


## AOGaussGrid Class
This is where a unique PSF is created for each star and placed into a 40 arcsecond array in preperation for creating the final image. **The user controls whether the high-order PSF is spatially variable or not through the keyword `fv_psf`.** The spatially variability of the PSF is a major characteristic of MCAO images. To create each star-specific PSF, the positional information in the `source` object is used to recover the closest four PSFs from a grid of 11 x 11 MAVIS Fourier PSFs sampling the science field of view. The four PSFs are then interpolated using a bilinear interpolation to create the final high-order PSF for each star. The high-order PSF is then convolved with the multivariate Gaussian for each star stored in the `source` object to capture the terms described in the previous section. The `AOGaussGrid` class returns two things, i) a 40 arcsecond array of stars represented by unique PSFs truncated at the AO control radius `ao_field` and ii) a 40 arcsecond array of multivariate Gaussians capturing the tip-tilt and static distortion information only, the `gauss_field`. If you have no interest in the information contained in the seeing wings, you could work with the `ao_field` only. 

### Relevant Command
```python
# Creating an image with a spatially variable high-order PSF
(ao_field, gauss_field) = AOGaussGrid(input_par, source, fv_psf=True).main() 
```

## SeeingGrid Class
The `SeeingGrid` class creates a complementary field to the `ao_field` by convolving the grid of multivariate Gaussians `gauss_field` created by `AOGaussGrid` with a single, large Fourier PSF. Before convolving the two, the information within the AO control radius of the large Fourier PSF is removed and the edges are ramped to avoid discontinuities. This is primarily a speed-saving measure to reduce the time required for convolusions, something that becomes significant in crowded fields. If you examine the `seeing_field` you'll see a stellar field composed only of seeing wings with dark centres. 

### Relevant Command
```python
# Creating a noise-free image
image = ao_field + seeing_field
```

## Creating the Final Image: Adding Detector Characteristics and Sky Background
The `add_all_noise` function handles adding noise and sky background to the image as well as converting from photons to ADU. Sky background is added first, system throughput is then accounted for then shot (Poisson) noise and read (Gaussian) noise is added to the noise-free image (<code>ao_field</code> + <code>seeing_field</code>). The final image (+noise) is converted first to electrons accounting for the detector quantum efficiency then to ADU assuming a detector gain and saturation point. Saturated stars are capped.

The final image (+noise) is larger than the actual MAVIS science field of view and must be trimmed down to simulate the 4k x 4k detector. The image is initially larger by `input_par.buffer` arcseconds to capture the effects of stray light from stars outside the MAVIS field. This is done in an effort to be more realistic. The `trim_image` function is then used to trim the final image to the correct size. It can also be used to trim a catalogue of the input source positions if desired (see the next section).

### Relevant Command
```python
# Creating the final image
final_image = add_all_noise(input_par, image, source.meta["exp_time"])
```

<div class="admonition tip">
<p class="admonition-title">Stacking Images</p>
<p>To simulate a stacked image simply loop this function, a unique noise profile will be added every time. Be wary of RAM though! These arrays are large. </p>
</div>

## *Optional Class:* InputCoo
This is an optional class that can be used to format an input catalogue of source positions (in pixels) to feed to a photometric software. This catalogue can be used to perform forced photometry for example. Thus far it has been tested on DAOPhot-IV (Stetson 1987, 1994). The final output catalogue is:

|Star| Flux | RA | Dec | CCD_Mapped_X | CCD_Mapped_PM_X | X | Static Dist X | CCD_Mapped_Y | CCD_Mapped_PM_Y | Y | Static Dist Y |
| -- | ---- | -- | --  |-- |-- |-- |-- |-- |-- |-- |-- |
| .. | photons | degrees | degrees | pixels |  pixels | arcseconds | pixels | pixels | pixels | arcseconds | pixels |

### Relevant Command
```python
# Creating a catalogue of the input source positions</p>
input_coo = InputCoo(input_par, source).main()
```





