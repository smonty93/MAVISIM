# Running MAVISIM 1.0: The Short Version

## Specifications
To create an image of a point source catalogue with the default settings, only a subset of settings must be specified. In the `input_parameter.py` file please edit the following:

- `path_to_mavisim`: full path to the location of the mavisim directory
- `path_to_data`: full path to the location of the mavisim data files (e.g. the PSF database)
- `input_file`: name of the input source catalogue (e.g. nbody file)
- `fv_psf_path`: name of the PSF database (assuming it's stored in the data directory)
- `filter`: the closest broad band filter to the simulation wavelength (recall Version 1.0 is monochromatic), choices are UBVRI
- `psf_wavelength`: PSF database wavelength
- `static_psf`: name of the preferred PSF to use for the static PSF case (no field variability)
- `tt_residual_map`: name of the tip-tilt residual map to use for the spatially variable tip-tilt case
- `tt_kernel`: FWHM in mas of a user-specified tip-tilt kernel

The other settings to specify are made when calling the `Source` and `AOGaussGrid` classes. These are as follows:

In the `Source` class:

- **`static_dist`** controls whether static distortion from the MAVIS optics is included (=`True`) in the simulation or not (=`False`). Default value is `False`.
- **`stat_amp`** controls the magnitude of the static distortion. Change this value to increase the distortion by a scalar in both the x- and y-direction. Default value is 1.0 (no amplification)
- **`tt_var`** controls whether the *spatial variability* of the tip-tilt residual error is included (=`True`) or not (=`False`). The choice of map is specified in the `input_parameters` file. Default value is `False`.
- **`user_tt`** specifies whether the user would like **fix the tip-tilt residual error** as a specific value described by a single Gaussian kernel. This option only works when the tip-tilt residual is *static* (`tt_var=False`). The user specifies the value of the residual in the `input_parameter` file as the `tt_kernel` parameter. Default value is `False`.
- **`tt_amp`** controls the magnitude of the tip-tilt residual error. Default value is 1.0 (no amplification)

In the `AOGaussGrid` class:

- **`fv_psf`** controls whether the high-order PSF is spatially variable or not. This is one of the main characteristics that makes MCAO images unique. 

## Commands to Create an Image
The following set of commands will create a MAVISIM image. These can also be found in the `/example` directory on the <a href="https://github.com/smonty93/MAVISIM" target="_blank">MAVISIM Github repository.</a>

<div class="admonition note">
<p class="admonition-title">Creating a MAVISIM image</p>
<p> 
	# Load the input source file <br>
	<code>>>> glob_clust = input_par.input_file </code> <br>
	<code>>>> exp_time = 30 # seconds </code>
</p>
<p> 
	# Create the Source object with default static distortion and spatially variable tip-tilt error<br>
	<code>>>> source = Source(glob_clust, exp_time, static_dist=True, tt_var=True).main()</code>
</p>
<p> 
	# Create the AOGrid and GaussGrid objects with a spatially variable high-order PSF<br>
	<code>>>> (ao_field, gauss_field) = AOGaussGrid(source, fv_psf=True).main()</code>
</p>
<p> 
	# Create a Noise-Free Image<br>
	<code>>> image = ao_field + seeing_field </code>
</p>
<p> 
	# Create the Final Image<br>
	<code>>> final_image = add_all_noise(image, source.meta["exp_time"])</code>
</p>
</div>

