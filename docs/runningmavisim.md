# Running MAVISIM: The Long Version
The current iteration of MAVISIM (v1.1.x) is mostly driven by the desire to simulate realistic images with astrometric statistics as close as possible to what MAVIS will have on-sky. This version is still limited to monochromatic images, and indeed the default database associated with MAVISIM is still exclusively at $\lambda=550$ nm.

## Parameter File
For a scientist to perform their own astrometric simulations using MAVISIM, they will be required to modify the `input_parameter.py` file. In the current iteration of MAVISIM, the only parameters that are likely to be modified are:

- `input_cat` : parsed input source catalogue (e.g. nbody file)
- `vib_term`/`cd_term` : uncorrected vibration and charge diffusion to be simulated,
- `static_distort` : parsed distortion field across field of view,
- `plate_scale`/`dynamic_amp` : plate scale and distortion scale factors (e.g., to simulate sensitivity of astrometry).
- `surf_bright` : sky background magnitude per acrsec^2

## `input_cat`
MAVISIM can **only** parse an input list of point sources at this time (stay tuned for MAVISIM $\geq$ 2.0 to test extended object science cases). The input file should contain the following columns. We have provided two example catalogues in the downloaded data.

| Column    | Units  | Description |
| --------- | ------ | ----------- |
| Star      | n/a    | Star ID integer. Helps track the stars across multiple epochs during proper motions studies|
| X         | arcsec | X-distance from the centre of the field set as (0", 0") |
| Y         | arcsec | Y-distance from the centre of the field set as (0", 0") |
| Flux      | phot/s | Used to scale the Gaussian representation of each star, multipled by the exposure time for final flux |
| RA        | deg    | Not essential (can set to 0) - Helps with deprojection effects |
| Dec       | deg    | Not essential (can set to 0) - Helps with deprojection effects if desired |
| PM_X      | mas/yr | Not essential (can set to 0) - can be projected or de-projected proper motion |
| PM_Y      | mas/yr | Not essential (can set to 0) - can be projected or de-projected proper motion |
