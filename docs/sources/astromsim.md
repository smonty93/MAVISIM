#


## AstromCalibSimGeneric
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L11)
```python 
AstromCalibSimGeneric(
   static_distort, *, pixel_size_as = 0.00736, pixel_size_mm = 0.01, dist_amp = 1.0,
   mask_scale = 1000.0/582, hole_position_std = 0.0, dx = 0.2, dy = 0.2, dx_meas = None,
   dy_meas = None, n_poly = 6, pin_pitch = 0.5, num_pin_x = 40
)
```


---
Generic class for astrometric calibration simulation.
See: AstromCalibSimAna for the analytical simulation
AstromCalibSimE2E for the end-to-end simulation


**Methods:**


### .input_dist
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L51)
```python
.input_dist(
   x, y
)
```

---
Evaluate interpolated input distortions at arbitrary coordinates.

`x` and `y` (in arcsec) can be anywhere in the science field, but must
be array-like and the same size.


**Args**

* **x**  : array-like float : field x-coordinates (arcsec)
* **y**  : array-like float : field y-coordinates (arcsec)

---
Returns
    out_y : array-like float : y-component of distortion at each coord

### .recovered_dist
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L109)
```python
.recovered_dist(
   x, y
)
```

---
Evaluate recovered/estimated input distortions at arbitrary coordinates.
This is the estimated distortion via the differential calibration method
based on the input static distortion.

`x` and `y` (in arcsec) can be anywhere in the science field, but must
be array-like and the same size.


**Args**

* **x**  : array-like float : field x-coordinates (arcsec)
* **y**  : array-like float : field y-coordinates (arcsec)

---
Returns
    out_y : array-like float : y-component of distortion at each coord

### .residual_dist
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L139)
```python
.residual_dist(
   x, y
)
```

---
Evaluate residual distortions at arbitrary coordinates.

`x` and `y` (in arcsec) can be anywhere in the science field, but must
be array-like and the same size.


**Args**

* **x**  : array-like float : field x-coordinates (arcsec)
* **y**  : array-like float : field y-coordinates (arcsec)

---
Returns
    out_y : array-like float : y-component of distortion at each coord

### ._hbvpoly
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L170)
```python
._hbvpoly(
   p, a, n_poly
)
```

---
Evaluate the homogenous bi-variate polynomial defined by
coefficients in a at position p.


**Arguments**

* **p**  : np.ndarray : position to evaluate polynomial at, (M,2)
* **a**  : np.ndarray : coefficients defining polynomial, (((N+1)(N+2))//2-1,)
* **N**  : int: maximum homogenous polynomial order to go to.


**Returns**

* **out**  : np.ndarray : evaluated polynomial, scalar or (M,1)


### ._hbvpoly_grad
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L193)
```python
._hbvpoly_grad(
   p, n_poly
)
```

---
Evaluate the gradient of the homogenous bi-variate polynomial
defined by coefficients in a at position p.


**Arguments**

* **p**  : np.ndarray : position to evaluate polynomial gradient at, (2,) or (M,2)
* **n_poly**  : int: maximum homogenous polynomial order to go to.


**Returns**

* **out**  : np.ndarray : evaluated polynomial gradient,


### ._make_pinhole_grid
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L225)
```python
._make_pinhole_grid(
   xshift = 0.0, yshift = 0.0, sigma = 0.0, grid = 'square', incl_dist = True,
   pins_per_side = 30, mask_scale = 1000.0/582, pin_pitch = 0.571,
   plate_scale = 0.00736/0.01, dist_func_degmm = None, seed = 1234
)
```

---
Generate arrays of x-y pinhole positions in pixels and arcsec.
Optionally pass a global shift in x/y, in mm, to shift pinhole grid
relative to distortion field. Can also provide an uncertainty in the hole
positions (also in mm), which is treated as Gaussian. Distortions are
included by default, but can be turned off to get "nominal" pinhole grid.



**Args**

* **xshift** (float, optional) : Shift amount in x-axis (arcsec). Defaults to 0..
* **yshift** (float, optional) : Shift amount in y-axis (arcsec). Defaults to 0..
* **sigma** (float, optional) : Standard deviation on pinhole positions. Defaults to 0..
* **incl_dist** (bool, optional) : Flag to include distortions in returned coordinates. Defaults to True.
* **pins_per_side** (int, optional) : Number of pins per side of the grid. Defaults to 30.
* **mask_scale** (float, optional) : arcsec/mm at the mask. Defaults to (1e3/582).
* **pin_pitch** (float, optional) : spacing between the pinholes (in mm). Defaults to 0.582.
* **plate_scale** (float, optional) : arcsec/mm at the sensor. Defaults to 7.36e-3/10e-3.
* **dist_func_degmm** (function, optional) : Function to obtain distortion.
    Takes argument in degrees in field, returns distortion in mm at sensor. Defaults to None.


**Returns**

* **ndarray**  : 2d array of x-y pinhole positions in arcsec.


----


## AstromCalibSimAna
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L376)
```python 
AstromCalibSimAna(
   *args, centroid_noise_std = 0.0, **kwargs
)
```



----


## AstromCalibSimE2E
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L404)
```python 
AstromCalibSimE2E(
   *args, pin_size = 0.01, pinhole_os = 4, pixel_os = 2, wavelength = 5.5e-07,
   pinhole_support_width = 128, noise_fun = None, centroid_win_rad = 0.2,
   centroid_threshold = 0.0, **kwargs
)
```




**Methods:**


### ._pinhole
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L519)
```python
._pinhole(
   size, x, y, radius
)
```

---
generate the sampled pinhole function


**Args**

* **size** (int) : size of the output array (in pixels)
* **x** (float) : x position of the pinhole (in pixels)
* **y** (float) : y position of the pinhole (in pixels)
* **radius** (float) : radius of the pinhole (in pixels)


**Returns**

* **ndarray**  : 2d array of the pinhole function

