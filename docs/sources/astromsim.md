#


## AstromCalibSim
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L4)
```python 
AstromCalibSim(
   static_distort, centroid_noise_std = 1e-05, hole_position_std = 0.01, dx = 0.2,
   dy = 0.2, n_poly = 6, pin_pitch = 0.5, num_pin_x = 30, mask_scale = 1000.0/582,
   pixel_size_as = 0.00736, pixel_size_mm = 0.01
)
```


---
Astrometric Calibration Simulator for MAVIS.

MAVIS will feature an astrometric calibration mask characterise the static
distortions present in the MAVIS system. This simulator allows one to quantify
the performance of that mask for a given distortion field, using the nominal
astrometric calibration identification algorithms intended for MAVIS.

The constructor takes in an astropy-parsed distortion field with the appropriate 
headers (see below), and performs a simulated astrometric calibration process. 
The public methods allow the evaluation of input, recovered, and residual 
distortion fields.

Static distortion file should have at least the following columns:
```    
  Field_x(deg)    Field_y(deg)    Predicted_x(mm)  Predicted_y(mm)    Real_x(mm)      Real_y(mm)
-4.16666667E-03 -4.16666667E-03   2.02613981E+01   2.02626354E+01   2.03494040E+01  2.03513749E+01
-4.16666667E-03 -4.07407407E-03   2.02613981E+01   1.98123546E+01   2.03496805E+01  1.98994423E+01
-4.16666667E-03 -3.98148148E-03   2.02613981E+01   1.93620738E+01   2.03499497E+01  1.94474891E+01
...
 4.07407407E-03  4.07407407E-03  -1.98111448E+01  -1.98123546E+01  -1.98997456E+01 -1.98980094E+01
```
and it should be parsed by astropy first, like:
```python
from astropy import ascii
static_distort = ascii(distort_file)
```


**Args**

* **static_distort**  : astropy table : table containing the distortions across the field
* **centroid_noise_std**  : float : standard deviation of Gaussian noise applied to centroids (arcsec).
* **hole_position_std**  : float : standard deviation of manufacturing error on hole positions (mm).
* **dx**  : float : shift applied in x direction for calibration process (mm).
* **dy**  : float : shift applied in y direction for calibration process (mm).
* **n_poly**  : int : maximum order of homogenous bivariate polynomial used to fit distortions.
* **pin_pitch**  : float : distance between pins at input plane (mm).
* **num_pin_x**  : int : number of pins across the x-dimension (sqrt(num_pins_total)).
* **mask_scale**  : float : arcsec/mm at input plane.
* **pixel_size_as**  : float : pixel size in arcsec of imager camera.
* **pixel_size_mm**  : float : physical pixel size in mm of imager camera.












**Methods:**


### .input_dist
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L76)
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
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L102)
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
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/astromsim.py/#L130)
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
