#


## Source
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/source.py/#L45)
```python 
Source(
   input_par, exp_time, static_dist = False, stat_amp = 1.0, tt_amp = 1.0,
   use_cov = False
)
```


---
The source object allows access to the parameters of each source (star) required
to compute the MAVISIM image.


**Args**

* **input_par**  : imported input parameter python file, e.g., `input_parameters.py`.
* **exp_time** (`float`) : exposure time in seconds to simulate.
* **static_dist** (`bool`, optional) : if `True`, add static distortion defined in `input_par`.
* **stat_amp** (`float`, optional) : scaling factor to apply to distortions.
* **tt_amp** (`float`, optional) : scaling factor to apply to tip-tilt blurring kernel width.


**Attributes**

* **exp_time** (`float`) : exposure time in seconds to simulate.
* **star** (`np.ndarray` of `int`) : unique ID of each star.
* **flux** (`np.ndarray` of `float`) : flux of each star.
* **ra** (`np.ndarray` of `float`) : RA of each star.
* **dec** (`np.ndarray` of `float`) : Dec of each star.
* **x_pos** (`np.ndarray` of `float`) : X position of each star (in arcsec).
* **x_pm** (`np.ndarray` of `float`) : X proper motion of each star (in pixels).
* **x_dist** (`np.ndarray` of `float`) : sub-pixel X position shift of each star (in pixels).
* **y_pos** (`np.ndarray` of `float`) : Y position of each star (in arcsec).
* **y_pm** (`np.ndarray` of `float`) : Y proper motion of each star (in mas/year).
* **y_dist** (`np.ndarray` of `float`) : sub-pixel Y position shift of each star (in pixels).
* **gauss_pos** (`np.ndarray` of `float`) : X/Y-position of each star (in arcsec).
* **gauss_cov** (`np.ndarray` of `float`) : covariance of Gaussian kernel to simulate tip-tilt blurring (in arcsec^2).
* **static_dist** (`np.ndarray` of `float`) : static distortion to apply to each source.



**Methods:**


### .build_source
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/source.py/#L91)
```python
.build_source()
```

---
From the data stored in the object, compute the source data as required


### .decimate
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/source.py/#L201)
```python
.decimate(
   nstar
)
```

---
Decimate the list of objects in the object (e.g., for faster simualtions).


**Args**

* **nstar** (`int`) : number of stars to reduce list down to. The resulting object will
only contain the first `nstar` stars.
