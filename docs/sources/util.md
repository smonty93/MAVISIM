#


### input_coo
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/util.py/#L19)
```python
.input_coo(
   input_par, source
)
```

---
Take input parameters and parse them as an astropy table for comparison with DAOphot.


**Args**

* **input_par**  : input parameters, e.g., from `input_parameters.py`.
* **source** (`Source` object) : the `Source` object containing all of the sources simulated.


**Returns**

trimmed_cat = an astropy table containing the TRUE input positions (static distortion,
sub-pixel positioning and proper motion taken into account) and additional information necessary to compare with the DAOPhot output.

----


### add_all_noise
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/util.py/#L71)
```python
.add_all_noise(
   input_par, image, exp_time
)
```

---
Adds all noise effects onto noiseless image.


**Args**

* **input_par**  : input parameters object.
* **image** (`np.ndarray`) : 2D noise-free image.
* **exp_time** (`float`) : total exposure time in seconds.


**Returns**

* **image_adu** (`np.ndarray`) : final image + noise in ADU


----


### make_static_dist_map
[source](https://github.com/smonty93/mavisim/blob/v1.1dev/mavisim/util.py/#L154)
```python
.make_static_dist_map(
   input_par
)
```

---
Make static distortion map from loaded distortion samples in `input_par`.


**Args**

* **input_par**  : input parameters object.


**Returns**

dist_x_func_degmm:
dist_y_func_degmm:
