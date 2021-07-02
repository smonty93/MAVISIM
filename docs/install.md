# Installation

MAVISIM is not currently uploaded to the Python Package Index and thus is not yet pip installable. To install MAVISIM 1.0 please download the main branch. Place the `mavisim` directory whereever you'd like and link to it using the `path_to_mavisim` parameter in the `input_parameters.py` file.

## Data Dependencies

To run MAVISIM 1.0 you'll need to download the `data` directory from <a href="http://www.mso.anu.edu.au/~montys/MAVISIM1/" target="_blank">this link.</a>. It contains the database of field variable Fourier PSFs at 550nm and other useful data including an example.

## External Dependencies

MAVISIM 1.0 was written for **Python 3** and relies on the following packages with the listed version being the tested version.

| Package     | Tested Version |
| -----------  | ----------- |
| Numpy      |   1.16.5    |
| Astropy      | 2.0.9    |
| Scipy | 1.2.1  |
