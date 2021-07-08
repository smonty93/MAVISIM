# Installation

MAVISIM is not currently uploaded to the Python Package Index and thus must be installed from source. The following should be sufficient in most cases to install MAVISIM 1.0 (currently the main branch):
```bash
git clone https://github.com/smonty93/MAVISIM.git
cd MAVISIM
pip install .
```

## Data Dependencies

To run MAVISIM 1.0 you'll need to download the `data` directory from <a href="http://www.mso.anu.edu.au/~montys/MAVISIM1/" target="_blank">this link.</a>. It contains the database of field variable Fourier PSFs at 550nm and other useful data including an example.

A helper script has been included to fetch this data. Assuming `wget` is available, the following should download the required data to `MAVISIM/data`:
```bash
./download_data.sh
```

## External Dependencies

MAVISIM 1.0 was written for **Python 3** and relies on the following packages with the listed version being the tested version.

| Package     | Tested Version |
| -----------  | ----------- |
| Numpy      |   1.16.5    |
| Astropy      | 2.0.9    |
| Scipy | 1.2.1  |
