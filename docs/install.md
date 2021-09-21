# Installation

## Current Version (MAVISIM 1.1)

MAVISIM is not currently uploaded to the Python Package Index and thus must be installed from source. The current version of MAVISIM can be installed via pip using the following command:

```bash
pip install "git+https://github.com/smonty93/MAVISIM.git"
```

## Older Version (MAVISIM 1.0)
The following should be sufficient in most cases to install MAVISIM 1.0:

```bash
pip install "git+https://github.com/smonty93/MAVISIM.git@v1.0"
```

[]: # Older code to install
[]: #```bash
[]: #git clone https://github.com/smonty93/MAVISIM.git
[]: #cd MAVISIM
[]: #pip install .
[]: #```

## Data Dependencies

## 

To run MAVISIM 1.0 you'll need to download the `data` directory from <a href="http://www.mso.anu.edu.au/~montys/MAVISIM1/" target="_blank">this link.</a>. It contains the database of field variable Fourier PSFs at 550nm and other useful data including an example.

A helper script has been included to fetch this data. Assuming `wget` is available, the following should download the required data to `MAVISIM/data`:
```bash
./download_data.sh
```

### External Dependencies

MAVISIM was written for **Python 3** and relies on the following packages with the listed version being the tested version.

| Package     | Tested Version |
| -----------  | ----------- |
| Numpy      |   1.16.5    |
| Astropy      | 2.0.9    |
| Scipy | 1.2.1  |
