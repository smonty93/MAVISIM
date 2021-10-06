# Installation

## Standard
For most users, the following should be sufficient to install MAVISIM and download the data dependencies:
```bash
pip install mavisim
git clone https://github.com/smonty93/mavisim
cd mavisim
./download_data.sh
```

In the root directory of the cloned git repository, an example Jupyter Notebook `mavisim.ipynb` can be found, which by default simulates the system defined in `input_parameters.py`. To run this notebook:
```bash
jupyter notebook mavisim.ipynb
```

## Explicit Version From PyPI

The latest stable version of MAVISIM can be directly installed via the Python Package Index:
```bash
pip install mavisim
```
For installing specific (e.g., legacy) versions, one can install from a previous `pip` release, e.g., to install `v1.0.x`:
```bash
pip install "mavisim>=1.0,<1.1" 
```

## Explicity Version From Github
For beta-releases/dev branches, it is possible to install MAVISIM directly from a github:
```bash
pip install "git+https://github.com/smonty93/MAVISIM@[branch or tag]"
```
e.g., for version `v1.0.2`:
```bash
pip install "git+https://github.com/smonty93/MAVISIM@v1.0.2"
```