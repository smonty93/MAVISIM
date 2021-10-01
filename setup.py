#!/usr/bin/env python

from distutils.core import setup
import os

thelibFolder = os.path.dirname(os.path.realpath(__file__))
requirementPath = thelibFolder + '/requirements.txt'
install_requires = []
if os.path.isfile(requirementPath):
    with open(requirementPath) as f:
        install_requires = f.read().splitlines()

setup(name='mavisim-gpu',
      version='1.1',
      description='MAVIS Image Simulator - with CUDA GPU acceleration',
      author='Stephanie Monty',
      author_email='Stephanie.Monty@anu.edu.au',
      url='https://www.github.com/smonty93/MAVISIM',
      install_requires=install_requires,
      package_dir={'mavisimgpu':'src/mavisim'},
      packages=['mavisimgpu'],
     )
