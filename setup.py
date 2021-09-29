#!/usr/bin/env python

from distutils.core import setup
import os

thelibFolder = os.path.dirname(os.path.realpath(__file__))
requirementPath = thelibFolder + '/requirements.txt'
install_requires = []
if os.path.isfile(requirementPath):
    with open(requirementPath) as f:
        install_requires = f.read().splitlines()

setup(name='MAVISIM',
      version='1.1',
      description='MAVIS Image Simulator',
      author='Stephanie Monty',
      author_email='Stephanie.Monty@anu.edu.au',
      url='https://www.github.com/smonty93/MAVISIM',
      install_requires=install_requires,
      package_dir={'mavisim':'src/mavisim1_1'},
      packages=['mavisim'],
     )
