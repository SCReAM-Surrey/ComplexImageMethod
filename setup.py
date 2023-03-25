#!/usr/bin/env python3

from setuptools import setup


setup(name='ComplexImageMethod',
      version='1.0',
      author='Orchisama Das',
      license='afl-3.0',
      package_dir={'ComplexImageMethod':'cim'},
      py_modules = [],
      scripts=[],
      description='Room acoustics simulation with ComplexImageMethod',
   	  long_description=open('README.md').read(),
      install_requires=[
       "numpy",
       "scipy",
       "mpmath",
       "matplotlib",
       "pyroomacoustics"
   ],
)