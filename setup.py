"""
Author: N. Abrate.

File: setup.py

Description: Installation file.
"""
from setuptools import setup, find_packages

setup(
   name='TEST',
   version='0.1.0',
   author='N. Abrate',
   author_email='nicolo.abrate@polito.it',
   url='https://nicolo_abrate@bitbucket.org/nicolo_abrate/test.git',
   package_name = ['TEST'],
   packages=find_packages(),
   license='LICENSE.md',
   description='A package for solving neutron transport equation in 1D plain geometry',
   long_description=open('README.md').read(),
   test_suite="tests",
   setup_requires=['pytest-runner'],
   tests_require=['pytest'],
   classifiers=[
        "Programming Language :: Python :: 3",
        "License :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
