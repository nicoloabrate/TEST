"""
Author: N. Abrate.

File: setup.py

Description: Installation file.
"""
from setuptools import setup, find_packages

requirements = "requirements.txt"

setup(
   name='TEST',
   version='1.0.0',
   author='N. Abrate',
   author_email='nicolo.abrate@polito.it',
   url='https://github.com/nicoloabrate/TEST.git',
   # package_name = ['TEST'],
   packages=find_packages(),
   license='LICENSE.md',
   description='TEST: solving neutron transport equation in 1D plain geometry',
   long_description=open('README.md').read(),
   test_suite="tests",
   setup_requires=['pytest-runner'],
   tests_require=['pytest'],
   include_package_data=True,
   classifiers=[
        "Programming Language :: Python :: 3",
        "License :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=open(requirements).read().splitlines(),
    python_requires='>=3.6',
)
