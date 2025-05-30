# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 13:30:11 2022

@author: danielgodinez
"""
from setuptools import setup, find_packages, Extension


setup(
    name="protoRT",
    version="0.9",
    author="Daniel Godines",
    author_email="danielgodinez123@gmail.com",
    description="Radiative transfer-based mass analysis of planetesimal formation simulations.",
    license='GPL-3.0',
    url = "https://github.com/Professor-G/protoRT",
    classifiers=[
		'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Programming Language :: Python :: 3', 
        'Operating System :: OS Independent',
],
    packages=find_packages('.'),
    install_requires = ['numpy','scipy','matplotlib', 'astropy'],
    python_requires='>=3.7',
    include_package_data=True,
    test_suite="nose.collector",
    package_data={
    '': ['all_grain_sizes.txt', 'all_wavelengths.txt', 'axis', 'density_cube_1.npy', 'density_cube_2.npy', 'kappa_abs.npy', 'kappa_sca.npy'],
},

)
