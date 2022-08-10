# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 13:30:11 2022

@author: danielgodinez
"""
from setuptools import setup, find_packages, Extension


setup(
    name="StreamingInstability_YJ14",
    version="0.1",
    author="Daniel Godines",
    author_email="danielgodinez123@gmail.com",
    description="Analysis of the planetesimal formation simulations by the streaming instability presented by Yang & Johansen (2014)",
    license='GPL-3.0',
    url = "https://github.com/Professor-G/StreamingInstability_YJ14",
    classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Developers',
		'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Programming Language :: Python :: 3',	   
],
    packages=find_packages('.'),
    install_requires = ['numpy','scipy','matplotlib', 'astropy'],
    python_requires='>=3.7,<4',
    include_package_data=True,
    test_suite="nose.collector",
    package_data={
    '': ['density_cube_1.npy', 'density_cube_2.npy', 'axis'],
},

)
