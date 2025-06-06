# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 13:30:11 2022

@author: danielgodinez
"""
from setuptools import setup, find_packages

setup(
    name="protoRT",
    version="1.0.2",  
    author="Daniel Godines",
    author_email="danielgodinez123@gmail.com",
    description="Radiative transfer-based mass analysis of planetesimal formation simulations.",
    license='GPL-3.0-or-later',
    url="https://github.com/Professor-G/protoRT",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Operating System :: OS Independent',
    ],
    packages=find_packages(),
    python_requires='>=3.9',
    install_requires=[
        "numpy",
        "scipy>=1.7.3",
        "matplotlib>=3.5.1",
        "astropy>=5.0.4"
    ],
    include_package_data=True,
    package_data={
        "protoRT": [
            "data/*.npy",
            "data/*.txt",
            "data/axis",
        ],
    },
)
