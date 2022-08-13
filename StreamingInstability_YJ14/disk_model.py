#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 09:25:13 2022

@author: daniel
"""
import numpy as  np 
import matplotlib.pyplot as plt  
import astropy.constants as const  


def calc_sigma_g(r, M_disk=0.2, r_c=30):
    """
    Calculates the gas surface density, as per the protoplanetary 
    disk model presented in Section 2 of Drazkowska et al (2022)
    See: https://arxiv.org/pdf/2101.01728.pdf

    Args:
        r (float, ndarray): The radius at which to 
        M_disk (float): Mass of the disk, in terms of solar mass. 
            Defaults to 0.2.
        r_c (float): Characteristic radius of the disk. This radius
            for compact disks is about 30 AU, and 300 AU for large disks.
            Defaults to 30.  

    Returns:
        Gas surface density.
    """

    M_disk *= const.M_sun.cgs.value
    r_c *= const.au.cgs.value

    sigma_g = ((M_disk/(2.*np.pi*r_c**2))*(r/r_c)**-1)*np.e**(-(r/r_c))

    return sigma_g

def calc_sigma_d(r, M_disk=0.2, r_c=30, Z=0.01):
    """
    Calculates the dust surface density, as per the protoplanetary 
    disk model presented in Section 2 of Drazkowska et al (2022)
    See: https://arxiv.org/pdf/2101.01728.pdf

    Args:
        r (float, ndarray): The radius at which to 
        M_disk (float): Mass of the disk, in terms of solar mass. 
            Defaults to 0.2.
        r_c (float): Characteristic radius of the disk. This radius
            for compact disks is about 30 AU, and 300 AU for large disks.
            Defaults to 30.  
        Z (float):  The global solids-to-gas ratio. Defaults to 0.01.

    Returns:
        Gas surface density.
    """

    sigma_g = calc_sigma_g(r=r, M_disk=M_disk, r_c=r_c)
    sigma_d = sigma_g * Z

    return sigma_d 


