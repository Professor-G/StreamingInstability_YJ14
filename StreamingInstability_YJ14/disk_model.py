#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 09:25:13 2022

@author: daniel
"""
import numpy as  np 
import matplotlib.pyplot as plt  
import astropy.constants as const  


def calc_sigma_g(r, r_c=30, M_disk=0.2):
    """
    Calculates the gas surface density, as per the protoplanetary 
    disk model presented in Section 2 of Drazkowska et al (2022)
    See: https://arxiv.org/pdf/2101.01728.pdf

    Args:
        r (float, ndarray): The radii of the disk at which to calculate
            the gas surface density. Must be in AU.
        r_c (float): Characteristic radius of the disk, in AU. This radius
            for compact disks is about 30 AU, and 300 AU for large disks.
            Defaults to 30.  
        M_disk (float): Mass of the disk, in terms of solar mass. 
            Defaults to 0.2.
        

    Returns:
        Gas surface density.
    """

    M_disk = M_disk*const.M_sun.cgs.value
    r_c = r_c*const.au.cgs.value
    r = r*const.au.cgs.value

    sigma_g = ((M_disk/(2.*np.pi*r_c**2))*(r/r_c)**-1)*np.e**(-(r/r_c))

    return sigma_g

def calc_sigma_d(r, r_c=30, M_disk=0.2, Z=0.01):
    """
    Calculates the dust surface density, as per the protoplanetary 
    disk model presented in Section 2 of Drazkowska et al (2022)
    See: https://arxiv.org/pdf/2101.01728.pdf

    Args:
        r (float, ndarray): The radii of the disk at which to calculate
            the gas surface density. Must be in AU. 
        r_c (float): Characteristic radius of the disk, in AU. This radius
            for compact disks is about 30 AU, and 300 AU for large disks.
            Defaults to 30.
        M_disk (float): Mass of the disk, in terms of solar mass. 
            Defaults to 0.2.  
        Z (float):  The global solids-to-gas ratio. Defaults to 0.01.

    Returns:
        Gas surface density.
    """

    sigma_g = calc_sigma_g(r, r_c, M_disk)
    sigma_d = sigma_g * Z

    return sigma_d 


def calc_stokes(grain_size, grain_rho, sigma_g):
    """
    Calculates the Stokes number according to the gas
    column density and the grain properties
    """

    return np.pi * grain_size * grain_rho / 2. / sigma_g


"""
#Disk model plot

import numpy as  np 
import matplotlib.pyplot as plt  
import astropy.constants as const

M_disk = 0.01 
r = np.arange(10,101,1)
compact, large = disk_model.calc_sigma_g(r, M_disk=M_disk), disk_model.calc_sigma_g(r,M_disk=M_disk,r_c=300)

  
plt.plot(r, compact, label='Compact Disk')
plt.plot(r, large, label='Large Disk')
plt.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
plt.legend(prop={'size':16})
plt.tick_params(labelsize=14, axis="both", which="both")
plt.xlabel('Radius (AU)', size=16)
plt.ylabel(r'$\Sigma_g \ (g/cm^3)$',size=16)
plt.title('Gas Column Density Profiles',size=18)
#plt.xscale('log')
#plt.yscale('log')
plt.show()


#Stokes number at different sigma_g and grain size

radii = np.array([10, 30, 100])
grain_size = np.array([1, 0.3, 0.1, 0.03])
sigma_g = disk_model.calc_sigma_g(radii, M_disk=M_disk, r_c=300)
stoke = calc_stokes(grain_size=grain_size, grain_rho=1, sigma_g=sigma_g)
"""



