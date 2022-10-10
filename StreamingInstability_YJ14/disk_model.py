#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 09:25:13 2022

@author: daniel
"""
import numpy as  np 
import matplotlib.pyplot as plt  
import astropy.constants as const  


def calc_sigma_g(r, r_c=300, M_disk=0.2):
    """
    Calculates the gas surface density, as per the protoplanetary 
    disk model presented in Section 2 of Drazkowska et al (2022)
    See: https://arxiv.org/pdf/2101.01728.pdf

    Args:
        r (float, ndarray): The radii of the disk at which to calculate
            the gas surface density. Must be in cgs units.
        r_c (float): Characteristic radius of the disk, in AU. This radius
            for compact disks is about 30 AU, and 300 AU for large disks.
            Defaults to 30 AU.  
        M_disk (float): Mass of the disk, in terms of solar mass. 
            Defaults to 0.2.
        
    Returns:
        Gas surface density at the specified r.
    """

    M_disk = M_disk*const.M_sun.cgs.value
    r_c = r_c*const.au.cgs.value

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

def calc_stokes(sigma_g, grain_size=0.1, grain_rho=1):
    """
    Calculates the Stokes number according to the gas
    column density and the grain properties
    
    Args:
        grain_size (float): Grain size, in cgs units. Defaults to 
            1 mm.
        grain_rho (float): Internal grain density in cgs units. Defaults to 1 g/cm2
        sigma_g (float): Gas column density, in cgs units.
        
    Returns:
        Stoke's number.
    """

    return np.pi * grain_size * grain_rho / 2. / sigma_g

def omega(M, r):
    """
    Keplerian frequency.

    Args:
        M (float): Mass of the central star, in cgs units.
        r (float): Distance from the central star, in cgs units.
    """

    return np.sqrt(const.G.cgs.value*M/r**3)

def toomreq(cs, M, r, sigma):
    """
    Toomre Q parameter (dimensionless)

    Args:
        cs (float): Sound speed, in cgs units.
        M (float): Mass of the central star, in cgs units.
        r (float): Distance from the central star, in cgs units.
        sigma (float): Gas column density, in cgs units.
        
    Returns:
        ToomreQ
    """

    return cs*omega(M,r) / np.pi / const.G.cgs.value / sigma

def T(r, T0=150, q=3./7):
    """
    Temperature model from Ida et al. 2016: https://arxiv.org/pdf/1604.01291.pdf

    Args:
        r (float): Distance from the central star, in cgs units.
        T0 (float): Temperature at r0, defaults to 150 K, which is
            the accepted standard.
        q (float): Power law index, defaults to 3/7, which is the 
            accepted standard.
    """

    return T0*(r / const.au.cgs.value)**(-q)

def aspect_ratio(r, M, mmol=2.3):
    """
    Aspect ratio, as a function of T 
    Equation 5: https://www.aanda.org/articles/aa/pdf/2015/03/aa24964-14.pdf

     Args:
        r (float): Distance from the central star, in cgs units.
        M (float): Mass of the central star, in cgs units.
        mmol (float): Mean molecular weight, defaults to 2.3, corresponding
            to 5 parts molecular H and 1 part He

    Returns:
        Aspect ratio.
    """

    return np.sqrt(T(r) * r * const.R.cgs.value / const.G.cgs.value / M / mmol )

def h(H, r):
    """
    Aspect ratio.

    Args:
        H (float): Scale height, in cgs units.
        r (float): Distance from the central star, in cgs units.

    Returns:
        Aspect ratio.
    """
    
    return H/r

def sound_speed(r, gamma=1.4, mmol=2.3):
    """
    Args:
        r (float): Distance from the central star, in parsecs
        gamma (float): Adiabatic index, defaults to 1.4.
        mmol (float): Mean molecular weight, defaults to 2.3, corresponding
            to 5 parts molecular H and 1 part He
    """
    
    Rgasmu = const.R.cgs.value / mmol
    cp = gamma*Rgasmu / (gamma-1)

    return np.sqrt(T(r)*cp*(gamma-1))

def G_tilde(Q):
    """
    G tilde self-gravity parameter.

    Args:
        Q (float): Toomre Q parameter

    """

    return np.sqrt(8/np.pi) / Q

def parameters(r, M, M_disk):
    """
    Calculates the following parameters
    according to the fiducial models

    Args:
        r (float): Distance from the central star, in cgs units.
        M (float): Mass of the central star, in cgs units.
        M_disk (float): Mass of the disk, in terms of solar mass. 
            Defaults to 0.2.

    Returns:    
        Dictionary containing the following model attributes: 
        sigma_g, cs, H, h, Q, G, beta, stoke

    """

    sigma = calc_sigma_g(r=r, M_disk=M_disk) 
    print("sigma_g: {}".format(sigma))

    cs = sound_speed(r)
    print("cs: {}".format(cs))

    temp = T(r)
    print("Temperature : {}".format(temp))
    
    H = cs/omega(M=M, r=r)
    print("Scale Height: {}".format(H))

    h1 = h(H, r=r)
    print("Aspect ratio 1: {}".format(h1))

    h2 = aspect_ratio(r=r, M=M)
    print("Aspect ratio 2: {}".format(h2))

    Q = toomreq(cs, M=M, r=r, sigma=sigma)
    print("Toomre Q: {}".format(Q))

    G = G_tilde(Q)
    print("G Tilde: {}".format(G))

    beta = -h1 * -2.28
    print("Beta: {}".format(beta))

    st = calc_stokes(sigma_g=sigma, grain_size=0.1)
    print("Stokes number: {}".format(st))

    return {'sigma_g': sigma, 'cs': cs, 'T':temp, 'H': H, 'h':h1, 'Q':Q, 'G_tilde':G, 'beta':beta, 'stoke':st}


######################################
r, M, M_disk = 30*const.au.cgs.value, const.M_sun.cgs.value, 0.2
print(parameters(r,M,M_disk))

"""
#Disk model plot

import numpy as  np 
import matplotlib.pyplot as plt  
import astropy.constants as const

M_disk = 0.2
r = np.arange(30,101,1)
compact, large = disk_model.calc_sigma_g(r, M_disk=M_disk), disk_model.calc_sigma_g(r,M_disk=M_disk,r_c=300)
temp = T(r*const.au.cgs.value)
  
plt.plot(r, compact, label='Compact Disk')


fig, axes = plt.subplots(nrows=2)
ax1, ax2 = axes
ax1.plot(r, large, c='k')
ax1.set_ylabel(r'$\Sigma_g \ [g \ cm^{-2}]$', fontsize=14)
ax1.set_ylim(0,31)
ax1.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
ax1.tick_params(labelsize=14, axis="both", which="both")
ax1.set_title('Disk Model',size=16)
ax1.set_xticklabels([])


ax2.plot(r, temp, c='k')
ax2.set_xlabel('Radius (AU)', size=14)
ax2.set_ylabel('T [K]', size=14)
ax2.tick_params(labelsize=14, axis="both", which="both")
ax2.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
ax2.set_ylim(20, 36)

plt.savefig('Disk_Model.png', dpi=300, bbox_inches='tight')



plt.plot(r, large, label='Large Disk')
plt.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
plt.legend(prop={'size':16})
plt.tick_params(labelsize=14, axis="both", which="both")
plt.xlabel('Radius (AU)', size=18)
plt.ylabel(r'$\Sigma_g \ (g \ cm^{-2})$', fontsize=18)
plt.title('Disk Model',size=20)
#plt.xscale('log')
#plt.yscale('log')
plt.savefig('Disk_Models', dpi=300, bbox_inches='tight')

plt.show()


#Stokes number at different sigma_g and grain size

radii = np.array([10, 30, 100])
grain_size = np.array([1, 0.3, 0.1, 0.03])
sigma_g = disk_model.calc_sigma_g(radii, M_disk=M_disk, r_c=300)
stoke = calc_stokes(grain_size=grain_size, grain_rho=1, sigma_g=sigma_g)
"""



