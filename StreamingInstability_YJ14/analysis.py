#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:01:11 2021

@author: daniel
"""
import warnings
import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as const
warnings.filterwarnings("ignore")
from pathlib import Path
import pkg_resources

from radiative_transport import calculate_tau

def load_cube():
    """
    Loads 256 x 256 x 256 density cube. Corresponds to one snapshot of a 100
    period simulation of the streaming instability in a protoplanetary disk,
    provided by Yang & Johansen (2014).

    Returns:
        The first output is the 3D array of the cube, the second output
        is the axis array. This array can be used as either axis (z,y,x).
    """
    resource_package = __name__
    resource_path = '/'.join(('data', 'density_cube_1.npy'))
    file = pkg_resources.resource_filename(resource_package, resource_path)
    density_cube_1 = np.load(file)

    resource_path = '/'.join(('data', 'density_cube_2.npy'))
    file = pkg_resources.resource_filename(resource_package, resource_path)
    density_cube_2 = np.load(file)

    resource_path = '/'.join(('data', 'axis'))
    file = pkg_resources.resource_filename(resource_package, resource_path)
    axis = np.loadtxt(file)

    return np.r_[density_cube_1, density_cube_2], axis


h,k,c,sb,au= 6.626e-27, 1.3807e-16, 2.997e10, 5.67e-5,  const.au.cgs.value
T = 30
kappa_0 = 2e-4
kappa = kappa_0 * T**2.1 #Opacity is function of T at low temperatures (See Table 2: https://iopscience.iop.org/article/10.1086/304514/pdf) & Figure 1 (https://arxiv.org/pdf/astro-ph/0308344.pdf)


rhop, axis = load_cube()

# Calculate optical depth
n = len(z) 
nu = np.logspace(11,15,n)

tau = calculate_tau(rhop, axis=axis, kappa=kappa, sigma=1000)

plt.contourf(axis, axis, np.log10(tau), np.linspace(-2,2,256))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Optical Depth')
plt.show()



#Calcualte outgoing flux
nx, ny = len(axis), len(axis)
dx, dy, dz = np.diff(axis)[0], np.diff(axis)[0], np.diff(axis)[0]
Nw = 1./(nx*ny)

box_mass_codeunits = np.sum(rhop)* dx * dy * dz 

H = 5*au 
Lx=np.abs(axis[0]-axis[-1])*H #Length of the box in cm
Ly=np.abs(axis[0]-axis[-1])*H #Length of the box in cm
area = Lx*Ly

src_fn = sb*T**4     #Estimate flux according to optical thickness

n = 10
sigma = np.linspace(1, 400, n)
filling_factor = np.zeros(n)
mass_excess = np.zeros(n)
for i in range(n):
    tau = calculate_tau(density=rhop, axis=axis, kappa=kappa, sigma=sigma[i])
    filling_factor[i] = len(np.where(tau > 1)[0]) * Nw
    #flux = calculate_flux(density=f.rhop, tau=tau, axis=z, T=T, kappa=kappa, sigma=sigma[i])

    #emissivity = srf_cn * kappa
    flux_approx = src_fn * (1-np.exp(-tau)) #If source fn is constant and region is optically thick (Eq. 5.120), assumes source function is constant from RT module notes)
    
    sigma_dust = np.mean(flux_approx) / (src_fn*kappa) #Sigma dust observer sees if optically thin assumption
    
    unit_mass = sigma[i] * H**2
    box_mass = box_mass_codeunits * unit_mass #Total mass in box (actual)
    observed_mass = sigma_dust*area  #Total mass of the dust in the box (observed)
    
    #tot_mass = sigma*dust_to_gas*area #Total mass in box (actual)
    mass_excess[i] = box_mass / observed_mass

    #arr = np.array([sigma_dust, mass, tot_mass, filling_factor, mass_excess])
    #np.savetxt(path+'s0_'+str(i)+'.txt', arr)


plt.plot(sigma, mass_excess, 'ro-')
plt.xlabel(r'$\Sigma_d \ (g / cm^2)$', size=20)
plt.ylabel('Mass Excess', size=20)
plt.show()


plt.plot(filling_factor, mass_excess, 'ro-')
plt.xlabel('Filling Factor', size=20)
plt.ylabel('Mass Excess', size=20)
plt.show()


plt.plot(sigma, filling_factor, 'ro-')
plt.xlabel(r'$\Sigma_d \ (g / cm^2)$', size=20)
plt.ylabel('Filling Factor', size=20)
plt.show()

print("Max flux :"+str(flux.max()))
print("Min flux :"+str(flux.min()))
print('Filling factor :'+str(len(np.where(tau > 1)[0]) / (tau.shape[0]*tau.shape[1])))
plt.contourf(x, y, flux, 256)
plt.title('Flux')
plt.colorbar()
plt.show()

