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

rhop = np.load("/Users/daniel/Desktop/density_cube.npy")

h,k,c,sb,au= 6.626e-27, 1.3807e-16, 2.997e10, 5.67e-5,  const.au.cgs.value

#r = 50 * au #position of shearing box
#h = 0.1 #scale ratio
#H = h * r #scale height (5 AU)

#x, y, z = f.x, f.y, f.z
#xp, yp, zp = fp.xp, fp.yp, fp.zp

#dust_to_gas = 1e-2
#sigma0 = 1000 #column density at 1 AU
#sigma = sigma0 * (r / au) ** (-1.)
#rho = sigma / H / np.sqrt(2*np.pi) #multiply to density to make unitless

T = 30
kappa_0 = 2e-4
kappa = kappa_0 * T**2.1 #Opacity is function of T at low temperatures (See Table 2: https://iopscience.iop.org/article/10.1086/304514/pdf) & Figure 1 (https://arxiv.org/pdf/astro-ph/0308344.pdf)


def blackbody(T, nu):
    """
    Planck's law, which describes the black body radiation 
    of a source in thermal equilibrium at a given temperature T.

    Args:
        T (float): Temperature of the source.
        nu (ndarray): Array of frequencies.

    Returns:
        Spectral radiance of the source. 
    """
    bb = 2*h*nu**3 / (c**2*(np.exp(h*nu / (k*T)) - 1))

    return bb

def calculate_tau(density, axis, kappa, sigma):
    """
    Integrates density cube along the first axis
    of the array. This function assumes a 3D cube with 
    a uniform grid.
    
    Args:
        density (ndarry): 3D array of densities.
        axis (ndarray): 1D array of the axis along which to integrate.
        kappa (float): Dust opacity coefficient.
        sigma (float): Column density of the dust
        
    Returns:
        2D array containing the integrated values along the third axis.
    """
    
    dz = np.diff(axis)[0]
    nx, ny = len(axis), len(axis)
    
    tau = np.zeros([ny, nx]) #python reverses order
           
    for i in range(nx):
        for j in range(ny):
            surface_density = np.trapz(density[:,j,i]) * dz * sigma
            tau[j, i] = surface_density * kappa
        
    return tau

def calculate_t(density, axis, kappa, sigma):
    """
    Optical depth with respect to z position along the column.
    Integrates from z to L. This is the optical depth as we move
    up the column, where as the optical_depth() function calculates
    along the entire column. 
    
    This function assumes a 3D cube with uniform grid.
    
    Args:
        density (ndarry): 3D array of densities.
        axis (ndarray): 1D array of the axis along which to integrate.
        kappa (float): Dust opacity coefficient 
        
    Returns:
        1D array containing the integrated values along the third axis.
    
    """
    
    dz = np.diff(axis)[0]
    nx, ny, nz = dim.nx, dim.ny, dim.nz
    
    t = np.zeros(nz)
    
    surface_density = 0
    for k in range(nz):        
        surface_density += density[k] * dz * sigma
        t[k] = surface_density * kappa
        
    return t   

def calculate_flux(density, tau, axis, T, kappa, sigma):
    """
    Calculate flux using sln for RT eqn (5.113)
    
    Args:
        density (ndarry): 3D array of densities.
        tau (ndarray): 2D array of optical depth.
        axis (ndarray): 1D array of the axis along which to integrate.
        T (float): Temperature.
        kappa (float): Dust opacity coefficient 
        
    Returns:
        2D array containing the integrated values along the third axis.
    """
    
    dx = np.diff(axis)[0]
    nx, ny, nz = dim.nx, dim.ny, dim.nz
    flux = np.zeros([ny, nx])
    
    for i in range(nx):
        print(i)
        for j in range(ny):
            bb = np.zeros(nz)
            rhod = density[:,j,i] #1D array
            t = calculate_t(rhod, axis=axis, kappa=kappa, sigma=sigma)
            index = np.where(rhod > 0)[0]
            bb[index] = sb*T**4
            flux[j, i] = np.trapz(bb*np.exp(-(tau[j,i]-t)), x=axis, dx=dx)
    
    return flux



# Calculate optical depth
n = len(z) 
nu = np.logspace(11,15,n)

tau = calculate_tau(bbb, axis=z, kappa=kappa, sigma=sigma)

plt.contourf(x, y, np.log10(tau), np.linspace(-2,2,256))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Optical Depth')
plt.show()



#Calcualte outgoing flux
#n = len(z) 
#nu = np.logspace(11,15,n)
nx, ny = dim.nx, dim.ny
Nw = 1./(nx*ny)
path = '/Users/daniel/Desktop/Planetary/StreamingInstability_YJ14/txt_files/'

box_mass_codeunits = np.sum(f.rhop)* f.dx * f.dy * f.dz 

H = 5*au 
Lx=np.abs(x[0]-x[-1])*H #Length of the box in cm
Ly=np.abs(y[0]-y[-1])*H #Length of the box in cm
area = Lx*Ly

src_fn = sb*T**4     #Estimate flux according to optical thickness

n = 10
sigma = np.linspace(1, 400, n)
filling_factor = np.zeros(n)
mass_excess = np.zeros(n)
for i in range(n):
    tau = calculate_tau(density=f.rhop, axis=z, kappa=kappa, sigma=sigma[i])
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

