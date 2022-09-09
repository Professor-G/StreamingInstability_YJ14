#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 9 06:08:45 2021

@author: daniel
"""
import numpy as np
import astropy.constants as const
import warnings 
warnings.filterwarnings("ignore")
from StreamingInstability_YJ14 import radiative_transfer 

####################
### These are standalone functions, for efficiency these
### have been integrated as methods in the shearing_box.density_cube() Class
####################

def calc_mass(data, axis, column_density, H):
    """
    Calculates the mass inside the density cube.
    
    Args:
        data (ndarray): 3D density cube.
        axis (ndarray): 1D array corresponding to the axis along which to integrate.
        column_density : Column density of the gas in g / cm^2.
        H: The scale height, this parameter must be in AU.

    Return:
        Entire mass inside the cube, in cgs units.
    """

    H = H*const.au.cgs.value
    Nx = Ny = Nz = len(axis)
    dx = dy = dz = np.diff(axis)[0]
    Lx = Ly = np.abs(axis[0] - axis[-1])*H 

    unit_sigma = column_density / np.sqrt(2*np.pi)
    area = Lx*Ly
    box_mass_codeunits = np.sum(data) * dx * dy * dz 
    unit_mass = unit_sigma * H**2
    mass = box_mass_codeunits * unit_mass 

    return mass 


def calc_mass_excess(data, axis, kappa, column_density, T, H, nu=230e9, mask_tau=False, threshold=1):
    """
    Calculates the mass_excess attribute.
    
    Args:
        data (ndarry): 3D density cube.
        axis (ndarray): 1D array of the axis along which to integrate.
        kappa (float): Dust opacity coefficient in cm^2 / g.
        column_density : Column density of the gas in g / cm^2.
        T (float): Temperature of the entire box, as this model asummes isothermal
            disk conditions.
        H: The scale height, this parameter must be in AU.
        nu (float): Frequency at which to calculate the flux, in Hz. 
            Defaults to 230 GHz, corresponding to 1mm emission.
        mask_tau (bool): If True the mass excess will be measured only in regions that are
            optically thick (tau > threshold). Defaults to False.
        threshold (float): If mask_tau is set to True, this paramater can be used
            to set the minimum threshold to use when creating the optical depth mask.
            Defaults to 1.

    Note:
        If mask_tau=True, the mass excess will be calculated only in optically thick regions, 
        which is not observationally possible as that would require resolution so high that 
        the filamentary structures that come about as a result of streaming instability could be 
        distinguished from one another. 

    Returns:
        Float.
    """

    mass = calc_mass(data=data, axis=axis, column_density=column_density, H=H)

    H = H*const.au.cgs.value
    Nx = Ny = Nz = len(axis)
    dx = dy = dz = np.diff(axis)[0]
    Lx = Ly = np.abs(axis[0] - axis[-1])*H 
    area = Lx*Ly

    tau = radiative_transfer.calc_tau(data=data, axis=axis, kappa=kappa, column_density=column_density)

    #Source function should be per frequency (1mm wavelength ~ 230GHz)
    src_fn = radiative_transfer.blackbody(nu=nu, T=T) 
    
    if mask_tau: 
        mask = np.argwhere(tau > threshold)
        ratio = len(mask)/(Nx*Ny)

        rhop_cropped = []
        for i in range(Nz):
            rhop_cropped.append(np.sum(data[i][mask]))

        unit_sigma = column_density / np.sqrt(2*np.pi)
        unit_mass = unit_sigma * H**2

        cropped_box_mass_codeunits = np.sum(rhop_cropped) * dx * dy * dz 
        cropped_mass = cropped_box_mass_codeunits * unit_mass 

        #If source fn is constant and region is optically thick (Eq. 5.120)
        flux_approx = src_fn * (1-np.exp(-tau[mask]))
    
        #Sigma dust observer sees if optically thin 
        sigma_dust = np.mean(flux_approx) / (src_fn*kappa)

        observed_mass = sigma_dust*area*ratio
        mass_excess = cropped_mass / observed_mass

        return mass_excess 
    
    flux_approx = src_fn * (1-np.exp(-tau))
    sigma_dust = np.mean(flux_approx) / (src_fn*kappa)

    observed_mass = sigma_dust*area  
    mass_excess = mass / observed_mass

    return mass_excess


def getmass(radius, npar, nx, rhopswarm, column_density):

    i = np.where(radius!=0)
    eps_dtog=0.03
    #column_density = 44.565547740180705  # g/cm2                                                                                                                                                         
    au = 1.49e13
    r = 20*au
    Hp = 0.1*r
    Lx = 4*Hp
    Ly = 16*Hp
    Lz = 2*Hp
    area=Lx*Ly
    mass_box = Lx*Ly*column_density
    dx = Lx/nx
    dy = Ly/nx
    dz = Lz/nx
    sigma = column_density
    rho0 = sigma / np.sqrt(2*np.pi) / Hp
    cell_volume = dx*dy*dz
    mp_code = eps_dtog * mass_box / npar
    mp = stats.mode(rhopswarm)[0]
    npclump = rhopswarm[i]/mp
    Mmars = 6.39e26
    Mearth=5.972e27
    tmp = np.log10(npclump)
    fac = mp_code                                                                                                                                                                                 
    ttmp = tmp + np.log10(fac)

    mass = 10**ttmp

    return np.sort(mass/Mearth)


#mass = getmass(fp.aps,npar,256,fp.rhopswarm)
#fp = pc.read_pvar(datadir=datadir,varfile=‘PVAR’+str(ivar))
#pdim=pc.read_pdim(datadir=datadir)
#npar = pdim.npar


