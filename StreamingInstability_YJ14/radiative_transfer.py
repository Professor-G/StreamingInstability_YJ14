#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 9 06:08:45 2021

@author: daniel
"""
import numpy as np
import astropy.constants as const
import warnings; warnings.filterwarnings("ignore")

####################
### These are standalone functions, for efficiency these
### have been integrated as methods in the shearing_box.density_cube() Class
####################

def blackbody(nu, T):
    """
    Planck's law, which describes the black body radiation 
    of a source in thermal equilibrium at a given temperature T.

    Args:
        nu (ndarray): Array of frequencies.
        T (float): Temperature of the entire box, as this model asummes isothermal
            disk conditions.

    Returns:
        Spectral radiance of the source. 
    """

    bb = 2*const.h.cgs.value*nu**3 / (const.c.cgs.value**2*(np.exp(const.h.cgs.value*nu / (const.k_B.cgs.value*T)) - 1))

    return bb

def calc_tau(data, axis, kappa, sigma, column_density):
    """
    Integrates density along the z-axis to compute the optical depth.

    Args:
        data (ndarry): 3D density cube.
        axis (ndarray): 1D array of the axis along which to integrate.
        kappa (float): Dust opacity coefficient in cm^2 / g.
        sigma (float): Dust scattering opacity coefficient in cm^2 / g.
        column_density (float): Column density of the gas in g / cm^2.
       
    Returns:
        2D array containing the optical depth values.
    """
    
    Nx = Ny = len(axis)
    dz = np.diff(axis)[0]
    unit_sigma = column_density / np.sqrt(2*np.pi)

    tau = np.zeros([Ny, Nx]) 
           
    for i in range(Nx):
        for j in range(Ny):
            surface_density = np.trapz(data[:,j,i]) * dz * unit_sigma
            chi = kappa + sigma 
            tau[j, i] = surface_density * chi 
        
    return tau

def calc_t(rhod, axis, kappa, sigma, column_density):
    """
    Optical depth with respect to z position along the column.
    Integrates from z to L. This is the optical depth as emission progresses
    up the column, where as the optical_depth() function calculates
    along the entire column. 
    
    Args:
        rhod (ndarray): 1D array of densities, along the given cell column.
        axis (ndarray): 1D array of the axis along which to integrate.
        kappa (float): Dust opacity coefficient in cm^2 / g.
        sigma (float): Dust scattering opacity coefficient in cm^2 / g.
        column_density (float): Column density of the gas in g / cm^2.

    Returns:
        1D array containing the cumulative optical depth along the column.
    """

    Nz = len(axis)
    dz = np.diff(axis)[0]
    unit_sigma = column_density / np.sqrt(2*np.pi)

    t = np.zeros(Nz)
    surface_density = 0

    for i in range(Nz):        
        surface_density += rhod[i] * dz * unit_sigma
        t[i] = surface_density * (kappa + sigma)
        
    return t 

def calc_flux(data, axis, kappa, sigma, T, column_density, nu=230e9):
    """
    Calculate outgoing flux using solution for RT eqn (5.113)
    
    Args:
        data (ndarry): 3D density cube.
        axis (ndarray): 1D array of the axis along which to integrate.
        kappa (float): Dust opacity coefficient in cm^2 / g.
        sigma (float): Dust scattering opacity coefficient in cm^2 / g.
        T (float): Temperature at which the evaluate the blackbody emission. 
        column_density (float): Column density of the gas in g / cm^2.
        nu (float): The frequency at which the evaluate the blackbody emission. 

    Returns:
        2D array containing the integrated values along the third axis.
    """

    Nx = Ny = Nz = len(axis)
    dz = np.diff(axis)[0]
    unit_sigma = column_density / np.sqrt(2*np.pi)
    
    tau = calc_tau(data=data, axis=axis, kappa=kappa, sigma=sigma, column_density=column_density)
    flux = np.zeros([Ny, Nx])

    #Flux assuming optically thick
    #src_fn = const.sigma_sb.cgs.value*T**4     

    for i in range(Nx):
        for j in range(Ny):
            bb = np.zeros(Nz)
            rhod = data[:,j,i] 
            t = calc_t(rhod=rhod, axis=axis, kappa=kappa, sigma=sigma, column_density=column_density)
            mask = (rhod > 0)
            chi = kappa + sigma 
            albedo = sigma / chi 
            epsilon = 1 - albedo 
            tau_d = 0
            denominator = np.exp(-np.sqrt(3*epsilon)*tau_d) + np.exp(np.sqrt(3*epsilon)*(tau_d-tau_d))
            numerator = np.exp(-np.sqrt(3*epsilon)*tau_d)*(np.sqrt(epsilon) - 1) - (np.sqrt(epsilon) + 1)
            J = 1 + (denominator/numerator)
            src_fn = albedo * J + (1 - albedo) * blackbody(nu, T)
            bb[mask] = src_fn
            flux[j, i] = np.trapz(bb*np.exp(-(tau[j,i]-t)), x=axis, dx=dz)

    return flux 















