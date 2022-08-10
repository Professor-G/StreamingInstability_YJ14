#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 9 06:08:45 2021

@author: daniel
"""

import astropy.constants as const
warnings.filterwarnings("ignore")
import numpy as np

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

def calc_tau(data, axis, kappa, column_density):
    """
    Integrates density along the z-axis to compute the optical depth.

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
    Returns:
        2D array containing the optical depth values.
    """
    
    Nx = Ny = len(axis)
    dz = np.diff(axis)[0]
    unit_sigma = column_density / np.sqrt(2*np.pi)

    tau = np.zeros([Ny, Nx]) #python reverses order
           
    for i in range(Nx):
        for j in range(Ny):
            surface_density = np.trapz(data[:,j,i]) * dz * unit_sigma
            tau[j, i] = surface_density * kappa
        
    return tau

def calc_t(rhod, axis, kappa, column_density):
    """
    Optical depth with respect to z position along the column.
    Integrates from z to L. This is the optical depth as emission progresses
    up the column, where as the optical_depth() function calculates
    along the entire column. 
    
    Args:
        rhod (ndarray): 1D array of densities, along the given cell column.

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
        t[i] = surface_density * kappa
        
    return t 

def calc_flux(data, axis, kappa, T, column_density):
    """
    Calculate outgoing flux using solution for RT eqn (5.113)
    
    Returns:
        2D array containing the integrated values along the third axis.
    """

    Nx = Ny = Nz = len(axis)
    dz = np.diff(axis)[0]
    unit_sigma = column_density / np.sqrt(2*np.pi)
    
    tau = calc_tau(data, axis, kappa, unit_sigma)
    flux = np.zeros([Ny, Nx])

    #Flux assuming optically thick
    src_fn = const.sigma_sb.cgs.value*T**4     

    for i in range(Nx):
        for j in range(Ny):
            bb = np.zeros(Nz)
            rhod = data[:,j,i] 
            t = calc_t(rhod, kappa, axis, unit_sigma)
            mask = (rhod > 0)
            bb[mask] = src_fn
            flux[j, i] = np.trapz(bb*np.exp(-(tau[j,i]-t)), x=axis, dx=dz)

    return flux 


