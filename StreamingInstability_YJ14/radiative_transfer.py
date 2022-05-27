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

h,k,c,sb,au= 6.626e-27, 1.3807e-16, 2.997e10, 5.67e-5,  const.au.cgs.value

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


