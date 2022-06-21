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

class density_cube:
    """
    Class for a density cube object. The class methods
    allow radiative transfer calculations along the z axis
    of the cube, enabling analysis given a range
    of conditions.

    Note:
        The class methods assume that the 3D cube is symmetrical and
        about each axis. 

    Args:
        density (ndarry): 3D density cube.
        axis (ndarray): 1D array of the axis along which to integrate.
        column_density : Column density of the gas (g / cm^2)
        T (float): Temperature of the entire box, as this model is isothermal.
            Defaults to 30.  
        kappa (float): Dust opacity coefficient (cm^2 / g)
        sigma (float): Column density of the gas 
        H (int): Scale height of the box, the value is multiplied by one AU.
            Defaults to 5.
    """
    
    def __init__(self, data=None, axis=None, column_density=100, T=30, kappa=1, H=5,
        Nz=None, Ny=None, Nx=None, dz=None, dy=None, dx=None, Lx=None, Ly=None):

        self.data = data
        self.axis = axis 
        self.column_density = column_density
        self.T = T 
        self.kappa = kappa
        self.H = H*const.au.cgs.value
        self.Lx = Lx 
        self.Ly = Ly
        self.Nz = Nz
        self.Ny = Ny 
        self.Nx = Nx 
        self.dz = dz 
        self.dy = dy 
        self.dx = dx 

        if self.data is None:
            print('No data input, automatically loading density cube.')
            if self.axis is None:
                self.data, self.axis = load_cube()
            else:
                self.data = load_cube()[0]

        self.unit_sigma = self.column_density / np.sqrt(2*np.pi)
        self.tau = None 
        self.flux = None
        self.mass = None
        self.mass_excess = None
        self.filling_factor = None 
        self.Nw = None
        self.area = None 
        self.mass = None

        self.configure()

    def configure(self):
        """
        Initializing parameters and creates flux, mass excess, and filling factor 
        attributes. If sigma, H, or dx/dy/dz attributes are updated, re-run this method 
        to re-configure the object.
        """
       #print('Initializing parameters and attributes...')
        if self.Nz is None:
            self.Nz = len(self.axis)
        if self.Ny is None:
            self.Ny = len(self.axis)
        if self.Nx is None:
            self.Nx = len(self.axis)

        self.Nw = 1./(self.Nx * self.Ny)

        if self.dz is None:
            self.dz = np.diff(self.axis)[0]
        if self.dy is None:
            self.dy = np.diff(self.axis)[0]
        if self.dx is None:
            self.dx = np.diff(self.axis)[0]

        if self.Lx is None:
            self.Lx = np.abs(self.axis[0] - self.axis[-1])*self.H 
        if self.Ly is None: 
            self.Ly = np.abs(self.axis[0] - self.axis[-1])*self.H 

        self.unit_sigma = self.column_density / np.sqrt(2*np.pi)
        #Observed mass of the dust in the box
        self.area = self.Lx * self.Ly
            
        #Code units
        box_mass_codeunits = np.sum(self.data)* self.dx * self.dy * self.dz 
        unit_mass = self.unit_sigma * self.H**2
        self.mass = box_mass_codeunits * unit_mass 
        self.calc_tau()
        #self.calc_flux()
        self.calc_mass_excess()
        self.calc_filling_factor()

        return 

    def blackbody(self, nu):
        """
        Planck's law, which describes the black body radiation 
        of a source in thermal equilibrium at a given temperature T.

        Args:
            nu (ndarray): Array of frequencies.

        Returns:
            Spectral radiance of the source. 
        """

        bb = 2*const.h.cgs.value*nu**3 / (const.c.cgs.value**2*(np.exp(const.h.cgs.value*nu / (const.k_B.cgs.value*self.T)) - 1))

        return bb

    def calc_tau(self):
        """
        Integrates density along the z-axis to compute the optical depth.
        
        Returns:
            2D array containing the optical depth values.
        """
        
        tau = np.zeros([self.Ny, self.Nx]) #python reverses order
               
        for i in range(self.Nx):
            for j in range(self.Ny):
                surface_density = np.trapz(self.data[:,j,i]) * self.dz * self.unit_sigma
                tau[j, i] = surface_density * self.kappa
            
        self.tau = tau

        return 

    def calc_t(self, rhod):
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
        
        t = np.zeros(self.Nz)
        surface_density = 0

        for i in range(self.Nz):        
            surface_density += rhod[i] * self.dz * self.unit_sigma
            t[i] = surface_density * self.kappa
            
        return t 

    def calc_flux(self):
        """
        Calculate outgoing flux using solution for RT eqn (5.113)
        
        Returns:
            2D array containing the integrated values along the third axis.
        """

        if self.tau is None:
            self.calc_tau(self)
    
        flux = np.zeros([self.Ny, self.Nx])
        #Estimate flux assuming optically thick
        src_fn = const.sigma_sb.cgs.value*self.T**4     

        for i in range(self.Nx):
            for j in range(self.Ny):
                bb = np.zeros(self.Nz)
                rhod = self.data[:,j,i] #1D array
                t = self.calc_t(rhod)
                mask = (rhod > 0)
                bb[mask] = src_fn
                flux[j, i] = np.trapz(bb*np.exp(-(self.tau[j,i]-t)), x=self.axis, dx=self.dx)
    
        self.flux = flux 

        return 
        
    def calc_mass_excess(self):
        """
        Calculates the mass_excess attributes.

        Returns:
            Float.
        """

        if self.tau is None:
            calc_tau(self)

        #Source function should be per frequency (1mm wavelength ~ 230GHz)
        src_fn_230 = self.blackbody(nu=230e9) 
        
        #If source fn is constant and region is optically thick (Eq. 5.120)
        flux_approx = src_fn_230 * (1-np.exp(-self.tau))
        
        #Sigma dust observer sees if optically thin 
        sigma_dust = np.mean(flux_approx) / (src_fn_230*self.kappa)

        self.observed_mass = sigma_dust*self.area  
    
        #Actual mass in box 
        self.mass_excess = self.mass / self.observed_mass

        return 
        
    def calc_filling_factor(self):
        """
        Calculates the filling factor attribute.

        Returns:
            Float.
        """

        if self.tau is None:
            self.calc_tau(self)

        self.filling_factor = len(np.where(self.tau >= 1)[0]) * self.Nw

        return 

        
    def plot_tau(self):
        """
        Plots the optical depth at the exit plane.
        """

        if self.tau is None:
            self.calc_tau(self)

        plt.contourf(self.axis, self.axis, np.log10(self.tau), np.linspace(-2,2,256))
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Optical Depth')
        plt.show()

    def plot_flux(self):
        """
        Plots the outgoing flux at the exit plane.
        """
        if self.flux is None:
            self.calc_flux(self)

        plt.contourf(self.axis, self.axis, self.flux, 256)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Flux')
        plt.colorbar()
        plt.show()

def plot_filling_factor(filling_factor, sigma):
    """
    Plots Gas Column Density vs Filling Factor

    Args:
        sigma = column density
    """


    plt.plot(sigma, filling_factor, 'ro-')
    plt.xscale('log')
    plt.xlabel(r'$\Sigma_g \ (g / cm^2)$', size=20)
    plt.ylabel('Filling Factor', size=20)
    plt.show()

def plot_mass_excess(column_density, mass_excess):
    """
    Plots Gas Column Density vs Mass Excess
    Column density of gas
    """

    plt.plot(column_density, mass_excess, 'ro-')
    plt.xscale('log')
    plt.xlabel(r'$\Sigma_g \ (g / cm^2)$', size=20)
    plt.ylabel('Mass Excess', size=20)
    plt.show()

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

    data, axis = np.r_[density_cube_1, density_cube_2], np.loadtxt(file)

    return data, axis


