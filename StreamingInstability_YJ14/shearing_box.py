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
from scipy.interpolate import UnivariateSpline

class density_cube:
    """
    Class for a density cube object. The class methods
    allow radiative transfer calculations along the z axis
    of the cube, enabling analysis given a range
    of conditions.

    Note:
        The class methods assume that the 3D cube is symmetrical
        about each axis. 

    Args:
        data (ndarry): 3D density cube, containing particle density.
        axis (ndarray): 1D array of the axis along which to integrate.
        column_density : Column density of the gas (g / cm^2). Defaults to 100.
        T (float): Temperature of the entire box, as this model is isothermal.
            Defaults to 30.  
        kappa (float): Dust opacity coefficient (cm^2 / g). Defaults to 1. 
        H (int): Scale height of the box, the value is multiplied by one AU,
            in cgs units. Defaults to 5.
    """
    
    def __init__(self, data=None, axis=None, column_density=100, T=30, kappa=1, H=5,
        Nz=None, Ny=None, Nx=None, dz=None, dy=None, dx=None, Lx=None, Ly=None,
        stoke, rho_grain):

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
        self.stoke = stoke 
        self.rho_grain = rho_grain
        try: 
            len(stoke)
            if len(stoke) != len(rho_grain)
                raise ValueError("If entering multiple stoke's numbers, the corresponding rho_grain paramater must be of same size!")
        except:
            pass 

        if self.data is None:
            print('No data input, automatically loading default density cube.')
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
        self.grain_size = None 


        self.configure()

    def configure(self, nu=230e9):
        """
        Initializing parameters and creates flux, mass excess, and filling factor 
        attributes. If sigma, H, or dx/dy/dz attributes are updated, re-run this method 
        to re-configure the object.

        Args:
            nu (float): Frequency at which to calculate the observed flux and thus mass excess.
        """

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
        self.area = self.Lx * self.Ly
        #Code units
        box_mass_codeunits = np.sum(self.data) * self.dx * self.dy * self.dz 
        unit_mass = self.unit_sigma * self.H**2
        self.mass = box_mass_codeunits * unit_mass 
        self.calc_tau()
        self.calc_mass_excess(nu=nu)
        self.calc_filling_factor()
        try:
            self.calc_grain_size()
        except:
            print('To calculate grain size input stoke and rho_grain parameters.')

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
            self.calc_tau()
    
        flux = np.zeros([self.Ny, self.Nx])
        #Estimate flux assuming optically thick
        src_fn = const.sigma_sb.cgs.value*self.T**4     

        for i in range(self.Nx):
            for j in range(self.Ny):
                bb = np.zeros(self.Nz)
                rhod = self.data[:,j,i] 
                t = self.calc_t(rhod)
                mask = (rhod > 0)
                bb[mask] = src_fn
                flux[j, i] = np.trapz(bb*np.exp(-(self.tau[j,i]-t)), x=self.axis, dx=self.dx)
    
        self.flux = flux 

        return 
        
    def calc_mass_excess(self, nu=230e9, mask_tau=False, threshold=1):
        """
        Calculates the mass_excess attribute.
        
        nu (float): Frequency at which to calculate the flux. Defaults to 1mm frequency.
        mask_tau (bool): If True the mass excess will be measured only in regions that are
            optically thick (tau > 1). Defaults to False.
        threshold (float): If mask_tau is set to True, this paramater can be used
            to set the minimum threshold to use when creating the optical depth mask.
            Defaults to 1.

        Returns:
            Float.
        """

        if self.tau is None:
            self.calc_tau()

        #Source function should be per frequency (1mm wavelength ~ 230GHz)
        src_fn_230 = self.blackbody(nu=nu) 
        
        if mask_tau:
            mask = np.argwhere(self.tau > threshold)
            ratio = len(mask)/(self.Nx*self.Ny)

            rhop_cropped = []
            for i in range(self.Nz):
                rhop_cropped.append(np.sum(self.data[i][mask]))

            cropped_box_mass_codeunits = np.sum(rhop_cropped) * self.dx * self.dy * self.dz 
            unit_mass = self.unit_sigma * self.H**2
            self.mass = cropped_box_mass_codeunits * unit_mass 

            #If source fn is constant and region is optically thick (Eq. 5.120)
            flux_approx = src_fn_230 * (1-np.exp(-self.tau[mask]))
        
            #Sigma dust observer sees if optically thin 
            sigma_dust = np.mean(flux_approx) / (src_fn_230*self.kappa)

            self.observed_mass = sigma_dust*self.area*ratio
            self.mass_excess = self.mass / self.observed_mass

            return 
        
        flux_approx = src_fn_230 * (1-np.exp(-self.tau))
        sigma_dust = np.mean(flux_approx) / (src_fn_230*self.kappa)

        self.observed_mass = sigma_dust*self.area  
        self.mass_excess = self.mass / self.observed_mass

        return 

    def calc_grain_size():
        """
        Calculates grain size given stokes number and 
        gas column density.

        st (float): Stoke's number
        rho_g (float): Internal grain density, approximately
            1 g/cm^3 for ices, and 3.5 g/cm^3 for silicates
        """

        if isinstance(self.stoke, np.ndarray) is False:
            self.grain_size = self.stoke * 2. * self.column_density / np.pi / self.rho_grain
        else:
            if isinstance(self.rho_grain, np.ndarray) is False:
                raise ValueError("If entering multiple stoke's numbers, the corresponding rho_grain paramater must be of same size!")
            self.grain_size = np.zeros(len(self.stoke))
            for grain in range(len(self.grain_size)):
                self.grain_size[grain] = self.stoke[grain] * 2. * self.column_density / np.pi / self.rho_grain[grain]
        
        return

    def extract_opacity():
        """
        Returns opacity according to grain size
        and wavelength
        """
        if self.grain_size is None:
            try:
                self.calc_grain_size()
            except:
                raise ValueError('Could not determine grain size, input stoke and rho_grain parameters and try again.')

        a, k_abs, k_sca = load_fig4_values()
        k_abs_fit = UnivariateSpline(a, k_abs, k=3)
        k_sca_fit = UnivariateSpline(a, k_sca, k=3)

        if isinstance(self.grain_size, np.ndarray) is False:
            if self.grain_size > a.max():
                raise ValueError('Maximum grain size supported is '+str(a.max()))
            if self.grain_size < a.min():
                raise ValueError('Minimum grain size supported is '+str(a.min()))
            self.kappa = k_abs_fit(self.grain_size) + k_sca_fit(self.grain_size)

        else:
            if self.grain_size.max() > a.max():
                raise ValueError('Maximum grain size supported is '+str(a.max()))
            if self.grain_size.min() < a.min():
                raise ValueError('Minimum grain size supported is '+str(a.min()))
            self.kappa = np.zeros(len(self.grain_size))
            for grain in self.grain_size:
                self.kappa[grain] = k_abs_fit(self.grain_size[grain]) + k_sca_fit(self.grain_size[grain])

    def calc_filling_factor(self):
        """
        Calculates the filling factor attribute.

        Returns:
            Float.
        """

        if self.tau is None:
            self.calc_tau()

        self.filling_factor = len(np.where(self.tau >= 1)[0]) * self.Nw

        return 
        
    def plot_tau(self):
        """
        Plots the optical depth at the exit plane.
        """

        if self.tau is None:
            self.calc_tau()

        plt.contourf(self.axis, self.axis, np.log10(self.tau), np.linspace(-2,2,256))
        plt.colorbar()
        plt.xlabel('x (H)', size=16)
        plt.ylabel('y (H)', size=16)
        plt.title('Optical Depth', size=18)
        plt.show()

    def plot_flux(self):
        """
        Plots the outgoing flux at the exit plane.
        """
        if self.flux is None:
            self.calc_flux()

        plt.contourf(self.axis, self.axis, self.flux, 256)
        plt.xlabel('x (H)', size=16)
        plt.ylabel('y (H)', size=16)
        plt.title('Flux', size=18)
        plt.colorbar()
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

def load_fig4_values():
    """
    Loads the opacity values taken from DSHARP
    """
    resource_package = __name__
    resource_path = '/'.join(('data', 'dsharp_fig4_values'))
    file = pkg_resources.resource_filename(resource_package, resource_path)
    values = np.loadtxt(file)
    a, k_abs, k_sca = values[:,0], values[:,1], values[:,2]

    order = np.array(a).argsort()
    a, k_abs, k_sca = a[order], k_abs[order], k_sca[order]

    return a, k_abs, k_sca

