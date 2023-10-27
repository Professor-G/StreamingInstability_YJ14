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
from scipy.interpolate import interp1d
from scipy import stats  

class density_cube:
    """
    Class for a density cube object. The class methods
    allow radiative transfer calculations along the z axis
    of the cube, enabling analysis given a range of conditions. 
    If no opacity information is input (kappa and/or sigma), the opacities 
    will be estimated using the DSHARP study, see: https://iopscience.iop.org/article/10.3847/2041-8213/aaf743/pdf

    Note:
        The class methods assume that the simulation is a 3D cube and thus symmetrical about each axis. 

    Args:
        data (ndarry): 3D density cube, containing particle density.
        axis (ndarray): 1D array of the axis along which to integrate.
        column_density : Column density of the gas in cgs units. Defaults to 100 [g / cm^2].
        T (float): Temperature of the entire box in cgs units, as the model is assumed to beisothermal.
            Defaults to 30 [K].
        H (float): The scale height of the box, in cgs units. Defaults to 5 au.
        kappa (float):  Dust absorption opacity coefficient in cgs units. If set to None, then the mm-wave
            dust absorption opacity will be calculated according to column density and grain size,
            using the DSHARP project. Defaults to None.
        sigma (float): Dust scattering opacity coefficient, in cgs units. If None then the mm-wave
            dust scattering opacity will be calculated according to column density and grain size,
            using the DSHARP code. Defaults to None.
        stoke (float): Stoke's number, either a float for monodisperse simulations or an array containing
            one value for each species in the case of polydisperse models.
        rho_grain (float): Internal grain density, in cgs units. A value of 1 [g / cm^3] is typical for ices, 
            and 3.5 [g / cm^3] for silicates. This input can either be a float for monodisperse simulations or 
            an array containing one value for each species in the case of polydisperse models. Must correspond
            with the input Stokes number(s).
        eps_dtog (float): Dust to gas ratio, defaults to 0.03. Only used to calculate the mass of protoplanets.
        npar (int): Number of particles in the simulation, defaults to one million. Used to calculate the mass of protoplanets
            as well as the multi-species weighted opacities. 
        aps (ndarray): pvar.aps, only used to calculate the mass of protoplanets.
        rhopswarm (ndarray): pvar.rhopswarm, only used to calculate the mass of protoplanets.
        include_scattering (bool): Whether scattering should be accounted for. If kappa is set to None, the DSHARP opacities will be applied. If 
            this parameter is set to True, the scattering opacity will be used in addition to the absorption opacities.
            Defaults to False, which will calculate the absorption opacities only.
        ipars (ndarray): The grain identifier array stored in the pvar files as the ipars attribute. Should be input as np.array(fp.ipars, dtype=int).
            Defaults to None. Is only used to compute the weighted opacity when the simulation is polydisperse.
        xp (ndarray): The x-position of the grains corresponding to the ipars array, saved as the xp attribute in the pvar files.
            Defaults to None. Is only used to compute the weighted opacity when the simulation is polydisperse.
        yp (ndarray): The x-position of the grains corresponding to the ipars array, saved as the xp attribute in the pvar files.
            Defaults to None. Is only used to compute the weighted opacity when the simulation is polydisperse.
        zp (ndarray): The z-position of the grains corresponding to the ipars array, saved as the zp attribute in the pvar files.
            Defaults to None. Is only used to compute the weighted opacity when the simulation is polydisperse.
    """
    
    def __init__(self, data=None, axis=None, column_density=100, T=30, H=5*const.au.cgs.value, kappa=None, sigma=None, q=None,
        stoke=0.3, rho_grain=1.0, eps_dtog=0.03, npar=1e6, aps=None, rhopswarm=None, init_var=None, 
        include_scattering=False, ipars=None, xp=None, yp=None, zp=None):

        self.data = data
        self.axis = axis 
        self.column_density = column_density
        self.T = T 
        self.H = H
        self.kappa = kappa 
        self.sigma = sigma 
        self.q = q 
        self.stoke = stoke 
        self.rho_grain = rho_grain
        self.eps_dtog = eps_dtog
        self.npar = npar
        self.aps = aps 
        self.rhopswarm = rhopswarm
        self.init_var = init_var
        self.include_scattering = include_scattering
        self.ipars = ipars  
        self.xp = xp  
        self.yp = yp   
        self.zp = zp

        if self.include_scattering:
            if self.kappa is not None:
                if self.sigma is None:
                    raise ValueError('The include_scattering paramater has been enabled but no scattering coefficient (sigma) was input!')
            if self.sigma is not None:
                if self.kappa is None:
                    raise ValueError('The include_scattering paramater has been enabled but no absorption coefficient (kappa) was input!')

        try:
            __ = len(stoke)
            if __ != len(rho_grain):
                raise ValueError("If entering multiple stoke's numbers, the corresponding rho_grain parameter must be of the same size!")
        except:
            pass

        if self.data is None:
            print('No data input, loading density cube from YJ14, orbit 100...')
            self.data, self.axis = load_cube()

        self.unit_sigma = None
        self.tau = None 
        self.flux = None
        self.mass = None
        self.mass_excess = None
        self.filling_factor = None 
        self.area = None 
        self.grain_size = None 
        self.proto_mass = None 

        #self.configure()

    def configure(self, nu=230e9):
        """
        Initializing parameters and creates flux, mass excess, and filling factor 
        attributes. If sigma, H, or dx/dy/dz attributes are updated, re-run this method 
        to re-configure the settings.

        Args:
            nu (float): Frequency at which to calculate the observed flux and thus mass excess.
        """

        self.Lz = self.Ly = self.Lx = np.abs(self.axis[0] - self.axis[-1])*self.H 
        self.dz = self.dy = self.dx = np.diff(self.axis)[0]
        self.Nz = self.Ny = self.Nx = len(self.axis)
        self.area = self.Lx * self.Ly

        self.unit_sigma = self.column_density / np.sqrt(2*np.pi)

        #Code units
        box_mass_codeunits = np.sum(self.data) if self.init_var is None else np.sum(self.init_var) 
        box_mass_codeunits = box_mass_codeunits * self.dx * self.dy * self.dz 
        unit_mass = self.unit_sigma * self.H**2
        self.mass = box_mass_codeunits * unit_mass 

        if self.kappa is None:
            try:
                self.calc_grain_size(); self.extract_opacity()
            except:
                raise ValueError('Cannot calculate kappa -- to calculate the appropriate grain size input the stoke and rho_grain parameters.')

        self.calc_mass_excess(nu=nu); self.calc_filling_factor()

        if self.aps is not None and self.rhopswarm is not None:
            self.get_proto_mass()

    def blackbody(self, nu=230e9):
        """
        Planck's law, which describes the black body radiation 
        of a source in thermal equilibrium at a given temperature T.

        Args:
            nu (float): Wavelength frequency. Defaults to 230e9 Hz (1 mm)

        Returns:
            Spectral radiance of the source. 
        """

        bb = 2*const.h.cgs.value*nu**3 / (const.c.cgs.value**2*(np.exp(const.h.cgs.value*nu / (const.k_B.cgs.value*self.T)) - 1))

        return bb

    def calc_tau(self):
        """
        Integrates density along the z-axis to compute the optical depth at every single column.
        Computes int_0^Lz k_eff * rhod * dz
        
        Returns:
            2D array containing the optical depth values.
        """

        self.unit_sigma = self.column_density / np.sqrt(2*np.pi)
  
        self.tau = np.zeros([self.Ny, self.Nx]) 
        
        # If polydisperse, create an array to store the number of dust grains in each (x,y) column on a per-species basis
        if isinstance(self.grain_size, np.ndarray):
            print(f"Detected {len(self.grain_size)} grain sizes, applying weighted opacity calculation...")
            species = np.ceil(self.ipars / (self.npar / len(self.grain_size))).astype(int) #Convert to numerical labels, first species is 1, second is 2, etc...
            
            #To store the number of each species located at each (x,y) column, shape (num_of_species, len(axis), len(axis))
            self.num_per_species = np.zeros((len(self.grain_size), self.Nz, self.Ny, self.Nx), dtype=int)

            #To store the weighted opacities at each cell
            self.effective_kappa = np.zeros([self.Nz, self.Ny, self.Nx]) 
            self.effective_sigma = np.zeros([self.Nz, self.Ny, self.Nx]) 

            #Calc the cell indices for each dust grain 
            x_indices = np.digitize(self.xp, self.axis) - 1  #Minus 1 due to 0-based indexing
            y_indices = np.digitize(self.yp, self.axis) - 1  
            z_indices = np.digitize(self.zp, self.axis) - 1

            #Iterate through each dust grain and assign it to a cell
            for i in range(len(self.xp)):
                x_idx = x_indices[i] - 1 #Minus 1 due to 0-based indexing
                y_idx = y_indices[i] - 1  
                z_idx = z_indices[i] - 1 

                #To ensure the particle is within a valid cell range, adjust particles on the outer boundary or beyond
                x_idx = 255 if x_idx >= 256 else x_idx
                y_idx = 255 if y_idx >= 256 else y_idx
                z_idx = 255 if z_idx >= 256 else z_idx

                #Add plus one to the appropriate species array and adjust for 0-based indexing since the species are numbered 1,2,3...
                self.num_per_species[species[i] - 1, z_idx, y_idx, x_idx] += 1 

        for i in range(self.Nx):
            for j in range(self.Ny):
                if isinstance(self.grain_size, np.ndarray) is False: 
                    #Monodisperse calculation
                    surface_density = np.trapz(self.data[:, j, i]) * self.dz * self.unit_sigma
                    self.tau[j, i] = surface_density * (self.kappa + self.sigma) if self.include_scattering else surface_density * self.kappa
                    
                else:
                    # If polydisperse, need to compute the cell-wise averaged opacities = (N1*κ1 + N2*κ2 + N3*κ3 + N4*κ4) / N
                    # Loop through all the z-cells in the particular (x,y) column and store in the self.effective_kappa
                    for k in range(self.Nz):
                        weighted_kappa, weighted_sigma = 0, 0
                        for species_type in range(len(self.grain_size)):
                            weighted_kappa += self.num_per_species[species_type][k, j, i] * self.kappa[species_type]   #N*k...
                            if self.include_scattering:
                                weighted_sigma += self.num_per_species[species_type][k, j, i] * self.sigma[species_type]

                        # Divide by total number of species in that (x,y,z) cell
                        if np.sum(self.num_per_species[:, k, j, i]) == 0:
                            #print(f"No particles present in the (x={i}, y={j}, z={k}) cell, setting opacity to 0...")
                            weighted_kappa = 0
                        else:
                            weighted_kappa /= np.sum(self.num_per_species[:, k, j, i]) #Divide by tot number of species in that (x,y,z) cell

                        self.effective_kappa[k, j, i] = weighted_kappa #Add to the weighted opacity for the entire column

                        if self.include_scattering:
                            if np.sum(self.num_per_species[:, k, j, i]) == 0: #Avoid division by zero 
                                weighted_sigma = 0
                            else:
                                weighted_sigma /= np.sum(self.num_per_species[:, k, j, i]) 

                            self.effective_sigma[k, j, i] = weighted_sigma
                   
                    if self.include_scattering:
                        self.tau[j, i] = np.trapz(self.data[:, j, i] * (self.effective_kappa[:, j, i] + self.effective_sigma[:, j, i])) * self.dz * self.unit_sigma
                    else:
                        self.tau[j, i] = np.trapz(self.data[:, j, i] * self.effective_kappa[:, j, i]) * self.dz * self.unit_sigma

        return 

    def calc_t(self, rhod, effective_kappa, effective_sigma):
        """
        Optical depth with respect to z position along a single column.
        Integrates from z to L. This is the optical depth as emission progresses
        up the column, where as the optical_depth() function calculates
        along the entire column. 
        Computes int_0^z k_eff * rhod * dz
        
        Args:
            rhod (ndarray): 1D array of densities, along the given cell column.
            effective_kappa (ndarray): 1D array of absorption opacities, along the given cell column corresponding to the input rhod.
            effective_sigma (ndarray): 1D array of scattering opacities, along the given cell column corresponding to the input rhod.

        Returns:
            1D array containing the cumulative optical depth along the column.
        """

        self.unit_sigma = self.column_density / np.sqrt(2*np.pi)

        t = np.zeros(self.Nz)

        for i in range(self.Nz):  
            if self.include_scattering:
                t[i] = np.trapz(rhod[:i] * (effective_kappa[:i] + effective_sigma[:i])) * self.dz * self.unit_sigma
            else:
                t[i] = np.trapz(rhod[:i] * effective_kappa[:i]) * self.dz * self.unit_sigma
            
        return t 

    def calc_flux(self):
        """
        Calculate outgoing flux using solution for RT eqn (5.113). Used only when polydisperse simulations are input!
    
        Returns:
            2D array containing the integrated values along the third axis.
        """

        #self.calc_tau()
        if self.tau is None:
            raise ValueError('No optical depth map exists! Run the calc_tau() class method first!')
    
        self.flux = np.zeros([self.Ny, self.Nx])

        #Estimate flux assuming optically thick emission
        #src_fn = const.sigma_sb.cgs.value*self.T**4     

        for i in range(self.Nx):
            for j in range(self.Ny):
                #bb = np.zeros(self.Nz)
                bb = self.src_fn[:, j, i]
                rhod = self.data[:, j, i] 
                kappa = self.effective_kappa[:, j, i]
                sigma = self.effective_sigma[:, j, i]
                t = self.calc_t(rhod, kappa, sigma)
                mask = (rhod == 0)
                bb[mask] = 0
                kappa[mask] = 0
                sigma[mask] = 0
                self.flux[j, i] = np.trapz(bb*np.exp(-(self.tau[j, i] - t)), x=self.axis, dx=self.dx)
    
        return 
        
    def calc_mass_excess(self, nu=230e9):
        """
        Calculates the mass_excess attribute.
        
        Args:
            nu (float): Frequency at which to calculate the flux. Defaults to 1mm frequency.
        
        Returns:
            Float.
        """
        
        self.unit_sigma = self.column_density / np.sqrt(2*np.pi)

        #Calculate mass in cgs units
        box_mass_codeunits = np.sum(self.data) if self.init_var is None else np.sum(self.init_var) 
        box_mass_codeunits = box_mass_codeunits * self.dx * self.dy * self.dz 
        unit_mass = self.unit_sigma * self.H**2 # H is in cgs
        self.mass = box_mass_codeunits * unit_mass # Mass is in grams

        self.calc_tau()

        #Source function should be per frequency (1mm wavelength ~ 230GHz)
        if self.include_scattering:
            # The scattering solution as approximated by Miyake & Nakagawa (1993)
            if isinstance(self.grain_size, np.ndarray) is False:
                albedo = self.sigma / (self.kappa + self.sigma)
            else:
                # The albedo is now 3-dimensional and thus so will the source fn!
                albedo = self.effective_sigma / (self.effective_kappa + self.effective_sigma)
                albedo[~np.isfinite(albedo)] = 0 #Replace all NaN

            epsilon = 1 - albedo
            mu = 1./np.sqrt(3.)
            tau_d = (2*mu) / 3
            tau_ = 0
            numerator = np.exp(-np.sqrt(3 * epsilon) * tau_) + np.exp(np.sqrt(3 * epsilon)*(tau_ - tau_d))
            denominator = (np.exp(-np.sqrt(3 * epsilon) * tau_d) * (1 - np.sqrt(epsilon))) + (np.sqrt(epsilon) + 1)
            J = self.blackbody(nu=nu) * (1 - (numerator / denominator))
            self.src_fn = albedo * J + (1 - albedo) * self.blackbody(nu=nu)
        else:
            self.src_fn = self.blackbody(nu=nu)
        
       # Compute the mass underestimation 
        if isinstance(self.grain_size, np.ndarray) is False:
            # If T is constant, this is what you get from integrating the RT equation with constant source fn and I(0)=0
            flux_approx = self.src_fn * (1 - np.exp(-self.tau))
            self.effective_kappa, self.effective_sigma = np.zeros((256, 256, 256)), np.zeros((256, 256, 256))
            self.effective_kappa[::] = self.kappa 
            self.effective_sigma[::] = self.sigma
            source_fn = np.zeros((256, 256, 256))
            source_fn[::] = self.src_fn
            self.src_fn = source_fn
            self.calc_flux()
            sigma_dust = np.mean(self.flux) / (np.mean(self.src_fn) * (self.kappa + self.sigma)) if self.include_scattering else np.mean(flux_approx) / (self.src_fn * self.kappa)
        else:
            # In this case T is still constant but the source function is now a function of position
            self.calc_flux()
            if self.include_scattering:
                assumed_opacity = (np.sum(self.kappa) + np.sum(self.sigma)) / len(self.grain_size)
                sigma_dust = np.mean(self.flux) / (np.mean(self.src_fn) * assumed_opacity)
            else:
                assumed_opacity = np.sum(self.kappa) / len(self.grain_size) #Assumes all grains in the simulation are equally distributed 
                sigma_dust = np.mean(self.flux) / (np.mean(self.src_fn) * assumed_opacity)

        self.observed_mass = sigma_dust * self.area  
        self.mass_excess = self.mass / self.observed_mass

        return 

    def calc_grain_size(self):
        """
        Calculates grain size given stokes number and 
        gas column density
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

    def extract_opacity(self):
        """
        Returns opacity according to grain size.
        """

        try:
            self.calc_grain_size()
        except:
            raise ValueError('Could not determine grain size, input stoke and rho_grain parameters and try again.')

        a, k_abs, k_sca = load_opacity_values(q=self.q) #Grain size, absorption opacity, and scattering opacity
        k_abs_fit, k_sca_fit = interp1d(a, k_abs), interp1d(a, k_sca)

        if isinstance(self.grain_size, np.ndarray) is False:
            if self.grain_size > a.max():
                raise ValueError('Maximum grain size supported is '+str(a.max())+' cm')
            if self.grain_size < a.min():
                raise ValueError('Minimum grain size supported is '+str(a.min())+' cm')

            if self.include_scattering:
                self.kappa, self.sigma = k_abs_fit(self.grain_size), k_sca_fit(self.grain_size)
            else:
                self.kappa = k_abs_fit(self.grain_size)
        else:
            if self.grain_size.max() > a.max():
                raise ValueError('Maximum grain size supported is '+str(a.max())+' cm')
            if self.grain_size.min() < a.min():
                raise ValueError('Minimum grain size supported is '+str(a.min())+' cm')

            self.kappa, self.sigma = np.zeros(len(self.grain_size)), np.zeros(len(self.grain_size))

            for grain in range(len(self.grain_size)):
                if self.include_scattering:
                    self.kappa[grain], self.sigma[grain] = k_abs_fit(self.grain_size[grain]), k_sca_fit(self.grain_size[grain])
                else:
                    self.kappa[grain] = k_abs_fit(self.grain_size[grain])
        
        return 

    def get_proto_mass(self):
        """
        Calculates the mass of the protoplanets

        Returns:
            Mass of the forming protoplanets, if no planetesimals the get_proto_mass attribute is zero 
        """

        mp_code = self.eps_dtog * self.mass / self.npar
        mp = stats.mode(self.rhopswarm)[0]

        index = np.where(self.aps != 0)[0]
        npclump = self.rhopswarm[index] / mp

        tmp = np.log10(npclump)
        ttmp = tmp + np.log10(mp_code)
        mass = 10**ttmp

        self.proto_mass = np.sort(mass)
        
        try:
            self.proto_mass.max()
        except ValueError:
            self.proto_mass = 0

        return

    def calc_filling_factor(self):
        """
        Calculates the filling factor attribute.

        Returns:
            Float.
        """

        #self.calc_tau()
        
        self.filling_factor = len(np.where(self.tau >= 1)[0]) / (self.Nx * self.Ny)

        return 
        
    """  
    def plot_tau(self, title='Optical Depth', savefig=False, filename='tau'):
        #Plots the optical depth at the exit plane.

        #self.calc_tau()

        plt.contourf(self.axis, self.axis, np.log10(self.tau), np.linspace(-2,2,256))
        plt.colorbar()
        plt.xlabel('x (H)', size=16)
        plt.ylabel('y (H)', size=16)
        plt.title(title, size=18)
        if savefig:
            plt.savefig(filename+'.png', bbox_inches='tight', dpi=300)
            plt.clf()
        else:
            plt.show()

    def plot_flux(self):
        #Plots the outgoing flux at the exit plane.

        #self.calc_flux()

        plt.contourf(self.axis, self.axis, self.flux, 256)
        plt.xlabel('x (H)', size=16)
        plt.ylabel('y (H)', size=16)
        plt.title('Flux', size=18)
        plt.colorbar()
        plt.show()
    """

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

def load_opacity_values(q=None):
    """
    Loads the opacity values taken from the DSHARP project

    Args:
        q (int, optional): The grain distribution power law index. If set to None the
            single grain opacities will be used. Can be 2.5 or 3.5, in which case the opacities
            will correspond to the maximum-grain size opacity from the distribution. Defaults to None. 

    Returns:
        Three arrays -- the grain size followed by the corresponding absorption and scattering opacities.
    """

    resource_package = __name__

    if q is None:
        resource_path = '/'.join(('data', 'a_single_opacities.txt'))
        file = pkg_resources.resource_filename(resource_package, resource_path)
        values = np.loadtxt(file)
        a, k_abs, k_sca = values[:,0], values[:,1], values[:,2]
    elif q == 2.5:
        resource_path = '/'.join(('data', 'a_max_opacities.txt'))
        file = pkg_resources.resource_filename(resource_package, resource_path)
        values = np.loadtxt(file)
        a, k_abs, k_sca = values[:,0], values[:,1], values[:,2]
    elif q == 3.5:
        resource_path = '/'.join(('data', 'a_max_opacities.txt'))
        file = pkg_resources.resource_filename(resource_package, resource_path)
        values = np.loadtxt(file)
        a, k_abs, k_sca = values[:,0], values[:,3], values[:,4]
    elif q == 4.0:
        resource_path = '/'.join(('data', 'a_max_opacities.txt'))
        file = pkg_resources.resource_filename(resource_package, resource_path)
        values = np.loadtxt(file)
        a, k_abs, k_sca = values[:,0], values[:,5], values[:,6]
    else:
        raise ValueError('Unsupported grain size distribution index, options are q=None for single grain size opacities or q=2.5, 3.5, or 4.0.')

    #order = np.array(a).argsort()
    #a, k_abs, k_sca = a[order], k_abs[order], k_sca[order]

    return a, k_abs, k_sca

