#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:01:11 2021

@author: daniel
"""
import warnings
warnings.filterwarnings("ignore")
from scipy import interpolate, stats
import astropy.constants as const
from pathlib import Path
import pkg_resources
import numpy as np
import math

class density_cube:
    """
    Class for a density cube object. The class methods
    allow radiative transfer calculations along the z axis
    of the cube, enabling analysis given a range of conditions. 
    If no opacity information is input (kappa and/or sigma), the opacities 
    will be estimated using coefficients from the DSHARP project, see: https://iopscience.iop.org/article/10.3847/2041-8213/aaf743/pdf

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
        grain_rho (float): Internal grain density, in cgs units. A value of 1 [g / cm^3] is typical for ices, 
            and 3.5 [g / cm^3] for silicates. This input can either be a float for monodisperse simulations or 
            an array containing one value for each species in the case of polydisperse models. Must correspond
            with the input Stokes number(s).
        eps_dtog (float): Dust to gas ratio, defaults to 0.03. Only used to calculate the mass of protoplanets.
        npar (int): Number of particles in the simulation, defaults to one million. Used to calculate the mass of protoplanets
            as well as the multi-species weighted opacities.
        code_rho (float): The midplane gas density in code units (rho0 param in the start.in file, &eos_init_pars). Defaults to 1.
        code_cs (float): The sound speed of the gas in code units (cs0 param in the start.in file, &eos_init_pars). Defaults to 2pi.
        code_omega (float): The Keplerian orbital timescale in code units (Omega param the start.in file, &density_init_pars). Defaults to 2pi.
        aps (ndarray): pvar.aps, only used to calculate the mass of protoplanets.
        rhopswarm (ndarray): pvar.rhopswarm, only used to calculate the mass of protoplanets.
        init_var (ndarray): 3D density cube of the first (initial) snapshot so that the initial mass of the cube is known and used
            for the mass underestimation. This is only applicable if the simulation includes self-gravity as the dust mass inside 
            the box will decrease as sink particles are formed and grow. Defaults to None.
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
        stoke=0.3, grain_rho=1.0, eps_dtog=0.03, npar=1e6, code_rho=1, code_cs=2*np.pi, code_omega=2*np.pi, aps=None, rhopswarm=None, init_var=None, 
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
        self.grain_rho = grain_rho
        self.eps_dtog = eps_dtog
        self.npar = npar
        self.code_rho = code_rho
        self.code_cs = code_cs 
        self.code_omega = code_omega 
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
                if self.sigma is None: raise ValueError('The include_scattering paramater has been enabled but no scattering coefficient (sigma) was input!')
            if self.sigma is not None:
                if self.kappa is None: raise ValueError('The include_scattering paramater has been enabled but no absorption coefficient (kappa) was input!')

        if self.data is None:
            print('No data input, loading density cube from YJ14, orbit 100...')
            self.data, self.axis = load_cube()

        self.tau = None 
        self.flux = None
        self.mass = None
        self.mass_excess = None
        self.filling_factor = None 
        self.area = None 
        self.grain_size = None 
        self.proto_mass = None 

    def configure(self, frequency=3e11):
        """
        Initializing parameters and creates flux, mass excess, and filling factor 
        attributes. If sigma, H, or dx/dy/dz attributes are updated, re-run this method 
        to re-configure the settings.

        Note:
            If running this method multiple times with different astrocentric parameters,
            ensure that the kappa and sigma parameters are reset to None if they were None 
            initially! Failure to do so would result in the grain size and thus the opacity 
            coefficients not properly updating. It is recommended to re-initialize the class
            if re-running with varying parameters.

        Args:
            frequency (float): Wavelength requency at which to calculate the observed flux and thus mass excess.
                Input must be Hz, defaults to 3e11 corresponding to the 1mm frequency.
        """

        # Dimensions of the box
        self.Lz = self.Ly = self.Lx = np.abs(self.axis[0] - self.axis[-1]) * self.H # Box length (assumes a cube!)
        self.dz = self.dy = self.dx = np.diff(self.axis)[0] # Cell length
        self.Nz = self.Ny = self.Nx = len(self.axis) # No. of cells
        self.area = self.Lx * self.Ly # Box area in cgs code units

        # Mass in code units
        box_mass_codeunits = np.sum(self.data) if self.init_var is None else np.sum(self.init_var) # init_var is the 0th snapshot and it's used for simulations with self-gravity (the initial mass)
        box_mass_codeunits = box_mass_codeunits * self.dx * self.dy * self.dz 
        
        # Convert the mass to cgs units
        self.unit_mass = (self.column_density * self.H**2) / math.sqrt(2 * math.pi) / (self.code_rho * (self.code_cs / self.code_omega)**3)
        self.mass = box_mass_codeunits * self.unit_mass 

        # To convert the dust surface density to cgs units (used when integrating to solve for tau and the RT equation)
        #self.unit_sigma = self.column_density / math.sqrt(2 * math.pi) / (self.code_rho * (self.code_cs / self.code_omega)**2)
        self.unit_sigma = self.column_density / math.sqrt(2 * math.pi) / (self.code_rho * (self.code_cs / self.code_omega))

        #self.unit_mass = self.column_density * self.H**2 / np.sqrt(2*np.pi) 
        #self.mass = box_mass_codeunits * self.unit_mass # Mass is in cgs
   
        # If opacity coefficients are not input calculate the grain sizes and extract the corresponding DSHARP opacities
        if self.kappa is None:
            try:
                self.calc_grain_size(); self.extract_opacity()
            except:
                raise ValueError('Cannot calculate kappa -- to calculate the appropriate grain size input the stoke and grain_rho parameters.')

        # Compute the mass underestimation 
        self.calc_mass_excess(frequency=frequency)

        # Compute the mass of the protoplanets if both the aps and rhopswarm arrays are input
        if self.aps is not None and self.rhopswarm is not None: self.get_proto_mass()

        return

    def blackbody(self, frequency=3e11):
        """
        Planck's law describing the black body radiation of a source in thermal equilibrium at a given temperature, T.
        
        Args:
            frequency (float): Fequency of the wavelength at which to calculate the observed flux and thus the mass excess.
                Input must be in units of Hz. Defaults to 3e11 Hz corresponding to the 1 mm wavelength.

        Returns:
            Spectral radiance of the source. 
        """

        return 2 * const.h.cgs.value * frequency**3 / (const.c.cgs.value**2 * (np.exp(const.h.cgs.value * frequency / (const.k_B.cgs.value * self.T)) - 1))

    def calc_tau(self):
        """
        Integrates density along the z-axis to compute the optical depth at every single column.
        
        Returns:
            2D array containing the optical depth values at each (x,y) position.
        """
  
        self.tau = np.zeros((self.Ny, self.Nx)) # To store the optical depth at each (x,y) column
        
        # If polydisperse, create a 'num_per_species' 4D array to store the number of dust grains in each 3D cell on a per-species basis
        if isinstance(self.grain_size, np.ndarray):
            # Convert the ipars array to numerical labels, first species is 1, second is 2, etc...
            species = np.ceil(self.ipars / (self.npar / len(self.grain_size))).astype(int) 
            num_species = len(np.unique(species))
            print(f"Detected {num_species} grain sizes, applying weighted opacity calculation...")

            #To store the number of each species located at each cell, shape (num_of_species, len(Nz), len(Ny), len(Nx))
            self.num_per_species = np.zeros((num_species, self.Nz, self.Ny, self.Nx), dtype=int)

            #To store the weighted opacity coefficients at each cell
            self.effective_kappa = np.zeros((self.Nz, self.Ny, self.Nx))
            self.effective_sigma = np.zeros((self.Nz, self.Ny, self.Nx))
            
            ###
            ### Compute the number of each species at each cell ###
            ###

            # This is a slow but very "fundamental" method as it manually loops through each individual (x,y,z) cell and counts the number of grains in the cell
            # NOTE: The cells are defined according to the grains within plus or minus half a cell from the specific axis position

            # The positions are assumed to be equally spaced out -- the same for all three axes
            axis_grid_x = axis_grid_y = axis_grid_z = self.axis

            # Half cell length
            half_cell_length = np.diff(self.axis)[0] / 2.

            for species_type in range(1, num_species+1):
                # Index the particular species
                index_species = np.where(species == species_type)[0]
                
                # Extract the (x,y,z) positions for all of these grains
                species_x, species_y, species_z = self.xp[index_species], self.yp[index_species], self.zp[index_species]
                
                # Start with the 0th x-axis position and find all the grains within plus or minus half a cell length of this 0th position, this is
                # the first cell. Loop through and go up one cell at a time. 
                for xx in range(self.Nx): 
                    print(f"Species {species_type}: {xx+1} of {self.Nx}")
                    # Find the grains in a given x-cell by setting minimum and maximum boundaries
                    xmin = axis_grid_x[xx] - half_cell_length # The left edge of the x-cell we are looking at
                    xmax = axis_grid_x[xx] + half_cell_length # The right edge of the x-cell we are looking at
                    
                    # Find all the dust grains that are within these x-positions (no need for leq or geq as very unlikely a grain would be at the exact bin edge)
                    index_x = np.where((species_x > xmin) & (species_x < xmax))[0]
                    
                    # Now loop through all the y-positions
                    for yy in range(self.Ny):
                        # Find the grains in a given y-cell by setting minimum and maximum boundaries
                        ymin = axis_grid_y[yy] - half_cell_length # The left edge of the y-cell we are looking at
                        ymax = axis_grid_y[yy] + half_cell_length # The right edge of the y-cell we are looking at
                        
                        # Find all the dust grains that are in the x-cell and in this newly defined y-position 
                        index_y = np.where((species_y[index_x] > ymin) & (species_y[index_x] < ymax))[0]
                        
                        # Finally loop through all the z-positions
                        for zz in range(self.Nz):
                            # Find the grains in a given z-cell by setting minimum and maximum boundaries
                            zmin = axis_grid_z[zz] - half_cell_length # The left edge of the z-cell we are looking at
                            zmax = axis_grid_z[zz] + half_cell_length # The right edge of the z-cell we are looking at
                            
                            # Find all the dust grains that are within both x and y-positions as defined above, and this newly defined z-position
                            index_z = np.where((species_z[index_x[index_y]] > zmin) & (species_z[index_x[index_y]] < zmax))[0]
                            
                            # Update num_per_species array with the count which is the total number of grains in index_z 
                            self.num_per_species[species_type - 1, zz, yy, xx] = len(index_z)

        ###
        ### Calculate the optical depth tau = kappa * rhop * dz ###
        ###
        for i in range(self.Nx):
            for j in range(self.Ny):

                # This is the Monodisperse case
                if isinstance(self.grain_size, np.ndarray) is False:
                    surface_density = np.trapz(self.data[:, j, i]) * self.dz * self.unit_sigma
                    self.tau[j, i] = surface_density * (self.kappa + self.sigma) if self.include_scattering else surface_density * self.kappa

                #This is the polydisperse case in which the cell-wise averaged opacities must be calculated (N1*κ1 + N2*κ2 + N3*κ3 + N4*κ4) / N
                else:
                    # The below code goes through all the z-cells in the particular (x,y) column to store in the self.effective_kappa array
                    for k in range(self.Nz):

                        # To count the averaged opacities at each (x,y,z) cell (equivalent to N*k)
                        weighted_kappa, weighted_sigma = 0, 0

                        # Add the number of grains in that cell for a given species times that species' opacity coefficient (stored in the self.kappa/self.sigma arrays)
                        for species_type in range(num_species):
                            weighted_kappa += self.num_per_species[species_type][k, j, i] * self.kappa[species_type] # This is the numerator (N1 * k1 + N2*k2 + ...)
                            weighted_sigma = weighted_sigma + self.num_per_species[species_type][k, j, i] * self.sigma[species_type] if self.include_scattering else 0
                        
                        # Divide by tot number of species in that (x,y,z) cell to compute the weighted mean and add to the 3D opacity array
                        self.effective_kappa[k, j, i] = weighted_kappa / np.sum(self.num_per_species[:, k, j, i]) if np.sum(self.num_per_species[:, k, j, i]) != 0 else 0 #Avoid division by zero 
                        if self.include_scattering:
                            self.effective_sigma[k, j, i] = weighted_sigma / np.sum(self.num_per_species[:, k, j, i]) if np.sum(self.num_per_species[:, k, j, i]) != 0 else 0 #Avoid division by zero 
                   
                    # Now that the opacities have been calculated for that column, integrate
                    self.tau[j, i] = np.trapz(self.data[:, j, i] * (self.effective_kappa[:, j, i] + self.effective_sigma[:, j, i])) * self.dz * self.unit_sigma
                    
        # Calculate the ratio of cells that are optically thick (tau >= 1)
        self.filling_factor = len(np.where(self.tau >= 1)[0]) / (self.Nx * self.Ny)

        return 

    def calc_t(self, rhod, effective_kappa, effective_sigma):
        """
        Optical depth with respect to z position along a single column.
        Integrates from z to L. This is the optical depth as emission progresses
        up the column, where as the optical_depth() function calculates
        along the entire column. 
        
        Args:
            rhod (ndarray): 1D array of dust densities, along the given column.
            effective_kappa (ndarray): 1D array of absorption opacities, along the given column corresponding to the input rhod.
            effective_sigma (ndarray): 1D array of scattering opacities, along the given column corresponding to the input rhod.

        Returns:
            1D array containing the cumulative optical depth along the column.
        """

        t = np.zeros(self.Nz) # To store the emission along the z-column

        # Integrate starting at the first cell of the column and move upward adding one cell at a time
        for i in range(self.Nz):  
            t[i] = np.trapz(rhod[:i] * (effective_kappa[:i] + effective_sigma[:i])) * self.dz * self.unit_sigma
            
        return t 

    def calc_flux(self):
        """
        Calculate outgoing flux using solution for RT eqn (5.113). Used only when polydisperse simulations are input!
    
        Returns:
            2D array containing the integrated values along the third axis.
        """

        if self.tau is None:
            raise ValueError('No optical depth map exists! Run the calc_tau() class method first!')
        
        # To store the outgoing flux at each (x,y) column
        self.flux = np.zeros((self.Ny, self.Nx))

        # Integrate each (x,y) column
        for i in range(self.Nx):
            for j in range(self.Ny):
                rhod = self.data[:, j, i] # The dust density in a particular column
                bb = self.src_fn[:, j, i] # The source function in a particular column
                kappa = self.effective_kappa[:, j, i] # The weighted absorption opacities in a particular column
                sigma = self.effective_sigma[:, j, i] # The weighted scattering opacities in a particular column (this is zero if include_scattering=False)
                
                # Mask where the particle density is zero along the column
                mask = (rhod == 0)

                # If density is zero then the source function and opacities should be zero as well
                bb[mask], kappa[mask], sigma[mask] = 0, 0, 0 

                # This is the optical depth as the emission progresses up the column 
                t = self.calc_t(rhod, kappa, sigma) 
                
                # Integrate to compute the output flux at a given (x,y) position
                self.flux[j, i] = np.trapz(bb * np.exp(-(self.tau[j, i] - t)), x=self.axis, dx=self.dx)
                #self.flux /= 2.0 # Does half the flux go up and half the flux go down? If so the calculated mass excess would double!
    
        return 
        
    def calc_mass_excess(self, frequency=3e11):
        """
        Calculates the mass_excess attribute.
        
        frequency (float): Wavelength requency at which to calculate the observed flux and thus mass excess.
            Input must be Hz, defaults to 3e11 corresponding to the 1mm frequency.

        Returns:
            Float.
        """

        # Calculate the optical depth map 
        self.calc_tau()

        if self.include_scattering:

            # The scattering solution of a thin slab as approximated by Miyake & Nakagawa (1993),
            # which has been used to compute the emergent intensity of protoplanetary disks including scattering.
            # This is applicable under the assumption that the disk temperature is constant
            # and that there are no incoming radiation fields at either the upper or lower disk surfaces.

            # Calculate the single scattering albedo
            if isinstance(self.grain_size, np.ndarray) is False:
                # Monodisperse
                albedo = self.sigma / (self.kappa + self.sigma)
            else:
                # Polydisperse -- the albedo and source function are now 3-dimensional!
                albedo = self.effective_sigma / (self.effective_kappa + self.effective_sigma)
                albedo[~np.isfinite(albedo)] = 0 #Replace all NaNs that come about when the denominator is zero

            # Similar format as Zhu. et al (2019) -- Section 2.1 (https://iopscience.iop.org/article/10.3847/2041-8213/ab1f8c/pdf)
            epsilon = 1.0 - albedo # For convinience 
            mu = 1.0 / np.sqrt(3.0) # The rays originate from the direction of cos(θ) = 1/sqrt(3) for all inclinations -- where θ is the angle between the intensity and the vertical direction
            tau_d = (2 * mu) / 3.0 # Total optical depth in the vertical direction? Or is this the specific depth according to the Eddington-Barbier relation?
            tau_ = 0.0 # Variable optical depth in the vertical direction? Or is this the optical depth at the surface of the slab, which is 0?

            # Same format as Eq. 8 of Zhu et al. (2019) -- (https://iopscience.iop.org/article/10.3847/2041-8213/ab1f8c/pdf)
            numerator = np.exp(-np.sqrt(3 * epsilon) * tau_) + np.exp(np.sqrt(3 * epsilon) * (tau_ - tau_d))
            denominator = (np.exp(-np.sqrt(3 * epsilon) * tau_d) * (1 - np.sqrt(epsilon))) + (np.sqrt(epsilon) + 1)
            J = self.blackbody(frequency=frequency) * (1 - (numerator / denominator))

            # With J known we can now solve for the source function (Eq. 6 of Zhu et al. (2019))
            self.src_fn = albedo * J + (1 - albedo) * self.blackbody(frequency=frequency)
        else:
            if isinstance(self.grain_size, np.ndarray) is False:
                # Monodisperse
                self.src_fn = self.blackbody(frequency=frequency)
            else:
                # Polydisperse -- make a 3D array in which the value in each cell is the source function
                # This will be used to integrate the RT equation (values where dust density is zero will be zeroed out)
                self.src_fn = np.zeros((self.Nz, self.Ny, self.Nx))
                self.src_fn[::] = self.blackbody(frequency=frequency)

       # Compute the mass underestimation 
        if isinstance(self.grain_size, np.ndarray) is False: # Monodisperse

            # Integrating the general RT equation with constant T and source function as well as I(0)=0 simplifies to the following
            self.flux = self.src_fn * (1 - np.exp(-self.tau))

            # Convolution theory -- take the mean of the output flux
            sigma_dust = np.mean(self.flux) / (self.src_fn * (self.kappa + self.sigma)) if self.include_scattering else np.mean(self.flux) / (self.src_fn * self.kappa)
        
        else: #Polydisperse -- in this case T is still constant but the source function is now a function of position
            
            # Need to integrate the general RT equation
            self.calc_flux() # Sets the flux attribute

            # The average opacity an observer would assume (assumes all grain sizes are equally distributed so scales with the inverse number of species)
            assumed_opacity = (np.sum(self.kappa) + np.sum(self.sigma)) / len(self.grain_size) if self.include_scattering else np.sum(self.kappa) / len(self.grain_size)
            
            # Under the assumption of optically thin emission, the observed flux scales with the column density of the dust, allowing us to analytically solve for Σd as 
            sigma_dust = np.mean(self.flux) / (np.mean(self.src_fn) * assumed_opacity)

        # The observed mass of the box can now be quantified as the product of Σd and the simulation area
        self.observed_mass = sigma_dust * self.area  

        # The mass underestimation, ratio of true box mass to the observed mass
        self.mass_excess = self.mass / self.observed_mass

        return 

    def calc_grain_size(self):
        """
        Calculates grain size given stokes number and gas column density
        """

        if isinstance(self.stoke, np.ndarray) is False:
            # Monodisperse
            self.grain_size = self.stoke * 2. * self.column_density / np.pi / self.grain_rho
        else:
            # Polydisperse
            if isinstance(self.grain_rho, np.ndarray) is False:
                raise ValueError("If entering multiple stoke's numbers, the corresponding grain_rho paramater must be a list/ndarray of same size!")
            
            # Calculate the grain sizes corresponding to each stokes number
            self.grain_size = np.zeros(len(self.stoke)) 

            for grain in range(len(self.grain_size)):
                self.grain_size[grain] = self.stoke[grain] * 2. * self.column_density / np.pi / self.grain_rho[grain]
                
        return

    def extract_opacity(self):
        """
        Returns the opacity coefficient(s) according to grain size.
        """

        try:
            self.calc_grain_size()
        except:
            raise ValueError('Could not determine grain size(s), input the stoke and grain_rho parameters and try again.')

        # Grain size, absorption opacity, and scattering opacity -- from DSHARP project
        a, k_abs, k_sca = load_opacity_values(q=self.q)

        # Interpolate 
        k_abs_fit, k_sca_fit = interpolate.interp1d(a, k_abs), interpolate.interp1d(a, k_sca)

        if isinstance(self.grain_size, np.ndarray) is False:
            # Monodisperse case
            if self.grain_size > a.max(): raise ValueError(f"Maximum grain size supported is {str(a.max())} cm.")
            if self.grain_size < a.min(): raise ValueError(f"Minimum grain size supported is {str(a.min())} cm.")

            # Extract the opacity coefficients for the grain size
            self.kappa = k_abs_fit(self.grain_size)
            self.sigma = k_sca_fit(self.grain_size) if self.include_scattering else None

        else:
            # Polydisperse case
            if self.grain_size.max() > a.max(): raise ValueError(f"Maximum grain size supported is {str(a.max())} cm.")
            if self.grain_size.min() < a.min(): raise ValueError(f"Minimum grain size supported is {str(a.min())} cm.")

            # Extract the opacity coefficients for all the grain sizes
            self.kappa, self.sigma = np.zeros(len(self.grain_size)), np.zeros(len(self.grain_size))

            for grain in range(len(self.grain_size)):
                self.kappa[grain] = k_abs_fit(self.grain_size[grain])
                self.sigma[grain] = k_sca_fit(self.grain_size[grain]) if self.include_scattering else 0

        return 

    def get_proto_mass(self):
        """
        Calculates the mass of the protoplanets.

        Returns:
            Mass of the planetesimals, stored in the proto_mass attribute. If no planetesimals are present the value is set to zero.
        """

        mp_code = self.eps_dtog * self.mass / self.npar
        mp = stats.mode(self.rhopswarm)[0]

        index = np.where(self.aps != 0)[0]
        npclump = self.rhopswarm[index] / mp

        tmp = np.log10(npclump)
        ttmp = tmp + np.log10(mp_code)
        mass = 10**ttmp

        self.proto_mass = np.sort(mass)
        
        # If there are no planetesimals set to 0
        try: 
            self.proto_mass.max()
        except ValueError:
            self.proto_mass = 0

        return


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

    Note:
        These coefficients correspond to the 1mm frequency!

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

    return a, k_abs, k_sca

