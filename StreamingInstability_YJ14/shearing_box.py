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

class density_cube:
    """
    Class for a density cube object. The class methods can be used to perform radiative transfer calculations 
    along the z-axis of the cube, enabling analysis given a range of physical disk conditions. 
    If no dust opacity coefficients are input (kappa and/or sigma), the opacities will be estimated 
    using coefficients from the DSHARP project (see: https://iopscience.iop.org/article/10.3847/2041-8213/aaf743/pdf)
    These DSHARP opacities are the 1mm wave opacities, therefore do not modify the frequency of analysis if
    these opacities are used (3e11 Hz corresponds to the 1mm wavelength)

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
        code_omega (float): The Keplerian orbital timescale in code units (Omega param the start.in file, &hydro_init_pars). Defaults to 2pi.
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
        xgrid (ndarray): x-coordinates of grid nodes. Saved as the x attribute in read_grid().
        ygrid (ndarray): y-coordinates of grid nodes. Saved as the y attribute in read_grid().
        zgrid (ndarray): z-coordinates of grid nodes. Saved as the z attribute in read_grid().
    """
    
    def __init__(self, data=None, axis=None, column_density=100, T=30, H=5*const.au.cgs.value, kappa=None, sigma=None, q=None,
        stoke=0.3, grain_rho=1.0, eps_dtog=0.03, npar=1e6, code_rho=1, code_cs=2*np.pi, code_omega=2*np.pi, aps=None, rhopswarm=None, init_var=None, 
        include_scattering=False, ipars=None, xp=None, yp=None, zp=None, xgrid=None, ygrid=None, zgrid=None):

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
        self.xgrid = xgrid
        self.ygrid = ygrid
        self.zgrid = zgrid

        if self.include_scattering:
            if self.kappa is not None:
                if self.sigma is None: raise ValueError('The include_scattering paramater has been enabled but no scattering coefficient (sigma) was input!')
            if self.sigma is not None:
                if self.kappa is None: raise ValueError('The include_scattering paramater has been enabled but no absorption coefficient (kappa) was input!')

        if self.data is None:
            print('No data input, loading density cube from YJ14, orbit 100...')
            self.data, self.axis = load_cube()

        self.tau = None 
        self.intensity = None
        self.mass = None
        self.mass_excess = None
        self.filling_factor = None 
        self.area = None 
        self.grain_size = None 
        self.proto_mass = None 

    def configure(self, frequency=3e11):
        """
        Initializing parameters and creates intensity, mass excess, and filling factor 
        attributes. If sigma, H, or dx/dy/dz attributes are updated, re-run this method 
        to re-configure the settings.

        Note:
            If running this method multiple times with different astrocentric parameters,
            ensure that the kappa and sigma parameters are reset to None if they were None 
            initially! Failure to do so will result in the grain size not updating and thus 
            the opacity coefficients will not be properly set. Likewise if include_scattering is
            changed, the presence of assigned kappa and sigma values will throw off the configuration. 
            It is recommended to re-initialize the class if re-running with varying parameters.

        Args:
            frequency (float): Wavelength requency at which to calculate the observed intensity and thus mass excess.
                Input must be Hz, defaults to 3e11 corresponding to the 1mm frequency.
        """

        # Dimensions of the box
        self.Lz = self.Ly = self.Lx = np.abs(self.axis[0] - self.axis[-1]) * self.H # Box length (assumes a cube!)
        self.dz = self.dy = self.dx = np.diff(self.axis)[0] # Cell length (axis is in units of scale height!)
        self.Nz = self.Ny = self.Nx = len(self.axis) # No. of cells 
        self.area = self.Lx * self.Ly # Box area in cgs code units

        # Mass of the box in code units
        box_mass_codeunits = np.sum(self.data) if self.init_var is None else np.sum(self.init_var) # init_var is the 0th snapshot and it's used for simulations with self-gravity (the initial mass)
        box_mass_codeunits = box_mass_codeunits * self.dx * self.dy * self.dz 
        
        # Convert the mass to cgs units
        self.unit_mass = (self.column_density * self.H**2) / np.sqrt(2 * np.pi) / (self.code_rho * (self.code_cs / self.code_omega)**3)
        self.mass = box_mass_codeunits * self.unit_mass 

        # To convert the dust surface density to cgs units (used when integrating the RT solution and to calculate tau)
        self.unit_sigma = self.column_density / np.sqrt(2 * np.pi) / (self.code_rho * (self.code_cs / self.code_omega))

        # Convert the density cube to cgs units
        self.unit_length = self.H
        self.unit_density = self.unit_mass / self.unit_length**3
        self.data *= self.unit_density

        self.dz *= self.unit_length # Convert to cgs units (for the vertical integrations)

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
            frequency (float): Fequency of the wavelength at which to calculate the observed intensity and thus the mass excess.
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
        
        # If polydisperse, create a 'density_per_species' 4D array to store the density of dust grains in each 3D cell on a per-species basis, used to calculate weighted opacities
        if isinstance(self.grain_size, np.ndarray):

            # To store the weighted opacity coefficients at each cell
            self.effective_kappa = np.zeros((self.Nz, self.Ny, self.Nx))
            self.effective_sigma = np.zeros((self.Nz, self.Ny, self.Nx))
            
            # Convert the ipars array to numerical labels, first species is 1, second is 2, etc...
            species = np.ceil(self.ipars / (self.npar / len(self.grain_size))).astype(int) 
            num_species = len(np.unique(species))
            print(f"Detected {num_species} grain sizes, applying weighted opacity calculation...")

            # To store the density of each species located at each cell, shape (num_of_species, len(Nz), len(Ny), len(Nx))
            self.density_per_species = np.zeros((num_species, self.Nz, self.Ny, self.Nx))

            for species_type in range(1, num_species+1):

                print(f"Converting grain positions to density field: {species_type} out of {num_species}")
                
                # Index the particular species
                index_species = np.where(species == species_type)[0]
                
                # Extract the (x,y,z) positions for all of these grains
                species_x, species_y, species_z = self.xp[index_species], self.yp[index_species], self.zp[index_species]
                
                # Convert positions of particles to a grid density field
                particle_density = particles_to_density(species_x, species_y, species_z, self.xgrid, self.ygrid, self.zgrid)

                # Update the density array
                self.density_per_species[species_type - 1] = particle_density

        ###
        ### Calculate the optical depth: tau = kappa * rhop * dz
        ###
        for i in range(self.Nx):
            for j in range(self.Ny):

                # This is the Monodisperse case, opacity coefficients are single values (self.kappa and self.sigma)
                if isinstance(self.grain_size, np.ndarray) is False:
                    surface_density = np.trapz(self.data[:, j, i]) * self.dz # * self.unit_sigma
                    self.tau[j, i] = surface_density * (self.kappa + self.sigma) if self.include_scattering else surface_density * self.kappa

                else:
                    # This is the polydisperse case in which the cell-wise averaged opacities must be calculated: (N1*κ1 + N2*κ2 + N3*κ3 + N4*κ4 + ...) / N
                    # The below code goes through all the z-cells in the particular (x,y) column and stores the opacities in self.effective_kappa and self.effective_sigma
                    for k in range(self.Nz):

                        # To count the opacities at each (x,y,z) cell (equivalent to N*k)
                        weighted_kappa, weighted_sigma = 0, 0

                        # Add the number of grains in that cell for a given species times that species' opacity coefficient (stored in self.kappa/self.sigma which in the polydisperse case are arrays)
                        for species_type in range(num_species):
                            weighted_kappa +=  self.density_per_species[species_type][k, j, i] * self.kappa[species_type] # This is the numerator (N1*k1 + N2*k2 + ...)
                            weighted_sigma = weighted_sigma + self.density_per_species[species_type][k, j, i] * self.sigma[species_type] if self.include_scattering else 0

                        # Divide by the total number of species in that (x,y,z) cell to compute the weighted mean and to the 3D opacity array
                        self.effective_kappa[k, j, i] = weighted_kappa / np.sum(self.density_per_species[:, k, j, i]) if np.sum(self.density_per_species[:, k, j, i]) != 0 else 0 # Avoid division by zero 
                        if self.include_scattering:
                            self.effective_sigma[k, j, i] = weighted_sigma / np.sum(self.density_per_species[:, k, j, i]) if np.sum(self.density_per_species[:, k, j, i]) != 0 else 0 # Avoid division by zero 
                   
                    # Now that the opacities have been calculated for that column, vertically integrate to find the opacity in that column
                    self.tau[j, i] = np.trapz(self.data[:, j, i] * (self.effective_kappa[:, j, i] + self.effective_sigma[:, j, i])) * self.dz # * self.unit_sigma
                    
        # Calculate the ratio of cells that are optically thick (tau >= 1)
        self.filling_factor = len(np.where(self.tau >= 1)[0]) / (self.Nx * self.Ny)

        return 

    def calc_t(self, rhod, effective_kappa, effective_sigma):
        """
        Optical depth with respect to z position along a single column.
        Integrates from 0 to z. This is the optical depth as emission progresses
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
            t[i] = np.trapz(rhod[:i] * (effective_kappa[:i] + effective_sigma[:i])) * self.dz
            
        return t 

    def calc_intensity(self):
        """
        Calculate outgoing intensity using solution for RT eqn (5.113). Used only when polydisperse simulations are input!
        
        Note:
            The 1st term of the general RT solution is the extinction of the original intensity, and the 2nd term is the emission at a 
            point t, extinguished in the path from t to tau. This implementation assumes that the first
            term is zero (no back-illumination).
        Returns:
            2D array containing the integrated values along the third axis.
        """

        if self.tau is None:
            raise ValueError('No optical depth map exists! Run the calc_tau() class method first!')
        
        # To store the outgoing intensity at each (x,y) column
        self.intensity = np.zeros((self.Ny, self.Nx))

        # Integrate each (x,y) column
        for i in range(self.Nx):

            print(f"Intensity calculation: {i+1} out of {self.Nx}")

            for j in range(self.Ny):

                rhod = self.data[:, j, i] # The dust density in a particular column
                #bb = self.src_fn[:, j, i] # The source function in a particular column
                kappa = self.effective_kappa[:, j, i] # The weighted absorption opacities in a particular column
                sigma = self.effective_sigma[:, j, i] # The weighted scattering opacities in a particular column (this is zero if include_scattering=False)

                # Mask where the particle density is zero along the column
                mask = (rhod == 0)

                # If density is zero then the source function and opacities should be zero as well
                kappa[mask], sigma[mask] = 0, 0

                # The directional average of the intensity from scattering solution (3D if polydisperse)
                if self.include_scattering:
                    if isinstance(self.grain_size, np.ndarray): #Polydisperse
                        J = self.J[:, j, i]
                        J[mask] = 0
                    else: #Monodisperse
                        J = self.J
                
                # This is the optical depth as the emission progresses up the column (0 to z integral)
                t = self.calc_t(rhod, kappa, sigma) 
                
                # The emissivity 
                emissivity = (kappa * self.B_nu) + (sigma * J) if self.include_scattering else kappa * self.B_nu

                # Integrate to compute the output intensity at a given (x,y) position
                #self.intensity[j, i] = np.trapz(bb * np.exp(-(self.tau[j, i] - t)), x=self.axis, dx=self.dx)
                self.intensity[j, i] = np.trapz(emissivity * rhod * np.exp(-(self.tau[j, i] - t))) * self.dz
                
        return 
        
    def calc_mass_excess(self, frequency=3e11):
        """
        Calculates the mass_excess attribute.
        
        frequency (float): Wavelength requency at which to calculate the observed intensity and thus mass excess.
            Input must be Hz, defaults to 3e11 corresponding to the 1mm frequency.

        Returns:
            Float.
        """

        # Calculate the optical depth map 
        self.calc_tau()

        # Compute the radiation assuming a Planckian black body
        self.B_nu = self.blackbody(frequency=frequency)

        ###
        ### Calculate the effective source function (which is just the Planckian if scattering is not accounted for)
        ###

        if self.include_scattering:

            # The scattering solution of a thin slab as approximated by Miyake & Nakagawa (1993),
            # which has been used to compute the emergent intensity of protoplanetary disks including scattering.
            # This is applicable under the assumption that the disk temperature is constant
            # and that there are no incoming radiation fields at either the upper or lower disk surfaces.

            # Calculate the single scattering albedo
            if isinstance(self.grain_size, np.ndarray) is False:
                # Monodisperse
                self.albedo = self.sigma / (self.kappa + self.sigma)
            else:
                # Polydisperse -- the albedo and source function are now 3-dimensional!
                self.albedo = self.effective_sigma / (self.effective_kappa + self.effective_sigma)
                self.albedo[~np.isfinite(self.albedo)] = 0 #Replace all NaNs that come about when the denominator is zero

            # Similar format as Zhu. et al (2019) -- Section 2.1 (https://iopscience.iop.org/article/10.3847/2041-8213/ab1f8c/pdf)
            epsilon = 1.0 - self.albedo # For convenience 
            mu = 1.0 / np.sqrt(3.0) # The rays originate from the direction of cos(θ) = 1/sqrt(3) for all inclinations -- where θ is the angle between the intensity and the vertical direction
            tau_d = (2 * mu) / 3.0 # Total optical depth in the vertical direction? Or is this the specific depth according to the Eddington-Barbier relation?
            tau_ = 0.0 # Variable optical depth in the vertical direction? Or is this the optical depth at the surface of the slab, which is 0?

            # Same format as Eq. 8 of Zhu et al. (2019) -- (https://iopscience.iop.org/article/10.3847/2041-8213/ab1f8c/pdf)
            numerator = np.exp(-np.sqrt(3 * epsilon) * tau_) + np.exp(np.sqrt(3 * epsilon) * (tau_ - tau_d))
            denominator = (np.exp(-np.sqrt(3 * epsilon) * tau_d) * (1 - np.sqrt(epsilon))) + (np.sqrt(epsilon) + 1)
            self.J = self.B_nu * (1 - (numerator / denominator))

            # With J known we can now solve for the source function (Eq. 6 of Zhu et al. (2019))
            if isinstance(self.grain_size, np.ndarray) is False: # Monodisperse case needs a cube for the RT integration
                self.src_fn = np.zeros((self.Nz, self.Ny, self.Nx))
                self.src_fn[::] = self.albedo * self.J + (1 - self.albedo) * self.B_nu
            else: # Polydisperse case is a cube already since albedo is 3D 
                self.src_fn = self.albedo * self.J + (1 - self.albedo) * self.B_nu
            
        else: #Absorption only case, source function is just the Planckian black body
            if isinstance(self.grain_size, np.ndarray) is False:
                # Monodisperse
                self.src_fn = self.B_nu
            else:
                # Polydisperse -- make a 3D array in which the value in each cell is the source function
                # This will be used to integrate the RT solution (values where dust density is zero will be zeroed out)
                self.src_fn = np.zeros((self.Nz, self.Ny, self.Nx))
                self.src_fn[::] = self.B_nu

        ###
        ### Compute the mass underestimation 
        ###

        # Monodisperse
        if isinstance(self.grain_size, np.ndarray) is False: 

            # In the absorption only case the, source function is independent of optical depth, therefore can use the following simplified form
            if self.include_scattering is False:
                
                # Integrating the general RT solution with constant T and source function as well as I(0)=0 simplifies to the following
                #self.intensity = self.B_nu * (1 - np.exp(-self.tau))

                # Set the opactities -- the opacity at every cell is the same (cells with zero dust density will be zeroed out during intensity calculation)
                self.effective_kappa, self.effective_sigma = np.zeros((self.Nz, self.Ny, self.Nx)), np.zeros((self.Nz, self.Ny, self.Nx))
                self.effective_kappa[::], self.effective_sigma[::] = self.kappa, 0 
                
                # Solve the RT solution to set the self.intensity parameter
                self.calc_intensity()

                # Under the assumption of optically thin emission, the observed intensity scales with the column density of the dust, allowing us to analytically solve for Σd as 
                self.sigma_dust = np.mean(self.intensity) / (self.B_nu * self.kappa) # Convolution theory -- take the mean of the output intensity
            
            # If including photon scattering we need to solve the general RT solution since the effective source function is now dependent on position
            else:
                
                # Set the opactities -- the opacity at every cell is the same (cells with zero dust density will be zeroed out during intensity calculation)
                self.effective_kappa, self.effective_sigma = np.zeros((self.Nz, self.Ny, self.Nx)), np.zeros((self.Nz, self.Ny, self.Nx))
                self.effective_kappa[::], self.effective_sigma[::] = self.kappa, self.sigma 

                # Solve the RT solution to set the self.intensity parameter
                self.calc_intensity()

                # Under the assumption of optically thin emission, the observed intensity scales with the column density of the dust, allowing us to analytically solve for Σd as 
                self.sigma_dust = np.mean(self.intensity) / (self.B_nu * (self.kappa + self.sigma)) # Convolution theory -- take the mean of the output intensity
        
        #Polydisperse
        else: 

            # Solve the RT solution to set the self.intensity parameter
            self.calc_intensity() 

            # The average opacity an observer would assume (assumes all grain sizes are equally distributed so scales with the inverse number of species)
            self.assumed_opacity = (np.sum(self.kappa) + np.sum(self.sigma)) / len(self.grain_size) if self.include_scattering else np.sum(self.kappa) / len(self.grain_size)
            
            # Under the assumption of optically thin emission, the observed intensity scales with the column density of the dust, allowing us to analytically solve for Σd as 
            self.sigma_dust = np.mean(self.intensity) / (self.B_nu * self.assumed_opacity) # Convolution theory -- take the mean of the output intensity and the source function

        # The observed mass of the box can now be quantified as the product of Σd and the domain area
        self.observed_mass = self.sigma_dust * self.area  

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
        These coefficients correspond to the 1mm wave opacities!

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

def particles_to_density(xp, yp, zp, x, y, z, rhop_swarm=6.30813659928, grid_func1='linear', grid_func2='linear', grid_func3='linear',
    mx=262, my=262, mz=262, nx=256, ny=256, nz=256, n1=3, n2=258, m1=3, m2=258, l1=3, l2=258, density=True):
    """
    Convert particle positions and weights to a density field on a grid.
    Author: Wladimir Lyra (adapted from Anders Johansen's IDL script) 

    Parameters:
        xp (array-like): x-coordinates of particles. Saved as the xp attribute in the pvar files.
        yp (array-like): y-coordinates of particles. Saved as the yp attribute in the pvar files.
        zp (array-like): z-coordinates of particles. Saved as the zp attribute in the pvar files.
        x (array-like): x-coordinates of grid nodes. Saved as the x attribute in read_grid().
        y (array-like): y-coordinates of grid nodes. Saved as the y attribute in read_grid().
        z (array-like): z-coordinates of grid nodes. Saved as the z attribute in read_grid().
        rhop_swarm (array-like or float): Density or weight of each particle, saved as the rhop_swarm attribute in read_param().
        grid_func1 (str, optional): Interpolation method for x-axis. Default is 'linear'. Saved as the grid_func[0] attribute in read_param().
        grid_func2 (str, optional): Interpolation method for y-axis. Default is 'linear'. Saved as the grid_func[1] attribute in read_param().
        grid_func3 (str, optional): Interpolation method for z-axis. Default is 'linear'. Saved as the grid_func[2] attribute in read_param().
        mx (int): Number of grid nodes in x-direction. Saved as the mx attribute in read_dim().
        my (int): Number of grid nodes in y-direction. Saved as the my attribute in read_dim().
        mz (int): Number of grid nodes in z-direction. Saved as the mz attribute in read_dim().
        nx (int): Number of interpolation points in x-direction. Saved as the nx attribute in read_dim().
        ny (int): Number of interpolation points in y-direction. Saved as the ny attribute in read_dim().
        nz (int): Number of interpolation points in z-direction. Saved as the nz attribute in read_dim().
        n1 (int): Lower index of the grid in z-direction. Saved as the n1 attribute in read_dim().
        n2 (int): Upper index of the grid in z-direction. Saved as the n2 attribute in read_dim().
        m1 (int): Lower index of the grid in y-direction. Saved as the m1 attribute in read_dim().
        m2 (int): Upper index of the grid in y-direction. Saved as the m2 attribute in read_dim().
        l1 (int): Lower index of the grid in x-direction. Saved as the l1 attribute in read_dim().
        l2 (int): Upper index of the grid in x-direction. Saved as the l2 attribute in read_dim().
        density (bool, optional): If True, compute density; if False, compute weight. Default is True.

    Returns:
        array-like: Density field on the specified grid.

    """

    dx, dy, dz = np.gradient(x), np.gradient(y), np.gradient(z) 
    dx1, dy1, dz1 = 1.0/dx, 1.0/dy, 1.0/dz
    dx2, dy2, dz2 = 1.0/dx**2, 1.0/dy**2, 1.0/dz**2

    if grid_func1 == 'linear': dx1_pt = dx1[0]
    if grid_func2 == 'linear': dy1_pt = dy1[0]
    if grid_func3 == 'linear': dz1_pt = dz1[0]

    nnp = np.zeros((mz, my, mx))

    for k in range(len(xp)):

        _xp_, _yp_, _zp_ = xp[k], yp[k], zp[k]

        ix0 = int(round((_xp_ - x[0]) * dx1_pt)) if grid_func1 == 'linear' else find_index_bisect(_xp_, x)
        if ix0 == l2 + 1: ix0 = ix0 - 1
        if ix0 == l1 - 1: ix0 = ix0 + 1
        dx_1, dx_2 = dx1[ix0], dx2[ix0]

        iy0 = int(round((_yp_ - y[0]) * dy1_pt)) if grid_func2 == 'linear' else find_index_bisect(_yp_, y)
        if iy0 == m2 + 1: iy0 = iy0 - 1
        if iy0 == m1 - 1: iy0 = iy0 + 1
        dy_1, dy_2 = dy1[iy0], dy2[iy0]

        iz0 = int(round((_zp_ - z[0]) * dz1_pt)) if grid_func3 == 'linear' else find_index_bisect(_zp_, z)
        if iz0 == n2 + 1: iz0 = iz0 - 1
        if iz0 == n1 - 1: iz0 = iz0 + 1
        dz_1, dz_2 = dz1[iz0], dz2[iz0]

        ixx0, ixx1 = ix0 - 1, ix0 + 1
        iyy0, iyy1 = iy0 - 1, iy0 + 1
        izz0, izz1 = iz0 - 1, iz0 + 1

        for ixx in np.arange(ixx0, ixx1 + 1):
            for iyy in np.arange(iyy0, iyy1 + 1):
                for izz in np.arange(izz0, izz1 + 1):

                    if ( ((ixx - ix0) == -1) or ((ixx - ix0) == +1) ):
                        weight_x = 1.125 - 1.5 * abs(_xp_ - x[ixx]) * dx_1 + 0.5 * abs(_xp_ - x[ixx])**2 * dx_2
                    else:
                        if nx != 1: weight_x = 0.75 - (_xp_ - x[ixx])**2 * dx_2

                    if ( ((iyy - iy0) == -1) or ((iyy - iy0) == +1) ):
                        weight_y = 1.125 - 1.5 * abs(_yp_ - y[iyy]) * dy_1 + 0.5 * abs(_yp_ - y[iyy])**2 * dy_2
                    else:
                        if ny != 1: weight_y = 0.75 - (_yp_ - y[iyy])**2 * dy_2

                    if ( ((izz - iz0) == -1) or ((izz - iz0) == +1) ):
                        weight_z = 1.125 - 1.5 * abs(_zp_ - z[izz]) * dz_1 + 0.5 * abs(_zp_ - z[izz])**2 * dz_2
                    else:
                        if nz != 1: weight_z = 0.75 - (_zp_ - z[izz])**2 * dz_2

                    if density:
                        weight = rhop_swarm if type(rhop_swarm) is float else rhop_swarm[k]
                    else:
                        weight = 1.0

                    if nx != 1: weight = weight * weight_x
                    if ny != 1: weight = weight * weight_y
                    if nz != 1: weight = weight * weight_z

                    nnp[izz, iyy, ixx] = nnp[izz, iyy, ixx] + weight

    return nnp[n1:n2 + 1, m1:m2 + 1, l1:l2 + 1]

def find_index_bisect(qpar, q):
    """
    Find the index of the element in the sorted list q that is closest to the given qpar using binary search.
    Author: Wladimir Lyra (adapted from Anders Johansen's IDL script) 

    Parameters:
        qpar (float): The value to search for.
        q (list of float): A sorted list of values to search within.

    Returns:
        int: The index of the element in the list q that is closest to qpar.
    """

    jl, ju = 0, len(q) - 1
    
    while (ju - jl) > 1:
        jm = (ju + jl) // 2

        if qpar > q[jm]:
            jl = jm
        else:
            ju = jm

    iq0 = jl if (qpar - q[jl] <= q[ju] - qpar) else ju

    return iq0
