#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:01:11 2021

@author: daniel
"""
from typing import Union, Optional
import warnings; warnings.filterwarnings("ignore")
from scipy import interpolate, stats
import astropy.constants as const
from pathlib import Path
import pkg_resources
import numpy as np

from StreamingInstability_YJ14 import compute_opacities

class density_cube:
    """
    Class for a density cube object. This class performs radiative transfer calculations along the z-axis of a 3D cube,
    enabling the analysis of protoplanetary disk properties under a range of physical conditions. If no dust opacity
    coefficients (kappa and/or sigma) are provided, the opacities will be estimated using the DSHARP project coefficients 
    (Birnstiel et al., 2018; https://iopscience.iop.org/article/10.3847/2041-8213/aaf743/pdf). 

    Note:
        This analysis assumes a cubic domain (symmetrical along each axis).

    Args:
        data (ndarray): Simulation Data. The 3D density data cube, containing dust particle density (rhop attribute). Defaults to None in 
            which case the particle densities for a single snapshop of a streaming instability simulation (Yang+14) will 
            be loaded, for testing purposes.
        axis (ndarray): Simulation Data. The 1D array of the axis along which to integrate (x, y, or z attribute -- cubic domain is assumed).
            This should be in units of gas scale height, the default when running SI simulations using the Pencil Code.
            Defaults to None in which case the z-axis from pre-saved model will be loaded, for testing purposes.
        code_rho (float): Data Normalization. Midplane gas density in code units (rho0 param in the start.in file, under &eos_init_pars), used to convert 
            data cube to cgs units. Defaults to 1.
        code_cs (float): Data Normalization. The sound speed of the gas in code units (cs0 param in the start.in file, under &eos_init_pars), used to convert 
            data cube to cgs units. Defaults to 2pi.
        code_omega (float): Data Normalization. The Keplerian orbital timescale in code units (Omega param the start.in file, under &hydro_init_pars), used to 
            convert data cube to cgs units. Defaults to 2pi.
        column_density (float): Disk Parameter. Column density of the gas in cgs units. Defaults to 100 (g / cm^2).
        T (float): Disk Parameter. Temperature of the entire box in cgs units, as we assume isothermality. Defaults to 30 (K).
        H (float): Disk Parameter. The scale height of the box, in cgs units. Defaults to 5 (au) ~ 7.5e13 (cm)
        stoke (float or ndarray): Dust Grain Parameter. Stokes number of the dust species present in the model. Must either be a float for 
            monodisperse simulations or an array containing one value for each species if the model is polydisperse.
        grain_rho (float): Dust Grain Parameter. Internal dust grain density, in cgs units. A value of 1 (g / cm^3) is typical for ices, 
            and 3.5 (g / cm^3) for silicates. This input can either be a float for monodisperse simulations or 
            an array containing one value for each species in the case of polydisperse models. Will be used to compute 
            the grain size according to the input Stokes number(s) gas column density. Defaults to None, which assumes
            a grain density of 1.675 (g / cm^3) for all grains, the density used in the DSHARP dust model.
        wavelength (float): Radiative Transfer Parameter. The wavelength at which to calculate the opacities, in cm. Value can range from 1e-5 to 10 cm.
            Defaults to 0.1 cm.        
        include_scattering (bool): Radiative Transfer Parameter. Whether the scattering opacities should be considered in the radiative transfer. Defaults to False,
            which calculates (if applicable) and employs only the absorption opacities. 
        kappa (float): Radiative Transfer Parameter. Dust absorption opacity coefficient in cgs units (cm^2 / g). Must either be a float for monodisperse 
            simulations or an array containing one value for each species if the model is polydisperse. If set to None, the
            dust absorption opacity will be calculated according to column density and grain size(s), using the DSHARP dust model
            from Birnstiel+18. Defaults to None.
        sigma (float): Radiative Transfer Parameter. Dust scattering opacity coefficient, in cgs units (cm^2 / g). Must either be a float for monodisperse 
            simulations or an array containing one value for each species if the model is polydisperse. If set to None, the
            dust scattering opacity will be calculated according to column density and grain size(s), using the DSHARP dust model
            from Birnstiel+18. Defaults to None.
        q (float): Radiative Transfer Parameter. The power law of the grain size distribution, n(a) scales with a^(-q) where a is the grain size. Defaults to 2.5 which
            is typical for protoplanetary disk environments (~ Class I/II).
        npar (int): Polydisperse Radiative Transfer Parameter. Number of particles in the simulation, used to calculate the density-weighted opacities if the model is polydisperse. 
            It is also used to calculate the mass of protoplanets if the model is self-gravitating (optional). Defaults to one million.
        ipars (ndarray, optional): Polydisperse Radiative Transfer Parameter. The grain identifier array stored in the pvar files as the ipars attribute. Should be input as np.array(fp.ipars, dtype=int).
            Defaults to None. It's only used to compute the weighted opacities when the simulation is polydisperse.
        xp (ndarray, optional): Polydisperse Radiative Transfer Parameter. The x-position of the grains corresponding to the ipars array, saved as the xp attribute in the pvar files.
            Defaults to None. It's only used to compute the weighted opacities when the simulation is polydisperse.
        yp (ndarray, optional): Polydisperse Radiative Transfer Parameter. The x-position of the grains corresponding to the ipars array, saved as the xp attribute in the pvar files.
            Defaults to None. It's only used to compute the weighted opacities when the simulation is polydisperse.
        zp (ndarray, optional): Polydisperse Radiative Transfer Parameter. The z-position of the grains corresponding to the ipars array, saved as the zp attribute in the pvar files.
            Defaults to None. It's only used to compute the weighted opacities when the simulation is polydisperse.
        xgrid (ndarray, optional): Polydisperse Radiative Transfer Parameter. The x-coordinates of grid nodes, used to compute the weighted opacities when the simulation is polydisperse. 
            Saved as the x attribute in read_grid(). Defaults to None.
        ygrid (ndarray, optional): Polydisperse Radiative Transfer Parameter. The y-coordinates of grid nodes, used to compute the weighted opacities when the simulation is polydisperse. 
            Saved as the y attribute in read_grid(). Defaults to None.
        zgrid (ndarray, optional): Polydisperse Radiative Transfer Parameter. The z-coordinates of grid nodes, used to compute the weighted opacities when the simulation is polydisperse. 
            Saved as the z attribute in read_grid(). Defaults to None.
        aps (ndarray, optional): Self-Gravity Parameter. The azimuthal particle positions (pvar.aps) from the Pencil Code simulation. Used for calculating the mass 
            distribution of protoplanets by integrating over azimuthal regions. Must be input to enable protomass calculation. Defaults to None.
        rhopswarm (ndarray, optional): Self-Gravity Parameter. The swarm density (pvar.rhopswarm) from the Pencil Code simulation. Must be input to enable protomass calculation. 
            Represents the local density of particle swarms, used to compute the total mass of dense clumps in the simulation domain. Defaults to None.
        eps_dtog (float, optional): Self-Gravity Parameter. Dust-to-gas ratio, used to calculate the mass of the protoplanets if the model is self-gravitating. Defaults to None.
        init_var (ndarray, optional): Self-Gravity Parameter. The 3D density cube of the first (initial) snapshot so that the initial mass of the cube is known and used
            for the mass excess calculations. This is only applicable if the simulation includes self-gravity as the dust mass inside the box will decrease
            as sink particles are formed and accrete. Defaults to None.
    
    Attributes:
        mass (float): Total mass of the box in cgs units.
        mass_excess (float): Ratio of true box mass to observed mass (underestimation factor).
        intensity (ndarray): Outgoing intensity integrated along the z-axis.
        tau (ndarray): Optical depth map integrated along the z-axis.
        filling_factor (float): Fraction of columns that are optically thick (self.tau > 1).
        grain_size (float or ndarray): Calculated grain size(s) based on Stokes number and grain density.
        proto_mass (float): Mass of protoplanets in the simulation, if self-gravity is enabled.
    
    Methods:
        configure():
            Initializes key parameters and attributes, including density normalization, box mass, and opacity calculations.
            Re-run this method if any disk or radiative transfer parameters are updated.
        blackbody():
            Calculates the blackbody spectral radiance of the source using Planck's law, given the cube's temperature.
        calc_tau():
            Computes the optical depth (self.tau) along the z-axis for all columns of the cube.
        calc_t(rhod, effective_kappa, effective_sigma):
            Computes the cumulative optical depth along a single column, integrating from the base to the given height, z.
        calc_intensity():
            Solves the radiative transfer equation to compute the outgoing intensity along the z-axis for all columns.
        calc_mass_excess():
            Calculates the ratio of the true box mass to the observed mass (mass excess), based on the computed intensity.
        calc_grain_size():
            Determines grain size(s) based on the Stokes number and gas column density.
        extract_opacity():
            Retrieves dust opacity coefficients (kappa, sigma) for the given grain size(s), based on the DSHARP model.
        get_proto_mass():
            Estimates the mass of protoplanets in self-gravitating simulations, based on the particle positions and densities.

    """

    # Type hints
    def __init__(
        self,
        data: Optional[np.ndarray] = None,
        axis: Optional[np.ndarray] = None,
        code_rho: float = 1.0,
        code_cs: float = 2 * np.pi,
        code_omega: float = 2 * np.pi,
        column_density: float = 100.0,
        T: float = 30.0,
        H: float = 5 * const.au.cgs.value,
        stoke: Union[float, np.ndarray] = 0.3,
        grain_rho: Union[float, np.ndarray] = 1.675,
        wavelength: float = 0.1,
        include_scattering: bool = False,
        kappa: Optional[Union[float, np.ndarray]] = None,
        sigma: Optional[Union[float, np.ndarray]] = None,
        q: float = 2.5,
        npar: int = int(1e6),
        ipars: Optional[np.ndarray] = None,
        xp: Optional[np.ndarray] = None,
        yp: Optional[np.ndarray] = None,
        zp: Optional[np.ndarray] = None,
        xgrid: Optional[np.ndarray] = None,
        ygrid: Optional[np.ndarray] = None,
        zgrid: Optional[np.ndarray] = None,
        aps: Optional[np.ndarray] = None,
        rhopswarm: Optional[np.ndarray] = None,
        eps_dtog: Optional[float] = None,
        init_var: Optional[np.ndarray] = None,
    ) -> None:

        """
        Initialize
        """

        # Simulation Data
        self.data: Optional[np.ndarray] = data
        self.axis: Optional[np.ndarray] = axis

        # Disk Parameters
        self.column_density: float = column_density
        self.T: float = T
        self.H: float = H

        # Dust Grain Parameters
        self.stoke: Union[float, np.ndarray] = stoke
        self.grain_rho: Union[float, np.ndarray] = grain_rho

        # Radiative Transfer Settings
        self.wavelength: float = wavelength
        self.include_scattering: bool = include_scattering
        self.kappa: Optional[Union[float, np.ndarray]] = kappa
        self.sigma: Optional[Union[float, np.ndarray]] = sigma
        self.q: float = q

        # Polydisperse Simulation Parameters
        self.npar: int = npar
        self.ipars: Optional[np.ndarray] = ipars
        self.xp: Optional[np.ndarray] = xp
        self.yp: Optional[np.ndarray] = yp
        self.zp: Optional[np.ndarray] = zp
        self.xgrid: Optional[np.ndarray] = xgrid
        self.ygrid: Optional[np.ndarray] = ygrid
        self.zgrid: Optional[np.ndarray] = zgrid

        # Data Normalization Parameters
        self.code_rho: float = code_rho
        self.code_cs: float = code_cs
        self.code_omega: float = code_omega

        # Self-Gravity Parameters
        self.aps: Optional[np.ndarray] = aps
        self.rhopswarm: Optional[np.ndarray] = rhopswarm
        self.eps_dtog: Optional[float] = eps_dtog
        self.init_var: Optional[np.ndarray] = init_var

        # Validations and Defaults
        if isinstance(self.grain_rho, np.ndarray) == False:
            if self.kappa is None and self.grain_rho != 1.675:
                raise ValueError(
                    f"The DSHARP dust model assumes a dust grain density of 1.675 g/cm^3! "
                    f"The input density is: {self.grain_rho} g/cm^3. "
                    f"Either change the dust grain density or input the appropriate opacities (kappa/sigma arguments)."
                )
        else:
            if self.kappa is None and np.any(self.grain_rho != 1.675):
                raise ValueError(
                    f"The DSHARP dust model assumes a dust grain density of 1.675 g/cm^3! "
                    f"The input density is: {self.grain_rho} g/cm^3. "
                    f"Either change the dust grain density or input the appropriate opacities (kappa/sigma arguments)."
                )

        if self.include_scattering:
            if self.kappa is not None and self.sigma is None:
                raise ValueError('The include_scattering parameter is enabled but no scattering coefficient (sigma) was provided!')
            if self.sigma is not None and self.kappa is None:
                raise ValueError('The include_scattering parameter is enabled but no absorption coefficient (kappa) was provided!')

        if self.data is None:
            print('No data input, loading density cube from Yang+14 paper, snapshot taken at orbit 100...')
            self.data, self.axis = load_cube()

        # Additional Attributes
        self.tau: Optional[np.ndarray] = None
        self.intensity: Optional[np.ndarray] = None
        self.mass: Optional[float] = None
        self.mass_excess: Optional[float] = None
        self.filling_factor: Optional[float] = None
        self.area: Optional[float] = None
        self.grain_size: Optional[Union[float, np.ndarray]] = None
        self.proto_mass: Optional[float] = None
        self.grain_size_bins: Optional[np.ndarray] = None

    def configure(self):
        """
        Initializes and configures parameters for the density cube object, including mass, optical depth, 
        and opacity coefficients. This method should be re-run if key parameters such as `sigma`, `H`, or 
        grid spacing (`dx/dy/dz`) are updated.

        Updates:
            - self.Lx, self.Ly, self.Lz: Dimensions of the cube (cm).
            - self.mass: Total mass of the cube (cgs units).
            - self.unit_mass, self.unit_density: Unit conversion factors (code to cgs units).
            - self.data: Density cube converted to cgs units.
            - self.frequency: Frequency corresponding to the analysis wavelength (Hz).
            - self.mass_excess: Ratio of true to observed mass.
            - self.tau: Optical depth (computed later).
            - self.proto_mass: Mass of protoplanets (if `aps` and `rhopswarm` are provided).

        Returns:
            None.
        """

        # Data configuration check
        if self.axis is None: raise ValueError("The `axis` array is not set, need to initialize `self.axis` before running the `configure()` method.")
        if self.data is None: raise ValueError("The density cube (`self.data`) is not set. Provide density data before running the `configure()` method.")

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

        # Wavelength frequency (in Hz) at which to calculate the output intensity and thus the mass excess
        self.frequency = const.c.cgs.value / self.wavelength

        # If opacity coefficients are not input calculate the grain sizes and extract the corresponding DSHARP opacities
        if self.kappa is None:
            try:
                self.calc_grain_size(); self.extract_opacity()
            except Exception as e:
                raise ValueError('Cannot calculate `kappa` -- to calculate the appropriate grain size input the `stoke` and `grain_rho` parameters. Error: ', e)

        # Compute the mass underestimation 
        self.calc_mass_excess()

        # Compute the mass of the protoplanets if both the aps and rhopswarm arrays are input
        if self.aps is not None and self.rhopswarm is not None: self.get_proto_mass()

        return

    def blackbody(self):
        """
        Computes the blackbody spectral radiance using Planck's law for a source in thermal equilibrium 
        at a temperature `T`. Absorption only should be

        Returns:
            None. Assigns the `B_nu` attribute, the spectral radiance in units of erg s⁻¹ cm⁻² Hz⁻¹ steradian⁻¹.
        """

        # Validate inputs
        if self.frequency is None or self.T is None: raise ValueError("Both `frequency` and `T` must be defined before calling the `blackbody()` method.")
        if self.T <= 0: raise ValueError("Temperature (`T`) must be positive!")

        #hf = const.h.cgs.value * self.frequency
        #kT = const.k_B.cgs.value * self.T
        #c2 = const.c.cgs.value**2

        self.B_nu = 2 * const.h.cgs.value * self.frequency**3 / (const.c.cgs.value**2 * (np.exp(const.h.cgs.value * self.frequency / (const.k_B.cgs.value * self.T)) - 1))
        #self.B_nu = 2 * hf**3 / (c2 * (np.exp(hf / kT) - 1))

        return

    def calc_tau(self):
        """
        Computes the optical depth (`self.tau`) for each column in the density cube by integrating along the z-axis.

        Notes:
            - `self.data` should contain the dust density in cgs units.
            - `self.dz` should be the grid spacing in cm.
            - For polydisperse cases, grain positions and grid coordinates must be defined.

        Returns:
            None. Assigns the `self.tau` and `filling_factor`` attributes.
        """

        # Validate
        if self.data is None or self.dz is None: raise ValueError("Density data (`self.data`) and grid spacing (`self.dz`) must be defined before calling `calc_tau`.")
        if self.kappa is None or (self.include_scattering and self.sigma is None): raise ValueError("Opacity coefficients (`kappa` and optionally `sigma`) must be defined.")
        if isinstance(self.grain_size, np.ndarray): #Polydisperse case
            if self.ipars is None or self.xp is None or self.xgrid is None:
                raise ValueError("Particle positions (`ipars`, `xp`, and grid coordinates) must be defined for polydisperse simulations.")

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
                    # This loop goes through all the z-cells in the particular (x,y) column and stores the opacities in self.effective_kappa and self.effective_sigma
                    for k in range(self.Nz):

                        # To count the opacities at each (x,y,z) cell (equivalent to N*k)
                        weighted_kappa, weighted_sigma = 0, 0
                        #denominator = 0
                        # Add the number of grains in that cell for a given species times that species' opacity coefficient (stored in self.kappa/self.sigma which in the polydisperse case are arrays)
                        for species_type in range(num_species):
                            # #
                            weighted_kappa += self.density_per_species[species_type][k, j, i] * self.kappa[species_type] # This is the numerator (N1*k1 + N2*k2 + ...)
                            weighted_sigma = weighted_sigma + self.density_per_species[species_type][k, j, i] * self.sigma[species_type] if self.include_scattering else 0
                            # #
                            #denominator += self.density_per_species[species_type][k, j, i]

                        # Divide by the total number of species in that (x,y,z) cell to compute the weighted mean and to the 3D opacity array
                        self.effective_kappa[k, j, i] = weighted_kappa / np.sum(self.density_per_species[:, k, j, i]) if np.sum(self.density_per_species[:, k, j, i]) != 0 else 0 # Avoid division by zero # Denominator is also np.sum(self.density_per_species[:, k, j, i])
                        if self.include_scattering:
                            self.effective_sigma[k, j, i] = weighted_sigma / np.sum(self.density_per_species[:, k, j, i]) if np.sum(self.density_per_species[:, k, j, i]) != 0 else 0 # Avoid division by zero 
                   
                    # Now that the opacities have been calculated for that column, vertically integrate to find the opacity in that column
                    self.tau[j, i] = np.trapz(self.data[:, j, i] * (self.effective_kappa[:, j, i] + self.effective_sigma[:, j, i])) * self.dz # * self.unit_sigma
                    
        # Calculate the ratio of cells that are optically thick (tau >= 1)
        self.filling_factor = len(np.where(self.tau >= 1)[0]) / (self.Nx * self.Ny)

        return 

    def calc_t(self, rhod, effective_kappa, effective_sigma):
        """
        Computes the cumulative optical depth along a single column, integrating from the base to a given height z.
        This is the optical depth as emission progresses up the column, where as the optical_depth() function calculates
        along the entire column. 
    
        Notes:
            - The grid spacing `self.dz` is assumed to be in cm.
            - All input arrays must have the same length as the number of vertical cells (`self.Nz`).

        Args:
            rhod (ndarray): 1D array of dust densities along the column (units: g/cm³).
            effective_kappa (ndarray): 1D array of absorption opacities along the column (units: cm²/g).
            effective_sigma (ndarray): 1D array of scattering opacities along the column (units: cm²/g).

        Returns:
            ndarray: 1D array of cumulative optical depth values at each height z (unitless).
        """

        # Validate inputs
        if rhod.ndim != 1 or effective_kappa.ndim != 1 or effective_sigma.ndim != 1: raise ValueError("All input arrays must be 1D.")
        if not (len(rhod) == len(effective_kappa) == len(effective_sigma) == self.Nz): raise ValueError("All input arrays must have the same length as the column height (Nz).")
        if self.dz <= 0: raise ValueError("Grid spacing `dz` must be positive.")

        # To store the emission along the z-column
        t = np.zeros(self.Nz) 

        # Integrate starting at the first cell of the column and move upward adding one cell at a time
        for i in range(self.Nz):
            t[i] = np.trapz(rhod[:i] * (effective_kappa[:i] + effective_sigma[:i])) * self.dz
            
        return t 

    def calc_intensity(self):
        """
        Calculates the outgoing intensity along the z-axis for all columns in the cube using the radiative transfer equation.
        
        This method implements the general solution for radiative transfer (RT) equation, assuming no back-illumination.
        It calculates the emissivity and integrates it along the column, taking into account optical depth and 
        scattering (if enabled). #Compute the bolometric flux
        
        Note:
            - The integration assumes the first term of the RT solution (extinction of original intensity) is zero.
            - This method is optimized for polydisperse simulations, but supports monodisperse cases as well.
            - Units of intensity are consistent with the input source function `self.B_nu`.

        Returns:
            None. Assigns the `intensity` attribute. 
        """

        # Validate 
        if self.tau is None: raise ValueError("No optical depth map exists! Run the `calc_tau()` method first.")
        if self.data is None or self.data.ndim != 3: raise ValueError("Density data `self.data` must be a 3D array (Nz, Ny, Nx).")
        if not (self.effective_kappa.shape == self.data.shape == self.effective_sigma.shape): raise ValueError("Opacity attributes (`effective_kappa`, `effective_sigma`) must match the shape of `self.data`.")

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
        
    def calc_mass_excess(self):
        """
        Calculates the mass underestimation (mass_excess) for the dust density cube.

        This method computes the optical depth map, evaluates the outgoing intensity 
        using radiative transfer, and calculates the observed dust mass based on the 
        assumptions of optically thin emission. It then compares the observed mass 
        with the true mass of the cube to estimate the mass excess factor.

        Returns:
            None. Assigns the mass underestimation factor, `self.mass_excess`.
        """

        # Validate
        if self.data is None: raise ValueError("Density data `self.data` is not initialized.")
        if self.mass is None: raise ValueError("The total mass of the box `self.mass` is not initialized. Run `configure()` first.")
        if self.kappa is None or (self.include_scattering and self.sigma is None): raise ValueError("Opacity coefficients (`kappa` and `sigma`) must be defined before calculating mass excess.")

        # Calculate the optical depth map 
        self.calc_tau()

        # Compute the radiation assuming a Planckian black body
        self.blackbody()

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
        ### Compute the mass excess / underestimation 
        ###

        # Monodisperse
        if isinstance(self.grain_size, np.ndarray) is False: 

            # In the absorption only case the, source function is independent of optical depth, therefore can use the following simplified form
            if self.include_scattering is False:
                
                # Set the opactities -- the opacity at every cell is the same (cells with zero dust density will be zeroed out during intensity calculation)
                self.effective_kappa, self.effective_sigma = np.zeros((self.Nz, self.Ny, self.Nx)), np.zeros((self.Nz, self.Ny, self.Nx))
                self.effective_kappa[::], self.effective_sigma[::] = self.kappa, 0 
                
                # Solve the RT solution to set the self.intensity parameter
                self.calc_intensity()
                
                # Integrating the general RT solution with constant T and source function as well as I(0)=0 simplifies to the following
                #self.intensity = self.B_nu * (1 - np.exp(-self.tau)) # Not used, rather solve the full RT equation using calc_intensity() method

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
                #self.sigma_dust = np.mean(self.intensity) / (self.B_nu * (self.kappa + self.sigma)) # Convolution theory -- take the mean of the output intensity
                self.sigma_dust = np.mean(self.intensity) / (self.B_nu * self.kappa) # Correction: Only the absorption opacity is considered when analytically approximating the observed dust mass
            
        #Polydisperse
        else: 

            # Solve the RT solution to set the self.intensity parameter
            self.calc_intensity() 

            # The average opacity an observer would assume (assumes all grain sizes are equally distributed so scales with the inverse number of species)
            #self.assumed_opacity = (np.sum(self.kappa) + np.sum(self.sigma)) / len(self.grain_size) if self.include_scattering else np.sum(self.kappa) / len(self.grain_size)
            self.assumed_opacity = np.sum(self.kappa) / len(self.grain_size) # Correction: Only the absorption opacity is considered when analytically approximating the observed dust mass
            
            # Under the assumption of optically thin emission, the observed intensity scales with the column density of the dust, allowing us to analytically solve for Σd as 
            self.sigma_dust = np.mean(self.intensity) / (self.B_nu * self.assumed_opacity) # Convolution theory -- take the mean of the output intensity and the source function

        # The observed mass of the box can now be quantified as the product of Σd and the domain area
        self.observed_mass = self.sigma_dust * self.area  

        # The mass underestimation, ratio of true box mass to the observed mass
        self.mass_excess = self.mass / self.observed_mass

        return 

    def calc_grain_size(self):
        """
        Calculates the grain size(s) based on the Stokes number and gas column density.

        This method computes the grain size(s) using the formula:

            a = St * (2 * Σ_g) / (π * ρ_g),

        where `a` is the grain size, `St` is the Stokes number, `Σ_g` is the gas column density,
        and `ρ_g` is the internal dust grain density.

        Returns:
            None. Assigns the `grain_size` attribute.
        """

        # Validate
        if self.stoke is None: raise ValueError("The `stoke` parameter must be provided to calculate grain sizes.")
        if self.column_density <= 0: raise ValueError("The gas column density (`column_density`) must be positive.")
        if self.grain_rho is None or np.any(np.array(self.grain_rho) <= 0): raise ValueError("The internal grain density (`grain_rho`) must be positive.")

        if isinstance(self.stoke, np.ndarray) is False:
            # Monodisperse
            self.grain_size = self.stoke * 2. * self.column_density / np.pi / self.grain_rho
        else:
            # Polydisperse
            if isinstance(self.grain_rho, np.ndarray) is False:
                raise ValueError("If entering multiple stoke's numbers, the corresponding grain_rho paramater must be a list/array of same size!")
            
            # Calculate the grain sizes corresponding to each stokes number
            self.grain_size = np.zeros(len(self.stoke)) 

            for grain in range(len(self.grain_size)):
                self.grain_size[grain] = self.stoke[grain] * 2. * self.column_density / np.pi / self.grain_rho[grain]
                
        return

    def extract_opacity(self):
        """
        Extracts the opacity coefficients (`kappa`, `sigma`) based on the grain size(s).

        This method computes the absorption (`kappa`) and scattering (`sigma`) opacity coefficients 
        for the dust grain sizes in the simulation. The coefficients are calculated using the 
        DSHARP model (Birnstiel et al., 2018).

        Returns:
            None. Assigns the `kappa`, `sigma`, and `grain_size_bins` attributes.
        """

        # Grain size, absorption opacity, and scattering opacity -- from DSHARP project. Note that the bin_approx determines whether the simulation is polydisperse.       
        self.kappa, self.sigma, self.grain_size_bins = compute_opacities.dsharp_model(self.q, self.wavelength, self.grain_size, bin_approx=isinstance(self.grain_size, np.ndarray))
       
        return 

    def get_proto_mass(self):
        """
        Calculates the mass of the protoplanets (planetesimals) based on dust-to-gas ratio,
        local swarm density, and azimuthal particle positions.

        The method computes the planetesimal masses using the formula:

            M_p = 10^(log10(N_clump) + log10(M_particle)),

        where `N_clump` is the number of particles in a clump (scaled by the swarm density),
        and `M_particle` is the particle mass derived from the dust-to-gas ratio.

        If no planetesimals are found, the `proto_mass` attribute is set to zero.

        Returns:
            None. Assigns the `proto_mass` attribute.
        """

        # Validate required inputs
        if self.eps_dtog is None: raise ValueError("The `eps_dtog` parameter (dust-to-gas ratio) must be provided.")
        if self.mass is None: raise ValueError("The `mass` parameter must be computed (via `configure()`) before calling this method.")
        if self.npar <= 0: raise ValueError("The number of particles (`npar`) must be a positive integer.")
        if self.rhopswarm is None: raise ValueError("The `rhopswarm` parameter (swarm density) must be provided.")
        if self.aps is None: raise ValueError("The `aps` parameter (azimuthal particle positions) must be provided.")

        # Calculate the mass of a single particle in code units
        mp_code = self.eps_dtog * self.mass / self.npar # Note: this is no longer in code units, as self.mass is coverted during the configuration
        mp = stats.mode(self.rhopswarm)[0]

        # Identify indices of particles in azimuthal clumps
        index = np.where(self.aps != 0)[0]

        # Calculate the number of particles in each clump
        npclump = self.rhopswarm[index] / mp

        tmp = np.log10(npclump)
        ttmp = tmp + np.log10(mp_code)
        mass = 10**ttmp # Total mass of each clump

        self.proto_mass = np.sort(mass)
        
        # If there are no planetesimals set to 0
        try: 
            self.proto_mass.max()
        except ValueError:
            self.proto_mass = 0

        return


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


def load_opacity_values_old(q=None):
    """
    Loads the opacity values taken from the DSHARP project.

    Note:
        These coefficients correspond to the 1mm wavelength opacities only! This function has been
        replaced with the compute_opacities module which allows for opacity calculations across
        a wide range of wavelengths.

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
