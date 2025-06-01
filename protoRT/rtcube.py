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

class RadiativeTransferCube:
    """
    Class for radiative transfer through a 3D dust density cube and for computing the mass excess.

    This class computes the optical depth and outgoing intensity through a
    protoplanetary disk model by performing radiative transfer along the z-axis
    of a 3D density cube. It supports both monodisperse and polydisperse dust
    distributions and can automatically compute dust opacities using the
    DSHARP model (Birnstiel et al. 2018) if they are not provided.

    It also computes the mass excess — defined as the ratio of the true dust mass in the cube
    to the mass inferred from the observed intensity under the assumption of optically thin emission.
    This metric quantifies how much disk mass may be underestimated in observations due to
    optically thick regions.

    Note
    ----
    - The current version of the code assumes that the cube is symmetric along all axes (i.e., cubic domain).
    - If no input data is provided, a default test dataset from a streaming
      instability simulation by Yang & Johansen (2014) is used. This is a snapshot at orbit 100.

    Parameters
    ----------
    data : ndarray, optional
        3D dust density cube (e.g., rhop field). Defaults to test data if None.
    axis : ndarray, optional
        1D coordinate array along the integration axis (must be `z`), in units of gas scale height.
    code_rho : float, optional
        Midplane gas density in code units (`rho0` in `start.in`). Default is 1.
    code_cs : float, optional
        Sound speed in code units (`cs0` in `start.in`). Default is 2π.
    code_omega : float, optional
        Orbital frequency in code units (`Omega` in `start.in`). Default is 2π.
    column_density : float, optional
        Gas column density in cgs units (g/cm²). Default is 100.
    T : float, optional
        Isothermal box temperature (K). Default is 30.
    H : float, optional
        Pressure scale height (cm). Default is 7.5e13 (≈ 5 AU).
    stoke : float or ndarray
        Stokes number(s) for the dust population.
    grain_rho : float or ndarray, optional
        Internal grain density (g/cm³). If None, defaults to 1.675, which is from the DSHARP dust model
        and must be used if no opacities are input as DSHARP opacities are used by default.
    wavelength : float, optional
        Wavelength at which opacities are computed (cm). Range: 1e-5 to 10. Default is 0.1.
    include_scattering : bool, optional
        If True, include scattering opacity in radiative transfer. Default is False.
    kappa : float or ndarray, optional
        Dust absorption opacity (cm²/g). If None, computed from DSHARP.
    sigma : float or ndarray, optional
        Dust scattering opacity (cm²/g). If None, computed from DSHARP.
    p : float, optional
        Power-law index for grain size distribution (n(a) ∝ a^{-p}). Default is 2.5.
    npar : int, optional
        Number of particles in the simulation. Only required for multi-species simulations. Default is 1,000,000.
    ipars, xp, yp, zp : ndarray, optional
        Particle species and positions. Only required for multi-species simulations.
    xgrid, ygrid, zgrid : ndarray, optional
        Grid coordinates used to bin particles into grid cells for opacity weighting. Only required for multi-species simulations.
    rhopswarm : ndarray, optional
        Local swarm densities used in density map conversion and proto-mass calculation. Only required for polydisperse and/or self-gravitating simulations.
    grid_func : {'linear'}, optional
        Grid interpolation scheme for each axis. Only `'linear'` is currently supported. Pencil Code stores these as attributes in read_param(). Default is 'linear'.
    num_grid_points : int
        Number of grid points in x, y, and z directions. Current code assumes cubic domain. Pencil code stores these as attributes in read_dim(). Default is 262.
    num_interp_points : int
        Number of interpolation points in x, y, and z. Current code assumes cubic domain. Controls smoothing kernel width. Pencil code stores these as attributes in read_dim(). Default is 256.
    index_limits_1 : int
        Grid index limits along x (l1), y (m1), and z (n1) directions for trimming ghost zones. Pencil code stores these as attributes in read_dim(). Default is 3.
    index_limits_2 : int
        Grid index limits along x (l2), y (m2), and z (n2) directions for trimming ghost zones. Pencil code stores these as attributes in read_dim(). Default is 258.
    aps : ndarray, optional
        Azimuthal positions of particles for proto-mass calculation. Only required for self-gravitating simulations. Must be input to enable protomass calculations.
    eps_dtog : float, optional
        Dust-to-gas ratio used to estimate proto-masses. Only required for self-gravitating simulations. Must be input to enable protomass calculations.
    init_var : ndarray, optional
        Initial density cube for mass excess calculation. Only required for self-gravitating simulations. In these cases the initial mass should be considered in the mass excess calculation.

    Attributes
    ----------
    axis : ndarray
        The 1D coordinate array corresponding to the integration axis (must be the z-axis), in units of gas scale height.
    data : ndarray
        The 3D dust density cube (rhop field), used for radiative transfer.
    filling_factor : float
        Fraction of columns with τ ≥ 1 (optically thick).
    grain_size : float or ndarray
        Derived grain size(s) based on Stokes number(s).
    H : float
        Gas pressure scale height of the box (cm).
    intensity : ndarray
        Outgoing intensity map computed along the z-axis.
    mass : float
        Total dust mass in the cube (g).
    mass_excess : float
        Ratio of true dust mass to inferred dust mass (from the optically thin approximation).
    proto_mass : float
        Total mass in gravitationally bound clumps (g), if simulation is self-gravitating.
    tau : ndarray
        Optical depth map computed along the z-axis.
    
    Methods
    -------
    configure()
        Initializes internal attributes and computes the opacity (if not provided), calculates the optical depth, 
        performs the radiative transfer, and then computes the mass excess.
    blackbody()
        Computes blackbody intensity from Planck's law at temperature T.
    calc_tau()
        Computes the 2D optical depth (τ) map.
    calc_t(rhod, kappa, sigma)
        Computes cumulative τ along one vertical column.
    calc_intensity()
        Solves the radiative transfer equation along the z-axis.
    calc_mass_excess()
        Calculates the dust mass underestimation factor.
    calc_grain_size()
        Converts Stokes number(s) into physical grain size(s).
    extract_opacity()
        Computes opacities using the DSHARP dust model.
    get_proto_mass()
        Estimates mass in bound clumps (protoplanets) using swarm density and azimuthal info.
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
        p: float = 2.5,
        npar: int = int(1e6),
        ipars: Optional[np.ndarray] = None,
        xp: Optional[np.ndarray] = None,
        yp: Optional[np.ndarray] = None,
        zp: Optional[np.ndarray] = None,
        xgrid: Optional[np.ndarray] = None,
        ygrid: Optional[np.ndarray] = None,
        zgrid: Optional[np.ndarray] = None,
        rhopswarm: Optional[np.ndarray] = None,
        grid_func: Optional[str] = 'linear',
        num_grid_points: Optional[int] = 262,
        num_interp_points: Optional[int] = 256,
        index_limits_1: Optional[int] = 3,
        index_limits_2: Optional[int] = 258,
        aps: Optional[np.ndarray] = None,
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
        self.p: float = p

        # Polydisperse Simulation Parameters
        self.npar: int = npar
        self.ipars: Optional[np.ndarray] = ipars
        self.xp: Optional[np.ndarray] = xp
        self.yp: Optional[np.ndarray] = yp
        self.zp: Optional[np.ndarray] = zp
        self.xgrid: Optional[np.ndarray] = xgrid
        self.ygrid: Optional[np.ndarray] = ygrid
        self.zgrid: Optional[np.ndarray] = zgrid
        self.rhopswarm: Optional[np.ndarray] = rhopswarm
        self.grid_func: Optional[str] = grid_func
        self.num_grid_points: Optional[int] = num_grid_points
        self.num_interp_points: Optional[int] = num_interp_points
        self.index_limits_1: Optional[int] = index_limits_1
        self.index_limits_2: Optional[int] = index_limits_2

        # Data Normalization Parameters
        self.code_rho: float = code_rho
        self.code_cs: float = code_cs
        self.code_omega: float = code_omega

        # Self-Gravity Parameters
        self.aps: Optional[np.ndarray] = aps
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
            print('No data input, loading density cube from Yang & Johansen (2014), snapshot taken at orbit 100...')
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
        Initialize and configure physical parameters for the density cube object.

        This method performs unit conversions, sets up spatial dimensions,
        computes the total mass of the cube in cgs units, and prepares
        quantities required for radiative transfer and mass excess analysis.

        It must be re-run if key inputs, such as the scale height `H`, the 
        absorption/scattering opacities (`kappa`, `sigma`), or grid spacing (`axis`),
        are updated.

        Raises
        ------
        ValueError
            If `self.data` or `self.axis` are not set prior to calling this method.

        Returns
        -------
        None
            This method updates the following internal attributes in place:

            - `self.Lx`, `self.Ly`, `self.Lz` : float
                Physical dimensions of the cube (cm).
            - `self.dx`, `self.dy`, `self.dz` : float
                Grid cell spacing (cm).
            - `self.mass` : float
                Total dust mass in the cube (g).
            - `self.unit_mass` : float
                Conversion factor from code mass units to CGS (g).
            - `self.unit_density` : float
                Conversion factor from code density units to CGS (g/cm³).
            - `self.unit_sigma` : float
                Conversion factor for dust surface density (g/cm²).
            - `self.data` : ndarray
                Density cube converted to CGS units (g/cm³).
            - `self.frequency` : float
                Frequency corresponding to the specified wavelength (Hz).
            - `self.mass_excess` : float
                Ratio of true dust mass to inferred mass assuming optically thin emission.
            - `self.tau` : ndarray
                Optical depth map (computed in later steps).
            - `self.proto_mass` : float
                Total protoplanet mass (if `aps` and `rhopswarm` are available).
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
        Compute the blackbody spectral radiance using Planck's law.

        This method calculates the specific intensity (spectral radiance) of a blackbody
        at temperature `T` and frequency `self.frequency`, assuming pure absorption
        (i.e., no scattering contribution). The result is stored in `self.B_nu`.

        Raises
        ------
        ValueError
            If `self.frequency` or `self.T` is not set, or if `T <= 0`.

        Returns
        -------
        None
            Updates the `self.B_nu` attribute with the spectral radiance
            in units of erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹.
        """

        # Validate inputs
        if self.frequency is None or self.T is None: raise ValueError("Both `frequency` and `T` must be defined before calling the `blackbody()` method.")
        if self.T <= 0: raise ValueError("Temperature (`T`) must be positive!")

        self.B_nu = 2 * const.h.cgs.value * self.frequency**3 / (const.c.cgs.value**2 * (np.exp(const.h.cgs.value * self.frequency / (const.k_B.cgs.value * self.T)) - 1))

        return

    def calc_tau(self):
        """
        Compute the optical depth map (`self.tau`) by integrating along the z-axis of the density cube.

        For monodisperse simulations, a uniform opacity coefficient is applied throughout.
        For polydisperse simulations, species-specific densities are reconstructed in 3D,
        and local weighted opacity coefficients are calculated at each cell using the particle data.

        This method updates `self.tau` for all (x, y) columns and computes the `self.filling_factor`,
        defined as the fraction of columns with optical depth τ ≥ 1.

        Notes
        -----
        - `self.data` must be in cgs units (g/cm³).
        - `self.dz` must be the vertical grid spacing in cm.
        - `self.kappa` (absorption opacity) is required.
        - If `include_scattering` is True, `self.sigma` (scattering opacity) must also be provided.
        - For polydisperse simulations:
            - `self.grain_size` must be an array.
            - `self.ipars`, `self.xp`, and `self.xgrid` must be defined to compute density per species.

        Raises
        ------
        ValueError
            If required attributes (`self.data`, `self.dz`, `self.kappa`, etc.) are not defined.
            If any necessary polydisperse inputs are missing.

        Returns
        -------
        None
            The following attributes are updated in place:

            - `self.tau` : ndarray
                2D optical depth map integrated along the z-axis for each (x, y) column.
            - `self.filling_factor` : float
                Fraction of columns with τ ≥ 1 (optically thick).
            - `self.effective_kappa` : ndarray
                3D array of locally averaged absorption opacities (polydisperse only).
            - `self.effective_sigma` : ndarray
                3D array of locally averaged scattering opacities (if `include_scattering=True`).
            - `self.density_per_species` : ndarray
                4D array of species-separated dust density grids (polydisperse only).
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
                
                # Index the local swarm density accordingly
                rhop_swarm = self.rhopswarm[index_species]

                # Convert positions of particles to a grid density field
                particle_density = particles_to_density(xp=species_x, yp=species_y, zp=species_z, 
                    x=self.xgrid, y=self.ygrid, z=self.zgrid,
                    rhop_swarm=rhop_swarm, grid_func1=self.grid_func, grid_func2=self.grid_func, grid_func3=self.grid_func, 
                    mx=self.num_grid_points, my=self.num_grid_points, mz=self.num_grid_points,
                    nx=self.num_interp_points, ny=self.num_interp_points, nz=self.num_interp_points, 
                    n1=self.index_limits_1, n2=self.index_limits_2, 
                    m1=self.index_limits_1, m2=self.index_limits_2, 
                    l1=self.index_limits_1, l2=self.index_limits_2, 
                    density=True)

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
        Compute the cumulative optical depth along a single vertical column.

        This method integrates the optical depth from the bottom of the cube (z = 0)
        upward to each height `z`, returning a 1D array of cumulative τ values.
        Unlike `calc_tau`, which computes total τ for an entire column, this returns
        the depth-dependent progression of optical depth as a function of height.

        Parameters
        ----------
        rhod : ndarray
            1D array of dust mass densities along the vertical column (g/cm³).
        effective_kappa : ndarray
            1D array of absorption opacity coefficients along the column (cm²/g).
        effective_sigma : ndarray
            1D array of scattering opacity coefficients along the column (cm²/g).

        Notes
        -----
        - All input arrays must be 1D and of length `self.Nz`.
        - The grid spacing `self.dz` is assumed to be in cm and must be positive.

        Raises
        ------
        ValueError
            If input arrays are not 1D or do not match the vertical grid size `self.Nz`,
            or if `self.dz` is not positive.

        Returns
        -------
        t : ndarray
            1D array of cumulative optical depth values at each vertical cell (unitless).
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
        Compute the outgoing intensity along the z-axis for all (x, y) columns.

        This method solves the radiative transfer equation using the emissivity
        term only, assuming no incident intensity at the base (i.e., no
        back-illumination). It integrates the column-wise emissivity weighted
        by the exponential attenuation factor exp(−τ) to determine the observed
        intensity at the top of the cube.

        Scattering is included if `self.include_scattering` is True, in which
        case the mean intensity `self.J` is also used in the emissivity term.
        Supports both monodisperse and polydisperse models.

        Notes
        -----
        - Requires `self.tau` to be precomputed using `calc_tau()`.
        - Assumes that `self.data` contains dust density in cgs units (g/cm³).
        - Assumes Planck function (`self.B_nu`) is constant along each column.
        - The returned intensity is in units of erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹.
        - If scattering is enabled, `self.effective_sigma` and `self.J` must be defined.

        Raises
        ------
        ValueError
            If required attributes such as `self.data`, `self.tau`, or opacity arrays are missing or have incorrect shapes.

        Returns
        -------
        None
            The following attribute is updated:

            - `self.intensity` : ndarray
                2D map of outgoing specific intensity at the top of the cube for each (x, y) position.
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
        Calculate the mass underestimation factor ("mass excess") for the dust density cube.

        This method solves the radiative transfer equation along each column, calculates
        the outgoing intensity, and derives the observed dust mass assuming optically thin
        emission. The ratio of the true dust mass to the inferred dust mass is stored
        as `self.mass_excess`, quantifying how much mass is missed due to optically
        thick regions or scattering.

        Scattering is handled via a two-stream approximation (Miyake & Nakagawa 1993; Zhu et al. 2019),
        which modifies the source function. Supports both monodisperse and polydisperse models.

        Notes
        -----
        - Requires `self.kappa` to be defined. If `include_scattering=True`, `self.sigma` is also required.
        - The cube must be configured first with `configure()`, which sets `self.mass` and grid info.
        - The outgoing intensity is computed using `calc_intensity()`, based on the full RT solution.
        - Assumes thermal equilibrium (LTE) and no external illumination (I₀ = 0).
        - For polydisperse simulations, the average opacity assumed in the optically thin approximation
          is taken to be the mean of the absorption opacities across species.

        Raises
        ------
        ValueError
            If required attributes such as density data, opacities, or mass are not initialized.

        Returns
        -------
        None
            The following attributes are updated:

            - `self.src_fn` : float or ndarray
                The effective source function used in radiative transfer integration.
            - `self.intensity` : ndarray
                2D array of emergent intensity for each (x, y) column.
            - `self.observed_mass` : float
                The dust mass inferred under the assumption of optically thin emission.
            - `self.mass_excess` : float
                The true mass divided by the observed mass. A value > 1 indicates underestimation.
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
        Compute the grain size(s) from the Stokes number and gas column density.

        This method uses the Stokes–Epstein drag relation to calculate the physical grain size `a`
        corresponding to a given Stokes number `St`, assuming all grains are in the Epstein regime.
        For monodisperse models, a single value is returned; for polydisperse models, an array of
        grain sizes is computed based on the per-species Stokes number and internal grain density.

        The grain size is given by:

            a = (St * 2 * Σ_g) / (π * ρ_g)

        where:
          - `a` is the grain size (cm),
          - `St` is the Stokes number (dimensionless),
          - `Σ_g` is the gas column density (g/cm²),
          - `ρ_g` is the internal grain density (g/cm³).

        Notes
        -----
        - Assumes Epstein drag regime is valid.
        - For polydisperse simulations, both `stoke` and `grain_rho` must be arrays of equal length.
        - Requires `self.column_density` > 0.

        Raises
        ------
        ValueError
            If `stoke` or `grain_rho` is not set, non-positive, or if their lengths mismatch in the polydisperse case.

        Returns
        -------
        None
            The computed grain sizes are stored in the following attribute:

            - `self.grain_size` : float or ndarray
                Computed grain size(s) in cm.
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
        Compute and assign dust opacity coefficients using the DSHARP model.

        This method calculates the absorption (`kappa`) and scattering (`sigma`) opacity
        coefficients based on the dust grain sizes in the simulation. Opacities are derived
        using the DSHARP dust opacity tables (Birnstiel et al. 2018), which assume a fixed
        internal grain composition and density.

        For polydisperse simulations (`grain_size` is an array), the method computes 
        opacity values using a binned approximation consistent with the multi-species model.

        Notes
        -----
        - The `grain_size` and `wavelength` attributes must be set prior to calling this method.
        - The method internally calls the compute_opacities module (`compute_opacities.dsharp_model()`).

        Raises
        ------
        ValueError
            If `grain_size` or `wavelength` is not set or invalid (implicitly by the underlying function in the compute_opacities module).

        Returns
        -------
        None
            Updates the following attributes in place:

            - `self.kappa` : float or ndarray
                Absorption opacity coefficient(s) (cm²/g).
            - `self.sigma` : float or ndarray
                Scattering opacity coefficient(s) (cm²/g).
            - `self.grain_size_bins` : ndarray
                Grain size bin edges used in the binned opacity approximation.
        """

        # Grain size, absorption opacity, and scattering opacity -- from DSHARP project. Note that the bin_approx determines whether the simulation is polydisperse.       
        self.kappa, self.sigma, self.grain_size_bins = compute_opacities.dsharp_model(self.p, self.wavelength, self.grain_size, bin_approx=isinstance(self.grain_size, np.ndarray))
       
        return 

    def get_proto_mass(self):
        """
        Estimate the masses of protoplanets (planetesimal clumps) from local particle densities.

        This method calculates the mass of planetesimal clumps based on the local swarm density
        (`rhopswarm`), azimuthal particle positions (`aps`), and the total box dust mass.
        The particle mass is inferred from the total number of particles and the global dust-to-gas ratio.

        Clumps are identified via non-zero azimuthal particle positions, and individual clump masses
        are computed as:

            M_p = 10 ** [log10(N_clump) + log10(M_particle)]

        where:
        - `N_clump` is the number of particles in a given clump (estimated from density),
        - `M_particle` is the dust mass per particle.

        Notes
        -----
        - Requires `configure()` to have been run (sets `self.mass`).
        - Assumes that `rhopswarm` is proportional to the number of particles per grid cell.
        - Returns clump masses in cgs units (g), sorted in ascending order.
        - If no clumps are found, `self.proto_mass` is set to 0.

        Raises
        ------
        ValueError
            If required attributes (`eps_dtog`, `mass`, `rhopswarm`, `aps`, `npar`) are not defined.

        Returns
        -------
        None
            The resulting clump masses are stored in the `self.proto_mass` attribute:

            - `self.proto_mass` : float or ndarray
                Sorted array of clump masses in grams, or 0 if no clumps are detected.
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


def particles_to_density(
    xp, yp, zp,
    x, y, z,
    rhop_swarm=6.30813659928,
    grid_func1='linear', grid_func2='linear', grid_func3='linear',
    mx=262, my=262, mz=262, nx=256, ny=256, nz=256,
    n1=3, n2=258, m1=3, m2=258, l1=3, l2=258,
    density=True):
    """
    Deposit particles onto a 3D grid to generate a dust density field.

    This function interpolates particle mass or count contributions onto a 3D mesh 
    using a quadratic weighting scheme adapted from an IDL routine by Anders Johansen.
    The method accounts for non-uniform grid spacing and allows interpolation weights 
    to be computed using either particle mass (`rhop_swarm`) or uniform weighting.

    Parameters
    ----------
    xp, yp, zp : array_like
        Particle positions in x, y, and z. Pencil Code stores these as attributes in the pvar files.
    x, y, z : array_like
        Grid node coordinates along x, y, and z axes, respectively. Pencil Code stores these as attributes in read_grid().
    rhop_swarm : float or array_like, optional
        Particle density or weight. If a float, all particles have the same weight;
        if an array, must match length of xp. Default is 6.30813659928. Pencil Code stores this as attribute in read_param().
    grid_func1, grid_func2, grid_func3 : {'linear'}, optional
        Grid interpolation scheme for each axis. Only `'linear'` is currently supported. Pencil Code stores these as attributes in read_param().
    mx, my, mz : int
        Number of grid points in x, y, and z directions. Pencil code stores these as attributes in read_dim().
    nx, ny, nz : int
        Number of interpolation points in x, y, and z. Controls smoothing kernel width. Pencil code stores these as attributes in read_dim().
    l1, l2, m1, m2, n1, n2 : int
        Grid index limits along x (l), y (m), and z (n) directions for trimming ghost zones. Pencil code stores these as attributes in read_dim().
    density : bool, optional
        If True, return a mass-weighted density field. If False, return particle count field. Default is True.

    Returns
    -------
    ndarray
        3D array of shape (n2 - n1 + 1, m2 - m1 + 1, l2 - l1 + 1), representing the deposited
        density or particle count on the grid.

    Notes
    -----
    - The interpolation uses a second-order accurate kernel.
    - The resulting array excludes ghost zones as defined by index limits.
    - Ensure `xp`, `yp`, `zp` fall within bounds of `x`, `y`, `z` to avoid edge errors.

    References
    ----------
    - W. Lyra, private communication (adapted from IDL code by A. Johansen)
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
    Find the index of the element in a sorted list closest to a target value using binary search.

    This function returns the index of the element in the sorted array `q` that is nearest 
    to the target value `qpar`, using an efficient bisection method.

    Parameters
    ----------
    qpar : float
        The target value for which the nearest index is sought.

    q : array_like
        A 1D sorted array of floats (in ascending order) to search within.

    Returns
    -------
    int
        Index of the element in `q` closest to `qpar`.

    Notes
    -----
    This method assumes that `q` is sorted in ascending order. If the array is not sorted,
    the result is undefined.

    References
    ----------
    - W. Lyra, private communication (adapted from IDL code by A. Johansen)
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
    Load a 256×256×256 dust density cube from a simulation snapshot.

    This function loads a 3D density cube corresponding to a single snapshot (orbit 100)
    of a 100-orbit streaming instability simulation in a protoplanetary disk,
    provided by Yang & Johansen (2014). The full cube is stored across two 
    `.npy` files and is concatenated along the z-axis.

    Returns
    -------
    data : ndarray
        A 3D NumPy array of shape (256, 256, 256) representing the dust 
        density cube.
    axis : ndarray
        A 1D NumPy array representing the spatial coordinates corresponding
        to each axis (z, y, x) in code units.
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
