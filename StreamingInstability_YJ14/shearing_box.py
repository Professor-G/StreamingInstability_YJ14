
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
        column_density : Column density of the gas (g / cm^2). Defaults to 100.
        T (float): Temperature of the entire box, as this model is isothermal.
            Defaults to 30.  
        H (float): Scale height of the box, the value is multiplied by one AU,
            in cgs units. Defaults to 5 AU.
        kappa (float):  Dust absorption opacity coefficient, if None then the mm-wave
            dust opacity will be calculated according to column density and grain size,
            using the DSHARP code. Defaults to None.
        sigma (float): Dust scattering opacity coefficient, if None then the mm-wave
            dust opacity will be calculated according to column density and grain size,
            using the DSHARP code. Defaults to None.
        stoke (float): Stoke's number, either a float or an ndarray
        rho_grain (float): Internal grain density, approximately
            1 g/cm^3 for ices, and 3.5 g/cm^3 for silicates. 
            Can either be a float or an ndarray. Must correspond
            with the stokes number.
        eps_dtog (float): Dust to gas ratio, defaults to 0.03. Only used
            to calculate the mass of protoplanets.
        npar (int): Number of particles in the simulation, defaults
            to one million. Only used to calculate the mass of protoplanets.
        aps (ndarray): pvar.aps, only used to calculate the mass of protoplanets.
        rhopswarm (ndarray): pvar.rhopswarm, only used to calculate the mass of protoplanets.
        include_scattering (bool): If kappa is set to None, the DSHARP opacities will be applied. If 
            this parameter is set to True, the scattering opacity will be used in addition to the absorption opacities.
            This is relevant in the context of SI as the mass is dominated by the largest grains. Defaults to False, which '
            will account for the absorption opacities only.
    """
    
    def __init__(self, data=None, axis=None, column_density=100, T=30, H=5, kappa=None, sigma=None,
        stoke=0.3, rho_grain=1.0, eps_dtog=0.03, npar=1e6, aps=None, rhopswarm=None, init_var=None, 
        include_scattering=False):

        self.data = data
        self.axis = axis 
        self.column_density = column_density
        self.T = T 
        self.H = H*const.au.cgs.value
        self.kappa = kappa 
        self.sigma = sigma 
        self.stoke = stoke 
        self.rho_grain = rho_grain
        self.eps_dtog = eps_dtog
        self.npar = npar
        self.aps = aps 
        self.rhopswarm = rhopswarm
        self.init_var = init_var
        self.include_scattering = include_scattering

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
                raise ValueError("If entering multiple stoke's numbers, the corresponding rho_grain parameter must be of same size!")
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

        self.Lz = None 
        self.Ly = None  
        self.Lx = None
        self.Nz = None
        self.Ny = None 
        self.Nx = None 
        self.dz = None 
        self.dy = None 
        self.dx = None 

        self.configure()

    def configure(self, nu=230e9):
        """
        Initializing parameters and creates flux, mass excess, and filling factor 
        attributes. If sigma, H, or dx/dy/dz attributes are updated, re-run this method 
        to re-configure the object.

        Args:
            nu (float): Frequency at which to calculate the observed flux and thus mass excess.
        """

        self.Lz = self.Ly = self.Lx = np.abs(self.axis[0] - self.axis[-1])*self.H 
        self.dz = self.dy = self.dx = np.diff(self.axis)[0]
        self.Nz = self.Ny = self.Nx = len(self.axis)
    
        self.unit_sigma = self.column_density / np.sqrt(2*np.pi)
        self.area = self.Lx * self.Ly

        #Code units
        box_mass_codeunits = np.sum(self.data) * self.dx * self.dy * self.dz 
        unit_mass = self.unit_sigma * self.H**2
        self.mass = box_mass_codeunits * unit_mass 

        if self.kappa is None:
            try:
                self.calc_grain_size(); self.extract_opacity()
            except:
                raise ValueError('Cannot calculate kappa -- to calculate the appropriate grain size input the stoke and rho_grain parameters.')

        self.calc_tau(); self.calc_mass_excess(nu=nu); self.calc_filling_factor()
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
        Integrates density along the z-axis to compute the optical depth.
        
        Returns:
            2D array containing the optical depth values.
        """

        self.unit_sigma = self.column_density / np.sqrt(2*np.pi)
  
        tau = np.zeros([self.Ny, self.Nx]) #python reverses order
               
        for i in range(self.Nx):
            for j in range(self.Ny):
                surface_density = np.trapz(self.data[:,j,i]) * self.dz * self.unit_sigma
                tau[j, i] = surface_density * (self.kappa + self.sigma) if self.include_scattering else surface_density * self.kappa
            
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

        self.unit_sigma = self.column_density / np.sqrt(2*np.pi)

        t = np.zeros(self.Nz)
        surface_density = 0

        for i in range(self.Nz):        
            surface_density += rhod[i] * self.dz * self.unit_sigma
            t[i] = surface_density * (self.kappa + self.sigma) if self.include_scattering else surface_density * self.kappa
            
        return t 

    def calc_flux(self, nu=230e9):
        """
        Calculate outgoing flux using solution for RT eqn (5.113). Not used!

        Returns:
            2D array containing the integrated values along the third axis.
        """

        self.calc_tau()
    
        flux = np.zeros([self.Ny, self.Nx])

        #Estimate flux assuming optically thick emission
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
            Defaults to tau=1.

        Returns:
            Float.
        """
        
        self.unit_sigma = self.column_density / np.sqrt(2*np.pi)

        #Code units
        box_mass_codeunits = np.sum(self.data) if self.init_var is None else np.sum(self.init_var) 
        box_mass_codeunits = box_mass_codeunits * self.dx * self.dy * self.dz 
        
        unit_mass = self.unit_sigma * self.H**2
        self.mass = box_mass_codeunits * unit_mass 

        self.calc_tau()

        #Source function should be per frequency (1mm wavelength ~ 230GHz)
        if self.include_scattering:
            albedo = self.sigma / (self.kappa + self.sigma)
            epsilon = 1 - albedo
            mu = 1./np.sqrt(3.)
            tau_d = 2*mu / 3
            tau = 0
            #numerator = np.exp(-np.sqrt(3*epsilon)*tau) + np.exp(np.sqrt(3*epsilon)*(tau-tau_d))
            #denominator = np.exp(-np.sqrt(3*epsilon)*tau_d)*(np.sqrt(epsilon) - 1) - (np.sqrt(epsilon) + 1)
            #J = self.blackbody(nu=nu)*(1 + (numerator/denominator))
            numerator = 2*np.sqrt(epsilon) * (np.exp(-np.sqrt(3*epsilon)*tau_d) - 1)
            denominator = np.exp(-np.sqrt(3*epsilon)*tau_d)*(np.sqrt(epsilon) - 1) - (np.sqrt(epsilon) + 1)
            J = self.blackbody(nu=nu) * (numerator / denominator)
            src_fn = albedo * J + (1 - albedo) * self.blackbody(nu=nu)
        else:
            src_fn = self.blackbody(nu=nu)
        
        if mask_tau: #If the filaments could actully be resolved! One day in the very distant future...
            mask = np.argwhere(self.tau > threshold)
            ratio = len(mask)/(self.Nx*self.Ny)

            rhop_cropped = []
            for i in range(self.Nz):
                rhop_cropped.append(np.sum(self.data[i][mask]))

            cropped_box_mass_codeunits = np.sum(rhop_cropped) * self.dx * self.dy * self.dz 
            unit_mass = self.unit_sigma * self.H**2
            self.mass = cropped_box_mass_codeunits * unit_mass 

            #If source fn is constant and region is optically thick (Eq. 5.120)
            flux_approx = src_fn * (1-np.exp(-self.tau[mask]))
        
            #Sigma dust observer sees if optically thin 
            sigma_dust = np.mean(flux_approx)/(src_fn*(self.kappa+self.sigma)) if self.include_scattering else np.mean(flux_approx)/(src_fn*self.kappa)
            
            self.observed_mass = sigma_dust*self.area*ratio
            self.mass_excess = self.mass / self.observed_mass

            return 
        
        flux_approx = src_fn * (1 - np.exp(-self.tau))
        sigma_dust = np.mean(flux_approx)/(src_fn*(self.kappa+self.sigma)) if self.include_scattering else np.mean(flux_approx)/(src_fn*self.kappa)
        
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
        
        #self.grain_size *= 2 #Make this the diameter
        
        return

    def extract_opacity(self):
        """
        Returns opacity according to grain size.
        """

        try:
            self.calc_grain_size()
        except:
            raise ValueError('Could not determine grain size, input stoke and rho_grain parameters and try again.')

        a, k_abs, k_sca = load_fig4_values() #Grain size, absorption opacity, scattering opacity
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

            self.kappa = np.zeros(len(self.grain_size))
            for grain in self.grain_size:
                if self.include_scattering:
                    self.kappa[grain], self.sigma[grain] = k_abs_fit(self.grain_size[grain]), k_sca_fit(self.grain_size[grain])
                else:
                    self.kappa = k_abs_fit(self.grain_size)
        return 

    def get_proto_mass(self):
        """
        Calculates the mass of the protoplanets

        Returns:
            Mass of the forming protoplanets
        """

        mp_code = self.eps_dtog * self.mass / self.npar
        mp = stats.mode(self.rhopswarm)[0]

        index = np.where(self.aps != 0)[0]
        npclump = self.rhopswarm[index] / mp

        tmp = np.log10(npclump)
        ttmp = tmp + np.log10(mp_code)
        mass = 10**ttmp

        self.proto_mass = np.sort(mass) 

        return

    def calc_filling_factor(self):
        """
        Calculates the filling factor attribute.

        Returns:
            Float.
        """

        self.calc_tau()
        
        self.filling_factor = len(np.where(self.tau >= 1)[0]) / (self.Nx * self.Ny)

        return 
        
    def plot_tau(self, title='Optical Depth', savefig=False, filename='tau'):
        """
        Plots the optical depth at the exit plane.
        """

        self.calc_tau()

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
        """
        Plots the outgoing flux at the exit plane.
        """

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
    Loads the opacity values taken from DSHARP, Figure 4
    """

    resource_package = __name__
    resource_path = '/'.join(('data', 'dsharp_fig4_values'))
    file = pkg_resources.resource_filename(resource_package, resource_path)
    values = np.loadtxt(file)
    a, k_abs, k_sca = values[:,0], values[:,1], values[:,2]

    order = np.array(a).argsort()
    a, k_abs, k_sca = a[order], k_abs[order], k_sca[order]

    return a, k_abs, k_sca

