
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 09:25:13 2022

@author: daniel
"""
import numpy as  np 
from pathlib import Path
import matplotlib.pyplot as plt  
from cycler import cycler
import astropy.constants as const  

class Model:
    """
    Creates the disk model.

    Disk model from Bai & Stone (2018): https://arxiv.org/pdf/1005.4982.pdf

    This model assumes a solar nebula with a vertically isothermal disk. 

    The disk is parameterized according to a set of power laws, see Youdin & Shu (2002)
    
    Args:
        r (float): Distance in cgs units. 
        b (float): Gas profile law index, defaults to 3/2.
        c (float): Temperature profile law index, defaults to 1/2.
        grain_size (float): Grain size, in cgs units. Only use if no stoke's number
            is input, as that will overwrite this value. Defaults to 1 mm.
        grain_rho (float): Internal grain density in cgs units. Defaults to 1 g/cm2.
        Z (float): The global solids-to-gas ratio. Defaults to 0.02.
        stoke (float): Stoke's number, if input then the grain sizes at each
            r will be computed. Defaults to None, in which case the stoke's
            number at each r will be calculated.
        
    Attributes:
        get_params: Calculates all disk parameters
    """

    def __init__(self, r, b=1.5, c=0.5, grain_size=0.1, grain_rho=1, Z=0.02, stoke=None):

        self.r = r 
        self.b = b
        self.c = c  
        self.grain_size = grain_size
        self.grain_rho = grain_rho
        self.Z = Z 
        self.stoke = stoke

        self.get_params(print_params=False)


    def get_params(self, print_params=True):
        """
        Calculates the disk parameters according to power law profiles.
        """

        self.calc_omega()
        self.calc_sigma_g() 
        self.calc_sigma_d()
        self.calc_T()
        self.calc_cs()
        self.calc_H()
        self.calc_h()
        self.calc_Q()
        self.calc_G()
        self.calc_vk()
        self.calc_big_pi()
        self.calc_beta()

        if self.stoke is None:
            if isinstance(self.grain_size, np.ndarray):
                raise ValueError('If stokes number is None, the grain size must contain only value!')
            self.calc_stokes()
            self.plot_stoke = True
        else:
            if isinstance(self.stoke, np.ndarray):
                raise ValueError('Only one stoke number can be input!')
            self.calc_grain_sizes()
            self.plot_stoke = False 

        if print_params:
            print('Sigma_g: {} g/cm2'.format(self.sigma_g))
            print('Sigma_d: {} g/cm2'.format(self.sigma_d))
            print('cs: {} cm/s'.format(self.cs))
            print('T : {} K'.format(self.T))
            print('H: {} AU'.format(self.H))
            print('h {}'.format(self.h))
            print('Q: {}'.format(self.Q))
            print('G: {}'.format(self.G))
            print('beta: {}'.format(self.beta))
            print('Stokes: {}'.format(self.stoke))
            print('Grain Radius: {}'.format(self.grain_size))

        return 

    def calc_sigma_g(self):
        """
        Gas surface mass density from the MMSN model (Hayashi 1981)
      
        Returns:
            Gas surface density profile [g/cm2]
        """

        self.sigma_g = 1700*(self.r/const.au.cgs.value)**-self.b

    def calc_sigma_d(self):
        """
        Calculates the dust surface density, according to the
        dust to gas ratio

        Returns:
            Dust surface density [g/cm2]
        """

        self.sigma_d = self.sigma_g * self.Z

    def calc_stokes(self):
        """
        Calculates the Stokes number according to the gas
        column density at each raddi and the corresponding
        grain properties

        Returns:
            Stoke's number.
        """

        self.stoke = np.pi * self.grain_size * self.grain_rho / 2. / self.sigma_g

    def calc_grain_sizes(self):
        """
        Calculates the grain sizes, to use if stoke's number is constant

        Returns:
            Grain size(s) [cm]
        """

        self.grain_size = 2 * self.stoke * self.sigma_g / np.pi / self.grain_rho

    def calc_cs(self):
        """
        Sound speed, assumes a mean molecular weight of 2.34

        Returns
            Sound speed in [cm/s]
        """

        #molecule_mass = self.mmol * const.m_p.cgs.value #Mass of hydrogen molecule
        #self.cs = np.sqrt(self.gamma * const.k_B.cgs.value * self.T / molecule_mass)
        self.cs = 1**(1./2)*(self.r/const.au.cgs.value)**(-self.c/2)
        self.cs = self.cs * 1e5 #km/s to cm/s

    def calc_omega(self):
        """
        Keplerian angular velocity.

        Returns:
            Keplerian velocity in [1/s]
        """

        #self.omega = np.sqrt(const.G.cgs.value * const.M_sun.cgs.value / self.r**3)
        self.omega = 2*np.pi*1**(1./2)*(self.r/const.au.cgs.value)**(-3./2)
        self.omega = self.omega / 3.156e7 #1/yr to 1/s conversion

    def calc_vk(self):
        """
        Rotational velocity

        Returns
            Velocity [cm/s]
        """

        self.vk = 30*1**(1./2)*(self.r/const.au.cgs.value)**(-1./2)
        self.vk = self.vk * 1e5 #km/s to cm/s

    def calc_Q(self):
        """
        Toomre Q parameter (dimensionless)

        Returns:
            ToomreQ
        """

        self.Q = self.cs * self.omega / np.pi / const.G.cgs.value / self.sigma_g

    def calc_T(self):
        """
        Radial temperature of the disk

        Args:
            f_t: Model parameter for the MSN model advocated by Hayashi (1981).
            c (float): Power law index, defaults to 3/2.

        Returns:
            Disk temperature T [K]
        """

        self.T = 280.0*(self.r/const.au.cgs.value)**-self.c

    def calc_H(self):
        """
        Gas scale height
    
        Returns:
            Gas scale height H [AU]
        """
        #self.H = self.cs / self.omega / const.au.cgs.value
        self.H = 3.4e-2*1**(1./2)*1**(-1./2)*(self.r/const.au.cgs.value)**((3-self.c)/2.)

    def calc_h(self):
        """
        Calculates the aspect ratio, h

        Returns:
            Aspect ratio.
        """
        
        self.h = self.H / (self.r/const.au.cgs.value)

    def calc_G(self):
        """
        G tilde self-gravity parameter.

        Returns:
            G Tilde parameter.
        """

        self.G = np.sqrt(8/np.pi) / self.Q

    def calc_big_pi(self):
        """
        Big pi calculated assuming solar parameters and MMSN
        power law indices. The z-dependence is ignored since z << Hg 

        Bai & Stone :https://arxiv.org/pdf/1005.4982.pdf
        """

        self.big_pi = 0.054*1**(1./2)*1**(-1./2)*(self.r/const.au.cgs.value)**(1./4)
        #self.big_pi = (3+2*self.b+self.c)/4. * self.cs/self.vk

    def calc_beta(self):
        """
        Calculates pressure gradient parameter
        """

        self.beta = 2.0 * self.big_pi 

    def plot(self, savefig=False, path=None, box_index=None):
        """
        Plots the model parameters.

        Args:
            save (bool): If True the figure will be saved only. Defaults to False.
        Returns:
            AxesPlot 
        """

        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(18,12))
        fig.suptitle("Protoplanetary Disk Model", fontsize=24, x=.5, y=.92)
        ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = axes

        ax1.plot(self.r/const.au.cgs.value, self.sigma_g, c='k')
        if box_index is not None:
            ax1.scatter(self.r[box_index]/const.au.cgs.value, self.sigma_g[box_index], marker='s', s=100)
        ax1.set_ylabel(r'$\Sigma_g \ [\rm g \ \rm cm^{-2}]$', fontsize=18)
        ax1.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax1.tick_params(labelsize=14, axis="both", which="both")
        ax1.set_xticklabels([])

        ax2.plot(self.r/const.au.cgs.value, self.T, c='k')
        if box_index is not None:
            ax2.scatter(self.r[box_index]/const.au.cgs.value, self.T[box_index], marker='s', s=100)
        ax2.set_ylabel('T [K]', size=18)
        ax2.tick_params(labelsize=14, axis="both", which="both")
        ax2.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax2.set_xticklabels([])

        ax3.plot(self.r/const.au.cgs.value, self.H/const.au.cgs.value, c='k')
        if box_index is not None:
            ax3.scatter(self.r[box_index]/const.au.cgs.value, self.H[box_index]/const.au.cgs.value, marker='s', s=100)
        ax3.set_ylabel('H [AU]', fontsize=18)
        ax3.tick_params(labelsize=14, axis="both", which="both")
        ax3.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax3.set_xticklabels([])
        """
        ax4.plot(self.r/const.au.cgs.value, self.sigma_d, c='k')
        ax4.set_ylabel(r'$\Sigma_d \ [g \ cm^{-2}]$', fontsize=18)
        ax4.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax4.tick_params(labelsize=14, axis="both", which="both")
        ax4.set_xticklabels([])
        
        ax5.plot(self.r/const.au.cgs.value, self.cs*1e-5, c='k')
        ax5.set_ylabel(r'$c_s \ [km \ s^{-1}]$', size=18)
        ax5.tick_params(labelsize=14, axis="both", which="both")
        ax5.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax5.set_xticklabels([])
        """
        ax4.plot(self.r/const.au.cgs.value, self.Q, c='k')
        if box_index is not None:
            ax4.scatter(self.r[box_index]/const.au.cgs.value, self.Q[box_index], marker='s', s=100)
        ax4.set_ylabel('Toomre Q', fontsize=18)
        ax4.tick_params(labelsize=14, axis="both", which="both")
        ax4.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax4.set_xticklabels([])

        if self.plot_stoke:
            ax5.plot(self.r/const.au.cgs.value, self.stoke, c='k')
            if box_index is not None:
                ax5.scatter(self.r[box_index]/const.au.cgs.value, self.stoke[box_index], marker='s', s=100)
            ax5.set_ylabel('Stoke', fontsize=18)
        else:
            ax5.plot(self.r/const.au.cgs.value, self.grain_size, c='k')
            if box_index is not None:
                ax5.scatter(self.r[box_index]/const.au.cgs.value, self.grain_size[box_index], marker='s', s=100)
            ax5.set_ylabel('Grain Radius [cm]', fontsize=18)
        ax5.tick_params(labelsize=14, axis="both", which="both")
        ax5.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax5.set_xlabel('Radius [AU]', size=18)

        ax6.plot(self.r/const.au.cgs.value, self.beta/2., c='k')
        if box_index is not None:
            ax6.scatter(self.r[box_index]/const.au.cgs.value, self.beta[box_index]/2., marker='s', s=100)
        ax6.set_ylabel(r'$\prod$', fontsize=18)
        ax6.set_xlabel('Radius [AU]', size=18)
        ax6.tick_params(labelsize=14, axis="both", which="both")
        ax6.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')

        if savefig:
            if path is None:
                path = str(Path.home())
            if path[-1] != '/':
                path+='/'
            _set_style_()
            plt.savefig(path+'Disk_Model.png', dpi=300, bbox_inches='tight')
            print('Figure saved in: {}'.format(path))
            plt.clf(); plt.style.use('default')
        else:
            plt.show()


class Model_2:
    """
    Creates the disk model.

    The gas surface density is defined in Section 2 of Drazkowska et al. (2022)
    See: https://arxiv.org/pdf/2101.01728.pdf

    The temperature model is from Ida et al. (2016) 
    See: https://arxiv.org/pdf/1604.01291.pdf

    Note that the temperature profile is appropriate for outer regions of the disk where viscous heating is negligible. 

    Note:
        All inputs must be in cgs units!

    Args:
        r (float, ndarray): The radii of the disk at which to calculate
            the gas surface density. Must be in cgs units.
        r_c (float): Characteristic radius of the disk, in cgs units. This radius
            for compact disks is about 30 AU, and 300 AU for large disks.
        M_star (float): Mass of the central star, in cgs units.
        M_disk (float): Mass of the disk, in cgs units. 
        grain_size (float): Grain size, in cgs units. Only use if no stoke's number
            is input, as that will overwrite this value. Defaults to 1 mm.
        grain_rho (float): Internal grain density in cgs units. Defaults to 1 g/cm2.
        T0 (float): Temperature at r=0, defaults to 150 K, which is an accepted standard.
        q (float): Power law index, defaults to 3/7, which is an accepted standard.
        gamma (float): Adiabatic index, defaults to 1.0 for an isothermal disk.
        mmol (float): Mean molecular weight, defaults to 2.3, corresponding
            to 5 parts molecular H and 1 part He1.
        Z (float): The global solids-to-gas ratio. Defaults to 0.01.
        st (float): Stoke's number, if input then the grain sizes at each
            r will be computed. Defaults to None. 
    
    Attributes:
        get_params: Calculates all disk parameters
    """
    def __init__(self, r, r_c, M_star, M_disk, grain_size=0.1, grain_rho=1,
        T0=150, q=3./7, gamma=1, mmol=2.3, Z=0.01, stoke=None):

        self.r = r 
        self.r_c = r_c
        self.M_star = M_star
        self.M_disk = M_disk
        self.grain_size = grain_size
        self.grain_rho = grain_rho
        self.T0 = T0
        self.q = q 
        self.gamma = gamma 
        self.mmol = mmol 
        self.Z = Z 
        self.stoke = stoke

        self.sigma_g = None
        self.sigma_d = None  
        self.cs = None 
        self.T = None 
        self.H = None 
        self.h = None 
        self.Q = None 
        self.G = None 
        self.beta = None 

        self.get_params(print_params=False)


    def get_params(self, print_params=True):
        """
        Calculates the disk parameters according to power law profiles.
        """

        self.calc_omega()
        self.calc_sigma_g() 
        self.calc_sigma_d()
        self.calc_T()
        self.calc_cs()
        self.calc_H()
        self.calc_h()
        self.calc_Q()
        self.calc_G()
        self.calc_beta()

        if self.stoke is None:
            if isinstance(self.grain_size, np.ndarray):
                raise ValueError('If stokes number is None, the grain size must contain only value!')
            self.calc_stokes()
            self.plot_stoke = True
        else:
            if isinstance(self.stoke, np.ndarray):
                raise ValueError('Only one stoke number can be input!')
            self.calc_grain_sizes()
            self.plot_stoke = False 

        if print_params:
            print('Sigma_g: {} g/cm2'.format(self.sigma_g))
            print('Sigma_d: {} g/cm2'.format(self.sigma_d))
            print('cs: {} cm/s'.format(self.cs))
            print('T : {} K'.format(self.T))
            print('H: {} AU'.format(self.H/const.au.cgs.value))
            print('h {}'.format(self.h))
            print('Q: {}'.format(self.Q))
            print('G: {}'.format(self.G))
            print('beta: {}'.format(self.beta))
            print('Stokes: {}'.format(self.stoke))
            print('Grain Radius: {}'.format(self.grain_size))

    def calc_sigma_g(self):
        """
        Calculates the gas surface density, as per the protoplanetary 
        disk model presented in Section 2 of Drazkowska et al (2022)
        See: https://arxiv.org/pdf/2101.01728.pdf

        Returns:
            Gas surface density at the specified r.
        """

        self.sigma_g = ((self.M_disk/(2.*np.pi*self.r_c**2))*(self.r/self.r_c)**-1)*np.e**(-(self.r/self.r_c))

    def calc_sigma_d(self):
        """
        Calculates the dust surface density, as per the protoplanetary 
        disk model presented in Section 2 of Drazkowska et al (2022)
        See: https://arxiv.org/pdf/2101.01728.pdf

        Returns:
            Gas surface density.
        """

        self.sigma_d = self.sigma_g * self.Z

    def calc_stokes(self):
        """
        Calculates the Stokes number according to the gas
        column density and the grain properties
        Returns:
            Stoke's number.
        """

        self.stoke = np.pi * self.grain_size * self.grain_rho / 2. / self.sigma_g

    def calc_grain_sizes(self):
        """
        Calculates the grain sizes, to use if stoke's number is constant
        
        Returns:
            Grain sizes, a, two times the grain radius.
        """

        self.grain_size = 2 * self.stoke * self.sigma_g / np.pi / self.grain_rho

    def calc_omega(self):
        """
        Keplerian angular velocity.
        """

        self.omega = np.sqrt(const.G.cgs.value*self.M_star/self.r**3)

    def calc_Q(self):
        """
        Toomre Q parameter (dimensionless)

        Returns:
            ToomreQ parameter.
        """

        self.Q = self.cs * self.omega / np.pi / const.G.cgs.value / self.sigma_g

    def calc_T(self):
        """
        Temperature model from Ida et al. 2016: https://arxiv.org/pdf/1604.01291.pdf
        
        Returns:
            Temperature as a function of distance, in cgs units.
        """

        self.T = self.T0*(self.r / const.au.cgs.value)**(-self.q)

    def calc_h(self):
        """
        Calculates the aspect ratio, h

        Returns:
            Aspect ratio.
        """
        
        self.h = self.H / self.r

    def calc_cs(self):
        """
        Calculates the sound speed.

        Returns
            cs [m/s]
        """

        molecule_mass = self.mmol * const.m_p.cgs.value #Mass of hydrogen molecule

        self.cs = np.sqrt(self.gamma * const.k_B.cgs.value * self.T / molecule_mass)

    def calc_H(self):
        """
        Calculates scale height
        """

        self.H = self.cs / self.omega

    def calc_G(self):
        """
        G tilde self-gravity parameter.

        Returns:
            Toomre Q parameter
        """

        self.G = np.sqrt(8/np.pi) / self.Q

    def calc_beta(self):
        """
        Calculates the beta parameter
        beta = -h * dlnp/dlnr
        dlnp/dlnr = dlnSigma/dlnr - 0.5*dlnT/dlnr + dlnOmega/dlnr = -1 -0.5*-3/7 - 3/2 = -2.2857
        
        Returns:
            beta
        """

        self.beta = -self.h * -2.0# -2.2857

    def plot(self, savefig=False, include_grid=False, path=None, box_index=None):
        """
        Plots the model parameters.

        Args:
            savefig (bool): If True the figure will be saved only. Defaults to False.
            include_grid (bool): Whether to include the plot gird. Defaults to False.
            path (str):
            box_index (int):

        Returns:
            AxesPlot 
        """
        plt.style.use('/Users/daniel/Documents/plot_style.txt')
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,8))
        fig.suptitle("Protoplanetary Disk Model", fontsize=24, x=.5, y=.97)
        ((ax1, ax2), (ax3, ax4)) = axes

        ax1.plot(self.r/const.au.cgs.value, self.sigma_g, c='k', label='Gas')#r'$\Sigma_g$')
        ax1.plot(self.r/const.au.cgs.value, self.sigma_d, c='k', linestyle='--', label='Dust')#r'$\Sigma_d$')
        if box_index is not None:
            ax1.scatter(self.r[box_index]/const.au.cgs.value, self.sigma_g[box_index], marker='s', s=100)
        ax1.set_ylabel(r'$\Sigma$ $[\rm g \ \rm cm^{-2}]$', fontsize=18)
        if include_grid:
            ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
            ax1.yaxis.set_major_locator(plt.MultipleLocator(1))
            ax1.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
            ax1.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
            ax1.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        ax1.tick_params(labelsize=14, axis="both", which="both")
        ax1.set_xlim((5,100)), ax1.set_ylim((1e-2,30))#, ax1.set_xlim((5,100))
        ax1.set_xticklabels([])
        ax1.legend(frameon=False, handlelength=1, loc='upper right', ncol=1, prop={'size': 18})
        ax1.set_yscale('log')#, ax1.set_xscale('log')

        ax2.plot(self.r/const.au.cgs.value, self.T, c='k')
        if box_index is not None:
            ax2.scatter(self.r[box_index]/const.au.cgs.value, self.T[box_index], marker='s', s=100)
        ax2.set_ylabel('T [K]', size=18)
        ax2.tick_params(labelsize=14, axis="both", which="both")
        if include_grid:
            ax2.xaxis.set_major_locator(plt.MultipleLocator(1))
            ax2.yaxis.set_major_locator(plt.MultipleLocator(1))
            ax2.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
            ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
            ax2.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        ax2.set_xticklabels([])
        ax2.set_xlim((5,100)), ax2.set_ylim((5,130))
        ax2.set_yscale('log')#, ax2.set_xscale('log')

        if self.plot_stoke:
            ax3.plot(self.r/const.au.cgs.value, self.stoke, c='k')
            if box_index is not None:
                ax3.scatter(self.r[box_index]/const.au.cgs.value, self.stoke[box_index], marker='s', s=100)
            ax3.set_ylabel('Stoke', fontsize=18)
        else:
            ax3.plot(self.r/const.au.cgs.value, self.grain_size, c='k')
            if box_index is not None:
                ax3.scatter(self.r[box_index]/const.au.cgs.value, self.grain_size[box_index], marker='s', s=100)
            ax3.set_ylabel('a [cm]', fontsize=18)
        ax3.tick_params(labelsize=14, axis="both", which="both")
        if include_grid:
            ax3.xaxis.set_major_locator(plt.MultipleLocator(1))
            ax3.yaxis.set_major_locator(plt.MultipleLocator(1))
            ax3.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
            ax3.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
            ax3.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        ax3.set_xlabel('Radius [AU]', size=18)
        ax3.set_xlim((5,100)), ax3.set_ylim((1e-2,0.4))
        ax3.set_yscale('log')#, ax3.set_xscale('log')

        ax4.plot(self.r/const.au.cgs.value, self.h, c='k')
        if box_index is not None:
            ax4.scatter(self.r[box_index]/const.au.cgs.value, self.h[box_index], marker='s', s=100)
        ax4.set_ylabel('H / r', fontsize=18)
        ax4.tick_params(labelsize=14, axis="both", which="both")
        if include_grid:
            ax4.xaxis.set_major_locator(plt.MultipleLocator(1))
            ax4.yaxis.set_major_locator(plt.MultipleLocator(1))
            ax4.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
            ax4.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
            ax4.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        ax4.set_xlabel('Radius [AU]', size=18)
        ax4.set_xlim((5,100)), ax4.set_ylim((0, 0.1)) #ax4.set_ylim((10,320))
        #ax4.set_yscale('log')#, ax4.set_xscale('log')

        if savefig:
            if path is None:
                path = str(Path.home())
            if path[-1] != '/':
                path+='/'
            _set_style_()
            plt.savefig(path+'Disk_Model.png', dpi=300, bbox_inches='tight')
            print('Figure saved in: {}'.format(path))
            plt.clf(); plt.style.use('default')
        else:
            plt.show()

        return 

def _set_style_():
    """
    Function to configure the matplotlib.pyplot style. This function is called before any images are saved,
    after which the style is reset to the default.
    """

    plt.rcParams["xtick.color"] = "323034"
    plt.rcParams["ytick.color"] = "323034"
    plt.rcParams["text.color"] = "323034"
    plt.rcParams["lines.markeredgecolor"] = "black"
    plt.rcParams["patch.facecolor"] = "#bc80bd"  # Replace with a valid color code
    plt.rcParams["patch.force_edgecolor"] = True
    plt.rcParams["patch.linewidth"] = 0.8
    plt.rcParams["scatter.edgecolors"] = "black"
    plt.rcParams["grid.color"] = "#b1afb5"  # Replace with a valid color code
    plt.rcParams["axes.titlesize"] = 16
    plt.rcParams["legend.title_fontsize"] = 12
    plt.rcParams["xtick.labelsize"] = 16
    plt.rcParams["ytick.labelsize"] = 16
    plt.rcParams["font.size"] = 15
    plt.rcParams["axes.prop_cycle"] = (cycler('color', ['#bc80bd', '#fb8072', '#b3de69', '#fdb462', '#fccde5', '#8dd3c7', '#ffed6f', '#bebada', '#80b1d3', '#ccebc5', '#d9d9d9']))  # Replace with valid color codes
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["font.family"] = "STIXGeneral"
    plt.rcParams["lines.linewidth"] = 2
    plt.rcParams["lines.markersize"] = 6
    plt.rcParams["legend.frameon"] = True
    plt.rcParams["legend.framealpha"] = 0.8
    plt.rcParams["legend.fontsize"] = 13
    plt.rcParams["legend.edgecolor"] = "black"
    plt.rcParams["legend.borderpad"] = 0.2
    plt.rcParams["legend.columnspacing"] = 1.5
    plt.rcParams["legend.labelspacing"] = 0.4
    plt.rcParams["text.usetex"] = False
    plt.rcParams["axes.labelsize"] = 17
    plt.rcParams["axes.titlelocation"] = "center"
    plt.rcParams["axes.formatter.use_mathtext"] = True
    plt.rcParams["axes.autolimit_mode"] = "round_numbers"
    plt.rcParams["axes.labelpad"] = 3
    plt.rcParams["axes.formatter.limits"] = (-4, 4)
    plt.rcParams["axes.labelcolor"] = "black"
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1
    plt.rcParams["axes.grid"] = False
    plt.rcParams["axes.spines.right"] = True
    plt.rcParams["axes.spines.left"] = True
    plt.rcParams["axes.spines.top"] = True
    plt.rcParams["figure.titlesize"] = 18
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams["figure.dpi"] = 300

    return
 
"""
def aspect_ratio(r, M, mmol=2.3):
    Aspect ratio, as a function of T 
    Equation 5: https://www.aanda.org/articles/aa/pdf/2015/03/aa24964-14.pdf
    return np.sqrt(T(r) * r * const.R.cgs.value / const.G.cgs.value / M / mmol )

def sound_speed(r, gamma=1.4, mmol=2.3):
    Rgasmu = const.R.cgs.value / mmol
    cp = gamma*Rgasmu / (gamma-1)
    return np.sqrt(T(r)*cp*(gamma-1))
"""

"""

#Disk model plot

import numpy as  np 
import matplotlib.pyplot as plt  
import astropy.constants as const
from StreamingInstability_YJ14 import disk_model

M_star, M_disk = const.M_sun.cgs.value, 0.02*const.M_sun.cgs.value
r, r_c = np.arange(5,101,1), 300
#r = 0.73
#r = 7.5
r, r_c = r*const.au.cgs.value, r_c*const.au.cgs.value

grain_rho = 1.675
stoke=0.05
Z = 0.02 
q = 3/7.
T0 = 150

q, T0 = 1., 600.


model = disk_model.Model_2(r, r_c, M_star, M_disk, grain_rho=grain_rho, Z=Z, stoke=stoke, q=q, T0=T0)


model.plot(savefig=True)

######

grain_rho = 1.675
Z = 0.03
q = 3/7.
T0 = 150

stoke = np.zeros((4,3))
sigma_g = np.zeros((4,3))
beta = np.zeros((4,3))
G = np.zeros((4,3))
i=0

for radius in [10, 30, 100]:
    r, r_c = radius*const.au.cgs.value, 300*const.au.cgs.value
    j=0
    for grain_size in [1, 0.3, 0.1, 0.03]:
        model = disk_model.Model_2(r, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=grain_size, Z=Z, stoke=None, q=q, T0=T0)
        stoke[j,i] = model.stoke
        sigma_g[j,i] = model.sigma_g
        beta[j,i] = model.beta
        G[j,i] = model.G
        j+=1
    i+=1

#Stokes number at different sigma_g and grain size

"""



