#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 09:25:13 2022

@author: daniel
"""

import numpy as  np 
from pathlib import Path
import matplotlib.pyplot as plt  
import astropy.constants as const  


class Model:
    """
    Creates the disk model.

    The gas surface density is defined in Section 2 of Drazkowska et al (2022)
    See: https://arxiv.org/pdf/2101.01728.pdf

    The temperature model from Ida et al. , see:2016: https://arxiv.org/pdf/1604.01291.pdf
    Note that the temperature profile is for outer regions of the disk only 
    where viscous heating is negligible. 

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

        return 

    def calc_sigma_g(self):
        """
        Calculates the gas surface density, as per the protoplanetary 
        disk model presented in Section 2 of Drazkowska et al (2022)
        See: https://arxiv.org/pdf/2101.01728.pdf

        Returns:
            Gas surface density at the specified r.
        """

        #M_disk = self.M_disk*const.M_sun.cgs.value
        #r_c = self.r_c*const.au.cgs.value

        self.sigma_g = ((self.M_disk/(2.*np.pi*self.r_c**2))*(self.r/self.r_c)**-1)*np.e**(-(self.r/self.r_c))

        return 

    def calc_sigma_d(self):
        """
        Calculates the dust surface density, as per the protoplanetary 
        disk model presented in Section 2 of Drazkowska et al (2022)
        See: https://arxiv.org/pdf/2101.01728.pdf

        Returns:
            Gas surface density.
        """

        self.sigma_d = self.sigma_g * self.Z

        return  

    def calc_stokes(self):
        """
        Calculates the Stokes number according to the gas
        column density and the grain properties

        Returns:
            Stoke's number.
        """

        self.stoke = np.pi * self.grain_size * self.grain_rho / 2. / self.sigma_g

        return

    def calc_grain_sizes(self):
        """
        Calculates the grain sizes, to use if stoke's number is constant

        Returns:
            Grain sizes, a.
        """

        self.grain_size = 2 * self.stoke * self.sigma_g / np.pi / self.grain_rho

    def calc_omega(self):
        """
        Keplerian angular velocity.

        """

        self.omega = np.sqrt(const.G.cgs.value*self.M_star/self.r**3)

        return

    def calc_Q(self):
        """
        Toomre Q parameter (dimensionless)

        Returns:
            ToomreQ
        """

        self.Q = self.cs * self.omega / np.pi / const.G.cgs.value / self.sigma_g

        return

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

        return 

    def calc_cs(self):
        """
        Calculates the sound speed

        Returns
            cs
        """

        molecule_mass = self.mmol * const.m_p.cgs.value #Mass of hydrogen molecule

        self.cs = np.sqrt(self.gamma * const.k_B.cgs.value * self.T / molecule_mass)

        return 

    def calc_H(self):
        """
        Calculates scale height
        """

        self.H = self.cs / self.omega

        return

    def calc_G(self):
        """
        G tilde self-gravity parameter.

        Returns:
            Toomre Q parameter
        """

        self.G = np.sqrt(8/np.pi) / self.Q

        return 

    def calc_beta(self):
        """
        Calculates the beta parameter
        beta = -h * dlnp/dlnr
        dlnp/dlnr = dlnSigma/dlnr - 0.5*dlnT/dlnr + dlnOmega/dlnr = -1 -0.5*-3/7 - 3/2 = -2.2857

        Returns:
            beta
        """

        self.beta = -self.h * -2.2857

        return 

    def plot(self, savefig=False, path=None):
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
        ax1.set_ylabel(r'$\Sigma_g \ [\rm g \ \rm cm^{-2}]$', fontsize=18)
        ax1.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax1.tick_params(labelsize=14, axis="both", which="both")
        ax1.set_xticklabels([])

        ax2.plot(self.r/const.au.cgs.value, self.T, c='k')
        ax2.set_ylabel('T [K]', size=18)
        ax2.tick_params(labelsize=14, axis="both", which="both")
        ax2.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax2.set_xticklabels([])

        ax3.plot(self.r/const.au.cgs.value, self.H/const.au.cgs.value, c='k')
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
        ax4.set_ylabel('Toomre Q', fontsize=18)
        ax4.tick_params(labelsize=14, axis="both", which="both")
        ax4.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax4.set_xticklabels([])

        if self.plot_stoke:
            ax5.plot(self.r/const.au.cgs.value, self.stoke, c='k')
            ax5.set_ylabel('Stoke', fontsize=18)
        else:
            ax5.plot(self.r/const.au.cgs.value, self.grain_size, c='k')
            ax5.set_ylabel('Grain Radius [cm]', fontsize=18)
        ax5.tick_params(labelsize=14, axis="both", which="both")
        ax5.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
        ax5.set_xlabel('Radius [AU]', size=18)

        ax6.plot(self.r/const.au.cgs.value, self.beta/2., c='k')
        ax6.set_ylabel(r'$\prod$', fontsize=18)
        ax6.set_xlabel('Radius [AU]', size=18)
        ax6.tick_params(labelsize=14, axis="both", which="both")
        ax6.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')

        if savefig:
            if path is None:
                path = str(Path.home())
            if path[-1] != '/':
                path+='/'
            plt.savefig(path+'Disk_Model.png', dpi=300, bbox_inches='tight')
            print('Figure saved in: {}'.format(path))
            plt.clf()
        else:
            plt.show()

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
r, r_c = np.arange(30,101,1), 300
r, r_c = r*const.au.cgs.value, r_c*const.au.cgs.value

model = disk_model.Model(r, r_c, M_star, M_disk, grain_size=0.1, grain_rho=1, Z=0.03, stoke=0.3)
model.plot(savefig=True)

stoke = np.zeros((4,3))
sigma_g = np.zeros((4,3))
beta = np.zeros((4,3))
G = np.zeros((4,3))
i=0
for radius in [10, 30, 100]:
    r, r_c = radius*const.au.cgs.value, 300*const.au.cgs.value
    j=0
    for grain_size in [1, 0.3, 0.1, 0.03]:
        model = disk_model.Model(r, r_c, M_star, M_disk, grain_size=grain_size, grain_rho=1, Z=0.02)
        stoke[j,i] = model.stoke
        sigma_g[j,i] = model.sigma_g
        beta[j,i] = model.beta
        G[j,i] = model.G
        j+=1
    i+=1

#Stokes number at different sigma_g and grain size

"""



