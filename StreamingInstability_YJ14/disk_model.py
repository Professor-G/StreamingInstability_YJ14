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
    Generates a 1D protoplanetary disk model (midplane).

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

        # These are the genereated disk parameters
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
        Calculates the disk parameters according.

        Args:
            print_params (bool): If set to True the parameters will be printed. Defaults to True.
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
            Grain size(s), a.
        """

        self.grain_size = 2 * self.stoke * self.sigma_g / np.pi / self.grain_rho

    def calc_omega(self):
        """
        Keplerian angular velocity.

        Returns:
            Omega
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

        Note:
            This temperature profile is only physical for the outer regions of the disk where viscous heating is negligible. 
        
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
        Calculates the gas scale height.

        Returns:
            The gas scale height, H.
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
        beta = h * dlnp/dlnr
        dlnp/dlnr = dlnSigma/dlnr - 0.5*dlnT/dlnr + dlnOmega/dlnr = -1 -0.5*-3/7 - 3/2 = -2.2857
        
        Returns:
            Beta parameter
        """

        self.beta = self.h * (-1.0 - (0.5*-self.q) - 1.5) #-2.0# -2.2857

    def plot(self, xticks=[5, 20, 40, 60, 80, 100], ax1_xlim=None, ax1_ylim=None, ax2_xlim=None, ax2_ylim=None,
        ax3_xlim=None, ax3_ylim=None, ax4_xlim=None, ax4_ylim=None, ax1_ylog=True, ax2_ylog=True, ax3_ylog=True, ax4_ylog=False, 
        include_grid=False, savefig=False, path=None, box_index=None, plot_vertical=False, title="Protoplanetary Disk Model"):
        """
        Plots four model parameters: gas and dust density profiles, temperature profile, grain size or stokes number, and aspect ratio.

        Args:
            xticks (list): The x-axis tick labels to display (the disk radii)
            ax1_xlim (tuple, optional): The x-limits of the first subplot (e.g. ax1_xlim=(5,100)). Defaults to None.
            ax1_ylim (tuple, optional): The y-limits of the first subplot. Defaults to None.
            ax2_xlim (tuple, optional): The x-limits of the second subplot. Defaults to None.
            ax2_ylim (tuple, optional): The y-limits of the second subplot. Defaults to None.
            ax3_xlim (tuple, optional): The x-limits of the third subplot. Defaults to None.
            ax3_ylim (tuple, optional): The y-limits of the third subplot. Defaults to None.
            ax4_xlim (tuple, optional): The x-limits of the fourth subplot. Defaults to None.
            ax4_ylim (tuple, optional): The y-limits of the fourth subplot. Defaults to None.
            ax1_ylog (bool): Whether the y-values of the first subplot should be log-scaled. Defaults to True.
            ax2_ylog (bool): Whether the y-values of the second subplot should be log-scaled. Defaults to True.
            ax3_ylog (bool): Whether the y-values of the third subplot should be log-scaled. Defaults to True.
            ax4_ylog (bool): Whether the y-values of the fourth subplot should be log-scaled. Defaults to False.
            include_grid (bool): Whether to include a plot grid. Defaults to False.
            savefig (bool): If True the figure will be saved and will not display. Defaults to False.
            path (str, optional): Path to where the figure should be saved (only applicable if savefig=True). Defaults to None,
                in which case the figure will be saved to the local home directory.
            box_index (int, optional): The location of the box within the disk model. Must be an integer index corresponding to the 
                location with respect to the radius of the disk (r attribute). Defaults to None.
            plot_vertical (bool): If set to True, the four subplots will be in vertical format. Defaults to False
                which returns a 2x2 display.
            title (str): The subtitle of the plot.
            
        Returns:
            AxesPlot 
        """

        if isinstance(self.r, np.ndarray):
            if len(self.r) <= 1:
                print('WARNING: The radius array (r) is too small for proper visualization!')
        else:
            print('WARNING: The radius array (r) is too small for proper visualization!')

        _set_style_() if savefig else plt.style.use('default')

        if plot_vertical:
            fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(5, 9), sharex=True)
            fig.suptitle(title, x=0.565, y=0.965)
            ax1, ax2, ax3, ax4 = axes
        else:
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,8))
            fig.suptitle(title, x=.5, y=.97)
            ((ax1, ax2), (ax3, ax4)) = axes

        ### ax1 ###

        ax1.plot(self.r/const.au.cgs.value, self.sigma_g, c='k', label='Gas')
        ax1.plot(self.r/const.au.cgs.value, self.sigma_d, c='k', linestyle='--', label='Dust')

        if box_index is not None: ax1.scatter(self.r[box_index]/const.au.cgs.value, self.sigma_g[box_index], marker='s', s=100)
        if include_grid: ax1.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        if ax1_xlim is not None: ax1.set_xlim(ax1_xlim) #(5, 100)
        if ax1_ylim is not None: ax1.set_ylim(ax1_ylim) #(1e-2, 30)
        if ax1_ylog: ax1.set_yscale('log')
        
        ax1.set_xticklabels([])
        ax1.set_ylabel(r'$\Sigma$ $[\rm g \ \rm cm^{-2}]$')
        ax1.legend(frameon=False, handlelength=1, loc='upper right', ncol=1)
        
        ### ax2 ###

        ax2.plot(self.r/const.au.cgs.value, self.T, c='k')

        if box_index is not None: ax2.scatter(self.r[box_index]/const.au.cgs.value, self.T[box_index], marker='s', s=100)
        if include_grid: ax2.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        if ax2_xlim is not None: ax1.set_xlim(ax2_xlim) #(5, 100)
        if ax2_ylim is not None: ax1.set_ylim(ax2_ylim) #(5, 250)
        if ax2_ylog: ax2.set_yscale('log')

        ax2.set_xticklabels([])
        ax2.set_ylabel('T [K]')
        
        ### ax3 ###

        if self.plot_stoke:

            ax3.plot(self.r/const.au.cgs.value, self.stoke, c='k')

            if box_index is not None: ax3.scatter(self.r[box_index]/const.au.cgs.value, self.stoke[box_index], marker='s', s=100)
            ax3.set_ylabel(r'$\tau_s$')

        else:

            ax3.plot(self.r/const.au.cgs.value, self.grain_size*10, c='k')
            
            if box_index is not None: ax3.scatter(self.r[box_index]/const.au.cgs.value, self.grain_size[box_index]*10, marker='s', s=100)
            ax3.set_ylabel('a [mm]')

        if include_grid: ax3.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        if ax3_xlim is not None: ax1.set_xlim(ax3_xlim) #(5, 100)
        if ax3_ylim is not None: ax1.set_ylim(ax3_ylim) #(5, 250)
        if ax3_ylog: ax3.set_yscale('log')

        if plot_vertical:
            ax3.set_xticklabels([])
        else:
            ax3.set_xticks(xticks); ax3.set_xticklabels(xticks)
            ax3.set_xlabel('Radius [au]')

        ### ax4 ###

        ax4.plot(self.r/const.au.cgs.value, self.h, c='k')

        if box_index is not None: ax4.scatter(self.r[box_index]/const.au.cgs.value, self.h[box_index], marker='s', s=100)
        if include_grid: ax4.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        if ax4_xlim is not None: ax1.set_xlim(ax4_xlim) #(5, 100)
        if ax4_ylim is not None: ax1.set_ylim(ax4_ylim) #(0.02, 0.08)
        if ax4_ylog: ax4.set_yscale('log')

        ax4.set_xticks(xticks); ax4.set_xticklabels(xticks)
        ax4.set_ylabel('H / r'); ax4.set_xlabel('Radius [au]')

        if savefig:
            path = str(Path.home()) if path is None else path
            path = path+'/' if path[-1] != '/' else path 
            plt.savefig(path+'Disk_Model.png', dpi=300, bbox_inches='tight')
            print('Figure saved in: {}'.format(path))
            plt.clf()
        else:
            plt.tight_layout()
            plt.show()

        return

def _set_style_():
    """Function to configure the matplotlib.pyplot style. This function is called before any images are saved, after which the style is reset to the default.
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
    #
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
import astropy.constants as const

plt.style.use('/Users/daniel/Downloads/plot_style.txt')

## Disk models for the simulations with self-gravity (monodisperse) ###

M_star = const.M_sun.cgs.value # Mast of the star 
r, r_c = np.arange(5,100.25,0.25), 300 # Radii which to model, and the characteristic radius of the disk (in [au])
r, r_c = r*const.au.cgs.value, r_c*const.au.cgs.value # Convert to cgs units 

grain_rho = 1.675 # Internal dust grain density (from DSHARP)
stoke = 0.314 # Stokes number of the grain
Z = 0.02 # Dust to gas ratio
q = 1.0 # Temperature power law index
T0 = 600 # Temperature at r = 1 au

M_disk = 0.01*const.M_sun.cgs.value
model_1a = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, Z=Z, stoke=stoke, q=q, T0=T0)

M_disk = 0.03*const.M_sun.cgs.value
model_1b = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, Z=Z, stoke=stoke, q=q, T0=T0)

M_disk = 0.05*const.M_sun.cgs.value
model_1c = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, Z=Z, stoke=stoke, q=q, T0=T0)

M_disk = 0.1*const.M_sun.cgs.value
model_1d = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, Z=Z, stoke=stoke, q=q, T0=T0)

## Disk model for the simulations without self-gravity (polydisperse) ###

grain_rho = 1.675 # Internal dust grain density (from DSHARP)
Z = 0.03    # Dust to gas ratio
q = 3/7. # Temperature power law index
T0 = 150 # Temperature at r = 1 au

M_disk = 0.01675 * const.M_sun.cgs.value

model_2a = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=1, Z=Z, stoke=None, q=q, T0=T0)
model_2b = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=0.3, Z=Z, stoke=None, q=q, T0=T0)
model_2c = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=0.1, Z=Z, stoke=None, q=q, T0=T0)
model_2d = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=0.03, Z=Z, stoke=None, q=q, T0=T0)

### Plot ###

fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(14, 12.5), sharex=True)
fig.suptitle("Protoplanetary Disk Models", x=0.51, y=0.975)

(ax1, ax5), (ax2, ax6), (ax3, ax7), (ax9, ax10), (ax4, ax8) = axes

# Plot the gas (red) on the right axis
ax1.plot(model_1a.r/const.au.cgs.value, model_1a.sigma_g, c='red', linestyle='-')#, label=r'$\rm M_{\rm disk} = 0.01 \, M_{\odot}$')
ax1.plot(model_1b.r/const.au.cgs.value, model_1b.sigma_g, c='red', linestyle='--')#, label=r'$\rm M_{\rm disk} = 0.03 \, M_{\odot}$')
ax1.plot(model_1c.r/const.au.cgs.value, model_1c.sigma_g, c='red', linestyle=':')#, label=r'$\rm M_{\rm disk} = 0.05 \, M_{\odot}$')
ax1.plot(model_1d.r/const.au.cgs.value, model_1d.sigma_g, c='red', linestyle='-.')#, label=r'$\rm M_{\rm disk} = 0.05 \, M_{\odot}$')

line_for_legend1, = ax1.plot([], [], c='black', linestyle='-', label=r'$\rm M_{\rm disk} = 0.01 \, M_{\odot}$')
line_for_legend2, = ax1.plot([], [], c='black', linestyle='--', label=r'$\rm M_{\rm disk} = 0.03 \, M_{\odot}$')
line_for_legend3, = ax1.plot([], [], c='black', linestyle=':', label=r'$\rm M_{\rm disk} = 0.05 \, M_{\odot}$')
line_for_legend4, = ax1.plot([], [], c='black', linestyle='-.', label=r'$\rm M_{\rm disk} = 0.10 \, M_{\odot}$')

# Plot the dust (blue) on the right axis
ax1_2 = ax1.twinx()
ax1_2.plot(model_1a.r/const.au.cgs.value, model_1a.sigma_d, c='blue', linestyle='-')#, label='Dust')
ax1_2.plot(model_1b.r/const.au.cgs.value, model_1b.sigma_d, c='blue', linestyle='--')#, label='Dust')
ax1_2.plot(model_1c.r/const.au.cgs.value, model_1c.sigma_d, c='blue', linestyle=':')#, label='Dust')
ax1_2.plot(model_1d.r/const.au.cgs.value, model_1d.sigma_d, c='blue', linestyle='-.')#, label='Dust')


ax1.set_ylabel(r'$\Sigma_g$ $[\rm g \ \rm cm^{-2}]$', color='red')
ax1.set_xlim((5, 100)); ax1.set_ylim((0.005, 1000))
ax1.set_xticklabels([])

ax1_2.set_ylabel(r'$\Sigma_d$ $[\rm g \ \rm cm^{-2}]$', color='blue')
ax1_2.set_xlim((5, 100)); ax1_2.set_ylim((0.005, 10000))
ax1_2.set_xticklabels([])
ax1.legend(handles=[line_for_legend1, line_for_legend2, line_for_legend3, line_for_legend4], frameon=False, handlelength=1.5, loc='upper center', ncol=2)
ax1.set_yscale('log'); ax1_2.set_yscale('log')
ax1.set_title('Without Self-Gravity')

ax2.plot(model_1a.r/const.au.cgs.value, model_1a.T, linestyle='-', c='k')
ax2.plot(model_1b.r/const.au.cgs.value, model_1b.T, linestyle='--', c='k')
ax2.plot(model_1c.r/const.au.cgs.value, model_1c.T, linestyle=':', c='k')
ax2.plot(model_1d.r/const.au.cgs.value, model_1d.T, linestyle='-.', c='k')
ax2.set_ylabel('T [K]')
ax2.set_xlim((5, 100)); ax2.set_ylim((5, 200))
ax2.set_xticklabels([])
ax2.set_yscale('log')

ax3.plot(model_1a.r/const.au.cgs.value, model_1a.grain_size*10, linestyle='-', c='k')#, label='St = 0.314')
ax3.plot(model_1b.r/const.au.cgs.value, model_1b.grain_size*10, linestyle='--', c='k')#, label='St = 0.314')
ax3.plot(model_1c.r/const.au.cgs.value, model_1c.grain_size*10, linestyle=':', c='k')#, label='St = 0.314')
ax3.plot(model_1d.r/const.au.cgs.value, model_1d.grain_size*10, linestyle='-.', c='k')#, label='St = 0.314')
ax3.text(0.985, 0.965, 'St = 0.314', transform=ax3.transAxes, ha='right', va='top', size=16)

ax3.set_ylabel('a [mm]')
ax3.set_xlim((5, 100)); ax3.set_ylim((0.3, 110))        
#ax3.legend(frameon=False, handlelength=1, loc='upper right', ncol=1)   
ax3.set_yscale('log')

ax4.plot(model_1a.r/const.au.cgs.value, model_1a.h, linestyle='-', c='k')
ax4.plot(model_1b.r/const.au.cgs.value, model_1b.h, linestyle='--', c='k')
ax4.plot(model_1c.r/const.au.cgs.value, model_1c.h, linestyle=':', c='k')
ax4.plot(model_1d.r/const.au.cgs.value, model_1d.h, linestyle='-.', c='k')

ax4.set_ylabel('H / r')
ax4.set_xlabel('Radius [au]')
ax4.set_xlim((5, 100)); ax4.set_ylim((0.02, 0.08))
xticks = [5, 20, 40, 60, 80, 100]
ax4.set_xticks(xticks)
ax4.set_xticklabels([5, 20, 40, 60, 80, 100])

ax9.plot(model_1a.r/const.au.cgs.value, model_1a.Q, c='k', linestyle='-')
ax9.plot(model_1b.r/const.au.cgs.value, model_1b.Q, c='k', linestyle='--')
ax9.plot(model_1c.r/const.au.cgs.value, model_1c.Q, c='k', linestyle=':')
ax9.plot(model_1d.r/const.au.cgs.value, model_1d.Q, c='k', linestyle='-.')

ax9.set_ylabel('Q')
ax9.set_yscale('log')
ax9.set_xlim((5, 100)); ax9.set_ylim((4, 600))    

ax5.plot(model_2a.r/const.au.cgs.value, model_2a.sigma_g, c='red', linestyle='-')#label='Gas')
ax5_2 = ax5.twinx()
ax5_2.plot(model_2a.r/const.au.cgs.value, model_2a.sigma_d, c='blue', linestyle='-')#, label='Dust')
line_for_legend, = ax5.plot([], [], c='black', linestyle='-', label=r'$\rm M_{\rm disk} = 0.01675 \, M_{\odot}$')
ax5.set_ylabel(r'$\Sigma_g$ $[\rm g \ \rm cm^{-2}]$', color='red')
ax5_2.set_ylabel(r'$\Sigma_d$ $[\rm g \ \rm cm^{-2}]$', color='blue')
ax5.set_xlim((5, 100)); ax5.set_ylim((1e-2, 30))
ax5_2.set_xlim((5, 100)); ax5_2.set_ylim((1e-2, 30))
ax5.set_xticklabels([]); ax5_2.set_xticklabels([])
ax5.legend(handles=[line_for_legend],frameon=False, handlelength=1.5, loc='upper center', ncol=1)
ax5.set_yscale('log'); ax5_2.set_yscale('log')
ax5.set_title('With Self-Gravity')

ax6.plot(model_2a.r/const.au.cgs.value, model_2a.T, c='k', linestyle='-')
ax6.set_ylabel('T [K]')
ax6.set_xlim((5, 100)); ax6.set_ylim((20, 80))
ax6.set_xticklabels([])

ax7.plot(model_2a.r/const.au.cgs.value, model_2a.stoke, c='green', linestyle='-', label='a = 1 cm')
ax7.plot(model_2c.r/const.au.cgs.value, model_2c.stoke, c='purple', linestyle='-', label='a = 1 mm')
ax7.plot(model_2b.r/const.au.cgs.value, model_2b.stoke, c='orange', linestyle='-', label='a = 3 mm')
ax7.plot(model_2d.r/const.au.cgs.value, model_2d.stoke, c='brown', linestyle='-', label='a = 0.3 mm')
ax7.legend(frameon=False, handlelength=1, loc='lower right', ncol=2)
ax7.set_ylabel('St')
ax7.set_xlim((5, 100)); ax7.set_ylim((0.001, 5))           
ax7.set_yscale('log')

ax8.plot(model_2a.r/const.au.cgs.value, model_2a.h, c='k', linestyle='-')
ax8.set_ylabel('H / r')
ax8.set_xlabel('Radius [au]')
ax8.set_xlim((5, 100))#; ax4.set_ylim((0.02, 0.08))
ax8.set_xticks([5, 20, 40, 60, 80, 100])
ax8.set_xticklabels([5, 20, 40, 60, 80, 100])

ax10.plot(model_2a.r/const.au.cgs.value, model_2a.Q, c='k', linestyle='-')
ax10.set_ylabel('Q')
ax10.set_yscale('log')
ax10.set_xlim((5, 100)); ax10.set_ylim((20, 400))    

plt.savefig('NewDisk_Model_.png', dpi=300, bbox_inches='tight')
    
######

import numpy as np
import matplotlib.pyplot as plt


single = np.loadtxt('/Users/daniel/Desktop/Folders/StreamingInstability_YJ14/StreamingInstability_YJ14/data/a_single_opacities.txt')
amax = np.loadtxt('/Users/daniel/Desktop/Folders/StreamingInstability_YJ14/StreamingInstability_YJ14/data/a_max_opacities.txt')

a, abs_single, sca_single = single[:,0], single[:,1], single[:,2]
abs_25, sca_25, abs_35, sca_35 = amax[:,1], amax[:,2], amax[:,3], amax[:,4]


plt.style.use('/Users/daniel/Documents/plot_style.txt')


plt.loglog(a, sca_25, color='blue', linestyle='--', label=r'Grain Size Averaged $\kappa_{sca}$')
plt.loglog(a, abs_25, color='blue', linestyle='-', label=r'Grain Size Averaged $\kappa_{abs}$')
#plt.loglog(a, sca_35, color='green', linestyle='--', label=r'Grain Size Averaged $\kappa_{sca}$, q=3.5')
#plt.loglog(a, abs_35, color='green', linestyle='-', label=r'Grain Size Averaged $\kappa_{abs}$, q=3.5')

plt.loglog(a, sca_single, color='r', linestyle='--', label=r'Single Grain Size $\kappa_{sca}$')
plt.loglog(a, abs_single, color='r', linestyle='-', label=r'Single Grain Size $\kappa_{abs}$')

plt.ylabel(r'$\kappa_{\rm 1mm} \ [\rm cm^2 \ \rm g^{-1}]$')
plt.xlabel('Grain Size [cm]')
plt.title('DSHARP 1mm Wave Opacities')
plt.xlim((0.03, 100))
plt.ylim((3e-3, 100))
plt.legend(ncol=1, loc='upper right', handlelength=1, frameon=False)
plt.savefig('/Users/daniel/Desktop/Folders/opacities_figure.png', dpi=300, bbox_inches='tight')
plt.show()




######

import numpy as  np 
from StreamingInstability_YJ14 import disk_model
import astropy.constants as const

M_star = 1*const.M_sun.cgs.value # Mass of the star in cgs units
M_disk = 0.01675 * const.M_sun.cgs.value # Mass of the disk in cgs units

r = [10, 30, 100] # The radii of the 3 simulations
#r = np.arange(5, 101, 1)
r = [radius * const.au.cgs.value for radius in r] # Convert to [cgs] units
r_c = 300 * const.au.cgs.value # The characteristic radius of the disk (in [cgs])

grain_rho = 1.675 # Internal dust grain density (from DSHARP)
Z = 0.03 # Dust to gas ratio
q = 3/7. # Temperature power law index
T0 = 150 # Temperature at r = 1 au

# Four grain sizes in the simulations (polydisperse), in cgs units
grain_sizes = np.array([1., 0.3, 0.1, 0.03]) 

grain_sizes = np.array([1.194, 0.3582, 0.1194, 0.03582]) 
grain_sizes = np.array([1.2, 0.36, 0.12, 0.036]) 

# To store the simulation parameters (4 rows, 3 columns)
# The rows are the dust grains from largest to smallest
# The columns are the increasing distances from the star

stoke = np.zeros((len(grain_sizes), len(r)))
beta = np.zeros((len(grain_sizes), len(r)))
G = np.zeros((len(grain_sizes), len(r)))

for i, radius in enumerate(r):
    for j, grain_size in enumerate(grain_sizes):
        model = disk_model.Model(radius, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=grain_size, Z=Z, stoke=None, q=q, T0=T0)
        stoke[j, i] = model.stoke
        beta[j, i] = model.beta
        G[j, i] = model.G



grain_rho = 1.675 # Internal dust grain density in cgs units (from DSHARP)
stoke = 0.314 # Stokes number of the grain (dimensionless)
Z = 0.02 # Dust to gas ratio (dimensionless)
q = 1.0 # Temperature profile power law index
T0 = 600 # Temperature at r = 1 au, in Kelvin


import numpy as  np 
from StreamingInstability_YJ14 import disk_model
import astropy.constants as const

M_star = const.M_sun.cgs.value # Mass of the star in cgs units
M_disk = const.M_sun.cgs.value * 0.02
grain_rho = 2
Z = 0.03
q = 3/7.
T0 = 150

stoke = np.zeros((4,3))
sigma_g = np.zeros((4,3))
beta = np.zeros((4,3))
G = np.zeros((4,3))
i=0

for radius in [10, 30, 100]:
    r, r_c = radius*const.au.cgs.value, 300*const.au.cgs.value #360
    j=0
    for grain_size in [1, 0.3, 0.1, 0.03]:
        model = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=grain_size, Z=Z, stoke=None, q=q, T0=T0)
        stoke[j,i] = model.stoke
        sigma_g[j,i] = model.sigma_g
        beta[j,i] = model.beta
        G[j,i] = model.G
        j+=1
    i+=1


#Stokes number at different sigma_g and grain size


grain_size = stoke * 2. * column_density / np.pi / rho_grain


model = disk_model.Model(r, b=1.5, c=0.5, grain_size=grain_size, grain_rho=grain_rho, Z=Z, stoke=None):


######
change the disk model figure to include the stokes numbers
include Q and maybe G as a separate plot 
need to disentangle the contribution of optically thick regions and planetesimals to the mass underestimation
would be good to have 4 separate simulations for one location (r=10 au)



"""
