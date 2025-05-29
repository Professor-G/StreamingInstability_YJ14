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
    Generate and visualize a 1D midplane model of a protoplanetary disk.

    This class defines a physically consistent disk structure, including gas and dust surface densities,
    temperature, and stability parameters. It enables proper unit conversion and parameter initialization 
    for use in simulations and radiative transfer calculations.
    
    Note
    ----------
    All units are assumed to be in CGS.

    Parameters
    ----------
    r : float or ndarray
        Disk radii where quantities are evaluated (cm).
    r_c : float
        Characteristic disk radius (cm), typically ~30 AU (compact) to 300 AU (large).
    M_star : float
        Mass of the central star (g).
    M_disk : float
        Total disk mass (g).
    grain_size : float, optional
        Grain size (cm); used only if `stoke` is None. Default is 0.1.
    grain_rho : float, optional
        Internal grain density (g/cm³). Default is 1.
    T0 : float, optional
        Reference temperature at 1 AU (K). Default is 150.
    q : float, optional
        Temperature power-law index. Default is 3/7.
    gamma : float, optional
        Adiabatic index. Default is 1 (isothermal).
    mmol : float, optional
        Mean molecular weight. Default is 2.3.
    Z : float, optional
        Global dust-to-gas mass ratio. Default is 0.01.
    stoke : float, optional
        Stokes number. If provided, grain size is derived from it.

    Attributes
    ----------
    r : float or ndarray
        Radial positions in the disk where quantities are evaluated (cm).
    r_c : float
        Characteristic radius of the disk (cm).
    M_star : float
        Mass of the central star (g).
    M_disk : float
        Total mass of the disk (g).
    grain_size : float or ndarray
        Grain size (cm); either input directly or derived from Stokes number.
    grain_rho : float
        Internal grain density (g/cm³).
    T0 : float
        Reference temperature at 1 AU (K).
    q : float
        Temperature power-law index.
    gamma : float
        Adiabatic index.
    mmol : float
        Mean molecular weight.
    Z : float
        Dust-to-gas mass ratio.
    stoke : float or ndarray
        Stokes number; either input or calculated from grain size.
    omega : ndarray
        Keplerian angular velocity (rad/s).
    sigma_g : ndarray
        Gas surface density profile (g/cm²).
    sigma_d : ndarray
        Dust surface density profile (g/cm²).
    T : ndarray
        Temperature profile (K).
    cs : ndarray
        Sound speed (cm/s).
    H : ndarray
        Gas pressure scale height (cm).
    h : ndarray
        Disk aspect ratio (H/r).
    Q : ndarray
        Toomre Q parameter.
    G : ndarray
        Self-gravity parameter (G-tilde).
    beta : ndarray
        Dimensionless radial pressure gradient parameter.
    plot_stoke : bool
        True if plotting Stokes number; False if plotting grain size.

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
        Compute all physical disk parameters and optionally print them.

        This method calculates the complete set of disk parameters including
        gas and dust surface densities, temperature, sound speed, disk thickness,
        and stability indicators like Toomre Q and the beta parameter. Depending
        on whether `stoke` is provided, it computes the corresponding grain size
        or Stokes number.

        Parameters
        ----------
        print_params : bool, optional
            If True, the computed parameters are printed to the console. Default is True.

        Returns
        -------
        None
            Updates instance attributes in-place.
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
        Calculate the gas surface density profile.

        Follows the exponential tapering model from Drazkowska et al. (2022).

        Returns
        -------
        None
            Sets `self.sigma_g` as a function of radius (g/cm²).
        """

        self.sigma_g = ((self.M_disk/(2.*np.pi*self.r_c**2))*(self.r/self.r_c)**-1)*np.e**(-(self.r/self.r_c))

    def calc_sigma_d(self):
        """
        Calculate the dust surface density profile.

        Computed as a constant fraction of the gas surface density using the
        global dust-to-gas mass ratio `Z`.

        Returns
        -------
        None
            Sets `self.sigma_d` as a function of radius (g/cm²).
        """

        self.sigma_d = self.sigma_g * self.Z

    def calc_stokes(self):
        """
        Calculate the Stokes number from grain properties and gas density.

        Assumes Epstein drag regime, valid for small particles in low-density gas.

        Returns
        -------
        None
            Sets `self.stoke` as a function of radius (dimensionless).
        """

        self.stoke = np.pi * self.grain_size * self.grain_rho / 2. / self.sigma_g

    def calc_grain_sizes(self):
        """
        Compute grain sizes from a fixed Stokes number.

        This inverts the Stokes number formula to recover grain size at each radius.

        Returns
        -------
        None
            Sets `self.grain_size` as a function of radius (cm).
        """

        self.grain_size = 2 * self.stoke * self.sigma_g / np.pi / self.grain_rho

    def calc_omega(self):
        """
        Compute Keplerian angular velocity at each radius.

        Returns
        -------
        None
            Sets `self.omega` in units of rad/s.
        """

        self.omega = np.sqrt(const.G.cgs.value*self.M_star/self.r**3)

    def calc_Q(self):
        """
        Compute the Toomre Q parameter for gravitational stability.

        Returns
        -------
        None
            Sets `self.Q` (dimensionless).
        """

        self.Q = self.cs * self.omega / np.pi / const.G.cgs.value / self.sigma_g

    def calc_T(self):
        """
        Compute the disk temperature profile from a power-law model.

        Follows Ida et al. (2016), applicable to outer disk regions where
        viscous heating is negligible.

        Returns
        -------
        None
            Sets `self.T` as a function of radius (K).
        """

        self.T = self.T0*(self.r / const.au.cgs.value)**(-self.q)

    def calc_h(self):
        """
        Compute the disk aspect ratio h = H / r.

        Returns
        -------
        None
            Sets `self.h` (dimensionless).
        """
        
        self.h = self.H / self.r

    def calc_cs(self):
        """
        Compute the sound speed based on local temperature and molecular weight (using hydrogen molecule mass)

        Returns
        -------
        None
            Sets `self.cs` as a function of radius (cm/s).
        """

        molecule_mass = self.mmol * const.m_p.cgs.value
        self.cs = np.sqrt(self.gamma * const.k_B.cgs.value * self.T / molecule_mass)

    def calc_H(self):
        """
        Compute the gas pressure scale height of the disk.

        Returns
        -------
        None
            Sets `self.H` as a function of radius (cm).
        """

        self.H = self.cs / self.omega

    def calc_G(self):
        """
        Compute the dimensionless self-gravity parameter G-tilde.

        Defined as G = sqrt(8 / pi) / Q.

        Returns
        -------
        None
            Sets `self.G` (dimensionless).
        """

        self.G = np.sqrt(8/np.pi) / self.Q

    def calc_beta(self):
        """
        Compute the radial pressure gradient parameter β.

        The dimensionless β parameter accounts for the radial variation in pressure,
        and is defined as:
            β = h × dlnp/dlnr

        Assuming power-law profiles for surface density, temperature, and Keplerian rotation,
        the gradient simplifies to:
            dlnp/dlnr ≈ −1 − 0.5 × q − 1.5

        For the default temperature index q = 3/7, this evaluates to:
            β ≈ h × (−2.2857)

        Returns
        -------
        None
            Updates the `self.beta` attribute (dimensionless).
        """

        self.beta = self.h * (-1.0 - (0.5*-self.q) - 1.5)

    def plot(self, xticks=[5, 20, 40, 60, 80, 100], ax1_xlim=None, ax1_ylim=None, ax2_xlim=None, ax2_ylim=None,
        ax3_xlim=None, ax3_ylim=None, ax4_xlim=None, ax4_ylim=None, ax1_ylog=True, ax2_ylog=True, ax3_ylog=True, ax4_ylog=False, 
        include_grid=False, savefig=False, path=None, box_index=None, plot_vertical=False, title="Protoplanetary Disk Model"):      
        """
        Plot the radial structure of the disk model. Only plots four key model parameters: gas and dust density profiles, temperature profile, grain size or stokes number, and aspect ratio.

        Parameters
        ----------
        xticks : list of float
            Tick positions on the x-axis (AU).
        ax{1..4}_xlim, ax{1..4}_ylim : tuple or None
            Axis limits for each subplot.
        ax{1..4}_ylog : bool
            Whether to use log scaling on the y-axis.
        include_grid : bool
            Add grid to all subplots.
        savefig : bool
            If True, save the figure instead of displaying.
        path : str or None
            Path to save the figure. Defaults to home directory.
        box_index : int or None
            Index of radius value to highlight.
        plot_vertical : bool
            Arrange plots vertically if True; otherwise in 2x2 grid.
        title : str
            Title of the plot.

        Returns
        -------
        None
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
        ax1.set_ylabel(r'$\Sigma$ $(\rm g \ \rm cm^{-2})$')
        ax1.legend(frameon=False, handlelength=1, loc='upper right', ncol=1)
        
        ### ax2 ###

        ax2.plot(self.r/const.au.cgs.value, self.T, c='k')

        if box_index is not None: ax2.scatter(self.r[box_index]/const.au.cgs.value, self.T[box_index], marker='s', s=100)
        if include_grid: ax2.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        if ax2_xlim is not None: ax1.set_xlim(ax2_xlim) #(5, 100)
        if ax2_ylim is not None: ax1.set_ylim(ax2_ylim) #(5, 250)
        if ax2_ylog: ax2.set_yscale('log')

        ax2.set_xticklabels([])
        ax2.set_ylabel('T (K)')
        
        ### ax3 ###

        if self.plot_stoke:

            ax3.plot(self.r/const.au.cgs.value, self.stoke, c='k')

            if box_index is not None: ax3.scatter(self.r[box_index]/const.au.cgs.value, self.stoke[box_index], marker='s', s=100)
            ax3.set_ylabel(r'$\tau_s$')

        else:

            ax3.plot(self.r/const.au.cgs.value, self.grain_size*10, c='k')
            
            if box_index is not None: ax3.scatter(self.r[box_index]/const.au.cgs.value, self.grain_size[box_index]*10, marker='s', s=100)
            ax3.set_ylabel('a (mm)')

        if include_grid: ax3.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        if ax3_xlim is not None: ax1.set_xlim(ax3_xlim) #(5, 100)
        if ax3_ylim is not None: ax1.set_ylim(ax3_ylim) #(5, 250)
        if ax3_ylog: ax3.set_yscale('log')

        if plot_vertical:
            ax3.set_xticklabels([])
        else:
            ax3.set_xticks(xticks); ax3.set_xticklabels(xticks)
            ax3.set_xlabel('Radius (au)')

        ### ax4 ###

        ax4.plot(self.r/const.au.cgs.value, self.h, c='k')

        if box_index is not None: ax4.scatter(self.r[box_index]/const.au.cgs.value, self.h[box_index], marker='s', s=100)
        if include_grid: ax4.grid(True, which='both', color='k', alpha=0.25, linewidth=1.5, linestyle='--')
        if ax4_xlim is not None: ax1.set_xlim(ax4_xlim) #(5, 100)
        if ax4_ylim is not None: ax1.set_ylim(ax4_ylim) #(0.02, 0.08)
        if ax4_ylog: ax4.set_yscale('log')

        ax4.set_xticks(xticks); ax4.set_xticklabels(xticks)
        ax4.set_ylabel('H / r'); ax4.set_xlabel('Radius (au)')

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
    """
    Configure matplotlib style settings for consistent plot appearance.

    This function updates `matplotlib.pyplot.rcParams` to define a custom
    plotting style, which is used before saving figures. 

    Returns
    -------
    None
        Applies changes to matplotlib's global `rcParams`.
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
 



