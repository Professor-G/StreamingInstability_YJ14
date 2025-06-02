.. protoRT documentation master file, created by
   sphinx-quickstart on Thu Mar 24 11:15:14 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to protoRT's documentation!
===============================
**protoRT** is a Python toolkit for radiative transfer and mass analysis in numerical simulations of planetesimal formation. It computes the optical depth, outgoing intensity, and the mass excess due to optically thick regions in dust-rich environments. This code was used to research the mass budget problem in protoplanetary disks, and accompanies the publication: Godines et al. 2025.

This documentation page details the analysis we conducted as well as examples on how to use the code to perform radiative transfer on numerical simulations of planet formation in the shearing box approximation. An in-depth overview of the code functionality is available in the `API Reference <https://protort.readthedocs.io/en/latest/autoapi/index.html>`_.

If you use this code for publication, we would appreciate citations to the paper. 

Questions? Please e-mail us: godines@nmsu.edu, wlyra@nmsu.edu   

Installation
==================

The current stable version can be installed via `pip`:

.. code-block:: bash

    pip install protoRT

You can also clone the development version:    

.. code-block:: bash

    git clone https://github.com/Professor-G/protoRT.git
    cd protoRT
    pip install .

Importing protoRT
==================

The code provides three main functionalities: Protoplanetary disk modeling, dust opacity calculations using `DSHARP <https://iopscience.iop.org/article/10.3847/2041-8213/aaf743>`_ opacities (supports mono and polydisperse models), and the main class `RadiativeTransferCube <https://protort.readthedocs.io/en/latest/autoapi/protoRT/rtcube/index.html#protoRT.rtcube.RadiativeTransferCube>`_ which conducts the radiative transfer, computing the optical depth, intensities, and resulting mass excess when optically thin emission is assumed.

To get started, a test dataset from a single-species streaming instability simulation without self-gravity, from `Yang and Johansen 2014 <https://iopscience.iop.org/article/10.1088/0004-637X/792/2/86>`_ is provided and automatically loaded if no data is input. This is a shearing box with a cubic domain, taken at orbit 100. The ``data`` will be loaded alongside the corresponding 1D coordinate ``axis`` array, in units of gas scale height.

.. code-block:: python

   from protoRT import rtcube

   cube = rtcube.RadiativeTransferCube()

.. figure:: _static/module_import.png
    :align: center
|
By default the class is initiated with a gas column density of 10 (g/cm²), a temperature of 30 (K), a gas scale height of 5 (au), and a stokes number of 0.3. The radiative transfer is also conducted at the 1 mm wavelength. There are a variety of simulation and model parameters to consider when using your own data, although note that a lot of these arguments are only required if analyzing multi-species simulations, with some others that are only necessary if the simulation is self-gravitating. Review the API documention for parameter details. 

The `configure <https://protort.readthedocs.io/en/latest/autoapi/protoRT/rtcube/index.html#id0>`_ class method will convert the simulation into physical units (cgs) and run through all the relevant calculations in order to perform the radiation transfer at the specified wavelength and compute the mass excess. When no opacities are input, the code will use the DSHARP opacities included in the `compute_opacities <https://protort.readthedocs.io/en/latest/autoapi/protoRT/compute_opacities/index.html>`_ module. When using the DSHARP opacities, ensure that the internal dust density of the dust grains (``grain_rho``) is the default value of 1.675 (g/cm³) to be consistent with the DSHARP dust model.

.. code-block:: python

   cube.configure()

Executing this method will automatically assign all of the relevant attributes including the optical depth at the output plane (``tau``), the corresponding intensity map (``intensity``), as well as the ``filling_factor``, ``mass_excess``, and the mass of each planetesimal in the simulation (``proto_mass``), if applicable.

.. code-block:: python

   import numpy as np
   import pylab as plt

   fig, axes = plt.subplots(1, 2, figsize=(12, 5))
   extent = [cube.axis[0], cube.axis[-1], cube.axis[0], cube.axis[-1]]

   # Plot optical depth with robust scaling
   med_tau = np.median(cube.tau)
   std_tau = np.median(np.abs(cube.tau - med_tau))
   vmin_tau, vmax_tau = med_tau - 3 * std_tau, med_tau + 10 * std_tau
   im0 = axes[0].imshow(cube.tau, vmin=vmin_tau, vmax=vmax_tau, cmap='viridis', extent=extent, origin='lower', aspect='equal')
   axes[0].set_title('Optical Depth', fontsize=18)
   axes[0].set_xlabel(r'$x/H$', fontsize=16)
   axes[0].set_ylabel(r'$y/H$', fontsize=16)
   cbar0 = plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
   cbar0.set_label(r'$\tau_{1.0\,\mathrm{mm}}^{\mathrm{eff}}$', fontsize=14)

   # Plot intensity with robust scaling
   med_int = np.median(cube.intensity)
   std_int = np.median(np.abs(cube.intensity - med_int))
   vmin_int, vmax_int = med_int - 3 * std_int, med_int + 10 * std_int
   im1 = axes[1].imshow(cube.intensity, vmin=vmin_int, vmax=vmax_int, cmap='inferno', extent=extent, origin='lower', aspect='equal')
   axes[1].set_title('Outgoing Intensity', fontsize=18)
   axes[1].set_xlabel(r'$x/H$', fontsize=16)
   axes[1].set_ylabel(r'$y/H$', fontsize=16)
   cbar1 = plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
   cbar1.set_label(r'$I_{1.0\,\mathrm{mm}} \, (\mathrm{erg}~\mathrm{s}^{-1}~\mathrm{cm}^{-2}~\mathrm{Hz}^{-1}~\mathrm{sr}^{-1})$', fontsize=14)

   fig.suptitle(f'RT Results (Absorption Only): Mass Excess = {float(cube.mass_excess):.3f}, Filling Factor = {float(cube.filling_factor):.3f}', fontsize=18)
   plt.tight_layout()
   plt.show()


.. figure:: _static/tau_intensity_abs_only.png
    :align: center
|
By default the scattering opacities are not considered in the radiative transfer, this is controlled via the ``include_scattering`` argument. The following shows the same analysis but with scattering. 

**IMPORTANT: The data cube is converted from code to physical units during the configuration, and is overwritten. As such, do not re-configure multiple times, instead re-instantiate the class object.**

.. code-block:: python

   cube = rtcube.RadiativeTransferCube(include_scattering=True)
   cube.configure() 

   # Plot
   fig, axes = plt.subplots(1, 2, figsize=(12, 5))
   extent = [cube.axis[0], cube.axis[-1], cube.axis[0], cube.axis[-1]]

   # Plot optical depth with robust scaling
   med_tau = np.median(cube.tau)
   std_tau = np.median(np.abs(cube.tau - med_tau))
   vmin_tau, vmax_tau = med_tau - 3 * std_tau, med_tau + 10 * std_tau
   im0 = axes[0].imshow(cube.tau, vmin=vmin_tau, vmax=vmax_tau, cmap='viridis', extent=extent, origin='lower', aspect='equal')
   axes[0].set_title('Optical Depth', fontsize=18)
   axes[0].set_xlabel(r'$x/H$', fontsize=16)
   axes[0].set_ylabel(r'$y/H$', fontsize=16)
   cbar0 = plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
   cbar0.set_label(r'$\tau_{1.0\,\mathrm{mm}}^{\mathrm{eff}}$', fontsize=14)

   # Plot intensity with robust scaling
   med_int = np.median(cube.intensity)
   std_int = np.median(np.abs(cube.intensity - med_int))
   vmin_int, vmax_int = med_int - 3 * std_int, med_int + 10 * std_int
   im1 = axes[1].imshow(cube.intensity, vmin=vmin_int, vmax=vmax_int, cmap='inferno', extent=extent, origin='lower', aspect='equal')
   axes[1].set_title('Outgoing Intensity', fontsize=18)
   axes[1].set_xlabel(r'$x/H$', fontsize=16)
   axes[1].set_ylabel(r'$y/H$', fontsize=16)
   cbar1 = plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
   cbar1.set_label(r'$I_{1.0\,\mathrm{mm}} \, (\mathrm{erg}~\mathrm{s}^{-1}~\mathrm{cm}^{-2}~\mathrm{Hz}^{-1}~\mathrm{sr}^{-1})$', fontsize=14)

   fig.suptitle(f'RT Results (Absorption + Scattering): Mass Excess = {float(cube.mass_excess):.3f}, Filling Factor = {float(cube.filling_factor):.3f}', fontsize=18)
   plt.tight_layout()
   plt.show()


.. figure:: _static/tau_intensity_abs_and_sca.png
    :align: center
|
In this example we can see that scattering effects further suppress the emergent intensity and increase the mass excess and filling factor by approximately 16\% and 13\%, respectively.

To learn how to use the code, please review the `following page <https://protort.readthedocs.io/en/latest/source/Godines%20et%20al%202025.html>`_ which details how the program was used in Godines et al. 2025, which employed multi-species simulations with self-gravity, at three locations in the disk. 

Pages
==================
.. toctree::
   :maxdepth: 1

   source/Godines et al 2025

Documentation
==================

Here is the documentation for all the modules:

.. toctree::
   :maxdepth: 1

   source/protoRT
