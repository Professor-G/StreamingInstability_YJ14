.. protoRT documentation master file, created by
   sphinx-quickstart on Thu Mar 24 11:15:14 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to protoRT's documentation!
===============================
**protoRT** is a Python toolkit for radiative transfer and mass analysis in numerical simulations of planetesimal formation. It computes the optical depth, outgoing intensity, and the mass excess due to optically thick regions in dust-rich environments. This code was used to research the mass budget problem in protoplanetary disks, and accompanies the publication: Godines et al. 2025.

This documentation page details the analysis we conducted as well as examples on how to use the code to perform radiative transfer on numerical simulations of planet formation in the shearing box approximation. An in-depth overview of the code functionality is available in the `API Reference <https://protort.readthedocs.io/en/latest/autoapi/protoRT/rtcube/index.html#protoRT.rtcube.RadiativeTransferCube>`_.

If you use this code for publication we would appreciate citations to the paper. 

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

The code provides three main functionalities: Disk model configuration, dust opacity calculations using DSHSARP opacities (supports mono and polydisperse models), and the main class `RadiativeTransferCube <https://protort.readthedocs.io/en/latest/autoapi/protoRT/rtcube/index.html#protoRT.rtcube.RadiativeTransferCube>`_ which conducts the radiative transfer, computing the optical depth, intensities, and resulting mass excess when optically thin emission is assumed.

To get started, a test dataset from a single-species streaming instability simulation without self-gravity, from Yang & Johansen (2014), is provided and automatically loaded if no data is input. This is a shearing box with cubic domain, taken orbit 100. It will be loaded alongside the corresponding 1D coordinate axis array, in units of gas scale height.

.. code-block:: python

   from protoRT import rtcube

   cube = rtcube.RadiativeTransferCube()

.. figure:: _static/module_import.png
    :align: center
|
By default the class is initiated with a gas column density of 30 (g/cm2), a temperature of 30 (K), a gas scale height of 5 (au), stokes number of 0.3 at a wavelength of 1 mm. There are a variety of simulation and model parameters to consider when using your own data, although note that some parameters are only required if analyzing multi-species simulations, and others that are required if the simulation is self-gravitating. please reivew the API documention for details. 

The `configure <https://protort.readthedocs.io/en/latest/autoapi/protoRT/rtcube/index.html#id0>`_ class method will run through all the relevant calculations in order to perform the radiation transfer at the specified wavelength and compute the mass excess. As no opacity parameters were input, the code will by default use the DSHARP opacities available in the `compute_opacities <https://protort.readthedocs.io/en/latest/autoapi/protoRT/compute_opacities/index.html>`_ module.

.. code-block:: python

   cube.configure()

This method will automatically assign all the relevant attributse including the optical depth at the output plane (``tau``), the corresponding intensity map (``intensity``), as well as the ``filling_factor``, the ``mass_excess``, and the mass of each planetesimal in the simulation, if applicable.

.. code-block:: python

   import pylab as plt

   plt.imshow(cube.tau)
   plt.imshow(cube.intensity)
   plt.title(f'Simulation results: Mass Excess: {cube.mass_excess}, Filling Factor: {cube.filling_factor}')


.. figure:: _static/module_import.png
    :align: center
|


Pages
==================
.. toctree::
   :maxdepth: 1

   source/Mass Budget Problem in Protoplanetary Disks
   source/Hydrodynamical Simulations of the Streaming Instability
   source/Godines et al 2025 - Figures

Documentation
==================

Here is the documentation for all the modules:

.. toctree::
   :maxdepth: 1

   source/protoRT
