.. _Godines_et_al_2025:


Godines et al. 2025
===========

Figure 1 - Disk Model
-----------

The following code shows how we set up our disk model. This model was used to configure the simulations, including the calculation of the stokes numbers, pressure gradient parameter and strength of the self-gravity at the three different disk locations.

.. code-block:: python

	import numpy as  np 
	import matplotlib.pyplot as plt  
	import astropy.constants as const
	from protoRT import disk_model

	try:
	    import scienceplots
	    plt.style.use('science')
	    plt.rcParams.update({'font.size': 26,})
	    plt.rcParams.update({'lines.linewidth': 1.5})
	except:
	    print('WARNING: Could not import scienceplots, please install via pip for proper figure formatting.')
	    plt.style.use('default')

	# Mass of the star 
	M_star = const.M_sun.cgs.value 

	# Mass of the protoplanetary disk
	M_disk = 0.02*const.M_sun.cgs.value 

	# Radii which to model, and the characteristic radius of the disk (in [au])
	r, r_c = np.arange(5,110.25,0.25), 300 

	# Convert to cgs units 
	r, r_c = r*const.au.cgs.value, r_c*const.au.cgs.value 


	## Disk model adopted in this work (Section 2.1 of the paper) ##
	## Used to convert our multi-species simulations with self-gravity from code units to physical units (cgs) ###

	# Internal dust grain density (from DSHARP dust model)
	grain_rho = 1.675

	# Dust to gas ratio (Equation 11)
	Z = 0.03 

	# Temperature power law index (Equation 2)
	q = 3/7. 

	# Temperature at r = 1 au (Equation 2)
	T0 = 150 

	# The sizes of the four dust grains in our simulations
	grain_sizes = np.array([1.194, 0.3582, 0.1194, 0.03582]) 

	# The Model class in the disk_model module only works with one grain size, so we define four unique models 
	# All other model parameteres are the same so these models only differ in their respective Stokes number
	model_2a = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=grain_sizes[0], Z=Z, stoke=None, q=q, T0=T0)
	model_2b = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=grain_sizes[1], Z=Z, stoke=None, q=q, T0=T0)
	model_2c = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=grain_sizes[2], Z=Z, stoke=None, q=q, T0=T0)
	model_2d = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, grain_size=grain_sizes[3], Z=Z, stoke=None, q=q, T0=T0)

	# Plot 
	fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(8, 20), sharex=True)
	fig.subplots_adjust(hspace=0.075) 

	# Plot 1: Surface Density
	axes[0].plot(model_2a.r / const.au.cgs.value, model_2a.sigma_g, c='k', linestyle='-', label='Gas')
	axes[0].plot(model_2a.r / const.au.cgs.value, model_2a.sigma_d, c='k', linestyle='--', label='Dust')
	axes[0].set_ylabel(r'$\Sigma$ (g cm$^{-2}$)')
	axes[0].set_yscale('log')
	axes[0].legend(frameon=False, loc='upper right', handlelength=0.75)
	axes[0].set_title('Protoplanetary Disk Model')

	# Plot 2: Temperature
	axes[1].plot(model_2a.r / const.au.cgs.value, model_2a.T, c='k', linestyle='-')
	axes[1].set_ylabel(r'T (K)')
	axes[1].set_ylim((20,80))

	# Plot 3: Stokes Number
	axes[2].plot(model_2a.r / const.au.cgs.value, model_2a.stoke, c='k', linestyle='-', label='12 mm')
	axes[2].plot(model_2c.r / const.au.cgs.value, model_2c.stoke, c='k', linestyle='--', label='1.2 mm')
	axes[2].plot(model_2b.r / const.au.cgs.value, model_2b.stoke, c='k', linestyle=':', label='3.6 mm')
	axes[2].plot(model_2d.r / const.au.cgs.value, model_2d.stoke, c='k', linestyle='-.', label='0.36 mm')
	axes[2].set_ylabel('St'); axes[2].set_yscale('log')
	axes[2].set_ylim((0.001, 5.25)) 

	# Add vertical lines
	axes[2].axvline(x=10, linestyle=':', linewidth=2.5, color='red', alpha=0.65)
	axes[2].axvline(x=30, linestyle=':', linewidth=2.5, color='red', alpha=0.65)
	axes[2].axvline(x=100, linestyle=':', linewidth=2.5, color='red', alpha=0.65)
	legend = axes[2].legend(loc='lower right', handlelength=0.75, ncol=2)

	# Plot 4: Scale Height
	axes[3].plot(model_2a.r / const.au.cgs.value, model_2a.h, c='k', linestyle='-')
	axes[3].set_ylabel(r'H/r')
	axes[3].set_ylim((0.04, 0.1))

	# Plot 5: Toomre Q Parameter
	axes[4].plot(model_2a.r / const.au.cgs.value, model_2a.Q, c='k', linestyle='-')
	axes[4].set_ylabel(r'$Q$'); axes[4].set_xlabel('Radius (au)')

	# X-axis formatting (hiding x-tick labels for all but the bottom plot)
	for ax in axes:
	    ax.set_xlim(5., 102.)
	    ax.label_outer() 

	# Add the vertical red dashed lines to denote location of our three simulations
	for i in range(5):
	    axes[i].axvline(x=10, linestyle=':', linewidth=2.5, color='red', alpha=0.65)
	    axes[i].axvline(x=30, linestyle=':', linewidth=2.5, color='red', alpha=0.65)
	    axes[i].axvline(x=100, linestyle=':', linewidth=2.5, color='red', alpha=0.65)

	# Add vertical text labels aligned with the lines (only lower plot)
	axes[4].text(10., 114.037, '10 au', color='red', rotation=90, verticalalignment='top', horizontalalignment='right')
	axes[4].text(30., 114.037, '30 au', color='red', rotation=90, verticalalignment='top', horizontalalignment='right')
	axes[4].text(100., 114.037, '100 au', color='red', rotation=90, verticalalignment='top', horizontalalignment='right')

	# Save
	plt.savefig('Disk_Model_SelfGravity_OneColumn.png', dpi=300, bbox_inches='tight')
	plt.show()


Figures 2 & 3 - Simulations
-----------

The following plot overlays the four species in our simulations, which evolve largely independently, on the strong clumping boundary for the streaming instability, as reported by `Lim et al 2025 <https://ui.adsabs.harvard.edu/abs/2025ApJ...981..160L/abstract>`_. 

.. code-block:: python

	import numpy as np
	import matplotlib.pyplot as plt

	try:
	    import scienceplots
	    plt.style.use('science')
	    plt.rcParams.update({'font.size': 26,})
	    plt.rcParams.update({'lines.linewidth': 1.5})
	except:
	    print('WARNING: Could not import scienceplots, please install via pip for proper figure formatting.')
	    plt.style.use('default')


	plt.figure(figsize=(8,8))

	# Simulation parameters, stokes numbers and pressure gradient parameter
	st10 = np.array([0.345, 0.103, 0.034, 0.0103])
	st30 = np.array([1.105, 0.331, 0.110, 0.033])
	st100 = np.array([4.651, 1.395, 0.465, 0.134])
	Pi = np.array([0.0545, 0.0745, 0.105])

	# The initial dust-to-gas ratio in our simulations (Equation 11)
	Total_Z = 0.03

	# In our simulations with four grain sizes and Z=0.03, the species act indepedent of one another (see Krapp et al. 2019)
	Collective_Z = 0.03 / 4.

	# The critical parameter adoped in this work (Z/Pi)
	ratio_collective = Collective_Z / (Pi)

	# Plot where the four species in each of the three simulations fall within this boundary
	plt.scatter(st10, [ratio_collective[0]]*4, marker='*', facecolor='#1f77b4', s=350, edgecolor='#1f77b4')
	plt.scatter(st30, [ratio_collective[1]]*4, marker='*', facecolor='#ff7f0e', s=350, edgecolor='#ff7f0e')
	plt.scatter(st100, [ratio_collective[2]]*4, marker='*', facecolor='#2ca02c', s=350, edgecolor='#2ca02c')

	# Adding horizontal dashed lines and text to denote where each simulation is in the disk
	plt.axhline(y = ratio_collective[0], linestyle='--', color='#1f77b4', alpha=0.5)
	plt.text(0.001, ratio_collective[0]+0.0017, "10 au", color='#1f77b4', fontweight="bold")

	plt.axhline(y = ratio_collective[1], linestyle='--', color='#ff7f0e', alpha=0.5)
	plt.text(0.001, ratio_collective[1]+0.0014, "30 au", color='#ff7f0e', fontweight="bold")

	plt.axhline(y = ratio_collective[2], linestyle='--', color='#2ca02c', alpha=0.5)
	plt.text(0.001, ratio_collective[2]+0.0008, "100 au", color='#2ca02c', fontweight="bold")

	# Now plot the Li+25 boundary
	# Stokes numbers to plot (x-axis)
	x = np.arange(-3, 0.74037, 0.01)
	St = 10 ** x

	# Boundary parameters
	Pi_ = 0.05 
	A = 0.10
	B = 0.07
	C = -2.36
	C = C - np.log10(Pi_) # To show Z / Pi. 

	# Plot the critical boundary
	Zcrit_array = (A * np.log10(St)**2) + (B * np.log10(St)) + C
	Zcrit_array = 10 ** Zcrit_array 
	plt.loglog(St, Zcrit_array, color='#d62728', label="Lim+25")

	# Label the clumping regions
	plt.title('Streaming Instability Strong Clumping Regime\nPolydisperse Simulations (4 Species)')
	plt.text(0.1, 0.155, "Strong Clumping", fontweight="bold")
	plt.text(0.00105, 0.155, "No Strong\nClumping", fontweight="bold")
	plt.xlabel("St"); plt.ylabel(r"$Z \ / \ \Pi$")
	plt.xlim(1e-3, 5.5)

	# Save
	plt.legend(ncol=1, loc='upper right', handlelength=1)
	plt.subplots_adjust(top=0.97, right=0.97, left=0.12, bottom=0.12)
	plt.savefig('SI_criteria_Independent.png', dpi=300, bbox_inches='tight')
	plt.show()



The following code shows the time evolution of the three simulations.


Figures 4 & 5 - Frequency-dependent Dust Opacities
-----------


Figure 6 & 7 - Optical Depth and Intensity Calculations
-----------


Figure 8 & 9 - Radiative Transfer Results
-----------

