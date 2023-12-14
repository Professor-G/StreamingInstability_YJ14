import numpy as np
import astropy.constants as const
from StreamingInstability_YJ14 import disk_model, shearing_box 

# Generate a disk model to extract disk parameters such as H, T and gas column density
M_star = const.M_sun.cgs.value # Mass of the star     
M_disk = 0.01*const.M_sun.cgs.value # Mass of the disk                                                       

r, r_c = np.arange(5,150,1), 300 # Radii which to model, and the characteristic radius of the disk (in [au])
r, r_c = r*const.au.cgs.value, r_c*const.au.cgs.value # Convert to cgs units

grain_rho = 1.675 # Internal dust grain density (from DSHARP)                                  
stoke = 0.314 # Stokes number of the grain                                                   
Z = 0.02 # Dust to gas ratio                                                                   
q = 1.0 # Temperature power law index                                                          
T0 = 600 # Temperature at r = 1 au                                                             
                                                                                               
disk_model = disk_model.Model(r, r_c, M_star, M_disk, grain_rho=grain_rho, Z=Z, stoke=stoke, q=q, T0=T0)                                                                                       

# Whether to include photon scattering as well as absorption 
include_scattering = True

# Simulation code units
code_omega = 2*np.pi
code_cs = 2*np.pi
code_rho = 1

filling_factors, mass_underestimations = [], []

for i in range(len(r)):
	#i = 44 # 49 au
	i = 144 # 149 au
	cube = shearing_box.density_cube(stoke=stoke, eps_dtog=Z, grain_rho=grain_rho, q=2.5,
		T=disk_model.T[i], H=disk_model.H[i], column_density=disk_model.sigma_g[i],
		code_rho=code_rho, code_cs=code_cs, code_omega=code_omega, include_scattering=include_scattering)
	cube.configure()
	print(f"{i+1} out of {len(r)} \nME: {cube.mass_excess}\nFF: {cube.filling_factor}")
	filling_factors.append(cube.filling_factor); mass_underestimations.append(cube.mass_excess)
	break





import matplotlib.pyplot as plt

fig, ax1 = plt.subplots()

# Plot the filling_factors on the left y-axis
ax1.plot(r / const.au.cgs.value, filling_factors, 'b-', label='Filling Factors')
ax1.set_xlabel('r [au]')
ax1.set_ylabel('Filling Factors', color='b')
ax1.tick_params('y', colors='b')

# Create a second y-axis on the right
ax2 = ax1.twinx()
ax2.plot(r / const.au.cgs.value, mass_underestimations, 'r-', label='Mass Underestimations')
ax2.set_ylabel('Mass Underestimations', color='r')
ax2.tick_params('y', colors='r')

# Add a vertical line where filling factors are zero
zero_fill_index = (np.array(filling_factors) == 0).nonzero()[0]
if zero_fill_index.size > 0:
    zero_fill_x = r[zero_fill_index[0]] / const.au.cgs.value
    ax1.axvline(x=zero_fill_x, color='g', linestyle='--')#, label='Zero Filling Factors')
    ax1.text(zero_fill_x, 0.005, 'Zero Filling Factors', color='g', ha='right', va='bottom', rotation='vertical')

# Display legends
#lines, labels = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax2.legend(lines + lines2, labels + labels2, loc='best')

ax1.set_ylim((0, 0.01))# ax1.set_ylim((0, 0.06))
ax2.set_ylim((1, 1.3)) # ax2.set_ylim((1, 1.8))
ax1.set_xlim((5, 300))

plt.title('CY14 Mono SI Simulation w/out Scattering (t/T = 100)')
plt.tight_layout()
plt.savefig('Monodisperse_case_without_scattering.png', dpi=900, bbox_inches='tight')
plt.show()





