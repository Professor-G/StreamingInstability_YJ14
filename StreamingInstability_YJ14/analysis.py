#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 02:53:37 2022

@author: daniel
"""
import numpy as np  
import matplotlib.pyplot as plt  
import pencil as pc   
from StreamingInstability_YJ14 import shearing_box

f = pc.read.var(trimall=True, ivar=var)

cube = shearing_box.density_cube(data=f.rhop,axis=f.x)
column_density = np.logspace(0, 2.7, 100)

freq=[3.44589e11, 2.30609e11, 9.99308e10] #0.87 mm, 1.3 mm, 3.0 mm
kappa = [1,1,1]

flux = []
for nu in freq: 
    print(nu)
    for j in [60]: #Column density iterations
        cube.column_density = j
        cube.configure(nu=nu)
        flux_bb = cube.blackbody(nu=nu) * (1-np.exp(-cube.tau))
        flux.append(np.mean(flux_bb))

spectral_index_alpha = np.log10(flux[0]/flux[1]) / np.log10(freq[0]/freq[1])
spectral_index_beta = np.log10(kappa[0]/kappa[1]) / np.log10(freq[0]/freq[1])

wavelength = [0.87, 1.3, 3]
plt.plot(wavelength, flux, 'ko--', label='ALMA')
plt.xlabel(r'$\lambda$ (mm)', fontweight='ultralight', size=18)
plt.ylabel(r'$F_\nu$', fontweight='ultralight', size=18)
plt.tick_params(labelsize=14, axis="both", which="both")
plt.grid(True, color='k', alpha=0.35, linewidth=1.5, linestyle='--')
plt.legend(prop={'size':16})
plt.rcParams["font.family"] = "Times New Roman"
plt.title('Spectral Indices', size=22)
#plt.text(0.93, 5.4e-13, r'$\alpha$ = ', fontsize=18)
#plt.text(0.93, 5.2e-13, r'$\beta$ = ', fontsize=18)
#plt.text(1.3, 2.979e-13, r'$\alpha$ = ', fontsize=18)
#plt.text(1.3, 2.779e-13, r'$\beta$ = ', fontsize=18)
#plt.text(3, 1e-13, r'$\alpha$ = ', fontsize=18)
#plt.text(3, 1.2e-13, r'$\beta$ = ', fontsize=18)
plt.show()



