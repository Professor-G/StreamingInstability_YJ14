import numpy as np  
from astropy.constants import au, G, M_sun

gamma=1.4
mmol=2.3
Rgas=8.314e7
Rgasmu=Rgas/mmol
cp=gamma*Rgasmu/(gamma-1)
cv=cp/gamma

rr=20
T0 = 500
Sigma0 = 1700
p = 1.5
q = 0.25
Sigma = Sigma0 * rr**(-p)
Temperature = T0  * rr**(-q)
r=rr*au.cgs.value
cs2=Temperature*cp*(gamma-1.)
cs=np.sqrt(cs2)
Omega = np.sqrt(G.cgs.value*M_sun.cgs.value/r**3)
ukep = Omega * r
H = cs/Omega
h = H/r
print('Temperature =',Temperature,' K')
print('Sound speed =',cs/1e5,' km/s')
print('Scale Height=',H,' cm')
print('Aspect Ratio=',H/r)
print('Keplerian velocity=',ukep/1e5,' km/s')
print('Column density =',Sigma,' g/cm2')
Q = cs*Omega/(np.pi*G.cgs.value*Sigma)
print('Toomre Q=',Q) 


Gtilde = 0.1
G=6.67384e-8
Msun=2e33
au=1.49597871e13
r = 20 * au
Omega = np.sqrt(G*Msun)/r**1.5
rho0 = Gtilde * Omega**2 / (4*np.pi*G)
Temperature = 150. * (r/au)**(-3./7)
cs=np.sqrt(Temperature*cp*(gamma-1.))
Omega = np.sqrt(G*Msun/r**3)
ukep = Omega * r
H = cs/Omega
Sigma0 = rho0  * np.sqrt(2*np.pi) * H
Q = cs*Omega/(np.pi*G*Sigma0)
print('Temperature =',Temperature,' K')
print('Sound speed =',cs/1e5,' km/s')
print('Scale Height=',H/au,' au')
print('Aspect Ratio=',H/r)
print('Column density =',Sigma0,' g/cm2')
print('Volume density =',rho0,' g/cm3')
print('Toomre Q=',Q) 





