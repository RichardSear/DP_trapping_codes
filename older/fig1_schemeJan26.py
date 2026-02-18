#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Python code to plot schematic of pipette flow field
Warren and Sear 2025/2026
'''

import numpy as np
from numpy import pi as π
import matplotlib.pyplot as plt
import final_flowfields
#
print(π)


drift = final_flowfields.drift_pipette_spherical

#
width=20
print('width of grid for streamplot in um',width)
#
fig, ax = plt.subplots(figsize=(4, 2.2))
# Set ticks to point inward for both axes
plt.tick_params(axis='both', which='both', direction='in')
Q = 1.0 # here value does not matter as just have streamlines
k=10.0
Γ = 1.0e-10
print('setting Γ to almost zero so drift field is just due to flow')
Ds = 1600.0
R1 = 1.0
α = 0.3
rcutoff=1.0

#
rstar = π*R1/α # where stokeslet and radial outflow match
v1 = Q / (π*R1**2) # flow speed (definition)
Pbyη = α*R1*v1 # from Secchi et al, should also = Q / r*
Pe = Pbyη / (4*π*Ds) # definition from Secchi et al
λ = Q / (4*π*Ds) # definition
λstar = λ / rstar # should be the same as Pe (salt)

kwargs = {"Q": Q, "k" : k,"Γ" : Γ,"rstar" : rstar,"v1" : v1, \
    "Pbyη" : Pbyη,"Pe" : Pe,"λ" : λ,"λstar" : λstar,"rcutoff" : rcutoff}

# 100 j means 100 points
x, z = np.mgrid[-width:width:100j, -width:width:100j]
print(x.shape)
ux=np.zeros((100,100))
uz=np.zeros((100,100))
r1=np.zeros(3)
salt=np.zeros((100,100))
for iz in range(0,100):
    for ix in range(0,100):
        r1[:]=x[iz,ix], 0.0, z[iz,ix]
        ux[iz,ix],dummy,uz[iz,ix]=drift(r1,Q,k,Γ,rstar,v1,Pbyη,Pe,λ,λstar,rcutoff)
#
        salt[iz,ix]=1.0/np.sqrt(z[iz,ix]**2+x[iz,ix]**2)
        salt[iz,ix]=np.log(salt[iz,ix])
        if(salt[iz,ix]>2): salt[iz,ix]
#        
ax.streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2,density=1.0)
ax.set_xlim(-width, width)
ax.set_ylim(-width, width)
# represent pipette
ax.plot([-width,0],[0.0,0.0],lw=4,c='k')
cf1=ax.contourf(z,x,salt,alpha=0.95, \
                 cmap='Greens')#,locator=ticker.LogLocator())
# Add labels and title
ax.set_aspect('equal')
#ax.tick_params(axis='both', which='major', labelsize=18)
ax.set_xlabel(r'$z$ $/\mathrm{\mu}$m',fontsize=14)
ax.set_ylabel(r'$x$ $/\mathrm{\mu}$m',fontsize=14)
ax.set_ylim([-10,10])

# Add a legend
#ax.legend(fontsize=14,loc='upper left')
plt.tight_layout()
plt.savefig('fig1.pdf')
plt.show()
