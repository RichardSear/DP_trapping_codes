#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Python code to plot schematic of pipette flow field
# Warren and Sear 2025/2026

import argparse
import numpy as np
from numpy import pi as π
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='figure 1 in manuscript')
parser.add_argument('-W', '--width', default=50.0, type=float, help='width of plot in um, default 20.0')
parser.add_argument('-H', '--height', default=25.0, type=float, help='height of plot in um, default 20.0')
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

Q = 1.0 # here value does not matter as just have streamlines
k = 10.0
Ds = 1610.0
R1 = 1.0
α = 0.3
rc = 1.0

rstar = π*R1/α # where stokeslet and radial outflow match
v1 = Q / (π*R1**2) # flow speed (definition)
Pbyη = α*R1*v1 # from Secchi et al, should also = Q / r*
Pe = Pbyη / (4*π*Ds) # definition from Secchi et al
λ = Q / (4*π*Ds) # definition
λstar = λ / rstar # should be the same as Pe (salt)

def flow_field(rvec):
    x, y, z = rvec[:] # z is normal distance = r cosθ
    ρ = np.sqrt(x**2 + y**2) # in-plane distance = r sinθ
    r = np.sqrt(x**2 + y**2 + z**2) # radial distance from origin
    cosθ, sinθ = z/r, ρ/r # polar angle, as cos and sin
    sinθ_cosφ, sinθ_sinφ = x/r, y/r # avoids dividing by ρ
    ur = Q/(4*π*r**2) + Pbyη*cosθ/(4*π*r) # radial flow velocity
    ux = ur*sinθ_cosφ - Pbyη*sinθ_cosφ*cosθ/(8*π*r) # radial components
    uy = ur*sinθ_sinφ - Pbyη*sinθ_sinφ*cosθ/(8*π*r) # avoiding dividing by ρ
    uz = ur*cosθ + Pbyη*sinθ**2/(8*π*r) # normal z-component of velocity
    return np.zeros_like(rvec) if r < rc else np.array((ux, uy, uz))

w, h = args.width, args.height

gen_lw, line_lw = 1.2, 1.2
fig, ax = plt.subplots(figsize=(4, 2.2))

# 100j means 100 points
x, z = np.mgrid[-h:h:100j, -w:w:50j]
nz, nx = x.shape

ux = np.zeros((nz, nx))
uz = np.zeros((nz, nx))
salt = np.zeros((nz, nx))

for iz in range(nz):
    for ix in range(nx):
        rvec = np.array([x[iz,ix], 0.0, z[iz,ix]])
        ux[iz,ix], _, uz[iz,ix] = flow_field(rvec)
        salt[iz,ix] = -np.log(np.sqrt(z[iz,ix]**2 + x[iz,ix]**2))

ax.streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2, density=0.5)
ax.contourf(z, x, salt, alpha=0.95, cmap='Greens')
ax.plot([-w, 0], [0, 0], lw=4, c='k') # represent pipette

ax.set_xlim(-w, w)
ax.set_ylim(-h, h)
ax.set_aspect('equal')

ax.tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=10)

#ax.set_xlabel('$A_\\mathrm{wall}$', fontsize=18)
#ax.set_ylabel('$\\gamma_s$', rotation=0, fontsize=18)
#ax.xaxis.set_label_coords(0.6, -0.07)
#ax.yaxis.set_label_coords(-0.1, 0.85)

for spine in ax.spines:
    ax.spines[spine].set_linewidth(gen_lw)

ax.set_xlabel(r'$z$ / µm', fontsize=10)
ax.set_ylabel(r'$x$ / µm', fontsize=10, labelpad=-2)

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
else:
    plt.show()
