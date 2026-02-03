#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Python code to plot schematic of pipette flow field
# Warren and Sear 2025/2026

import argparse
import numpy as np
import matplotlib.pyplot as plt
from models import Model

parser = argparse.ArgumentParser(description='figure 1 in manuscript')
parser.add_argument('-W', '--width', default=50.0, type=float, help='width of plot in um, default 20.0')
parser.add_argument('-H', '--height', default=25.0, type=float, help='height of plot in um, default 20.0')
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

pipette = Model("pipette")
pipette.Gamma = 0
pipette.update()

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
        ux[iz,ix], _, uz[iz,ix] = pipette.flow_field(rvec)
        salt[iz,ix] = -np.log(np.sqrt(z[iz,ix]**2 + x[iz,ix]**2))

ax.streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2, density=0.5) # flow field
ax.contourf(z, x, salt, alpha=0.95, cmap='Greens') # salt concentration
ax.plot([-w, 0], [0, 0], lw=4, c='k') # represent pipette

ax.set_xlim(-w, w)
ax.set_ylim(-h, h)
ax.set_aspect('equal')

ax.tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=10)

for spine in ax.spines:
    ax.spines[spine].set_linewidth(gen_lw)

ax.set_xlabel(r'$z$ / µm', fontsize=10)
ax.set_ylabel(r'$x$ / µm', fontsize=10, labelpad=-2)

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
else:
    plt.show()
