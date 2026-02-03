#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Python code to plot schematic of pipette flow field
# Warren and Sear 2025/2026

import argparse
import numpy as np
import matplotlib.pyplot as plt
from models import Model

parser = argparse.ArgumentParser(description='figure 1 in manuscript')
parser.add_argument('-W', '--width', default=50.0, type=float, help='half width of plot in um, default 50.0')
parser.add_argument('-H', '--height', default=25.0, type=float, help='half height of plot in um, default 25.0')
parser.add_argument('-v', '--verbose', action='count', default=0, help='increasing verbosity')
parser.add_argument('-n', '--no-show', action='store_true', help="don't display the figure")
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

pipette = Model("pipette")
pipette.update(Γ=0) # turn off DP contribution for the flow field

if args.verbose or args.no_show:
    print(pipette.info)

w, h = args.width, args.height

tick_fs, label_fs = 10, 12
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
        ux[iz,ix], _, uz[iz,ix] = pipette.drift(rvec)
        salt[iz,ix] = -np.log(np.sqrt(z[iz,ix]**2 + x[iz,ix]**2))

ax.streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2, density=0.5) # flow field
ax.contourf(z, x, salt, alpha=0.95, cmap='Greens') # salt concentration
ax.plot([-w, 0], [0, 0], lw=4, c='k') # represent pipette

ax.set_xlim(-w, w)
ax.set_ylim(-h, h)
ax.set_aspect('equal')

ax.tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)

for spine in ax.spines:
    ax.spines[spine].set_linewidth(gen_lw)

ax.set_xlabel(r'$z$ / µm', fontsize=label_fs)
ax.set_ylabel(r'$x$ / µm', fontsize=label_fs)#, labelpad=-2)

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
elif not args.no_show:
    plt.show()
