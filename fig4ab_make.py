#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Python code to plot schematic of pipette flow field
# Warren and Sear 2025/2026

import argparse
import numpy as np
import matplotlib.pyplot as plt
from models import Model

parser = argparse.ArgumentParser(description='figure 4ab in manuscript')
parser.add_argument('-W', '--width', default=4.0, type=float, help='half width of plot in um, default 4')
parser.add_argument('-H', '--height', default=4.0, type=float, help='half width of plot in um, default 4')
parser.add_argument('-P', '--pars', default='0,1;150,0.8', help='pair of Γ, Q values, default 0,1;150,0.8')
parser.add_argument("-v", "--verbose", action="count", default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

w = args.width

tick_fs, label_fs = 12, 14
gen_lw, line_lw = 1.2, 1.2
xticks = [0, 2, 4, 6, 8]
yticks = [-4, -2, 0, 2, 4]

fig, ax = plt.subplots(1, 2, figsize=(6, 3.2), sharex=True, sharey=True)

pore = Model("pore")

pars = [eval(f'[{x}]') for x in args.pars.split(';')] # split on ';', then convert to lists

x, z = np.mgrid[-w:w:100j, 0:2*w:100j]
nz, nx = x.shape
ux = np.zeros((nz, nx))
uz = np.zeros((nz, nx))

for i, (Γ, Q) in enumerate(pars):

    pore.update(Q=Q, Γ=Γ)

    if args.verbose:
        print(pore.info)

    for iz in range(nz):
        for ix in range(nx):
            rvec = np.array([x[iz,ix], 0.0, z[iz,ix]])
            ux[iz,ix], _, uz[iz,ix] = pore.drift(rvec)

    ax[i].streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2, density=1.2) # drift streamlines
#    ax[i].streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2, density=1.5)

    if pore.fixed_points is not None:
        z1, z2 = pore.fixed_points
        # ax[i].scatter(z1, 0, s=80, color='red', lw=4, marker='o', zorder=99)
        ax[i].scatter(z2, 0, s=120, color='magenta', lw=3, marker='+', zorder=99)

    ax[i].plot([0, 0], [pore.R1, w], lw=6, c='k') # represent pore ..
    ax[i].plot([0, 0], [-pore.R1, -w], lw=6, c='k') # .. other side

    ax[i].set_xlim(0, 2*w)
    ax[i].set_ylim(-w, w)
    ax[i].set_aspect('equal')

    ax[i].tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)

    for spine in ax[i].spines:
        ax[i].spines[spine].set_linewidth(gen_lw)

    ax[i].set_xlabel(r'$z$ / µm', fontsize=label_fs)
    ax[i].set_xticks(xticks)
    

ax[0].set_yticks(yticks)
ax[0].set_ylabel(r'$x$ / µm', fontsize=label_fs, labelpad=0)

for i, label in enumerate(['a', 'b']):
    ax[i].annotate(f'({label})', (6.45, 3.1), fontsize=label_fs, backgroundcolor='w')
    #ax[i].text(0.8*w, 0.8*w, label, fontsize=label_fs, bbox=dict(facecolor='w', alpha=1.0, edgecolor='w'))

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
elif not args.verbose:
    plt.show()
