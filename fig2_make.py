#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Python code to plot schematic of pipette flow field
# Warren and Sear 2025/2026

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from models import Model

parser = argparse.ArgumentParser(description='figure 1 in manuscript')
parser.add_argument('-W', '--width', default=100.0, type=float, help='width of plot in um, default 2100')
parser.add_argument('-Q', '--Qvals', default='10,100', help='pair of Q values to use in pL/s, default 10,100')
parser.add_argument("-v", "--verbose", action="count", default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

w = args.width

tick_fs, label_fs = 12, 14
gen_lw, line_lw = 1.2, 1.2
ticks = [-100, -50, 0, 50, 100]

fig, ax = plt.subplots(1, 2, figsize=(6, 3.2), sharex=True, sharey=True)

pipette = Model("pipette")

def drift(s, y):
    rvec = np.array([y[0], 0, y[1]])
    dxdt, _, dzdt = pipette.drift(rvec)
    dsdt = np.sqrt(dxdt**2 + dzdt**2)
    dxds, dzds = dxdt/dsdt, dzdt/dsdt
    return np.array([dxds, dzds])

x, z = np.mgrid[-w:w:100j, -w:w:100j]
nz, nx = x.shape
ux = np.zeros((nz, nx))
uz = np.zeros((nz, nx))

for i, Q in enumerate(eval(f'[{args.Qvals}]')):

    pipette.update(Q=Q)

    if args.verbose:
        print(pipette.info)

    if pipette.fixed_points is not None:
        z1, z2 = pipette.fixed_points
        sol = solve_ivp(lambda s, y: -drift(s, y), [0, 100], [0.01, z2], dense_output=True)
        s_max = np.max(sol.t)
        s = np.linspace(0, s_max, 20)
        y = sol.sol(s)
        x_sep, z_sep = y[0], y[1]
        x_sep = np.concatenate((-x_sep[::-1], x_sep))
        z_sep = np.concatenate((z_sep[::-1], z_sep))

    for iz in range(nz):
        for ix in range(nx):
            rvec = np.array([x[iz,ix], 0.0, z[iz,ix]])
            ux[iz,ix], _, uz[iz,ix] = pipette.drift(rvec)
        
    ax[i].streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2,density=1.0) # drift streamlines

    if pipette.fixed_points is not None: # add markers for fixed points and separatrix
        ax[i].scatter(z1, 0, s=120, color='magenta', lw=3, marker='+', zorder=99)
        ax[i].scatter(z2, 0, s=80, color='red', lw=4, marker='o', zorder=99)
        ax[i].plot(z_sep, x_sep, lw=3, label='sep', color='red', zorder=19, ls='dashed')

    ax[i].plot([-w, 0], [0, 0], lw=4, c='k') # represent pipette

    ax[i].set_xlim(-w, w)
    ax[i].set_ylim(-w, w)
    ax[i].set_aspect('equal')

    ax[i].tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)

    for spine in ax[i].spines:
        ax[i].spines[spine].set_linewidth(gen_lw)

    ax[i].set_xlabel(r'$z$ / µm', fontsize=label_fs)
    ax[i].set_xticks(ticks)
    

ax[0].set_yticks(ticks)
ax[0].set_ylabel(r'$x$ / µm', fontsize=label_fs, labelpad=-10)

for i, label in enumerate(['(a)', '(b)']):
    ax[i].annotate(label, (-92, 77), fontsize=label_fs, backgroundcolor='w')

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
else:
    plt.show()
