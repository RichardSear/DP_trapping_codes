#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Sampson flow and drift field for wall pore case
# Warren and Sear 2025/2026

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from models import Model

parser = argparse.ArgumentParser(description='figure 4 in manuscript')
parser.add_argument('-W', '--width', default=4.0, type=float, help='half width of plot in um, default 4')
parser.add_argument('-H', '--height', default=4.0, type=float, help='half width of plot in um, default 4')
parser.add_argument('-Q', '--Qvals', default='0.15,0.20', help='pair of Q values to use in pL/s, default 0.18,0.20')
parser.add_argument('-k', default=35, type=float, help='salt ratio k, default 35')
parser.add_argument('--dpi', default=72, type=int, help='resolution (dpi) for image output, default (for pdf) 72')
parser.add_argument("-v", "--verbose", action="count", default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

Qvals = eval(f'[{args.Qvals}]')

w = args.width

pore = Model("pore")

def drift(s, y):
    rvec = np.array([y[0], 0, y[1]])
    dxdt, _, dzdt = pore.drift(rvec)
    dsdt = np.sqrt(dxdt**2 + dzdt**2)
    dxds, dzds = dxdt/dsdt, dzdt/dsdt
    return np.array([dxds, dzds])

x, z = np.mgrid[-w:w:100j, 0:2*w:100j]
nz, nx = x.shape
ux = np.zeros((nz, nx))
uz = np.zeros((nz, nx))

tick_fs, label_fs = 12, 14
gen_lw, line_lw = 1.2, 1.2
xticks = [0, 2, 4, 6, 8]
yticks = [-4, -2, 0, 2, 4]

fig, ax = plt.subplots(1, 2, figsize=(6, 3.2), sharex=True, sharey=True, dpi=args.dpi)

for i, Q in enumerate(Qvals):

    pore.update(Q=Q, k=args.k)

    if args.verbose:
        print(pore.info)

    for iz in range(nz):
        for ix in range(nx):
            rvec = np.array([x[iz,ix], 0.0, z[iz,ix]])
            ux[iz,ix], _, uz[iz,ix] = pore.drift(rvec)

    ax[i].streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2, density=1.2) # drift streamlines

    if pore.fixed_points is not None:

        z1, z2 = pore.fixed_points

        sol = solve_ivp(lambda s, y: -drift(s, y), [0, 20], [0.01, z1], dense_output=True)
        s_max = np.max(sol.t)
        s = np.linspace(0, s_max, 50)
        y = sol.sol(s)
        x_sep, z_sep = y[0], y[1]
        x_sep = np.concatenate((-x_sep[::-1], x_sep))
        z_sep = np.concatenate((z_sep[::-1], z_sep))

        ax[i].scatter(z1, 0, s=80, color='tab:red', lw=4, marker='o', zorder=9)
        ax[i].scatter(z2, 0, s=120, color='tab:orange', lw=3, marker='+', zorder=9)
        ax[i].plot(z_sep, x_sep, lw=3, label='sep', color='tab:red', zorder=19, ls='dashed')

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

bbox = dict(boxstyle='round', fc='w', ls='') #lw=gen_lw)

for i, label in enumerate(['(a)', '(b)']):
    ax[i].annotate(label, (6.45, 3.1), fontsize=label_fs, bbox=bbox, zorder=10)
    #ax[i].text(0.8*w, 0.8*w, label, fontsize=label_fs, bbox=dict(facecolor='w', alpha=1.0, edgecolor='w'))

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
elif not args.verbose:
    plt.show()
