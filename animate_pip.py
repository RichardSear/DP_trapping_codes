#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Animate pipette drift field and fixed points for a range of Q values
# Warren and Sear 2025/2026

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from models import Model

parser = argparse.ArgumentParser(description='figure 2 in manuscript')
parser.add_argument('-W', '--width', default=100.0, type=float, help='half width of plot in um, default 100')
parser.add_argument('-S', '--shift', default=50.0, type=float, help='shift right in um, default 50')
parser.add_argument('-Q', '--Qvals', default='5,15,5', help='range of Q values to use in pL/s, default 5,15,5')
parser.add_argument('-g', '--geom', action='store_true', help='use geomspace rather than linspace')
parser.add_argument('--dpi', default=150, type=int, help='resolution (dpi) for image output, default 150')
parser.add_argument('-v', '--verbose', action='count', default=0)
parser.add_argument('-t', '--template', default='frames/f%05d.png', help='template for frames, default frames/f%05d.png')
args = parser.parse_args()

Qvals = eval(f'[{args.Qvals}]')

w, s = args.width, args.shift
w1, w2 = -w+s, w+s

pipette = Model("pipette")

def drift(s, y):
    rvec = np.array([y[0], 0, y[1]])
    dxdt, _, dzdt = pipette.drift(rvec)
    dsdt = np.sqrt(dxdt**2 + dzdt**2)
    dxds, dzds = dxdt/dsdt, dzdt/dsdt
    return np.array([dxds, dzds])

x, z = np.mgrid[-w:w:100j, w1:w2:100j]
nz, nx = x.shape
ux = np.zeros((nz, nx))
uz = np.zeros((nz, nx))

tick_fs, label_fs = 12, 14
gen_lw, line_lw = 1.2, 1.2
xticks = [-50, 0, 50, 100, 150]
yticks = [-100, -50, 0, 50, 100]

fig, ax = plt.subplots(1, 1, figsize=(6, 6), dpi=args.dpi)

np_space = np.geomspace if args.geom else np.linspace

Qvals = np_space(*eval(args.Qvals))

for i, Q in enumerate(Qvals):

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
        
    ax.streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2,density=1.0) # drift streamlines

    if pipette.fixed_points is not None: # add markers for fixed points and separatrix
        ax.scatter(z1, 0, s=120, color='tab:orange', lw=3, marker='+', zorder=99)
        ax.scatter(z2, 0, s=80, color='tab:red', lw=4, marker='o', zorder=99)
        ax.plot(z_sep, x_sep, lw=3, label='sep', color='tab:red', zorder=19, ls='dashed')

    ax.plot([-w, 0], [0, 0], lw=4, c='k') # represent pipette

    ax.set_xlim(w1, w2)
    ax.set_ylim(-w, w)
    ax.set_aspect('equal')

    ax.tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)

    for spine in ax.spines:
        ax.spines[spine].set_linewidth(gen_lw)

    ax.set_xlabel(r'$z$ / µm', fontsize=label_fs)
    ax.set_xticks(xticks)

    ax.set_yticks(yticks)
    ax.set_ylabel(r'$x$ / µm', fontsize=label_fs, labelpad=-10)

    bbox = dict(boxstyle='round', fc='w', ls='') # lw=gen_lw)
    label = f'Q = {Q:5.2f} pL/s'
    ax.annotate(label, (-42, 77), fontsize=label_fs, bbox=bbox)

    frame = args.template % i
    plt.savefig(frame, bbox_inches='tight', pad_inches=0.05)
    print(f'frame ({label}) saved to {frame}')

    plt.cla()

print('To make a movie run something like:')
print(f'ffmpeg -framerate 15 -i {args.template} movie.mp4')

