#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Pipette drift field and fixed points for two Q values
# Warren and Sear 2025/2026

# Layout here massively helped by Claude Sonnet 5.5

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import blended_transform_factory
from scipy.integrate import solve_ivp
from models import Model

parser = argparse.ArgumentParser(description='figure 2 in manuscript')
parser.add_argument('-W', '--width', default=100.0, type=float, help='half width of plot in um, default 100')
parser.add_argument('-S', '--shift', default=50.0, type=float, help='shift right in um, default 50')
parser.add_argument('-Q', '--Qvals', default='10,100', help='pair of Q values to use in pL/s, default 10,100')
parser.add_argument('--dpi', default=72, type=int, help='resolution (dpi) for image output, default (for pdf) 72')
parser.add_argument("-v", "--verbose", action="count", default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
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

lw, ms = 2, 8
tick_fs, label_fs, legend_fs = 12, 14, 12
gen_lw, line_lw = 1.2, 1.2
xticks = [-50, 0, 50, 100, 150]
yticks = [-100, -50, 0, 50, 100]

# fig, ax = plt.subplots(1, 2, figsize=(6, 3.2), sharex=True, sharey=True, dpi=args.dpi)

fig = plt.figure(figsize=(6, 6), dpi=args.dpi)
gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 0.6], wspace=0.15, hspace=0.15)

# Top row: two square plots

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)

# Bottom row: one rectangular plot spanning both columns
ax3 = fig.add_subplot(gs[1, :])

ax = [ax1, ax2, ax3]

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
        
    ax[i].streamplot(z, x, uz, ux, linewidth=1.5, arrowsize=2,density=1.0) # drift streamlines

    if pipette.fixed_points is not None: # add markers for fixed points and separatrix
        ax[i].scatter(z1, 0, s=120, color='tab:orange', lw=3, marker='x', zorder=99)
        ax[i].scatter(z2, 0, s=80, color='tab:red', lw=4, marker='o', zorder=99)
        ax[i].plot(z_sep, x_sep, lw=3, label='sep', color='tab:red', zorder=19, ls='dashed')

    ax[i].plot([-w, 0], [0, 0], lw=4, c='k') # represent pipette

    ax[i].set_xlim(w1, w2)
    ax[i].set_ylim(-w, w)
    ax[i].set_aspect('equal')

    ax[i].set_xlabel(r'$z$ / µm', fontsize=label_fs)
    ax[i].set_xticks(xticks)

for ax in ax1, ax2, ax3:
    ax.minorticks_off()
    ax.tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)
    for spine in ax.spines:
        ax.spines[spine].set_linewidth(gen_lw)

ax1.set_yticks(yticks)
ax1.set_ylabel(r'$x$ / µm', fontsize=label_fs, labelpad=-10)

ax2.tick_params(labelleft=False)

bbox = dict(boxstyle='round', fc='w', ls='') # lw=gen_lw)

ax1.annotate('(a)', (-42, 77), fontsize=label_fs, bbox=bbox)
ax2.annotate('(b)', (-42, 77), fontsize=label_fs, bbox=bbox)

ax3.annotate('(c)', (30, 1.4), fontsize=label_fs)

# Now make the polar axis velocity plot

r = np.geomspace(0.01, 150, 160)
r = r[r > 0]

k = pipette.k
Γ = pipette.Γ
π = np.pi
Ds = pipette.Ds
α = pipette.α
Rt = pipette.R1
rstar = pipette.rstar

Q = 10 * 1e3
vt = Q/(π*Rt**2)
Pbyeta = α*Rt*vt
λ = Q/(4*π*Ds)
vr = -Γ*k*λ/(r*(r+k*λ)) + Q/(4*π*r**2) + Pbyeta/(4*π*r)
ax3.plot(r, vr, c='tab:blue', lw=lw, label='Q = 10 pL s$^{-1}$')
marker_sf = 0.8
ax3.scatter(z1, 0, s=marker_sf*120, color='tab:orange', lw=3, marker='x', zorder=99)
ax3.scatter(z2, 0, s=marker_sf*80, color='tab:red', lw=4, marker='o', zorder=99)

Qc = pipette.Qcrit
vt = Qc/(π*Rt**2)
Pbyeta = α*Rt*vt
λ = Qc/(4*π*Ds)
vr = -Γ*k*λ/(r*(r+k*λ)) + Qc/(4*π*r**2) + Pbyeta/(4*π*r)
ax3.plot(r, vr, lw=lw, c='tab:green', label='Q = Q$_c$ $\\simeq$ ' + f'{1e-3*Qc:0.2f}' + ' pL s$^{-1}$')
# the quadratic for the roots is z² − (kΓbyD − kλ* − 1)z + kλ* = 0 where z is in units of r*

kλ = k*Qc/(4*π*Ds) # this is now a scalar
zc = 0.5*rstar*(Γ*k/Ds-kλ/rstar-1) # bifurcation point solves 2z − (kΓbyD − kλ* − 1) = 0
ax3.scatter(zc, 0, s=marker_sf*80,  color='tab:brown', lw=4, marker='o', zorder=99)

Q = 100 * 1e3
vt = Q/(π*Rt**2)
Pbyeta = α*Rt*vt
λ = Q/(4*π*Ds)
vr = -Γ*k*λ/(r*(r+k*λ)) + Q/(4*π*r**2) + Pbyeta/(4*π*r)
ax3.plot(r, 1e-2*vr, lw=lw, c='tab:red', label='Q = 100 pL s$^{-1}$')
ax3.annotate('× 10$^{-2}$', (30, 0.5), color='tab:red', size=label_fs)

ax3.axhline(0, color='k', ls=':', lw=lw)

ax3.set_xlim(0, 140)
ax3.set_xticks(np.arange(0, 150, 20))

ax3.set_ylim(-1, 2)
ax3.set_yticks([-1, 0, 1, 2])

# ax3.legend()
legend = ax3.legend(title_fontsize=legend_fs, fontsize=legend_fs, labelspacing=0.5, frameon=False)

ax3.set_xlabel('$r$ / µm', size=label_fs)
ax3.set_ylabel('$u_r$ ($\\theta=0$) / µm s$^{-1}$', fontsize=label_fs, labelpad=-10)

fig.subplots_adjust(left=0.08, right=0.95, top=0.93, bottom=0.08)

# align left hand labels 

label_x = 0.04   # figure fraction from the left edge

for ax in ax1, ax3:
    tr = blended_transform_factory(fig.transFigure, ax.transAxes)
    ax.yaxis.set_label_coords(label_x, 0.58, transform=tr) # offset here shifts y-axis label upwards a bit
    
if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
elif not args.verbose:
    plt.show()
