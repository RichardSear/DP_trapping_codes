#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# BD data for wall pore case
# Warren and Sear 2025/2026

import argparse
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from numpy import pi as π
from models import Model

parser = argparse.ArgumentParser(description='figure 5 in manuscript')
parser.add_argument('datafile', help='input data spreadsheet, *.ods, *.xlsx')
parser.add_argument('-D', '--Dpvals', default='10,50', help='set of Dp values to use, default 10,50')
parser.add_argument('-Q', '--Qrange', default='1e-3,1e2', help='Q range in pL/s, default 1e-4,1e2')
parser.add_argument('-s', '--saddle', action='store_true', help='show the saddle point also')
parser.add_argument('-e', '--epsilon', default=1e-6, type=float, help='nearness to Qcrit, default 1e-6')
parser.add_argument('-f', '--frac', default=0.7, type=float, help='fraction of Qcrit, default 0.7')
parser.add_argument('-n', '--npt', default=80, type=int, help='number of points, default 80')
parser.add_argument('--dpi', default=72, type=int, help='resolution (dpi) for image output, default (for pdf) 72')
parser.add_argument('-j', '--justify', action='store_true', help='attempt to right-justify labels in legend')
parser.add_argument('-v', '--verbose', action='count', default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

Q1, Q2 = np.array(eval(f'[{args.Qrange}]'))

Dpvals = np.array(eval(f'[{args.Dpvals}]'), dtype=float)
data = dict([(Dp, pd.read_excel(args.datafile, sheet_name=f'Dp={Dp}')) for Dp in Dpvals])

pore = Model('pore')

if args.verbose:
    print(pore.info)

# pick up these parameter values from the model
k, Γ, Ds = pore.k, pore.Γ, pore.Ds
ΓkbyDs, R1, rc = pore.ΓkbyDs, pore.R1, pore.rc

Qc, Qx = pore.Qcrit, pore.Qcrit/args.frac
Qa = np.geomspace(Qc, Qx, args.npt)
Qb = np.geomspace(Qx, 1e3*Q2, args.npt) # convert Q2 to um^3/sec
Q = np.concatenate([Qa, Qb[1:]])

# the quadratic for the roots is (1/2)(kΓ/Ds−3) z² − 3kλ z + (1/2)(kΓ/Ds) R1² = 0

kλ = k*Q/(4*π*Ds) # this will be an array, as also the things below
Δ = 9*kλ**2 - R1**2*ΓkbyDs*(ΓkbyDs - 3) # the discriminant
z1 = (3*kλ - np.sqrt(Δ)) / (ΓkbyDs - 3) # lower root (saddle)
z2 = (3*kλ + np.sqrt(Δ)) / (ΓkbyDs - 3) # upper root (stable fixed point)
kλ = k*Qc/(4*π*Ds) # this is a scalar
zc = 3*kλ/(ΓkbyDs-3) # bifurcation point, solves (kΓ/Ds−3) z − 3kλ = 0

lw, ms = 2, 8
gen_lw, line_lw = 1.2, 1.2
tick_fs, label_fs, legend_fs = 12, 14, 12
umsqpersec = r'µm$^2\,$s$^{-1}$' # ensure commonality between legend and axis label

fig, ax = plt.subplots(figsize=(6, 3.2), dpi=args.dpi)
renderer = fig.canvas.get_renderer() # used below to right-justify legend labels

ylims = 0.3, 3e3

Qa = 1e3*np.array([Q1, Q2])
za = 3*k*Qa/(2*π*(Γ*k-3*Ds))

ax.loglog(1e-3*Qa, za, color='tab:orange', ls='--', lw=lw, zorder=4) # orange, stable fixed point (approximation)

if args.saddle:
    ax.loglog(1e-3*Q, z2, color='tab:orange', lw=lw, zorder=5) # orange, stable fixed point
    ax.loglog(1e-3*Q, z1, color='tab:red',lw=lw, zorder=4) # red, saddle point
    ax.loglog(1e-3*Qc, zc, 'o', color='tab:brown', ms=ms, zorder=6) # bifurcation, black citcle

symbol = ['o', 's', 'D']
color = [f'tab:{c}' for c in ['green', 'blue', 'purple']]

for i, Dp in enumerate(Dpvals):
    df = data[Dp]
    ax.plot(df.Q, df.RMSD, symbol[i], color=color[i], label=f'{Dp}')
    ax.errorbar(df.Q, df.RMSD, 2*df.std_err, fmt='.', color=color[i], capsize=3, capthick=2)

ax.set_xscale('log')
ax.set_yscale('log')

legend = ax.legend(loc='lower left', bbox_to_anchor=(0.05, 0.15),
                   title='$D_p$ / {units}'.format(units=umsqpersec), frameon=False, markerscale=1.3,
                   title_fontsize=legend_fs, fontsize=legend_fs, labelspacing=0.5)

# The following right-justifies the legend texts, from
# https://stackoverflow.com/questions/7936034/text-alignment-in-a-matplotlib-legend
# Doesn't work properly when plot saved as PDF, only as PNG with dpi specified in
# the subplots() call; see http://github.com/matplotlib/matplotlib/issues/15497
# A fix for PDF output is to specify the dpi as 72.

if args.justify:
    legend_txts = legend.get_texts()
    w_max = max([txt.get_window_extent(renderer).width for txt in legend_txts])
    for txt in legend_txts:
        txt.set_ha('right')  # ha is alias for horizontalalignment
        Δw = w_max - txt.get_window_extent().width
        txt.set_position((Δw, 0))

ax.set_xlim(Q1, Q2)
xticks = [1e-3, 1e-2, 0.1, 1, 10, 100]
xlabels = ['$10^{-3}$', '$10^{-2}$', '0.1', '1', '10', '$10^2$']
ax.set_xticks(xticks, labels=xlabels)
ax.set_xlabel(r'Q / pL$\,$s$^{-1}$', fontsize=label_fs)

ax.set_ylim(*ylims)
ax.set_yticks([1, 10, 100, 1e3], labels=['1', '10', '$10^{2}$', '$10^{3}$'])
ax.set_ylabel('RMSD / µm', fontsize=label_fs)

ax.minorticks_off()
ax.tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)
for spine in ax.spines:
    ax.spines[spine].set_linewidth(gen_lw)

for tick in ax.xaxis.get_majorticklabels():
    tick.set_verticalalignment('bottom') # force the tick label alignment to the bottom ..

ax.tick_params(axis='x', which='major', pad=20) # .. which then needs padding out

plt.tight_layout()

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
else:
    plt.show()
