#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Plot curves as a parameter is varied
# Warren and Sear 2025/2026

import argparse
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from numpy import log as ln
from numpy import pi as π
from models import Model

parser = argparse.ArgumentParser(description='figure 3 in manuscript')
parser.add_argument('datafile', help='input data spreadsheet, *.ods, *.xlsx')
parser.add_argument('-Q', '--Qrange', default='1e-3,1e2', help='Q range in pL/s, default 1e-3,1e2')
parser.add_argument('-e', '--epsilon', default=1e-6, type=float, help='nearness to Qcrit, default 1e-6')
parser.add_argument('-f', '--frac', default=0.7, type=float, help='fraction of Qcrit, default 0.7')
parser.add_argument('-n', '--npt', default=80, type=int, help='number of points, default 80')
parser.add_argument('--dpi', default=72, type=int, help='resolution (dpi) for image output, default (for pdf) 72')
parser.add_argument('-v', '--verbose', action='count', default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

Q1, Q2 = np.array(eval(f'[{args.Qrange}]'))

codes = {'pore10k': '$10^4$', 'pore100k': '$10^5$', 'pore1M': '$10^6$'}
data = dict([(k, pd.read_excel(args.datafile, sheet_name=f'code={k}')) for k in codes])

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

fig, ax = plt.subplots(figsize=(6, 4), dpi=args.dpi)

ylims = 0.5, 1e3

ax.loglog(1e-3*Q, z1, color='tab:red',lw=lw, zorder=4) # red, saddle point
ax.loglog(1e-3*Q, z2, color='tab:orange', lw=lw, zorder=4) # orange, stable fixed point
ax.loglog(1e-3*Qc, zc, 'o', color='tab:brown', ms=ms, zorder=6) # bifurcation, black citcle

symbol = ['o', 's', 'D']
color = [f'tab:{c}' for c in ['red', 'green', 'blue']]

for i, k in enumerate(codes):
    df = data[k]
    ax.plot(df.Q, df.RMSD, symbol[i], color=color[i], label=codes[k])
    ax.errorbar(df.Q, df.RMSD, 2*df.std_err, fmt='.', color=color[i], capsize=3, capthick=2)

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend(loc='lower left', bbox_to_anchor=(0.05, 0.2),
          title='# steps', frameon=False, markerscale=1.3,
          title_fontsize=legend_fs, fontsize=legend_fs, labelspacing=0.5)

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
