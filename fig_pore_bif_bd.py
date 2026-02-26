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
parser.add_argument('datafile', help='input data spreadsheet proforma, *.ods, *.xlsx')
parser.add_argument('-D', '--Dpvals', default='1,2,5', help='set of Dp values to use, default 1,2,5')
parser.add_argument('-Q', '--Qrange', default='1e-3,1e2', help='Q range in pL/s, default 1e-3,1e2')
parser.add_argument('-y', '--ylims', default='0.3,5e3', help='y limits, default 0.3,5e3')
parser.add_argument('-k', '--kvals', default='35,30,25', help='salt ratios k, default 35,30,25')
parser.add_argument('-a', '--asymp', action='store_true', help='show the asymptotic fixed point')
parser.add_argument('-e', '--epsilon', default=1e-6, type=float, help='nearness to Qcrit, default 1e-6')
parser.add_argument('-f', '--frac', default=0.7, type=float, help='fraction of Qcrit, default 0.7')
parser.add_argument('-n', '--npt', default=80, type=int, help='number of points, default 80')
parser.add_argument('--dpi', default=72, type=int, help='resolution (dpi) for image output, default (for pdf) 72')
parser.add_argument('-v', '--verbose', action='count', default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

Q1, Q2 = np.array(eval(f'[{args.Qrange}]'))
Dpvals = np.array(eval(f'[{args.Dpvals}]'), dtype=float)
kvals = eval(f'[{args.kvals}]')

pore = Model('pore')

lw, ms = 2, 8
gen_lw, line_lw = 1.2, 1.2
tick_fs, label_fs, legend_fs = 12, 14, 12
umsqpersec = r'µm$^2\,$s$^{-1}$' # ensure commonality between legend and axis label

symbol = ['o', 's', 'D']
color = [f'tab:{c}' for c in ['green', 'blue', 'purple']]

ylims = np.array(eval(f'[{args.ylims}]'))

fig, ax = plt.subplots(3, 1, figsize=(6, 8.2), sharex=True, dpi=args.dpi)

for i, k in enumerate(kvals):

    datafile = args.datafile.format(k=k)
    
    data = dict([(Dp, pd.read_excel(datafile, sheet_name=f'Dp={Dp}')) for Dp in Dpvals])

    pore.update(k=k)

    if args.verbose:
        print(pore.info)

    # pick up these parameter values from the model
    k, Γ, Ds = pore.k, pore.Γ, pore.Ds
    ΓkbyDs, R1, rc = pore.ΓkbyDs, pore.R1, pore.rc

    # the quadratic for the roots is (1/2)(kΓ/Ds−3) z² − 3kλ z + (1/2)(kΓ/Ds) R1² = 0

    if ΓkbyDs > 3:
        Qc, Qx = pore.Qcrit, pore.Qcrit/args.frac
        Qa = np.geomspace(Qc, Qx, args.npt)
        Qb = np.geomspace(Qx, 1e3*Q2, args.npt) # convert Q2 to um^3/sec
        Q = np.concatenate([Qa, Qb[1:]])
        kλ = k*Q/(4*π*Ds) # this will be an array, as also the things below
        Δ = 9*kλ**2 - R1**2*ΓkbyDs*(ΓkbyDs - 3) # the discriminant
        z1 = (3*kλ - np.sqrt(Δ)) / (ΓkbyDs - 3) # lower root (saddle)
        z2 = (3*kλ + np.sqrt(Δ)) / (ΓkbyDs - 3) # upper root (stable fixed point)
        kλ = k*Qc/(4*π*Ds) # this is a scalar
        zc = 3*kλ/(ΓkbyDs-3) # bifurcation point, solves (kΓ/Ds−3) z − 3kλ = 0
    else:
        Q = np.geomspace(1e3*Q1, 1e3*Q2, args.npt) # convert Q2 to um^3/sec
        kλ = k*Q/(4*π*Ds) # this will be an array, as also the things below
        Δ = 9*kλ**2 - R1**2*ΓkbyDs*(ΓkbyDs - 3) # the discriminant
        z1 = (3*kλ - np.sqrt(Δ)) / (ΓkbyDs - 3) # upper root (saddle)


    Qa = 1e3*np.array([Q1, Q2])
    za = 3*k*Qa/(2*π*(Γ*k-3*Ds))

    if ΓkbyDs > 3:
        ax[i].loglog(1e-3*Q, z2, color='tab:orange', lw=lw, zorder=5) # orange, stable fixed point
        ax[i].loglog(1e-3*Q, z1, color='tab:red', lw=lw, zorder=4) # red, saddle point
        ax[i].loglog(1e-3*Qc, zc, 'o', color='tab:brown', ms=ms, zorder=6) # bifurcation, black citcle
    else:
        ax[i].loglog(1e-3*Q, z1, color='tab:red', lw=lw, zorder=4) # red, saddle point

    if args.asymp:
        ax[i].loglog(1e-3*Qa, za, color='tab:orange', ls='--', lw=lw, zorder=4) # orange, stable fp (approximation)


    for j, Dp in enumerate(Dpvals):
        df = data[Dp]
        ax[i].plot(df.Q, df.RMSD, symbol[j], color=color[j], label=f'{Dp}')
        ax[i].errorbar(df.Q, df.RMSD, 2*df.std_err, fmt='.', color=color[j], capsize=3, capthick=2)

    ax[i].set_xscale('log')
    ax[i].set_yscale('log')

    ax[i].legend(loc='lower left', bbox_to_anchor=(0.70, 0.03),
                 title='$D_p$ / {units}'.format(units=umsqpersec), frameon=False, markerscale=1.3,
                 title_fontsize=legend_fs, fontsize=legend_fs, labelspacing=0.5)

    ax[i].set_ylim(*ylims)
    yticks = np.array([1, 10, 100, 1e3, 1e4])
    ylabels = np.array(['1', '10', '$10^2$', '$10^3$', '$10^4$'])
    ycut = (yticks < ylims[0]) | (yticks > ylims[1])
    ax[i].set_yticks(yticks[~ycut], labels=list(ylabels[~ycut]))
    ax[i].set_ylabel('RMSD / µm', fontsize=label_fs)

    ax[i].minorticks_off()
    ax[i].tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)
    for spine in ax[i].spines:
        ax[i].spines[spine].set_linewidth(gen_lw)


    ax[i].set_xlim(Q1, Q2)
    xticks = np.array([1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100])
    xlabels = np.array(['$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '0.1', '1', '10', '$10^2$'])
    xcut = (xticks < Q1) | (xticks > Q2)
    
ax[2].set_xticks(xticks[~xcut], labels=list(xlabels[~xcut]))
ax[2].set_xlabel(r'Q / pL$\,$s$^{-1}$', fontsize=label_fs)

for tick in ax[2].xaxis.get_majorticklabels():
    tick.set_verticalalignment('bottom') # force the tick label alignment to the bottom ..

ax[2].tick_params(axis='x', which='major', pad=20) # .. which then needs padding out

for i, label in enumerate(['(a)', '(b)', '(c)']):
    ax[i].annotate(label, (1.8e-3, 1.3e3), fontsize=label_fs, zorder=10)

ax[0].annotate(f'k = {kvals[0]}', (3e-3, 10), fontsize=label_fs, zorder=10)
ax[1].annotate(f'k = {kvals[1]}', (4e-3, 13), fontsize=label_fs, zorder=10)
ax[2].annotate(f'k = {kvals[2]}', (5e-3, 10), fontsize=label_fs, zorder=10)

plt.tight_layout()

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
else:
    plt.show()
