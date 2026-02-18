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
parser.add_argument('-D', '--Dpvals', default='2,20,50', help='set of Dp values to use, default 2,20,50')
parser.add_argument('-Q', '--Qrange', default='1e-4,1e2', help='Q range in pL/s, default 1e-4,1e2')
parser.add_argument('-e', '--epsilon', default=1e-6, type=float, help='nearness to Qcrit, default 1e-6')
parser.add_argument('-f', '--frac', default=0.7, type=float, help='fraction of Qcrit, default 0.7')
parser.add_argument('-t', '--t-final', default=3600, type=float, help='duration, default 3600 sec')
parser.add_argument('-n', '--npt', default=80, type=int, help='number of points, default 80')
parser.add_argument('--dpi', default=72, type=int, help='resolution (dpi) for image output, default (for pdf) 72')
parser.add_argument('-j', '--justify', action='store_true', help='attempt to right-justify labels in legend')
parser.add_argument('-v', '--verbose', action='count', default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

Dpvals = [eval(x.split('=')[1]) for x in pd.ExcelFile(args.datafile).sheet_names]
data = dict([(Dp, pd.read_excel(args.datafile, sheet_name=f'Dp={Dp}')) for Dp in Dpvals])

pip = Model('pipette')

if args.verbose:
    print(pip.info)

Q1, Q2 = 1e3 * np.array(eval(f'[{args.Qrange}]')) # convert to um^3/sec
Qx, Qc = args.frac * pip.Qcrit, pip.Qcrit
Qa = np.geomspace(Q1, Qx, args.npt)
Qb = Qc - np.geomspace(Qc-Qx, args.epsilon, args.npt)
Q = np.concatenate([Qa, Qb[1:]])

# pick up these parameter values from the model
k, Γ, Ds = pip.k, pip.Γ, pip.Ds,
rstar, α, R1, rc = pip.rstar, pip.α, pip.R1, pip.rc

# the quadratic for the roots is z² − (kΓbyD − kλ* − 1)z + kλ* = 0 where z is in units of r*

kλ = k*Q/(4*π*Ds) # this will be an array, as also the things below
b = Γ*k/Ds - kλ/rstar - 1 # the equation is r^2 − br + c = 0
Δ = b**2 - 4*kλ/rstar # discriminant of above
z1 = 0.5*rstar*(b - np.sqrt(Δ)) # lower root (stable fixed point)
z2 = 0.5*rstar*(b + np.sqrt(Δ)) # upper root (saddle)
kλ = k*Qc/(4*π*Ds) # this is now a scalar
zc = 0.5*rstar*(Γ*k/Ds-kλ/rstar-1) # bifurcation point solves 2z − (kΓbyD − kλ* − 1) = 0

# The drift field uz = Γ d(ln c)/dz + Q/4πz² + P/4πηz
# This integrates to action = - Γ ln(c) + Q/4πz - P/4πη ln(z)
# Note the sign comes from integrating -uz

def S(z, Q):
    v1 = Q / (π*R1**2) # flow speed (definition)
    Pbyη = α*R1*v1 # from Secchi et al
    kλ = k*Q/(4*π*Ds)
    S = - Γ*ln(kλ/z + 1) + Q/(4*π*z) - Pbyη*ln(z)/(4*π)
    return S

ΔS = np.array([(S(z2, Q) - S(max(z1, rc), Q)) for z1, z2, Q in zip(z1, z2, Q)])

qc_ser = pd.Series([1e-3 * np.interp(10*Dp, ΔS[Q<1e3], Q[Q<1e3], right=np.nan) for Dp in Dpvals],
                   index=Dpvals).dropna()

if args.verbose:
    print(qc_ser)

lw, ms = 2, 8
gen_lw, line_lw = 1.2, 1.2
tick_fs, label_fs, legend_fs = 12, 14, 12
umsqpersec = r'µm$^2\,$s$^{-1}$' # ensure commonality between legend and axis label

gs_kw = {'height_ratios': [2, 1]}
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8), sharex=True, dpi=args.dpi, gridspec_kw=gs_kw)
renderer = fig.canvas.get_renderer() # used below to right-justify legend labels

ylims1 = 1, 3e3
ylims2 = 0.1, 1e3

symbol = ['o', 's', 'D', '<', '^', '>', 'v']
color = [f'tab:{c}' for c in ['red', 'orange', 'olive', 'green', 'blue', 'purple', 'pink']]

for i, Dp in enumerate(Dpvals):
    df = data[Dp]
    freed = np.sqrt(6*Dp*args.t_final)
    ax1.axhline(freed, ls=':', color=color[i], zorder=3)
    ax1.plot(df.Q, df.RMSD, symbol[i], color=color[i], label=f'{Dp}')
    ax1.errorbar(df.Q, df.RMSD, 2*df.std_err, fmt='.', color=color[i], capsize=3, capthick=2)

for ax, ylims in (ax1, ylims1), (ax2, ylims2):
    for i in range(1, qc_ser.size):
        Qc1, Qc2 = qc_ser.iloc[i-1], qc_ser.iloc[i]
        ax.fill_betweenx(ylims, [Qc1]*2, [Qc2]*2, color=color[i-1], alpha=0.2)
    Qc1, Qc2 = qc_ser.iloc[i], 1e-3*Qc
    ax.fill_betweenx(ylims, [Qc1]*2, [Qc2]*2, color=color[i], alpha=0.2)

for ax in ax1, ax2:
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-3*Q1, 1e-3*Q2)

legend = ax1.legend(loc='lower left', bbox_to_anchor=(0.0, 0.03),
                    title='$D_p$ / {units}'.format(units=umsqpersec),
                    frameon=False, #markerscale=1.3,
                    title_fontsize=legend_fs, fontsize=legend_fs, labelspacing=0.5)

if args.justify:
    legend_txts = legend.get_texts()
    w_max = max([txt.get_window_extent(renderer).width for txt in legend_txts])
    for txt in legend_txts:
        txt.set_ha('right')  # ha is alias for horizontalalignment
        Δw = w_max - txt.get_window_extent().width
        txt.set_position((Δw, 0))

ax1.set_ylim(*ylims1)
ax1.set_yticks([1, 10, 100, 1e3], labels=['1', '10', '$10^{2}$', '$10^{3}$'])
ax1.set_ylabel('RMSD / µm', fontsize=label_fs)

for i, Dp in enumerate(qc_ser.index):
    ax2.axhline(10*Dp, ls=':', color=color[i], lw=lw)

ax2.loglog(1e-3*Q, ΔS, color='peru', lw=lw)
ax2.set_ylim(*ylims2)
ax2.set_yticks([0.1, 10, 1e3], labels=['0.1', '10', r'10$^3$'])

xticks = [1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100]
xlabels = ['$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '0.1', '1', '10', '$10^2$']
ax2.set_xticks(xticks, labels=xlabels)

ax2.set_xlabel(r'Q / pL$\,$s$^{-1}$', fontsize=label_fs)
ax2.set_ylabel(r'$\Delta\mathfrak{{S}}$ / {units}'.format(units=umsqpersec),
               fontsize=label_fs) # double {{..}} to escape the braces in the format string

for ax in ax1, ax2:
    ax.minorticks_off()
    ax.tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)
    for spine in ax.spines:
        ax.spines[spine].set_linewidth(gen_lw)

for tick in ax2.xaxis.get_majorticklabels():
    tick.set_verticalalignment('bottom') # force the tick label alignment to the bottom ..

ax2.tick_params(axis='x', which='major', pad=20) # .. which then needs padding out

plt.tight_layout()

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
else:
    plt.show()
