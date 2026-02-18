#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Plot individual endpoints for trajectories
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
parser.add_argument('dataset', help='input raw data spreadsheet, *.dat.gz')
parser.add_argument('-Q', '--Qrange', default='1e-4,1e2', help='Q range in pL/s, default 1e-4,1e2')
parser.add_argument('--Dp', default=2.0, type=float, help='particle diffusion coeff, default 2.0 um^2/s')
parser.add_argument('--dpi', default=72, type=int, help='resolution (dpi) for image output, default (for pdf) 72')
parser.add_argument('-j', '--justify', action='store_true', help='attempt to right-justify labels in legend')
parser.add_argument('-n', '--ntraj', default=20, type=int, help='number of trajectories per block, default 20')
parser.add_argument('-v', '--verbose', action='count', default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

Q1, Q2 = np.array(eval(f'[{args.Qrange}]'))

schema= {'k':float, 'Γ':float, 'Ds':float, 'Dp':float, 'R1':float, 
         'α':float, 'Q':float, 'rc':float, 't_final':float, 
         'ntrial':int, 'nsuccess':int, 't':float, 'Δt_final':float, 'Δr2':float, 
         'traj':int, 'block':int, 'ntraj':int, 'nblock':int, 'code':str}

df = pd.read_csv(args.dataset, sep='\t', names=schema.keys(), dtype=schema)

df['Δr'] = np.sqrt(df.Δr2)
ntrial_max = df.ntrial.max()
df['ntrial_frac'] = df.ntrial / ntrial_max

Dp = args.Dp
dfx = df[(df.Dp == Dp) & (df.traj < args.ntraj) & (df.Q > Q1) & (df.Q < Q2)]
t_final_max = dfx.t_final.max()
Δr_free = np.sqrt(6*Dp*t_final_max)

lw, ms = 2, 8
gen_lw, line_lw = 1.2, 1.2
tick_fs, label_fs, legend_fs = 12, 14, 12

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8), sharex=True, dpi=args.dpi)
renderer = fig.canvas.get_renderer() # used below to right-justify legend labels

kwargs = {'log_scale': True, 'native_scale': True, 'zorder': 2}
sb.stripplot(data=dfx, x='Q', y='Δr', hue='ntrial_frac', palette='autumn_r', ax=ax1, **kwargs)
sb.stripplot(data=dfx, x='Q', y='Δr', hue='t', palette='winter_r', ax=ax2, **kwargs)

kwargs = {'frameon': False, 'markerscale': 1.3,
          'title_fontsize': legend_fs, 'fontsize': legend_fs, 'labelspacing': 0.4}
ax1.legend(title='# steps / $10^5$', **kwargs)
legend = ax2.legend(title='time / s', **kwargs)

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

for ax in ax1, ax2:
    ax.axhline(Δr_free, ls=':', color='k', zorder=3)

ax2.set_xlim([Q1, Q2])
ax2.set_xticks([1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100],
               labels=['$10^{-4}$', '$10^{-3}$', r'$10^{-2}$', '0.1', '1', '10', '$10^2$'])
ax2.set_xlabel(r'Q / pL$\,$s$^{-1}$', fontsize=label_fs)

for ax in ax1, ax2:
    ax.minorticks_off()
    ax.tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)
    for spine in ax.spines:
        ax.spines[spine].set_linewidth(gen_lw)

for ax in ax1, ax2:
    ax.set_ylim(0.03, 3000)
    ax.set_yticks([0.1, 1, 10, 100, 1e3], labels=['0.1', '1', '10', '10$^2$', '10$^3$'])
    ax.set_ylabel('$\Delta r$ / µm', fontsize=label_fs)

for tick in ax2.xaxis.get_majorticklabels():
    tick.set_verticalalignment('bottom') # force the tick label alignment to the bottom ..

ax2.tick_params(axis='x', which='major', pad=20) # .. which then needs padding out

ax1.annotate('(a)', (20, 0.1), fontsize=label_fs)
ax2.annotate('(b)', (20, 0.1), fontsize=label_fs)

umsqpersec = r'µm$^2\,$s$^{-1}$' # ensure commonality between legend and axis label

for ax in ax1, ax2:
    ax.annotate('$D_p$ = {Dp}$\,${units}'.format(Dp=Dp, units=umsqpersec), (0.06, 20), fontsize=legend_fs)

plt.tight_layout()

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
elif not args.verbose:
    plt.show()
