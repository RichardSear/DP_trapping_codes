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

def range_str(v, vals): # convert a list of values to a singleton, several values, or a range
    s = ', '.join([str(x) for x in vals]) if len(vals) < 10 else '--'.join([str(f(vals)) for f in [min, max]])
    return v, '  '+s, f'{len(vals):10}'

parser = argparse.ArgumentParser(description='figure 3 in manuscript')
parser.add_argument('dataset', help='input raw data spreadsheet, *.dat.gz')
parser.add_argument('--Dp', default=2.0, type=float, help='particle diffusion coeff, default 2.0 um^2/s')
parser.add_argument('-n', '--ntraj', default=20, type=int, help='number of trajectories per block, default 20')
parser.add_argument('-v', '--verbose', action='count', default=0)
parser.add_argument('-d', '--describe', action='store_true', help='print a summary of the columns in the raw data')
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

schema= {'k':float, 'Γ':float, 'Ds':float, 'Dp':float, 'R1':float, 
         'α':float, 'Q':float, 'rc':float, 't_final':float, 
         'ntrial':int, 'nsuccess':int, 't':float, 'Δt_final':float, 'Δr2':float, 
         'traj':int, 'block':int, 'ntraj':int, 'nblock':int, 'code':str}

df = pd.read_csv(args.dataset, sep='\t', names=schema.keys(), dtype=schema)

if args.describe:
    dff = pd.DataFrame([range_str(v, df[v].unique()) for v in df.columns], columns=['column', 'range', 'count'])
    header_row = pd.DataFrame(index=[-1], columns=dff.columns)
    dff = pd.concat([header_row, dff])
    dff.loc[-1] = dff.columns
    print('Dataset', args.dataset)
    print('\n'.join(dff.to_string(justify='left', index=False).split('\n')[1:]))
    exit()

df['Δr'] = np.sqrt(df.Δr2)
ntrial_max = df.ntrial.max()
t_final_max = df.t_final.max()
df['ntrial_frac'] = df.ntrial / ntrial_max
df['t_frac'] = df.t / t_final_max

lw, ms = 2, 8
gen_lw, line_lw = 1.2, 1.2

tick_fs, label_fs, legend_fs = 12, 14, 12
xticks = [-50, 0, 50, 100, 150]
yticks = [-100, -50, 0, 50, 100]

umsqpersec = r'µm$^2\,$s$^{-1}$' # ensure commonality between legend and axis label

Dp = args.Dp
dfx = df[(df.Dp == Dp) & (df.traj < 20) & (df.Q > 1e-4) & (df.Q < 1e2)]
freed = np.sqrt(6*Dp*t_final_max)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 10), sharex=True, gridspec_kw={'height_ratios':[1,1]})

sp1 = sb.stripplot(data=dfx, x='Q', y='Δr', log_scale=True, hue='ntrial_frac', palette='autumn_r', 
                   native_scale=True, ax=ax1)

sp1.legend(title='# steps / $10^5$')

sp2 = sb.stripplot(data=dfx, x='Q', y='Δr', log_scale=True, hue='t', palette='winter_r', 
                   native_scale=True, ax=ax2)

sp2.legend(title='$t$ / s')

ax1.axhline(freed, ls='--', color='k', zorder=3)
ax2.set_xlim([1e-4, 1e2])
ax2.set_xlabel('Q')

plt.show()
exit()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 5), sharex=True, gridspec_kw={'height_ratios':[1.5,1]})

ylims = [1, 1e4]
ax1.fill_betweenx(ylims, [1.5e-2]*2, [10.0]*2, color='darkcyan', alpha=0.2)

ax1.set_xlim(1e-3*Q1, 1e-3*Q2)
ax1.set_ylim(*ylims)

ax1.set_yticks([1, 10, 100, 1e3, 1e4],
               labels=['1', '10', r'$10^{2}$', r'$10^{3}$', r'$10^{4}$'])

if data is not None:
    ax1.legend(title=r'$D_p$ / {units}'.format(units=umsqpersec), frameon=False, markerscale=1.3,
               title_fontsize=legend_fs, fontsize=legend_fs, labelspacing=0.3)

ax1.set_ylabel('RMSD / µm', fontsize=label_fs)

ax2.set_xlim(1e-3*Q1, 1e-3*Q2)
ax2.set_ylim(0.1, 1e4)
ax2.set_yticks([0.1, 10, 1e3], labels=['0.1', '10', r'10$^3$'])
ax2.set_xticks([1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100],
               labels=[r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', '0.1', '1', '10', '100'])

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

ax1.annotate('(a)', (2e-4, 3e3), fontsize=label_fs)
ax2.annotate('(b)', (2e-4, 1e3), fontsize=label_fs)

plt.tight_layout()

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
elif not args.verbose:
    plt.show()
