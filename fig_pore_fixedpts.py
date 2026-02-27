#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Sampson flow and drift field for wall pore case
# Warren and Sear 2025/2026

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from models import Model

#plt.rcParams['text.usetex'] = True

parser = argparse.ArgumentParser(description='figure 4 in manuscript')
parser.add_argument('--dpi', default=72, type=int, help='resolution (dpi) for image output, default (for pdf) 72')
parser.add_argument("-v", "--verbose", action="count", default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

lw, ms = 2, 8
tick_fs, label_fs = 12, 14
gen_lw, line_lw = 1.2, 1.2

fig, ax = plt.subplots(figsize=(6, 4), dpi=args.dpi)

# The roots solve (g−3)z²/2 − 3kλz + g/2 = 0 where g = Γ k / D_salt
# and length units are 'b' (aperture radius).  This quadratic solves
# to kλ = [(g−3)z² + g] / (6z).  Alternatively for g < 3, the workable
# soln is z = [√(9kλ² + g(3−g)) − 3kλ] / (3−g).  Exactly at g = 3, one
# has z kλ = 1/2.  For g > 3 critical is where (g−3)z − 3kλ = 0, so
# that the critical value zc = 3kλ / (g−3), where 9kλ² + g(3−g) = 0
# solves to kλ(crit) = 1/3 √(g(g−3)), hence zc = √(g/(g−3)).

def f(z):
    return ((g-3)*z**2 + g)/(6*z)

def h(kλ):
    return (np.sqrt(9*kλ**2 + g*(3-g)) - 3*kλ) / (3-g)

for g in [3.01, 3.1, 3.5]:
    zc = np.sqrt(g/(g-3))
    z1 = np.geomspace(0.1, zc, 80)
    z2 = np.geomspace(zc, 100, 80)
    ax.plot(f(z1), z1, '-', lw=lw, c='tab:red')
    ax.plot(f(z2), z2, '-', lw=lw, c='tab:orange')
    ax.plot(f(zc), zc, 'o', ms=ms, c='tab:brown')

z = np.array([0.1, 100])
plt.plot(1/(2*z), z, '--', lw=lw, c='tab:red')

kλ = np.geomspace(0.01, 1, 80)
for g in [2.5, 2.9, 2.99]:
    plt.plot(kλ, h(kλ), '-', lw=lw, c='tab:red')

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(1e-2, 10)
xticks = [1e-2, 0.1, 1, 10]
xlabels = ['0.01', '0.1', '1', '10']
ax.set_xticks(xticks, labels=xlabels)
ax.set_xlabel(r'$kQ$ / $4\pi b D_{\mathrm{s}}$', fontsize=label_fs)

ax.set_ylim(1, 1e2)
ax.set_yticks([1, 10, 100], labels=['1', '10', '$10^{2}$'])
ax.set_ylabel('$z$ / $b$', fontsize=label_fs)

ax.minorticks_off()
ax.tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)
for spine in ax.spines:
    ax.spines[spine].set_linewidth(gen_lw)

ax.text(1.2, 6, r'$$\frac{\Gamma k}{D_{\mathrm{s}}}$$', usetex=True, fontsize=16)#1.2*label_fs)
ax.text(2.0, 6, '= 3.5', fontsize=label_fs)
ax.text(0.4, 14, '3.1', fontsize=label_fs)
ax.text(0.10, 30, '3.01', fontsize=label_fs)
ax.text(0.015, 42, '3 (exact)', fontsize=label_fs)
ax.text(0.015, 8, '2.99', fontsize=label_fs)
ax.text(0.015*1.05, 3.2, '2.9', fontsize=label_fs)
ax.text(0.015*1.05, 1.4, '2.5', fontsize=label_fs)

plt.tight_layout()

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
elif not args.verbose:
    plt.show()
