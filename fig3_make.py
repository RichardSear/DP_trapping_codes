#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Python code to plot schematic of pipette flow field
# Warren and Sear 2025/2026

import argparse
import numpy as np
import matplotlib.pyplot as plt
from numpy import log as ln
from numpy import pi as π
from models import Model

parser = argparse.ArgumentParser(description='figure 3b in manuscript')
parser.add_argument('-Q', '--Qrange', default='1e-4,1e2', help='Q range in pL/s, default 1e-4,1e2,80')
parser.add_argument('-e', '--epsilon', default=1e-6, type=float, help='nearness to Qcrit, default 1e-6')
parser.add_argument('-f', '--frac', default=0.7, type=float, help='fraction of Qcrit, default 0.7')
parser.add_argument('-n', '--npt', default=80, type=int, help='number of points, default 80')
parser.add_argument("-v", "--verbose", action="count", default=0)
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

tick_fs, label_fs = 12, 14
gen_lw, line_lw = 1.2, 1.2
xticks = [-50, 0, 50, 100, 150]
yticks = [-100, -50, 0, 50, 100]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 5), sharex=True,
                               gridspec_kw={'height_ratios':[1.5,1]})

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

# the quadratic for the roots is z^2 − (kΓbyD − kλ* − 1)z + kλ* = 0 where z is in units of r*

kλ = k*Q/(4*π*Ds) # this will be an array, as also the things below

b = Γ*k/Ds - kλ/rstar - 1 # the equation is r^2 − br + c = 0
Δ = b**2 - 4*kλ/rstar # discriminant of above

z1 = 0.5*rstar*(b - np.sqrt(Δ))
z2 = 0.5*rstar*(b + np.sqrt(Δ))

# The double root is where the discriminant vanishes at Q = Qc
kλ = k*Qc/(4*π*Ds) # this is now a scalar
b = Γ*k/Ds - kλ/rstar - 1 # the equation is r^2 − br + c = 0 ; 2*r − b = 0
zc = 0.5*rstar*b

ax1.loglog(1e-3*Q, z1, 'm-')
ax1.loglog(1e-3*Q, z2, 'c-')
ax1.loglog(1e-3*Qc, zc, 'ok')

ax1.set_xlim(1e-3*Q1, 1e-3*Q2)
ax1.set_ylim(1, 1e4)

ax1.set_yticks([1, 10, 100, 1e3, 1e4],
               labels=['1', '10', r'$10^{2}$', r'$10^{3}$', r'$10^{4}$'])

ax1.set_ylabel('RMSD / µm', fontsize=label_fs)

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

ax2.loglog(1e-3*Q, ΔS, 'r-')

ax2.set_xlim(1e-3*Q1, 1e-3*Q2)
ax2.set_ylim(0.1, 1e4)
ax2.set_yticks([0.1, 10, 1e3], labels=['0.1', '10', r'10$^3$'])
ax2.set_xticks([1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100],
               labels=[r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', '0.1', '1', '10', '100'])

ax2.set_xlabel(r'Q / pL$\,$s$^{-1}$', fontsize=label_fs)
ax2.set_ylabel(r'$\Delta\mathfrak{S}$ / µm$^2\,$s$^{-1}$', fontsize=label_fs)

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
