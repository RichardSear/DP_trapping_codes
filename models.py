#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Python code to implement pipette and hole models
# Warren and Sear 2025/2026

import numpy as np
import pandas as pd
from numpy import pi as π

class Model:

    def __init__(self, injector): # initialise base parameters, units are um and seconds
        self.Q = 1.0 # injection rate
        self.Γ = 150 # drift coefficient, here for DNA in LiCl
        self.k = 200 # in um^3/sec ; note conversion 1 pL/sec = 10^3 um^3/sec
        self.Ds = 1610 # for NaCl
        self.R1 = 1.0 # pipette radius, or pore radius
        self.alpha = 0.3 # from Secchi at al
        self.rc = 1.0 # default cut off
        self.drift = self.pipette_drift if injector == 'pipette' else None
        self.refresh()

    def update(self, Q=1e-3, Γ=150, k=200, Ds=1610, R1=1.0, α=0.3, rc=1.0): # add more here as required
        self.Q = 1e3*Q # convert to um^3/sec
        self.Γ = Γ
        self.k = k
        self.Ds = Ds
        self.R1 = R1
        self.alpha = α
        self.rc = rc
        self.refresh()

    def refresh(self): # calculate derived quantities
        self.rstar = π*self.R1/self.alpha # where stokeslet and radial outflow match
        self.v1 = self.Q / (π*self.R1**2) # flow speed (definition)
        self.Pbyη = self.alpha*self.R1*self.v1 # from Secchi et al, should also = Q / r*
        self.Pe = self.Pbyη / (4*π*self.Ds) # definition from Secchi et al
        self.λ = self.Q / (4*π*self.Ds) # definition
        self.λstar = self.λ / self.rstar # should be the same as Pe (salt)
        self.kλ = self.k * self.λ
        self.kλΓ = self.k * self.λ * self.Γ
        self.Qcrit = 4*π*self.Ds*self.rstar*(np.sqrt(self.Γ/self.Ds)-np.sqrt(1/self.k))**2 # critical upper bound on Q
        # the quadratic for the roots is z^2 − (kΓbyD − kλ* − 1)z + kλ* = 0 where z is in units of r*
        b = self.k*self.Γ/self.Ds - self.k*self.λstar - 1 # the equation is r^2 − br + c = 0
        Δ = b**2 - 4*self.k*self.λstar # discriminant of above
        if self.Q > 0 and Δ > 0 and b > 0: # condition for roots to exist
            self.fixed_points = np.array((0.5*self.rstar*(b-np.sqrt(Δ)), 0.5*self.rstar*(b+np.sqrt(Δ))))
        else:
            self.fixed_points = None
        self.info = self.report()

    def report(self):
        um, umpersec, um2persec, um3persec, none = 'µm', 'µm/s', 'µm²/s', 'µm³/s', ''
        names = ['Q', 'Qcrit', 'Γ', 'k', 'Ds', 'R1', 'α', 'rc',
                 'r*', 'v1', 'P/η', 'Pe', 'λ', 'λ*', 'kλ', 'kλ*']
        values = [1e-3*self.Q, 1e-3*self.Qcrit, self.Γ, self.k, self.Ds, self.R1, self.alpha, self.rc,
                  self.rstar, 1e-3*self.v1, self.Pbyη, self.Pe, self.λ, self.λstar, self.kλ, self.k*self.λstar]
        units = ['pL/s', 'pL/s', um2persec, none, um2persec, um, none, um,
                 none, 'mm/s', um2persec, none, um, none, um, none]
        if self.fixed_points is not None:
            names.extend(['z1', 'z2'])
            values.extend(list(self.fixed_points))
            units.extend([um, um])
        table_d = {'name':pd.Series(names, dtype=str),
                   'value':pd.Series(values, dtype=float),
                   'unit':pd.Series(units, dtype=str)}
        table = pd.DataFrame(table_d).set_index('name')
        table.value = table.value.apply(lambda x: round(x, 3))
        pd.set_option('display.max_columns', None)
        return '\n'.join(table.to_string().split('\n')[2:]) # lose the first two lines

    def pipette_drift(self, rvec):
        x, y, z = rvec[:] # z is normal distance = r cosθ
        ρ = np.sqrt(x**2 + y**2) # in-plane distance = r sinθ
        r = np.sqrt(x**2 + y**2 + z**2) # radial distance from origin
        cosθ, sinθ = z/r, ρ/r # polar angle, as cos and sin
        sinθ_cosφ, sinθ_sinφ = x/r, y/r # avoids dividing by ρ
        ur = self.Q/(4*π*r**2) + self.Pbyη*cosθ/(4*π*r) - self.kλΓ/(r*(r+self.kλ)) # radial drift velocity
        ux = ur*sinθ_cosφ - self.Pbyη*sinθ_cosφ*cosθ/(8*π*r) # radial components
        uy = ur*sinθ_sinφ - self.Pbyη*sinθ_sinφ*cosθ/(8*π*r) # avoiding dividing by ρ
        uz = ur*cosθ + self.Pbyη*sinθ**2/(8*π*r) # normal z-component of velocity
        return np.zeros_like(rvec) if r < self.rc else np.array((ux, uy, uz))
