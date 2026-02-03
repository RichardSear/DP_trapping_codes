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
        self.flow_field = self.pipette_flow_field if injector == 'pipette' else None
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
        self.info = self.report()

    def report(self):
        um, umpersec, um2persec, um3persec, none = 'µm', 'µm/s', 'µm²/s', 'µm³/s', ''
        names = ['Q', 'Q', 'Γ', 'k', 'Ds', 'R1', 'α', 'rc',
                 'r*', 'v1', 'P/η', 'Pe', 'λ', 'λ*', 'kλ', 'kλΓ']
        values = [self.Q, 1e-3*self.Q, self.Γ, self.k, self.Ds, self.R1, self.alpha, self.rc,
                  self.rstar, 1e-3*self.v1, self.Pbyη, self.Pe, self.λ, self.λstar, self.kλ, self.kλΓ]
        units = [um3persec, 'pL/s', um2persec, none, um2persec, um, none, um,
                 none, 'mm/s', um2persec, none, um, none, um, um3persec]
        table_d = {'name':pd.Series(names, dtype=str),
                   'value':pd.Series(values, dtype=float),
                   'unit':pd.Series(units, dtype=str)}
        table = pd.DataFrame(table_d).set_index('name')
        table.value = table.value.apply(lambda x: round(x, 3))
        pd.set_option('display.max_columns', None)
        return '\n'.join(table.to_string().split('\n')[2:]) # lose the first two lines

    def pipette_flow_field(self, rvec):
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
