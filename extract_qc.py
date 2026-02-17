#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Analyse raw Brownian dynamics simulation output and compile to a
# spreadsheet.  For example:

# ./vardp_analyse.py data/vardp.dat.gz -o vardp_results.ods

import gzip
import argparse
import numpy as np
import pandas as pd
from numpy import exp
from numpy import log as ln

def range_str(v, vals): # convert a list of values to a singleton, several values, or a range
    s = ', '.join([str(x) for x in vals]) if len(vals) < 10 else '--'.join([str(f(vals)) for f in [min, max]])
    return v, '  '+s, f'{len(vals):10}'

def getQc(df, Dp, cutoff=0.1): # clunky scheme to extract threshold
    ser = df.loc[Dp]['rat']
    if np.any(ser < cutoff):
        Qc2 = ser[ser < cutoff].dropna().index[0]
        i = ser.index.get_loc(Qc2)
        Qc1 = ser[ser.index < Qc2].index[-1]
        v1, v2 = ser[Qc1], ser[Qc2] # linearly interpolate to cutoff
        lnQc1, lnQc2 = ln(Qc1), ln(Qc2)
        m = (v2-v1)/(lnQc2-lnQc1)
        lnQc = lnQc1 + (cutoff-v1)/m
        Qc = exp(lnQc)
    else:
        Qc, Qc1, Qc2 = np.nan, np.nan, np.nan
    return [Dp, Qc, Qc1, Qc2]

parser = argparse.ArgumentParser(description='compile raw BD data to a spreadsheet')
parser.add_argument('dataset', help='raw input data file, eg *.dat.gz')
parser.add_argument('-d', '--describe', action='store_true', help='print a summary of the columns in the raw data')
parser.add_argument('-c', '--cut-off', default=0.1, type=float, help='cut off for identifying threshold, default 0.1')
#parser.add_argument('-o', '--output', help='output compiled data to a spreadsheet, eg .ods, .xlsx')
args = parser.parse_args()

schema= {'k':float, 'Γ':float, 'Ds':float, 'Dp':float, 'R1':float, 
         'α':float, 'Q':float, 'rc':float, 't_final':float, 
         'ntrial':int, 'nsuccess':int, 't':float, 'Δt_final':float, 'Δr2':float, 
         'traj':int, 'block':int, 'ntraj':int, 'nblock':int, 'code':str}

with gzip.open(args.dataset, 'rt') as f:
    first_line = f.readline()

if len(first_line.split('\t')) < len(schema): # wrangle dataset type, pipette or wall pore
    del schema['α'] # if wall pore then there is no α column

df = pd.read_csv(args.dataset, sep='\t', names=schema.keys(), dtype=schema)
df.sort_values(['Dp', 'Q', 'traj'], inplace=True)

if args.describe:
    df2 = pd.DataFrame([range_str(col, df[col].unique()) for col in df.columns], columns=['column', 'range', 'count'])
    header_row = pd.DataFrame(index=[-1], columns=df2.columns)
    df2 = pd.concat([header_row, df2])
    df2.loc[-1] = df2.columns
    print('Dataset', args.dataset, 'contains', df.shape[0], 'records')
    print('\n'.join(df2.to_string(justify='left', index=False).split('\n')[1:]))
    exit()

df['rat'] = df.Δr2 / (6*df.Dp*df.t_final) # compared to free diffusion
df2 = df[['Dp', 'Q', 'rat']].groupby(['Dp', 'Q']).mean()

Dpvals = df.Dp.unique()
res = np.array([getQc(df2, Dp, cutoff=args.cut_off) for Dp in Dpvals])
df3 = pd.DataFrame(res, columns=['Dp', 'Qc', 'Qc1', 'Qc2'])

print(df3)
