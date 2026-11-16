"""Axis Bragg weight (2D junction set, window q/2pi in [0.25,2] sqrt T) vs T for each side rule.
-> results/bragg2.json"""
import json, numpy as np, qdiv
from pisot import axis_top
from diffraction import vertices
from oned import theta
TAU = qdiv.TAU; RHO = theta(3, 2)[0]
CASES = {'golden': TAU ** -2, 'golden-side': ('side', TAU ** -2), 'golden-avoid': ('avoid', TAU ** -2, 'low'),
         '(5,1)': RHO ** -5, '(5,1)-avoid': ('avoid', RHO ** -5, 'low'), '(3,2)': RHO ** -3, '(3,2)-avoid': ('avoid', RHO ** -3, 'low'),
         'r0.3': 0.3, 'uniform': 'uniform'}
out = []
for nm, sp in CASES.items():
    for T in (4095, 16383, 65535, 262143, 1048575):
        I, k = axis_top(vertices(qdiv.run(T, ('largest', sp, 'perp'), seed=0)), T)
        out.append(dict(rule=nm, T=T, Imax=float(I), k=float(k))); print(nm, T, round(I, 4), k, flush=True)
        json.dump(out, open('results/bragg2.json', 'w'), indent=1)
