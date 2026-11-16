"""Finite-thickness crossover: m_snap(rho_bar) for the anti-aligned splits (see geo.py).
rho_g = largest density with m_snap = 1. -> results/snap.json"""
import json, numpy as np, qdiv
from mech3 import RULES
from geo import m_snap
out = []
for nm in ('golden-avoid', '(5,1)-avoid', '(3,2)-avoid'):
    for T in (1000, 10000, 100000, 1000000):
        g = qdiv.run(T, RULES[nm]); Lt = np.linalg.norm(g['B'] - g['A'], axis=1).sum()
        rs = np.round(np.arange(0.1, 1.01, 0.05), 2)
        ms = [int(m_snap(g, r / Lt)) for r in rs]
        out.append(dict(rule=nm, T=T, rho=rs.tolist(), m_snap=ms)); print(nm, T, list(zip(rs.tolist(), ms)), flush=True)
json.dump(out, open('results/snap.json', 'w'), indent=1)
