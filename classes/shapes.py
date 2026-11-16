"""Finite local complexity test: number of distinct leaf shapes (w, h up to the 90-degree
rotation, in units 1/sqrt(T+1)) and of distinct leaf areas vs T. -> results/shapes.json"""
import json, numpy as np, qdiv
from diffraction import RULES
out = []
for nm in ('regular', 'golden', 'plastic', 'r0.3', 'golden-side'):
    for T in (300, 1000, 3000, 10000, 30000, 100000, 300000):
        g = qdiv.run(T, RULES[nm], seed=0)
        wh = np.sort(g['leaf_wh'], axis=1) * np.sqrt(T + 1)
        ns = len(set(map(tuple, np.round(wh, 7)))); na = len(set(np.round(g['leaf_area'] * (T + 1), 7)))
        out.append(dict(rule=nm, T=T, shapes=ns, areas=na)); print(nm, T, ns, na, flush=True)
json.dump(out, open('results/shapes.json', 'w'), indent=1)
