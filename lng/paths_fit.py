"""Fits of <mu>, <ln mu>, <d> vs ln T from results/paths_scaling.csv -> results/paths_fit.txt"""
import csv, collections, math, os
import numpy as np
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
d = collections.defaultdict(list)
for r in csv.DictReader(l for l in open(os.path.join(OUT, 'paths_scaling.csv')) if not l.startswith('#')):
    d[r['model']].append(r)
lines = []
for m, rs in d.items():
    Ts = sorted(set(int(r['T']) for r in rs))
    av = lambda key, T: np.mean([float(r[key]) for r in rs if int(r['T']) == T])
    x = np.log(Ts)
    for key in ('meand', 'mu', 'lnmu'):
        y = np.array([av(key, T) for T in Ts])
        p1 = np.polyfit(x, y, 1); p2 = np.polyfit(x, y, 2)
        lines.append(f"{m:7s} {key:6s}: " + ' '.join(f"{v:.2f}" for v in y) +
                     f" | linear in lnT slope={p1[0]:.3f} | quad coef={p2[0]:.3f}")
    c2 = np.array([av('c2', T) for T in Ts]); mu = np.array([av('mu', T) for T in Ts])
    lines.append(f"{m:7s} <mu>/(c2/e): " + ' '.join(f"{v:.2f}" for v in mu / (c2 / math.e)))
lines.insert(0, 'T = ' + ' '.join(str(T) for T in Ts))
open(os.path.join(OUT, 'paths_fit.txt'), 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
