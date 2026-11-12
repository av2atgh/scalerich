"""Finite-K prediction of the fitness (length x age) model: k_i = 2 + Poisson(4 l_i (sqrt T - sqrt i)/c_L),
with the actual birth lengths l_i, vs the measured k_i (same realisations). Output: results/tail2.txt"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
T = 1000000; Ks = (16, 32, 64, 128); rng = np.random.default_rng(5)
M, F = [], []
for s in (11, 12, 13):
    o = srcells.run(T, seed=s); h = o['hit']; L = o['length']
    k = np.bincount(h[h >= 0], minlength=T) + 2
    cL = L.sum() / np.sqrt(T)            # measured total-length constant of this realisation
    i = np.arange(T)
    mu = 4 * L / cL * (np.sqrt(T) - np.sqrt(i + 1))
    kf = 2 + rng.poisson(mu)
    M.append([K * K * (k >= K).mean() for K in Ks]); F.append([K * K * (kf >= K).mean() for K in Ks])
    print(s, cL, M[-1], F[-1], flush=True)
M, F = np.array(M), np.array(F)
with open(os.path.join(OUT, 'tail2.txt'), 'w') as f:
    for j, K in enumerate(Ks):
        f.write(f"K={K}: measured {M[:,j].mean():.3f}+-{M[:,j].std(ddof=1)/np.sqrt(3):.3f}  fitness-model (Poisson, actual lengths) {F[:,j].mean():.3f}+-{F[:,j].std(ddof=1)/np.sqrt(3):.3f}  BA {6*K/(K+1):.3f}\n")
print(open(os.path.join(OUT, 'tail2.txt')).read())
