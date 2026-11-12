"""Tail amplitude K^2 P(k>=K) at T=1e6 over 5 realisations, vs 4<x^2> (length/fitness) and 6K/(K+1) (BA).
Output: results/tail.txt"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
T = 1000000; Ks = (16, 32, 64, 128); A = []; X = []
for s in (11, 12, 13, 14, 15):
    o = srcells.run(T, seed=s); h = o['hit']
    k = np.bincount(h[h >= 0], minlength=T) + 2
    i = np.arange(T); x = o['length'] * np.sqrt(i + 1) / 1.163
    X.append(np.mean(x[i > 1000] ** 2)); A.append([K * K * (k >= K).mean() for K in Ks])
    print(s, X[-1], A[-1], flush=True)
A = np.array(A); X = np.array(X)
with open(os.path.join(OUT, 'tail.txt'), 'w') as f:
    f.write(f"T=1e6, 5 realisations. <x^2> = {X.mean():.3f} +- {X.std(ddof=1)/np.sqrt(5):.3f} -> 4<x^2> = {4*X.mean():.2f}\n")
    for j, K in enumerate(Ks):
        f.write(f"K={K}: K^2P(k>=K) = {A[:,j].mean():.3f} +- {A[:,j].std(ddof=1)/np.sqrt(5):.3f}   BA 6K/(K+1) = {6*K/(K+1):.3f}\n")
    pooled = A[:, 1:].mean(1)
    f.write(f"pooled K=32..128 per realisation mean {pooled.mean():.3f} +- {pooled.std(ddof=1)/np.sqrt(5):.3f}; BA mean over same K {np.mean([6*K/(K+1) for K in Ks[1:]]):.3f}\n")
print(open(os.path.join(OUT, 'tail.txt')).read())
