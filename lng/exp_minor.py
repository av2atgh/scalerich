"""Referee minor comments: (9) MLE bias on an exact BA m=2 sample and uncertainties;
(12) direct per-hit triangle probability p(k); boundary hits.
Output: results/minor.txt"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, graphstats as gs
from exp_topology import edges_from_hits
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
f = open(os.path.join(OUT, 'minor.txt'), 'w')
def P(s): print(s, flush=True); f.write(s + '\n'); f.flush()
rng = np.random.default_rng(0)
# exact p_k sample: P(k>=K) = 6/[K(K+1)] -> inverse transform on continuous approx then floor check
K = np.arange(2, 10**6); cdf = 1 - 6 / ((K + 1) * (K + 2))
for n in (10**5, 3 * 10**5):
    u = rng.random(n); k = K[np.searchsorted(cdf, u)]
    P(f"exact BA m=2 sample n={n}: " + "  ".join(
        f"MLE(k>={km})={gs.gamma_mle(k, km)[0]:.3f}+-{(gs.gamma_mle(k, km)[0]-1)/np.sqrt(gs.gamma_mle(k, km)[1]):.3f}"
        for km in (6, 10, 20, 40)))
T = 100000
for s in (1, 2):
    o = srcells.run(T, seed=s); h = o['hit']
    nb = (h < 0).sum(); P(f"SR seed {s}: boundary ends {nb} of {2*T} ({nb/(2*T):.4f}); square boundary is not a node")
    k = np.bincount(h[h >= 0], minlength=T) + 2
    P("   MLE: " + "  ".join(f"k>={km}: {gs.gamma_mle(k, km)[0]:.3f}+-{(gs.gamma_mle(k, km)[0]-1)/np.sqrt(gs.gamma_mle(k, km)[1]):.3f}" for km in (6, 10, 20, 40)))
    # per-hit triangle probability: at time t newcomer hits a and b; the hit on a creates a
    # triangle at a iff b is already a neighbour of a.
    nbr = [set() for _ in range(T)]; deg = np.full(T, 0)
    rows = []
    for t in range(T):
        a, b = h[t]
        if a >= 0 and b >= 0:
            adj = b in nbr[a]
            rows.append((deg[a], adj)); rows.append((deg[b], adj))
        for x in (a, b):
            if x >= 0: nbr[x].add(t); nbr[t].add(x); deg[x] += 1
        deg[t] = sum(1 for x in (a, b) if x >= 0)
    r = np.array(rows, float)
    P(f"   per-hit triangle prob overall {r[:,1].mean():.3f}")
    for lo, hi in ((2, 3), (3, 5), (5, 10), (10, 20), (20, 50), (50, 10**6)):
        m = (r[:, 0] >= lo) & (r[:, 0] < hi)
        P(f"     degree at hit in [{lo},{hi}): n={int(m.sum())} p={r[m,1].mean():.3f}")
f.close()
