"""Clustering mechanics: triangles per node e_i vs k_i (SR, T=1e5), slope = P(hit closes
a triangle at i); C(k) prediction 2[p0+p(k-2)]/(k(k-1)); per-hit triangle probability."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, graphstats as gs
from exp_topology import edges_from_hits
T = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
o = srcells.run(T, seed=11)
h = o['hit']; E = edges_from_hits(h); A = gs.adjacency(E, T); k = gs.degrees(A)
C, tri = gs.local_clustering(A, k)
born = np.arange(T)
# interior nodes (both ends on ligaments) only
inter = (h >= 0).all(1)
# birth triangle: hits dual-adjacent
nb = inter
p0 = (np.asarray(A[h[nb, 0], h[nb, 1]]).ravel() > 0).mean()
# per-hit: for each new node t with hits a,b (both >=0), triangle at a iff a~b
x = k[inter] - 2; y = tri[inter]
m = x > 0
p = np.sum(x[m] * (y[m] - p0)) / np.sum(x[m] ** 2)   # regression through (0,p0)
pp = np.polyfit(x[m & (x > 5)], y[m & (x > 5)], 1)
print(f"T={T} p0(birth triangle)={p0:.4f}  slope p (through p0)={p:.4f}  free fit (k>7): slope={pp[0]:.4f} icpt={pp[1]:.3f}")
print("k  C_meas  C_pred=2[p0+p(k-2)]/(k(k-1))  n")
for lo, hi in [(2, 3), (3, 4), (4, 5), (5, 7), (7, 10), (10, 15), (15, 25), (25, 40), (40, 80), (80, 200), (200, 10**6)]:
    s = inter & (k >= lo) & (k < hi)
    if s.sum() < 3: continue
    kk = k[s]
    pred = np.mean(2 * (p0 + p * (kk - 2)) / (kk * (kk - 1)))
    print(f"[{lo},{hi})  {C[s].mean():.4f}  {pred:.4f}  {s.sum()}")
# total-length law: L(t) = sum of lengths ~ cL sqrt(t); k_i - 2 vs (4 l_i / cL)(sqrt T - sqrt i)
L = np.cumsum(o['length']); tt = np.arange(1, T + 1)
cL = np.polyfit(np.sqrt(tt[T // 10:]), L[T // 10:], 1)
print("L(t) = a sqrt(t) + b:", cL)
l0 = np.polyfit(np.log(tt[1000:]), np.log(o['length'][1000:] + 1e-300), 1)
mean_l = [o['length'][a:2 * a].mean() * np.sqrt(1.5 * a) for a in (1000, 10000, 40000)]
print("<l_t> sqrt(t) at t~1.5e3,1.5e4,6e4:", mean_l)
pred = 4 * o['length'] / cL[0] * (np.sqrt(T) - np.sqrt(tt))
s = inter & (tt > 100) & (tt < 0.5 * T)
print("k-2 vs pred: corr", np.corrcoef(k[s] - 2, pred[s])[0, 1], " mean ratio", (k[s] - 2).sum() / pred[s].sum())
# x = l_i / <l_i> at fixed birth time: <x^2>
xs = []
for a in (100, 1000, 10000):
    seg = o['length'][a:2 * a] * np.sqrt(np.arange(a, 2 * a) + 1)
    xs.append((seg ** 2).mean() / seg.mean() ** 2)
print("<x^2> (cohorts 1e2,1e3,1e4):", xs, " prediction 4<x^2>/K^2 vs BA 6/K^2")
