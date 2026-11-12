"""Do high-degree ligaments border large voids (cells)?
At the end of a run (cells kept), for each ligament i: the cells it bounds,
their areas (mean, largest, total) against degree k_i. Areas in units of the
mean cell area. Usage: python exp_void_degree.py T seed. Output appended to
results/void_degree.txt"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells
T = int(sys.argv[1]); seed = int(sys.argv[2])
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'void_degree.txt')
f = open(OUT, 'a')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); f.write(s + '\n'); f.flush()
o = srcells.run(T, seed=seed, keep_cells=True)
T = o['T']; h = o['hit']
tg = h.ravel(); tg = tg[tg >= 0]
k = np.bincount(tg, minlength=T) + 2
A = np.array(o['cellA']); live = A > 0
Abar = A[live].mean()
nb = [[] for _ in range(T)]
for Lc, Ac in zip(o['cellL'], A):
    if Ac <= 0: continue
    for lab in set(Lc):
        if lab >= 0: nb[lab].append(Ac / Abar)
mean_a = np.array([np.mean(x) if x else np.nan for x in nb])
max_a = np.array([np.max(x) if x else np.nan for x in nb])
tot_a = np.array([np.sum(x) if x else np.nan for x in nb])
ncell = np.array([len(x) for x in nb])
# interior ligaments only (both ends on ligaments) to avoid wall effects
inter = (h[:, 0] >= 0) & (h[:, 1] >= 0) & (ncell > 0)
P(f"== T={T} seed={seed}  cells={live.sum()}  (areas in units of mean cell area)")
P(f"{'k bin':>10} {'n':>7} {'cells/k':>8} {'mean adj':>9} {'median adj':>10} {'largest adj':>11} {'total adj':>10}")
for lo, hi in [(2,3),(3,4),(4,6),(6,10),(10,20),(20,40),(40,80),(80,160),(160,10**9)]:
    s = inter & (k >= lo) & (k < hi)
    if s.sum() < 10: continue
    P(f"{f'[{lo},{hi})':>10} {s.sum():7d} {np.mean(ncell[s]/k[s]):8.2f} {np.mean(mean_a[s]):9.2f} "
      f"{np.median(mean_a[s]):10.2f} {np.mean(max_a[s]):11.2f} {np.mean(tot_a[s]):10.1f}")
P(f"Spearman(k, mean adjacent area) = "
  f"{np.corrcoef(np.argsort(np.argsort(k[inter])), np.argsort(np.argsort(mean_a[inter])))[0,1]:.3f}; "
  f"Spearman(k, largest adjacent area) = "
  f"{np.corrcoef(np.argsort(np.argsort(k[inter])), np.argsort(np.argsort(max_a[inter])))[0,1]:.3f}")
# area-weighted: share of all void area bordered, by degree class (what a random nucleation point sees)
