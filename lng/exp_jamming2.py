"""Second-round referee item 2: finite-T jamming thresholds and the 2D depletion recursion.
Criterion: jammed = free area exhausted (total free-cell area <= 1e-14, the simulator stops) before step T.
Depletion recursion (2D analogue of the 3D one): l_t = c(t) (F/t)^{1/2} =>  dF/dt = -lambda0 c(t) F^{1/2} t^{-alpha-1/2}
=> sqrt(F_inf) = 1 - (lambda0/2) sum_t c(t) t^{-alpha-1/2}, jam iff lambda0 >= 2 / sum_t c(t) t^{-alpha-1/2}.
c(t) = <l_t> sqrt(t) from zero-thickness runs (F = 1); tail beyond the runs with c = 1.163.
Output: results/jamming2.txt"""
import sys, os
import numpy as np
from scipy.special import zeta
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
f = open(os.path.join(OUT, 'jamming2.txt'), 'w')
def P(s): print(s, flush=True); f.write(s + '\n'); f.flush()
T0 = 100000
Ls = np.mean([srcells.run(T0, seed=s)['length'] for s in range(4)], axis=0)
t = np.arange(1, T0 + 1); c = Ls * np.sqrt(t)
for alpha in (0.6, 0.75, 1.0, 1.5):
    S = np.sum(c * t ** (-alpha - 0.5)) + 1.163 * (zeta(alpha + 0.5) - np.sum(t ** (-alpha - 0.5)))
    P(f"alpha={alpha}: depletion recursion lambda0* = {2/S:.3f} (c(t) measured, T=inf); with c=1.163 throughout: {2/(1.163*zeta(alpha+0.5)):.3f}; Both et al.: {1/(0.946*zeta(2*alpha)):.3f}")
for alpha in (0.6, 0.75):
    for T in (20000, 100000):
        lo, hi = 0.1, 0.6
        for _ in range(9):
            mid = np.sqrt(lo * hi)
            j = np.mean([srcells.run(T, alpha=alpha, lam0=mid, seed=s)['T'] < T for s in (1, 2)])
            if j >= 0.5: hi = mid
            else: lo = mid
        P(f"alpha={alpha} T={T}: simulated lambda0* = {np.sqrt(lo*hi):.3f}")
f.close()
