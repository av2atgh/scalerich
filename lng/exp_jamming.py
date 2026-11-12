"""Referee comment 3: lambda0 branch of the jamming boundary for alpha > 1/2.
For each alpha, bisect lambda0 on whether the free area is exhausted before T steps
(finite-width simulator srcells: nucleation in free area, ligaments stop at surfaces).
Compare with the overlap-free estimate 1/sum_t <l_t> t^-alpha and with (l0 zeta(2 alpha))^-1 of Both et al.
Output: results/jamming.txt"""
import sys, os
import numpy as np
from scipy.special import zeta
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
f = open(os.path.join(OUT, 'jamming.txt'), 'w')
def P(s): print(s, flush=True); f.write(s + '\n'); f.flush()
# <l_t> from zero-width runs
T0 = 20000
Ls = np.mean([srcells.run(T0, seed=s)['length'] for s in range(4)], axis=0)
T = 20000
for alpha in (0.6, 0.75, 1.0, 1.5):
    t = np.arange(1, T0 + 1)
    s_free = 1 / np.sum(Ls * t ** -alpha)
    paper = 1 / (0.946 * zeta(2 * alpha)) if 2 * alpha > 1 else float('nan')
    lo, hi = 0.2 * s_free, 3.0 * s_free
    def jams(l0):
        r = [srcells.run(T, alpha=alpha, lam0=l0, seed=s)['T'] < T for s in (1, 2)]
        return np.mean(r)
    for _ in range(9):
        mid = np.sqrt(lo * hi)
        if jams(mid) >= 0.5: hi = mid
        else: lo = mid
    P(f"alpha={alpha}: lambda0* (sim, T={T}) = {np.sqrt(lo*hi):.3f}  overlap-free 1/sum<l_t>t^-a = {s_free:.3f}  "
      f"1/(c_l zeta(a+1/2)) = {1/(1.163*zeta(alpha+0.5)):.3f}  paper 1/(l0 zeta(2a)) = {paper:.3f}")
f.close()
