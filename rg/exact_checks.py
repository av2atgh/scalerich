"""Checks of the exact consequences of the BST mapping, and the density /
jamming bookkeeping (Goal 1, 2c).

1. Cell areas at step T are uniform spacings (Dirichlet(1,...,1) on T+1 cells):
   sum A^2 = 2/(T+2) exactly in mean; P(TA > x) -> e^-x.
2. Height of the splitting tree ~ 4.311 ln T (Devroye 1986), saturation level
   ~ 0.373 ln T; depth of a uniform leaf 2 ln T.
3. <l_t>: l_1 and c_l = lim l_t sqrt t; zero-width density bookkeeping
   rho_inf(lam0, alpha) = lam0 S(alpha), S = sum_t <l_t> t^-alpha (alpha > 1/2),
   so the zero-width jamming threshold is lam0* = 1/S(alpha); compare with the
   paper's (l0 zeta(2 alpha))^-1 and with (c_l zeta(alpha + 1/2))^-1.
Writes results/exact_checks.json
"""
import json
import numpy as np
from scipy.special import zeta
import srtree

out = {}
T = 20000; NS = 20
s2 = []; tail = []; H = []; sat = []; Ls = np.zeros(T)
for s in range(NS):
    g = srtree.run(T, 1000 + s)
    a = g['cell_area'][g['leaf']]
    s2.append(np.sum(a ** 2)); tail.append([np.mean(T * a > x) for x in (0.5, 1, 2, 4)])
    dep = g['cell_depth'][g['leaf']]
    H.append(dep.max())
    # saturation level: deepest level with all 2^d nodes present = (min leaf depth)
    sat.append(dep.min())
    Ls += np.linalg.norm(g['B'] - g['A'], axis=1)
Ls /= NS
out['sumA2'] = dict(measured=float(np.mean(s2)), exact=2 / (T + 2))
out['P(TA>x)'] = dict(x=[0.5, 1, 2, 4], measured=np.mean(tail, 0).round(4).tolist(),
                      exp=np.exp(-np.array([0.5, 1, 2, 4])).round(4).tolist())
out['height'] = dict(mean=float(np.mean(H)), per_lnT=float(np.mean(H) / np.log(T)), devroye=4.311,
                     saturation=float(np.mean(sat)), sat_per_lnT=float(np.mean(sat) / np.log(T)), devroye_sat=0.373)
t = np.arange(1, T + 1)
cl = float(np.mean(Ls[T // 2:] * np.sqrt(t[T // 2:])))
out['lengths'] = dict(l1=float(Ls[0]), c_l=cl)
S = {}
for a in (0.6, 0.75, 1.0, 1.5):
    Sa = float(np.sum(Ls * t ** -a) + cl * (zeta(a + 0.5) - np.sum(t ** (-a - 0.5))))
    S[a] = dict(S=Sa, lam0_star=1 / Sa, paper=1 / (0.946 * zeta(2 * a)), cl_zeta=1 / (cl * zeta(a + 0.5)))
out['jamming_lam0'] = S
print(json.dumps(out, indent=1))
json.dump(out, open('results/exact_checks.json', 'w'), indent=1)
