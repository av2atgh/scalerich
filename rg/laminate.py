"""Hierarchical (diamond-like) surrogate: the SR splitting tree read as a random
sequential laminate (Goal 2a).

Bond moving in continuum form: the chord of cell c (direction n, thickness lam,
length l, cell area A) is smeared into its cell as a rank-one medium
(lam l / A) n n^T, and the two children are taken as a laminate with
interface tangent n (scale separation between generations, as on a
hierarchical lattice).  The recursion is exact on that surrogate:

  a_i = m.S_i.m, b_i = n.S_i.m, c_i = n.S_i.n   (m = normal)
  a* = <1/a>^-1, b* = a* <b/a>, c* = <c - b^2/a> + a* <b/a>^2
  S_c = laminate(S_0, S_1; f_0 = r, f_1 = 1 - r) + (lam l / A) n n^T,
  leaves S = 0.

Root tensor -> Phi_lam = tr S_root / tr S_aff, S_aff = sum lam l n n^T, to be
compared with the network values Phi (conductance.py / flow.py).
Writes results/laminate.json
"""
import json
import numpy as np
import srtree


def laminate(S0, S1, f0, n):
    m = np.array([-n[1], n[0]])
    f = (f0, 1 - f0)
    a = [m @ S @ m for S in (S0, S1)]; b = [n @ S @ m for S in (S0, S1)]; c = [n @ S @ n for S in (S0, S1)]
    if min(a) <= 1e-300:
        astar = 0.0; bstar = 0.0
        cstar = sum(fi * (ci - (bi * bi / ai if ai > 1e-300 else 0.0)) for fi, ai, bi, ci in zip(f, a, b, c))
    else:
        astar = 1 / sum(fi / ai for fi, ai in zip(f, a))
        ba = sum(fi * bi / ai for fi, ai, bi in zip(f, a, b))
        bstar = astar * ba
        cstar = sum(fi * (ci - bi * bi / ai) for fi, ai, bi, ci in zip(f, a, b, c)) + astar * ba * ba
    N = np.c_[n, m]                      # columns n, m
    Sl = np.array([[cstar, bstar], [bstar, astar]])
    return N @ Sl @ N.T


def run(g, lam):
    T = g['T']
    S = {}
    for t in range(T - 1, -1, -1):
        c = g['parent'][t]; k0, k1 = g['kids'][t]
        d = g['B'][t] - g['A'][t]; l = np.linalg.norm(d); n = d / l
        S0 = S.pop(k0, np.zeros((2, 2))); S1 = S.pop(k1, np.zeros((2, 2)))
        S[c] = laminate(S0, S1, g['r'][t], n) + lam[t] * l / g['area'][t] * np.outer(n, n)
    d = g['B'] - g['A']; l = np.linalg.norm(d, axis=1)
    Saff = np.einsum('t,ti,tj->ij', lam / l, d, d)
    return S[0], Saff


if __name__ == '__main__':
    out = {}
    for T in (3000, 30000):
        for a in (0.0, 0.2, 0.4, 0.5, 0.6, 1.0):
            ph = []
            for s in range(4):
                g = srtree.run(T, 500 + s)
                Sr, Sa = run(g, srtree.thickness(T, 1.0, a))
                ph.append(np.trace(Sr) / np.trace(Sa))
            out[f'{T}_{a}'] = dict(Phi_lam=float(np.mean(ph)), sd=float(np.std(ph)))
            print(T, a, out[f'{T}_{a}'], flush=True)
    json.dump(out, open('results/laminate.json', 'w'), indent=1)
