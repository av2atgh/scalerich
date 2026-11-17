"""Exact Tc(B) of the periodic square wire network (de Gennes-Alexander): at nodes,
sum_j exp(i gamma_ij) psi_j = 4 cos(ka) psi_i, i.e. the Hofstadter-Harper problem. At f = p/q:
psi_{m+1} + psi_{m-1} + 2 cos(2 pi f m + kx) psi_m = eps psi_m, psi_{m+q} = exp(i ky) psi_m.
Lowest k from the largest eps over Bloch phases; k^2 A = (ka)^2, A = a^2.
Same f grid as the torus runs (f = n/4096). -> results/supercond_torusA.json, key periodic_T1023"""
import json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.linalg import eigh_tridiagonal
from math import gcd
T = 1023; D = 4 * (T + 1)
def emax(p, q):
    best = -9
    for kx in (0.0, np.pi / q):
        for ky in (0.0, np.pi):
            m = np.arange(q) - q // 2; d = 2 * np.cos(2 * np.pi * p / q * m + kx)   # cosine maximum mid-chain
            if q == 1: best = max(best, d[0] + 2 * np.cos(ky)); continue
            if q == 2: H = np.diag(d).astype(complex); H[0, 1] = H[1, 0] = 1 + np.exp(1j * ky); best = max(best, np.linalg.eigvalsh(H).max()); continue
            H = sp.diags([d, np.ones(q - 1), np.ones(q - 1)], [0, 1, -1], format='lil').astype(complex)
            H[q - 1, 0] = np.exp(1j * ky); H[0, q - 1] = np.exp(-1j * ky)
            if q <= 400: e = np.linalg.eigvalsh(H.toarray()).max()
            else:   # top state localised in m (Landau level): the corner hop is irrelevant -> open chain
                e = eigh_tridiagonal(d, np.ones(q - 1), eigvals_only=True, select='i', select_range=(q - 1, q - 1))[0]
            best = max(best, float(np.real(e)))
    return min(4.0, best)
ns = np.arange(0, 3 * D + 1, int(round(0.01 * D)))
k2 = []
for n in ns:
    g = gcd(n, D); p, q = n // g, D // g
    k2.append(float(np.arccos(emax(p % q, q) / 4) ** 2))
out = json.load(open('results/supercond_torusA.json'))
out['periodic_T1023'] = dict(rule='periodic', T=T, f=(ns / D).tolist(), k2=[x * (T + 1) for x in k2],
                             pr=[None] * len(ns), Acell=1.0 / (T + 1), nodes=None, method='Harper exact')
json.dump(out, open('results/supercond_torusA.json', 'w'))
for x in (0.25, 0.5, 1.0, 1.5): i = np.argmin(abs(ns / D - x)); print(round(ns[i] / D, 4), round(k2[i], 4))
print('pi/4 squared', (np.pi / 4) ** 2)
