"""Exact pinned ceiling for axis-parallel divisions in the framed mirror tiling:
C11_pin = E lam (m_h + 1), C22_pin = E lam (m_v + 1), C66_pin = 0 (checked against FE in mech2/
ceiling check). With rho = lam (L_tot + 2) / 2 per unit area of the 2x2 supercell... computed
directly: n_x,max^2 = rho / C11. Counts only -> T up to 1e6. -> results/ceiling.json"""
import json, numpy as np, qdiv
from oned import theta
TAU = qdiv.TAU
def span_by_dir(g, tol=1e-9):
    A, B = g['A'], g['B']; res = {}
    for ax in (0, 1):
        m = np.abs(A[:, ax] - B[:, ax]) < tol
        c = np.round(A[m, ax], 9); lo = np.minimum(A[m, 1 - ax], B[m, 1 - ax]); hi = np.maximum(A[m, 1 - ax], B[m, 1 - ax])
        o = np.lexsort((lo, c)); c, lo, hi = c[o], lo[o], hi[o]; n = 0; i = 0
        while i < len(c):
            j = i; s, e = lo[i], hi[i]
            while j + 1 < len(c) and c[j + 1] == c[i]:
                j += 1
                if lo[j] <= e + tol: e = max(e, hi[j])
                else: n += e - s > 1 - 1e-6; s, e = lo[j], hi[j]
            n += e - s > 1 - 1e-6; i = j + 1
        res[ax] = n
    return res
RULES = {'golden': ('largest', 1 / TAU ** 2, 'perp'), 'golden-side': ('largest', ('side', 1 / TAU ** 2), 'perp'),
         'golden-avoid': ('largest', ('avoid', 1 / TAU ** 2, 'low'), 'perp'), 'plastic51-avoid': ('largest', ('avoid', theta(3, 2)[0] ** -5, 'low'), 'perp'),
         'r0.3': ('largest', 0.3, 'perp'), 'uniform': ('largest', 'uniform', 'perp')}
if __name__ == '__main__':
    out = []
    for nm, rule in RULES.items():
        for T in (1000, 4000, 16000, 64000, 256000, 1000000):
            for seed in ((0, 1, 2) if nm in ('golden-side', 'uniform') else (0,)):
                g = qdiv.run(T, rule, seed=seed); L = np.linalg.norm(g['B'] - g['A'], axis=1).sum()
                sb = span_by_dir(g)
                # supercell 2x2: mass = lam * (4 L + 4 frame lengths of 2 counted once... ) per area 4
                rho_over_lam = (4 * L + 8) / 4.0            # tile ligaments x4 images, frame: 4 lines of length 2
                n2x = rho_over_lam / (sb[1] + 1); n2y = rho_over_lam / (sb[0] + 1)
                out.append(dict(rule=nm, T=T, seed=seed, m=int(sb[0] + sb[1]), m_h=int(sb[1]), m_v=int(sb[0]), L=float(L),
                                nmax_x=float(np.sqrt(n2x)), nmax_y=float(np.sqrt(n2y))))
                print(nm, T, seed, sb, round(L, 1), round(np.sqrt(n2x), 2), round(np.sqrt(n2y), 2), flush=True)
            json.dump(out, open('results/ceiling.json', 'w'), indent=1)
