"""Long-wavelength index of the split-law family: framed periodic tiling (beyond/periodic.py,
the model of the regular article), alpha = 0, uniform thickness at areal density rho_t.
Coincident junctions (X-crossings) are merged into one rigid joint; zero-length spans dropped.
Also: maximal straight lines (collinear ligaments joined through crossings) -> number m_span
of lines spanning the square (frame excluded) and the length-weighted mean line length.
-> results/mech_T1000.json"""
import sys, json, time, numpy as np
sys.path[:0] = ['../beyond', '../rg']
import qdiv, divisions, periodic, waves, srtree, frame
from diffraction import RULES

def build(g, lam, frame_lam, bend=1.0):
    nodes, seg, sl, bn = srtree.network(g, frame=True)
    lam_ext = np.r_[lam, np.full(4, frame_lam)]
    X, S, Ls = periodic.mirror_supercell(nodes, seg, lam_ext[sl])
    keep = S[:, 0] != S[:, 1]
    S, Ls = S[keep], Ls[keep]
    key = np.sort(S, axis=1); _, first = np.unique(key, axis=0, return_index=True)
    k2 = np.zeros(len(S), bool); k2[first] = True
    on2 = lambda ax: (np.abs(X[S[:, 0], ax] - 2) < 1e-9) & (np.abs(X[S[:, 1], ax] - 2) < 1e-9)
    k2 &= ~on2(0) & ~on2(1); S, Ls = S[k2], Ls[k2]
    K = frame.assemble(X, S, Ls, bend=bend)[0]
    master, shift = periodic.pairing(X); M, L = periodic.mass_vector(X, S, Ls)
    return dict(X=X, S=S, lam=Ls, K=K, M=M, master=master, shift=shift, rho=float((Ls * L).sum() / 4))

def lines(g, tol=1e-9):
    """maximal straight lines of axis-parallel ligaments; returns lengths (frame excluded)."""
    A, B = g['A'], g['B']; out = []
    for ax in (0, 1):                     # ax = 0: vertical ligaments (constant x)
        m = np.abs(A[:, ax] - B[:, ax]) < tol
        c = np.round(A[m, ax], 9); lo = np.minimum(A[m, 1 - ax], B[m, 1 - ax]); hi = np.maximum(A[m, 1 - ax], B[m, 1 - ax])
        for x in np.unique(c):
            k = c == x; iv = sorted(zip(lo[k], hi[k])); s, e = iv[0]
            for a, b in iv[1:]:
                if a <= e + tol: e = max(e, b)
                else: out.append(e - s); s, e = a, b
            out.append(e - s)
    return np.array(out)

if __name__ == '__main__':
    T = 1000; out = []; t0 = time.time()
    names = ['regular', 'golden', 'plastic', 'r0.3', 'golden-side', 'largest-uniform', 'balanced', 'SR']
    for name in names:
        for seed in (0, 1, 2, 3):
            g = qdiv.run(T, RULES[name], seed=seed) if name in RULES else divisions.run(T, 100 + seed, name)
            L = np.linalg.norm(g['B'] - g['A'], axis=1)
            ln = lines(g) if name in RULES else L
            for rho_t in (0.1, 0.2, 0.4):
                lam = np.full(T, rho_t / L.sum())
                for bend in (1.0, 0.0):
                    sc = build(g, lam, lam[0], bend=bend)
                    try:
                        C = periodic.periodic_C(sc)
                    except Exception as e:
                        print('fail', name, seed, rho_t, bend, e); continue
                    th, n = waves.christoffel_n(C, sc['rho'])
                    o = dict(rule=name, seed=seed, rho_t=rho_t, bend=bend, rho=sc['rho'], n_mean=float(n.mean()),
                             n_x=float(n[0]), eta=float(n.min() / n.max()),
                             m_span=int((ln > 1 - 1e-6).sum()), line_w=float((ln ** 2).sum() / ln.sum()))
                    out.append(o)
                    print(name, seed, rho_t, bend, 'n=%.2f eta=%.2f m=%d lw=%.3f' % (o['n_mean'], o['eta'], o['m_span'], o['line_w']), round(time.time() - t0), flush=True)
            if name in ('regular', 'golden', 'plastic', 'r0.3'): break   # deterministic
    json.dump(out, open('results/mech_T1000.json', 'w'), indent=1)
