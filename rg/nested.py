"""Exact chygraph contraction of the SR conductance problem along the splitting
tree (Goal 1/4).

Complex = cell c.  Its interior = segments of the chords created inside c
(its descendants).  Its atoms = endpoints of those chords lying on c's
boundary.  Message = Dirichlet-to-Neumann (Schur complement) matrix of the
interior network onto the atoms.  Parent = union of the two children's
messages + the splitting chord's segments; eliminate the junctions on the
chord (endpoints of inner stems hosted by it).  At the root apply the
electrodes.  This is exact; it is nested dissection with the chords as
separators.  Measures the width w(c) = #atoms against n(c) = #chords inside c.
Writes results/nested.json
"""
import json, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import srtree


def contract(g, lam):
    T = g['T']
    nodes, seg, sl, bn = srtree.network(g)
    Ls = np.linalg.norm(nodes[seg[:, 0]] - nodes[seg[:, 1]], axis=1)
    gs = lam[sl] / Ls
    segs_of = [[] for _ in range(T)]
    for k, t in enumerate(sl):
        segs_of[t].append(k)
    host = np.empty(2 * T, int); host[0::2] = g['hostA']; host[1::2] = g['hostB']
    ncell = len(g['cell_area'])
    msg = [None] * ncell            # (atoms array, M)
    nin = np.zeros(ncell, int)
    widths = []
    flops = 0
    for t in range(T - 1, -1, -1):
        c = g['parent'][t]; k0, k1 = g['kids'][t]
        parts = [m for m in (msg[k0], msg[k1]) if m is not None]
        nin[c] = 1 + nin[k0] + nin[k1]
        ids = set()
        for at, _ in parts:
            ids.update(at.tolist())
        for k in segs_of[t]:
            ids.update(seg[k].tolist())
        ids = np.array(sorted(ids)); pos = {v: i for i, v in enumerate(ids)}
        n = len(ids); M = np.zeros((n, n))
        for at, Mk in parts:
            ix = np.array([pos[v] for v in at])
            M[np.ix_(ix, ix)] += Mk
        for k in segs_of[t]:
            a, b = pos[seg[k, 0]], pos[seg[k, 1]]; w = gs[k]
            M[a, a] += w; M[b, b] += w; M[a, b] -= w; M[b, a] -= w
        elim = np.array([host[v] == t for v in ids])   # junctions on chord t
        keep = ~elim
        if elim.any():
            Aii = M[np.ix_(elim, elim)]; Aib = M[np.ix_(elim, keep)]
            M = M[np.ix_(keep, keep)] - Aib.T @ np.linalg.solve(Aii, Aib)
            flops += elim.sum() ** 3 + elim.sum() ** 2 * keep.sum()
        msg[c] = (ids[keep], M)
        widths.append((int(nin[c]), int(keep.sum()), float(g['cell_area'][c]), int(g['cell_depth'][c])))
        msg[k0] = msg[k1] = None
    at, M = msg[0]
    hs = np.array([host[v] for v in at])
    V = np.full(len(at), np.nan); V[hs == -4] = 0.0; V[hs == -2] = 1.0
    fx = ~np.isnan(V); fr = ~fx
    if fr.any():
        V[fr] = np.linalg.solve(M[np.ix_(fr, fr)], -M[np.ix_(fr, fx)] @ V[fx])
    P = float(V @ M @ V)
    return P, np.array(widths), flops


def direct(g, lam):
    from conductance import solve
    nodes, seg, sl, bn = srtree.network(g)
    Ls = np.linalg.norm(nodes[seg[:, 0]] - nodes[seg[:, 1]], axis=1)
    V, P = solve(nodes, seg, lam[sl] / Ls, bn[-4], bn[-2])
    return float(P.sum())


if __name__ == '__main__':
    out = {}
    for T in (300, 1000, 3000, 10000):
        t0 = time.time()
        g = srtree.run(T, 7)
        lam = srtree.thickness(T, 1.0, 0.3)
        P, W, fl = contract(g, lam)
        Pd = direct(g, lam)
        # w vs n
        n, w = W[:, 0], W[:, 1]
        sel = n >= 10
        b = np.polyfit(np.log(n[sel]), np.log(w[sel] + 1), 1)
        root_kids = W[np.argsort(-n)[:3]]
        out[T] = dict(P_contract=P, P_direct=Pd, rel_err=abs(P - Pd) / Pd, w_max=int(w.max()),
                      w_root_children=root_kids[:, 1].tolist(), n_root_children=root_kids[:, 0].tolist(),
                      slope_logw_logn=float(b[0]), mean_w=float(w.mean()), flops=float(fl),
                      w_by_n={str(k): float(w[(n >= k) & (n < 3 * k)].mean()) for k in (1, 3, 10, 30, 100, 300, 1000, 3000) if ((n >= k) & (n < 3 * k)).any()},
                      time=time.time() - t0)
        print(T, out[T], flush=True)
        json.dump(out, open('results/nested.json', 'w'), indent=1)
