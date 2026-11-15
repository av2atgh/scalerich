"""Planar-graph analysis shared by tissues (label images) and venation (vein graphs).

Input: vertices (y, x); edges as (u, v, pts) with pts an (m, 2) polyline from u to v
(may be just the endpoints), optional per-edge width.
Output: junction classes (T / Y / X / ...), 'lines' (edges merged through collinear
pairs at T and X), dual graph of lines (stem line -> through line at T), line widths,
corners per face (when faces are given).
"""
import numpy as np


def arm_dir(pts, r1, r2):
    """direction of a polyline leaving pts[0], from points at arc distance r1..r2."""
    d = np.hypot(*(pts - pts[0]).T)
    sel = (d >= r1) & (d <= r2)
    if sel.sum() == 0:
        sel = d > 0
        if sel.sum() == 0:
            return None
        j = np.argmax(d >= min(d.max(), r2)) if (d >= r2).any() else len(d) - 1
        v = pts[j] - pts[0]
    else:
        v = pts[sel].mean(0) - pts[0]
    n = np.linalg.norm(v)
    return v / n if n > 0 else None


def classify(V, E, r1=3.0, r2=10.0, tol=25.0):
    """E: list of (u, v, pts). Returns junction records and per-edge ids."""
    arms = {}
    for e, (u, v, pts) in enumerate(E):
        pts = np.asarray(pts, float)
        if len(pts) < 2:
            pts = np.array([V[u], V[v]], float)
        a = arm_dir(pts, r1, r2); b = arm_dir(pts[::-1], r1, r2)
        if a is not None:
            arms.setdefault(u, []).append((e, a))
        if b is not None:
            arms.setdefault(v, []).append((e, b))
    cos_t = np.cos(np.radians(180 - tol))
    recs = []
    for j, a in arms.items():
        m = len(a)
        if m < 3:
            continue
        Vd = np.array([x[1] for x in a]); ids = [x[0] for x in a]
        Cm = Vd @ Vd.T
        pairs = sorted((Cm[p, q], p, q) for p in range(m) for q in range(p + 1, m))
        used = set(); through = []
        for c, p, q in pairs:
            if c > cos_t:
                break
            if p in used or q in used:
                continue
            used |= {p, q}; through.append((ids[p], ids[q]))
        stems = [ids[p] for p in range(m) if p not in used]
        # angles: largest angle between arms (180 = perfect T), for diagnostics
        amax = float(np.degrees(np.arccos(np.clip(Cm[np.triu_indices(m, 1)].min(), -1, 1))))
        if m == 3:
            kind = 'T' if through else 'Y'
        elif m == 4:
            kind = 'X' if len(through) == 2 else ('T+' if through else 'Y+')
        else:
            kind = str(m)
        recs.append(dict(j=j, m=m, kind=kind, through=through, stems=stems, amax=amax))
    return recs


def lines(nE, recs, width=None, length=None):
    par = np.arange(nE)

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for r in recs:
        for p, q in r['through']:
            a, b = f(p), f(q)
            if a != b:
                par[a] = b
    root = np.array([f(x) for x in range(nE)])
    _, lid = np.unique(root, return_inverse=True)
    nL = lid.max() + 1 if nE else 0
    D = []
    for r in recs:
        if not r['through']:
            continue
        host = lid[r['through'][0][0]]
        for s in r['stems']:
            if lid[s] != host:
                D.append((lid[s], host))
    D = np.array(D, int).reshape(-1, 2)
    out = dict(lid=lid, nL=nL, D=D)
    if length is not None:
        out['L'] = np.bincount(lid, weights=length, minlength=nL)
        if width is not None:
            out['W'] = np.bincount(lid, weights=width * length, minlength=nL) / np.maximum(out['L'], 1e-12)
    return out


def dual_stats(nL, D, L=None, W=None):
    k = np.bincount(D.ravel(), minlength=nL) if len(D) else np.zeros(nL, int)
    adj = [set() for _ in range(nL)]
    for a, b in D:
        if a != b:
            adj[a].add(b); adj[b].add(a)
    kk = np.array([len(s) for s in adj])
    C = np.zeros(nL)
    for i in range(nL):
        n = list(adj[i]); m = len(n)
        if m >= 2:
            t = sum(1 for x in range(m) for y in range(x + 1, m) if n[y] in adj[n[x]])
            C[i] = 2 * t / (m * (m - 1))
    out = dict(nL=int(nL), kmean=float(k.mean()) if nL else 0, kmax=int(k.max()) if nL else 0,
               C=float(C[kk >= 2].mean()) if (kk >= 2).any() else float('nan'))
    for q in (2, 4, 8, 16, 32):
        out[f'P(k>={q})'] = float((k >= q).mean()) if nL else 0
    if (k >= 4).sum() >= 10:
        x = k[k >= 4].astype(float); out['gamma4'] = float(1 + len(x) / np.sum(np.log(x / 3.5)))
    ck = []
    for lo, hi in [(2, 3), (3, 5), (5, 8), (8, 13), (13, 21), (21, 1000)]:
        s = (kk >= lo) & (kk < hi)
        if s.sum() >= 5:
            ck.append((lo, int(s.sum()), round(float(C[s].mean()), 3)))
    out['Ck'] = ck
    if L is not None and nL > 10:
        s = (L > 0)
        out['corr_k_L'] = float(np.corrcoef(k[s], L[s])[0, 1])
        if W is not None:
            out['corr_k_W'] = float(np.corrcoef(k[s], W[s])[0, 1])
            m = s & (k > 0) & (W > 0)
            if m.sum() > 20:
                X = np.c_[np.ones(m.sum()), np.log(L[m])]
                rk = np.log(k[m]) - X @ np.linalg.lstsq(X, np.log(k[m]), rcond=None)[0]
                rw = np.log(W[m]) - X @ np.linalg.lstsq(X, np.log(W[m]), rcond=None)[0]
                out['partial_k_W_given_L'] = float(np.corrcoef(rk, rw)[0, 1])
    return out, k


def kinds(recs):
    c = {}
    for r in recs:
        c[r['kind']] = c.get(r['kind'], 0) + 1
    n3 = c.get('T', 0) + c.get('Y', 0)
    return c, (c.get('T', 0) / n3 if n3 else float('nan'))


# ------------------------------------------------------------ label image -> graph
def label_graph(Lab, merge=2):
    """2D label image (0 = background) -> vertices where >=3 labels meet (2x2 windows),
    walls (a, b) as polylines between vertices. Returns V, E [(u, v, pts, (a, b))], cells."""
    from scipy import ndimage as ndi
    a = Lab[:-1, :-1]; b = Lab[1:, :-1]; c = Lab[:-1, 1:]; d = Lab[1:, 1:]
    W = np.stack([a, b, c, d], -1)
    srt = np.sort(W, -1)
    nd = 1 + (srt[..., 1] != srt[..., 0]) + (srt[..., 2] != srt[..., 1]) + (srt[..., 3] != srt[..., 2])
    vmask = nd >= 3
    vl, nv = ndi.label(ndi.binary_dilation(vmask, iterations=merge) & vmask | vmask, structure=np.ones((3, 3)))
    vl, nv = ndi.label(ndi.binary_dilation(vmask, iterations=merge), structure=np.ones((3, 3)))
    vl = vl * 1  # vertex cluster id map on the window grid
    Vc = np.array(ndi.center_of_mass(vmask, vl, np.arange(1, nv + 1))) + 0.5
    # labels meeting at each vertex
    vy, vx = np.nonzero(vmask)
    vlabs = {}
    for y, x in zip(vy, vx):
        vlabs.setdefault(vl[y, x] - 1, set()).update(W[y, x].tolist())
    # wall pixels: windows with exactly 2 labels, outside vertex clusters
    wmask = (nd == 2) & (vl == 0)
    lo = srt[..., 0]; hi = srt[..., 3]
    key = lo.astype(np.int64) * (Lab.max() + 1) + hi
    wy, wx = np.nonzero(wmask)
    kk = key[wy, wx]
    order = np.argsort(kk); wy, wx, kk = wy[order], wx[order], kk[order]
    uk, st = np.unique(kk, return_index=True); en = np.r_[st[1:], len(kk)]
    # vertices for each label pair
    pv = {}
    for v, labs in vlabs.items():
        labs = sorted(labs)
        for i in range(len(labs)):
            for j in range(i + 1, len(labs)):
                pv.setdefault(labs[i] * (Lab.max() + 1) + labs[j], []).append(v)
    E = []
    for k_, s, e in zip(uk, st, en):
        vs = pv.get(k_, [])
        if len(vs) != 2:
            continue                      # wall with != 2 end vertices (closed / split): skip
        u, v = vs
        P = np.c_[wy[s:e], wx[s:e]] + 0.5
        # order wall pixels by distance from u
        du = np.hypot(*(P - Vc[u]).T)
        P = P[np.argsort(du)]
        pts = np.vstack([Vc[u], P, Vc[v]])
        a_, b_ = divmod(int(k_), int(Lab.max() + 1))
        E.append((u, v, pts, (a_, b_)))
    return Vc, E, vlabs


def corners_per_cell(Vc, E, recs, vlabs, tol=25.0):
    """geometric corners of each cell: vertices on its boundary that are not the flat
    (through) side of a T for that cell."""
    rec_by = {r['j']: r for r in recs}
    cells = {}
    for v, labs in vlabs.items():
        for c in labs:
            if c == 0:
                continue
            cells.setdefault(c, []).append(v)
    E_lab = [e[3] for e in E]
    out = {}
    for c, vs in cells.items():
        n = 0
        for v in vs:
            r = rec_by.get(v)
            if r is None or r['m'] != 3 or not r['through']:
                n += 1; continue
            p, q = r['through'][0]
            # flat for c iff c is on both through walls and not on the stem
            if c in E_lab[p] and c in E_lab[q] and c not in E_lab[r['stems'][0]]:
                continue
            n += 1
        out[c] = (n, len(vs))
    return out
