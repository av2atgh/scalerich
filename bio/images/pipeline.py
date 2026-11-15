"""Fibre-network extraction for planar cytoskeleton images.

image -> ridge filter -> binary -> skeleton -> pruned skeleton graph
      -> junction classes (T / Y / X / higher) from arm directions
      -> fibres (branches merged through collinear pairs at T and X)
      -> dual graph (fibres = nodes, T-contacts stem -> through = edges)

Run with bio/.venv/bin/python (scikit-image).
"""
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import sato, threshold_otsu, gaussian
from skimage.morphology import skeletonize, remove_small_objects, dilation, disk

K8 = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])


def segment(img, sigmas=(1, 2, 3), mask=None, thr_scale=1.0, min_obj=30):
    img = img.astype(float)
    img = (img - np.percentile(img, 1)) / (np.percentile(img, 99.8) - np.percentile(img, 1) + 1e-12)
    img = np.clip(img, 0, 1)
    R = sato(img, sigmas=sigmas, black_ridges=False)
    vals = R[mask] if mask is not None else R.ravel()
    th = threshold_otsu(vals) * thr_scale
    B = R > th
    if mask is not None:
        B &= mask
    B = remove_small_objects(B, max_size=min_obj)
    return img, R, B


def _nbr(S):
    return ndi.convolve(S.astype(np.uint8), K8, mode='constant') * S


def prune(S, lmin):
    """remove terminal branches shorter than lmin (iterated twice)."""
    S = S.copy()
    for _ in range(2):
        nb = _nbr(S)
        J = (nb >= 3)
        lab, n = ndi.label(S & ~dilation(J, disk(1)), structure=np.ones((3, 3)))
        if n == 0:
            break
        sizes = ndi.sum(np.ones_like(lab), lab, index=np.arange(1, n + 1))
        ends = (nb == 1)
        has_end = ndi.maximum(ends, lab, index=np.arange(1, n + 1)) > 0
        Jd = dilation(J, disk(2))
        touchJ = ndi.maximum(Jd, lab, index=np.arange(1, n + 1)) > 0
        kill = np.where(has_end & touchJ & (sizes < lmin))[0] + 1
        S[np.isin(lab, kill)] = False
        S = skeletonize(S)
    return S


def graph(S, r1=6, r2=18, merge_len=8):
    """junction clusters (merged across branches shorter than merge_len),
    branches, arm directions fitted on branch pixels at distance r1..r2."""
    nb = _nbr(S)
    Jpix = nb >= 3
    Jlab, nJ = ndi.label(dilation(Jpix, disk(1)) & S, structure=np.ones((3, 3)))
    Bmask = S & (Jlab == 0)
    Blab, nB = ndi.label(Bmask, structure=np.ones((3, 3)))
    Blen = ndi.sum(np.ones_like(Blab), Blab, np.arange(1, nB + 1))
    JlabD = ndi.grey_dilation(Jlab, size=(3, 3))
    ys, xs = np.nonzero(Bmask & (JlabD > 0))
    pairs = set(zip(Blab[ys, xs] - 1, JlabD[ys, xs] - 1))
    bj = {}
    for b, j in pairs:
        bj.setdefault(b, set()).add(j)
    # merge junctions joined by short branches
    par = np.arange(nJ)

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    internal = set()
    for b, js in bj.items():
        js = list(js)
        if len(js) == 2 and Blen[b] < merge_len:
            a, c = f(js[0]), f(js[1])
            if a != c:
                par[a] = c
            internal.add(b)
    root = np.array([f(x) for x in range(nJ)])
    ur, jm = np.unique(root, return_inverse=True)
    nJm = len(ur)
    Jl2 = np.zeros_like(Jlab); Jl2[Jlab > 0] = jm[Jlab[Jlab > 0] - 1] + 1
    # internal branch pixels belong to their merged junction
    for b in internal:
        Jl2[Blab == b + 1] = jm[list(bj[b])[0]] + 1
    Jc = np.array(ndi.center_of_mass(np.ones_like(Jl2), Jl2, np.arange(1, nJm + 1))).reshape(-1, 2)
    by, bx = np.nonzero(Bmask)
    bl = Blab[by, bx] - 1
    order = np.argsort(bl)
    by, bx, bl = by[order], bx[order], bl[order]
    starts = np.searchsorted(bl, np.arange(nB)); stops = np.searchsorted(bl, np.arange(nB), 'right')
    arms = {j: [] for j in range(nJm)}
    for b, js in bj.items():
        if b in internal:
            continue
        for j0 in set(jm[list(js)]):
            py, px = by[starts[b]:stops[b]], bx[starts[b]:stops[b]]
            c = Jc[j0]
            d = np.hypot(py - c[0], px - c[1])
            sel = (d >= r1) & (d <= r2)
            if sel.sum() < 3:
                sel = d <= np.sort(d)[min(len(d) - 1, 8)]
            v = np.array([py[sel].mean() - c[0], px[sel].mean() - c[1]])
            nv = np.linalg.norm(v)
            if nv > 0:
                arms[j0].append((b, v / nv))
    return dict(Jc=Jc, nJ=nJm, Blab=Blab, nB=nB, Blen=Blen, arms=arms, Jlab=Jl2)


def classify(G, tol=25.0):
    """T: one arm pair within tol of collinear (3 arms); Y otherwise.
    X: 4 arms forming two collinear pairs. Returns junction records."""
    recs = []
    cos_t = np.cos(np.radians(180 - tol))
    for j, a in G['arms'].items():
        m = len(a)
        if m < 3:
            continue
        V = np.array([x[1] for x in a]); ids = [x[0] for x in a]
        Cm = V @ V.T
        pairs = sorted(((Cm[p, q], p, q) for p in range(m) for q in range(p + 1, m)))
        used = set(); through = []
        for c, p, q in pairs:                 # most anti-parallel first
            if c > cos_t:
                break
            if p in used or q in used:
                continue
            used |= {p, q}; through.append((ids[p], ids[q]))
        stems = [ids[p] for p in range(m) if p not in used]
        if m == 3:
            kind = 'T' if len(through) == 1 else 'Y'
        elif m == 4:
            kind = 'X' if len(through) == 2 else ('T+' if len(through) == 1 else 'Y+')
        else:
            kind = f'{m}'
        recs.append(dict(j=j, m=m, kind=kind, through=through, stems=stems))
    return recs


def fibres(G, recs):
    """union branches along through-pairs; dual edges stem-fibre -> through-fibre."""
    par = np.arange(G['nB'])

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for r in recs:
        for p, q in r['through']:
            a, b = f(p), f(q)
            if a != b:
                par[a] = b
    root = np.array([f(x) for x in range(G['nB'])])
    uniq, fid = np.unique(root, return_inverse=True)
    nF = len(uniq)
    flen = np.bincount(fid, weights=G['Blen'], minlength=nF)
    E = []; nonedge = 0
    for r in recs:
        if not r['through']:
            nonedge += 1
            continue
        host = fid[r['through'][0][0]]
        for s in r['stems']:
            if fid[s] != host:
                E.append((fid[s], host))
        if r['kind'] == 'X' or r['kind'] == 'T+':
            # crossing: second through pair is a fibre crossing the first (no T-contact)
            pass
    return fid, nF, flen, np.array(E, dtype=int).reshape(-1, 2)


def fibre_profile(fid, G, img, dist):
    """mean intensity and width (2 x distance transform) per fibre."""
    Blab = G['Blab']; m = Blab > 0
    f_of_pix = fid[Blab[m] - 1]
    nF = fid.max() + 1
    cnt = np.bincount(f_of_pix, minlength=nF)
    I = np.bincount(f_of_pix, weights=img[m], minlength=nF) / np.maximum(cnt, 1)
    W = np.bincount(f_of_pix, weights=2 * dist[m], minlength=nF) / np.maximum(cnt, 1)
    return I, W


def analyse(img, mask=None, sigmas=(1, 2, 3), lmin=8, r1=6, r2=18, merge_len=8, tol=25.0, thr_scale=1.0):
    imgn, R, B = segment(img, sigmas, mask, thr_scale)
    S = skeletonize(B)
    S = prune(S, lmin)
    G = graph(S, r1, r2, merge_len)
    recs = classify(G, tol)
    fid, nF, flen, E = fibres(G, recs)
    dist = ndi.distance_transform_edt(B)
    I, W = fibre_profile(fid, G, imgn, dist)
    kinds = {}
    for r in recs:
        kinds[r['kind']] = kinds.get(r['kind'], 0) + 1
    k = np.bincount(E.ravel(), minlength=nF) if len(E) else np.zeros(nF, int)
    return dict(B=B, S=S, G=G, recs=recs, kinds=kinds, fid=fid, nF=nF, flen=flen, E=E,
                k=k, I=I, W=W)


def summary(res):
    k = res['k']; kinds = res['kinds']; n3 = sum(v for kk, v in kinds.items())
    out = dict(junctions=n3, kinds=kinds,
               T_frac=kinds.get('T', 0) / max(1, kinds.get('T', 0) + kinds.get('Y', 0)),
               X_per_junction=kinds.get('X', 0) / max(1, n3),
               nF=int(res['nF']), kmean=float(k.mean()), kmax=int(k.max()) if len(k) else 0)
    for q in (2, 4, 8, 16, 32):
        out[f'P(k>={q})'] = float((k >= q).mean())
    if (k >= 4).sum() > 10:
        x = k[k >= 4].astype(float)
        out['gamma4'] = float(1 + len(x) / np.sum(np.log(x / 3.5)))
    s = res['flen'] > 0
    if s.sum() > 10:
        out['corr_k_len'] = float(np.corrcoef(k[s], res['flen'][s])[0, 1])
        out['corr_k_width'] = float(np.corrcoef(k[s], res['W'][s])[0, 1])
        out['corr_k_int'] = float(np.corrcoef(k[s], res['I'][s])[0, 1])
    return out


# ---------------------------------------------------------------- synthetic
def render(A, B, width, M=1024, psf=1.5, noise=0.05, seed=0):
    """draw segments A->B (unit coords) with given widths (px) into an M x M image."""
    rng = np.random.default_rng(seed)
    img = np.zeros((M, M))
    A = A * (M - 1); B = B * (M - 1)
    for a, b, w in zip(A, B, width):
        lo = np.floor(np.minimum(a, b) - 3 * w - 2).astype(int).clip(0, M - 1)
        hi = np.ceil(np.maximum(a, b) + 3 * w + 2).astype(int).clip(0, M - 1)
        yy, xx = np.mgrid[lo[1]:hi[1] + 1, lo[0]:hi[0] + 1]
        d = b - a; L2 = d @ d + 1e-12
        t = (((xx - a[0]) * d[0] + (yy - a[1]) * d[1]) / L2).clip(0, 1)
        dist = np.hypot(xx - (a[0] + t * d[0]), yy - (a[1] + t * d[1]))
        v = np.clip(w / 2 + 0.5 - dist, 0, 1)
        sub = img[lo[1]:hi[1] + 1, lo[0]:hi[0] + 1]
        np.maximum(sub, v, out=sub)
    img = gaussian(img, psf)
    img = img / img.max() + noise * rng.standard_normal(img.shape)
    return img
