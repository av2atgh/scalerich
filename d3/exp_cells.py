"""Goals 1-2: volume splitting, BST, cell/plate bookkeeping, plate areas.

usage: python3 exp_cells.py run1 [run2 ...]  -> results/cells.txt
"""
import sys, os
import numpy as np
from common import load, RES

lines = []


def out(s):
    print(s, flush=True); lines.append(s)


def first_plate_area(n=400_000, seed=0):
    """MC: mean area/perimeter of the section of the unit cube by an isotropic
    plane through a uniform point (the 3D analogue of ell_0 = 0.946)."""
    rng = np.random.default_rng(seed)
    V = np.array([[x, y, z] for x in (0, 1) for y in (0, 1) for z in (0, 1)], float)
    Ed = [(i, j) for i in range(8) for j in range(i + 1, 8) if np.abs(V[i] - V[j]).sum() == 1]
    nvec = rng.normal(size=(n, 3)); nvec /= np.linalg.norm(nvec, axis=1)[:, None]
    p = rng.random((n, 3)); c = (nvec * p).sum(1)
    d = V @ nvec.T - c  # 8 x n
    pts = np.full((12, n, 3), np.nan)
    for e, (i, j) in enumerate(Ed):
        s = d[i] / (d[i] - d[j])
        cr = (d[i] * d[j]) < 0
        q = V[i][None, :] + s[:, None] * (V[j] - V[i])[None, :]
        pts[e][cr] = q[cr]
    cen = np.nanmean(pts, axis=0)
    e1 = np.nan_to_num(pts - cen[None]); e1 = np.where(np.isnan(pts), 0, pts - cen[None])
    # basis in the plane
    a = np.where(np.abs(nvec[:, :1]) < 0.9, np.array([[1, 0, 0]]), np.array([[0, 1, 0]]))
    b1 = np.cross(nvec, a); b1 /= np.linalg.norm(b1, axis=1)[:, None]
    b2 = np.cross(nvec, b1)
    X = (e1 * b1[None]).sum(2); Y = (e1 * b2[None]).sum(2)
    ang = np.where(np.isnan(pts[..., 0]), 99, np.arctan2(Y, X))
    o = np.argsort(ang, axis=0)
    X = np.take_along_axis(X, o, 0); Y = np.take_along_axis(Y, o, 0)
    valid = np.take_along_axis(~np.isnan(pts[..., 0]), o, 0)
    cnt = valid.sum(0)
    A = np.zeros(n); P = np.zeros(n)
    for i in range(12):
        j = np.where(i + 1 < cnt, i + 1, 0)
        xj = np.take_along_axis(X, j[None], 0)[0]; yj = np.take_along_axis(Y, j[None], 0)[0]
        use = i < cnt
        A += np.where(use, X[i] * yj - xj * Y[i], 0)
        P += np.where(use, np.hypot(xj - X[i], yj - Y[i]), 0)
    A = np.abs(A) / 2
    return A.mean(), A.std() / np.sqrt(n), P.mean(), cnt.mean()


for name in sys.argv[1:]:
    S, E, C = load(name)
    T = len(S)
    g = S[:, 0] > 0
    t = S[g, 0]; dep = S[g, 1]; m = S[g, 3]; nc = S[g, 4]; a = S[g, 5]; p = S[g, 6]; u = S[g, 8]
    out(f'== {name}: T={T}, rejected steps={(~g).sum()}')
    # --- 1. uniform split and BST
    h, _ = np.histogram(u, bins=10, range=(0, 1))
    ks = np.max(np.abs(np.sort(u) - (np.arange(1, len(u) + 1) - 0.5) / len(u)))
    out(f'[split] daughter volume fraction: mean {u.mean():.4f} var {u.var():.5f} (1/12={1/12:.5f}); '
        f'decile counts/expected: {np.round(h / (len(u) / 10), 3).tolist()}; KS D={ks:.5f} (n={len(u)})')
    sel = t >= 1000
    lt = np.log(t[sel])
    bins = np.linspace(lt.min(), lt.max(), 30)
    ib = np.digitize(lt, bins)
    xs, mu, va = [], [], []
    for b in np.unique(ib):
        s = ib == b
        if s.sum() > 200:
            xs.append(lt[s].mean()); mu.append(dep[sel][s].mean()); va.append(dep[sel][s].var())
    out(f'[BST] depth of split at t: slope d<depth>/dln t = {np.polyfit(xs, mu, 1)[0]:.4f} (2), '
        f'd var/dln t = {np.polyfit(xs, va, 1)[0]:.4f} (2); <depth>(T)-2lnT = {mu[-1] - 2 * xs[-1]:.3f}')
    v = C['v']
    x = v * len(v)
    out(f'[volumes] cells={len(v)} sum V={v.sum():.6f}; P(NV>x) at x=0.5,1,2,3,5: '
        + ' '.join(f'{(x > q).mean():.4f}/{np.exp(-q):.4f}' for q in (0.5, 1, 2, 3, 5))
        + f'; sum V^2 = {np.sum(v**2):.4e} vs 2/(T+2) = {2/(T+2):.4e}')
    # --- 2. bookkeeping
    nf = C['nf']; nvt = C['nv']
    Vc = nvt / 3; Ec = nvt / 2
    mb = m.mean()
    out(f'[cells] <F>={nf.mean():.4f} pred (6+sum(m+2))/(T+1)={(6 + np.sum(m + 2)) / (T + 1):.4f} = 2+<m> -> {2 + mb:.4f}; '
        f'<V>={Vc.mean():.4f} pred (8+2 sum m)/(T+1)={(8 + 2 * m.sum()) / (T + 1):.4f}; '
        f'<E>={Ec.mean():.4f} pred (12+3 sum m)/(T+1)={(12 + 3 * m.sum()) / (T + 1):.4f}')
    fr = np.bincount(nf, minlength=14)[4:14] / len(nf)
    out('[cells] number fractions F=4..13: ' + ' '.join(f'{q:.4f}' for q in fr)
        + f'; volume-weighted <F>={np.sum(nf * v) / v.sum():.4f}; cells touching cube: {C["tc"].mean():.4f}')
    out(f'[m] <m> all steps={mb:.4f}; interior plates (no cube face) <m>={m[nc == 0].mean():.4f} '
        f'(fraction interior {np.mean(nc == 0):.4f})')
    fm = np.bincount(m[nc == 0].astype(int), minlength=12)[3:12] / np.sum(nc == 0)
    out('[m] interior m=3..11 fractions: ' + ' '.join(f'{q:.4f}' for q in fm))
    for lo, hi in ((1e2, 1e3), (1e3, 1e4), (1e4, 1e5), (1e5, 1e6), (1e6, 1e7)):
        s = (t >= lo) & (t < hi) & (nc == 0)
        if s.sum() > 100:
            out(f'[m(t)] t in [{lo:.0e},{hi:.0e}): interior <m>={m[s].mean():.4f} +- {m[s].std() / np.sqrt(s.sum()):.4f}  '
                f'parent <F>={S[g, 2][s].mean():.4f}  <a t^2/3>={np.mean(a[s] * t[s]**(2/3)):.4f}  '
                f'<p t^1/3>={np.mean(p[s] * t[s]**(1/3)):.4f}')
    # area law
    A = np.cumsum(a)
    for q in (1e3, 1e4, 1e5, 1e6, T):
        i = int(q) - 1
        if i < len(A):
            out(f'[A(t)] t={int(q)}: A={A[i]:.3f}  A/t^(1/3)={A[i] / q**(1/3):.4f}  L_perim={np.sum(p[:i+1]) / q**(2/3):.4f} t^(2/3)')
    s = t > 1e4
    ca = np.polyfit(np.log(t[s]), np.log(a[s]), 1)[0]
    out(f'[a_t] fit log a_t vs log t (t>1e4): slope {ca:.4f} (-2/3); <a_1>={a[0]:.4f} (single draw)')
    # plate area distribution of all plates (power law a^-5/2 expected)
    aT = np.mean(a[t > T / 2])
    for f in (10, 30, 100):
        xx = a[a >= f * aT]
        out(f'[P(a)] MLE exponent of plate areas a >= {f} <a_T> (n={len(xx)}): {1 + len(xx) / np.sum(np.log(xx / (f * aT))):.4f} (5/2)')
    # area relative to parent: section area / V^(2/3)
    r = a / S[g, 7] ** (2 / 3)
    out(f'[shape] <a/V_parent^(2/3)>={r.mean():.4f}, <p/V^(1/3)>={np.mean(p / S[g, 7]**(1/3)):.4f}, '
        f'<m> for parent F=4,5,6,7,8: ' + ' '.join(f'{m[S[g, 2] == f].mean():.3f}' for f in (4, 5, 6, 7, 8)))

A1, eA1, P1, n1 = first_plate_area()
out(f'[a_1] MC (4e5 draws) first plate in the unit cube: <a_1>={A1:.4f}+-{eA1:.4f}, <p_1>={P1:.4f}, <m_1>={n1:.4f}')
open(os.path.join(RES, 'cells.txt'), 'w').write('\n'.join(lines) + '\n')
