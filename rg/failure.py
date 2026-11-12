"""Quasistatic buckling cascade on beam frames: a cheap proxy for the
compaction-band question (Goal 2b/3).

Uniaxial compression along x (left ux=0, right ux=-eps, top/bottom free).
Segment s buckles when its compressive axial force exceeds Euler's load
P_c = pi^2 EI / L^2 (E=1, I=lam^3/12; L_eff = max(L, lam)).  A buckled
segment keeps its bending stiffness, its axial stiffness drops to 1e-3.
Extremal dynamics: raise eps to the next threshold, buckle that segment,
re-solve at fixed eps, buckle every segment now over threshold (avalanche),
repeat.  Stop after NF buckled segments.

Systems (rho = 0.2, T = 1000 ligaments, or a square grid of equal density):
  SR alpha = 0, 0.3, 0.6;  SR alpha = 0.6 with thicknesses shuffled among
  ligaments (same P(lam), no age-thickness correlation);  square grid with
  5% thickness noise.
Localisation metrics on the first NF buckled segments (midpoints):
  band  = max fraction inside an x-strip of width 0.1 (uniform: ~0.1-0.2)
  nnr   = mean nearest-neighbour distance / that of an equal-size random
          sample of segment midpoints (clustered < 1)
  aval  = largest avalanche, mean avalanche size
Writes results/failure.json
"""
import json, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import srtree, frame


def grid(n, lam, rng, noise=0.05):
    xs = (np.arange(n) + 0.5) / n
    pts = {}; nodes = []; seg = []; segl = []
    def nid(x, y):
        k = (round(x, 9), round(y, 9))
        if k not in pts:
            pts[k] = len(nodes); nodes.append([x, y])
        return pts[k]
    lid = 0
    for y in xs:
        xx = np.r_[0, xs, 1]
        for a, b in zip(xx[:-1], xx[1:]):
            seg.append((nid(a, y), nid(b, y))); segl.append(lid)
        lid += 1
    for x in xs:
        yy = np.r_[0, xs, 1]
        for a, b in zip(yy[:-1], yy[1:]):
            seg.append((nid(x, a), nid(x, b))); segl.append(lid)
        lid += 1
    nodes = np.array(nodes); seg = np.array(seg); segl = np.array(segl)
    laml = lam * (1 + noise * rng.standard_normal(lid))
    bn = {-4: [i for i, p in enumerate(nodes) if p[0] == 0], -2: [i for i, p in enumerate(nodes) if p[0] == 1],
          -1: [i for i, p in enumerate(nodes) if p[1] == 0], -3: [i for i, p in enumerate(nodes) if p[1] == 1]}
    return nodes, seg, laml[segl], bn


def cascade(nodes, seg, lam_s, bn, NF, rng):
    N = len(nodes)
    d = nodes[seg[:, 1]] - nodes[seg[:, 0]]
    L = np.linalg.norm(d, axis=1); c = d[:, 0] / L; s = d[:, 1] / L
    Le = np.maximum(L, lam_s)
    EA = lam_s.copy(); EI = lam_s ** 3 / 12
    Pc = np.pi ** 2 * EI / Le ** 2
    dof = np.c_[3 * seg[:, 0], 3 * seg[:, 0] + 1, 3 * seg[:, 0] + 2,
                3 * seg[:, 1], 3 * seg[:, 1] + 1, 3 * seg[:, 1] + 2]
    I = np.repeat(dof, 6, axis=1).ravel(); J = np.tile(dof, (1, 6)).ravel()
    lo, hi = np.array(bn[-4]), np.array(bn[-2])
    fixedmask = np.zeros(3 * N, bool); fixedmask[3 * lo] = True; fixedmask[3 * hi] = True; fixedmask[3 * lo[0] + 1] = True
    free = ~fixedmask
    ub = np.zeros(3 * N); ub[3 * hi] = -1.0

    def solve(EA):
        kg, R, kl = frame.element_k(Le, c, s, EA, EI)
        K = sp.coo_matrix((kg.ravel(), (I, J)), shape=(3 * N, 3 * N)).tocsr()
        Kff = K[free][:, free].tocsc()
        Kff = Kff + 1e-14 * abs(Kff.diagonal()).max() * sp.eye(Kff.shape[0], format='csc')
        u = ub.copy()
        u[free] = spl.spsolve(Kff, -(K[free][:, fixedmask] @ ub[fixedmask]))
        ul = np.einsum('nij,nj->ni', R, u[dof])
        Nax = EA * (ul[:, 3] - ul[:, 0]) / Le      # tension > 0 at unit eps
        E = float(u @ (K @ u))
        return Nax, E

    buck = np.zeros(len(seg), bool)
    Nax, E0 = solve(EA)
    eps = 0.0; events = []; aval = []; E = E0
    while buck.sum() < NF:
        comp = np.where((~buck) & (Nax < 0))[0]
        if len(comp) == 0:
            break
        ec = Pc[comp] / (-Nax[comp])
        k = comp[np.argmin(ec)]; eps = float(ec.min())
        size = 0
        newly = [k]
        while newly:
            for q in newly:
                buck[q] = True; EA[q] = 1e-3 * lam_s[q]
                events.append((eps, q)); size += 1
            Nax, E = solve(EA)
            over = np.where((~buck) & (Nax < 0) & (-Nax * eps >= Pc))[0]
            newly = list(over)
            if buck.sum() >= NF:
                break
        aval.append(size)
    mid = 0.5 * (nodes[seg[:, 0]] + nodes[seg[:, 1]])
    fx = np.array([q for _, q in events[:NF]])
    P = mid[fx]
    xs = np.sort(P[:, 0])
    band = max(np.searchsorted(xs, x0 + 0.1) - np.searchsorted(xs, x0) for x0 in np.linspace(0, 0.9, 91)) / len(xs)
    def nnd(Q):
        D = np.linalg.norm(Q[:, None] - Q[None], axis=2); np.fill_diagonal(D, np.inf)
        return D.min(1).mean()
    rnd = np.mean([nnd(mid[rng.choice(len(seg), len(fx), replace=False)]) for _ in range(5)])
    # spread of buckling thresholds at eps=0 (initial), coefficient of variation of log critical strain
    return dict(band=float(band), nnr=float(nnd(P) / rnd), aval_max=int(max(aval)), aval_mean=float(np.mean(aval)),
                eps_first=float(events[0][0]), eps_last=float(events[min(NF, len(events)) - 1][0]),
                E_drop=float(E / E0), nseg=int(len(seg)))


def sr_system(T, a, rho, seed, shuffle=False):
    g = srtree.run(T, seed)
    l = np.linalg.norm(g['B'] - g['A'], axis=1)
    lam = srtree.thickness(T, 1.0, a)
    if shuffle:
        lam = np.random.default_rng(seed + 99).permutation(lam)
    lam *= rho / np.sum(lam * l)
    nodes, seg, sl, bn = srtree.network(g)
    return nodes, seg, lam[sl], bn


if __name__ == '__main__':
    T = 1000; rho = 0.2; NF = 150; NS = 3
    out = {}
    systems = [('SR a=0', lambda s: sr_system(T, 0.0, rho, 600 + s)),
               ('SR a=0.3', lambda s: sr_system(T, 0.3, rho, 600 + s)),
               ('SR a=0.6', lambda s: sr_system(T, 0.6, rho, 600 + s)),
               ('SR a=0.6 shuffled', lambda s: sr_system(T, 0.6, rho, 600 + s, True)),
               ('square 5% noise', lambda s: grid(25, rho / 50, np.random.default_rng(700 + s)))]
    for name, mk in systems:
        rs = []
        for s in range(NS):
            t0 = time.time()
            try:
                r = cascade(*mk(s), NF, np.random.default_rng(s))
            except Exception as e:
                print('FAILED', name, s, e); continue
            rs.append(r)
            print(name, s, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, 'time %.1f' % (time.time() - t0), flush=True)
        out[name] = dict(runs=rs, mean={k: float(np.mean([r[k] for r in rs])) for k in rs[0]} if rs else {})
        json.dump(out, open('results/failure.json', 'w'), indent=1)
    for name, v in out.items():
        print(name, {k: round(x, 3) for k, x in v['mean'].items()})
