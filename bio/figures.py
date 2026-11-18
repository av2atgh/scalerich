"""Figures for bio/manuscript.tex. Run from bio/: .venv/bin/python -I figures.py
Reads saved results in bio/results; the HLD null (fig 2b) is recomputed with
images/planar.py + beyond/divisions.py and cached in figs/cache_hld_null.npz."""
import os, sys, json, glob, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
B = os.path.dirname(os.path.abspath(__file__))
R = f'{B}/results'; F = f'{B}/figs'
os.makedirs(F, exist_ok=True)
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'lines.linewidth': 1.2})
C1, C2, C3, C4 = '#1f5fa8', '#c0392b', '#2e8b57', '#7f7f7f'


def panel(ax, s):
    ax.text(-0.2, 1.04, s, transform=ax.transAxes, fontweight='bold', fontsize=9)


# ------------------------------------------------------------------ fig 1: models
def fig_theory():  # -> figs/fig2.pdf
    sc = json.load(open(f'{R}/scan_N20000.json'))
    Pis = sorted(sc, key=float)
    x = np.array([float(p) for p in Pis])
    g = np.array([np.mean([o['gamma10'] for o in sc[p]]) for p in Pis])
    C = np.array([np.mean([o['C'] for o in sc[p]]) for p in Pis])
    km = np.array([np.mean([o['kmax'] for o in sc[p]]) for p in Pis])
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.2))
    a = ax[0]
    a.semilogx(x, g, 'o-', color=C1, label=r'$\gamma$ ($k\geq10$)')
    a.set_xlabel(r'$\Pi$'); a.set_ylabel(r'$\gamma$', color=C1); a.set_ylim(2.5, 10)
    a2 = a.twinx(); a2.semilogx(x, C, 's--', color=C2); a2.set_ylabel(r'$\langle C\rangle$', color=C2)
    a2.set_ylim(0, 0.45); panel(a, '(a)')
    # crossover n* from cohort mean degree relative to Pi -> 0
    N = 20000; edges = np.unique(np.logspace(0, np.log10(N), 30).astype(int))
    cen = np.sqrt(edges[:-1] * edges[1:])

    def cohort(P):
        ks = np.mean([np.load(f) for f in sorted(glob.glob(f'{R}/k_N{N}_Pi{P}_s*.npy'))], 0)
        return np.array([ks[lo:hi].mean() for lo, hi in zip(edges[:-1], edges[1:])])
    ref = cohort('0.001')
    a = ax[1]; ns = []
    for P, c in zip(['1000', '10000', '1e+06'], [C1, C2, C3]):
        r = cohort(P) / ref
        a.semilogx(cen, r, '-', color=c, label=rf'$\Pi=${float(P):.0e}'.replace('e+0', 'e').replace('e+', 'e'))
        j = np.where(r >= 0.5)[0]; ns.append((float(P), cen[j[0]]))
    a.axhline(0.5, color=C4, lw=0.6, ls=':')
    a.set_xlabel('birth rank $i$'); a.set_ylabel(r'$\bar k(i)/\bar k_{\Pi\to0}(i)$'); a.set_ylim(0, 1.6)
    a.legend(frameon=False, loc='upper left', fontsize=6); panel(a, '(b)')
    ins = a.inset_axes([0.66, 0.12, 0.3, 0.3])
    ns = np.array(ns); ins.loglog(ns[:, 0], ns[:, 1], 'ko', ms=3)
    xx = np.logspace(3, 6, 10); ins.loglog(xx, 0.04 * xx ** (2 / 3), 'k-', lw=0.7)
    ins.set_xlabel(r'$\Pi$', fontsize=6, labelpad=0); ins.set_ylabel(r'$n^*$', fontsize=6, labelpad=0)
    ins.tick_params(labelsize=5)
    a = ax[2]
    q = [2, 4, 8, 16, 32, 64]
    sys.path.insert(0, f'{B}/model')
    from pargrow import run
    sr = np.mean([[(np.bincount(r['E'].ravel(), minlength=1500) >= v).mean() for v in q]
                  for r in (run(1500, 1e-3, seed=s) for s in range(4))], 0)   # SR reference, N = 1500
    print('SR N=1500 ccdf', np.round(sr, 4))
    a.loglog(q, sr, 'k-', label='no turnover')
    for s, c in zip([0, 1, 2], [C2, C1, C3]):
        fn = f'{R}/turnover_nss1500_stab{s}.json'
        t = json.load(open(fn))['0.01']
        cc = np.mean([[v for _, v in o['ccdf']] for o in t], 0)
        qq = [k for k, _ in t[0]['ccdf']]
        m = cc > 0
        a.loglog(np.array(qq)[m], cc[m], 'o--', color=c, ms=3, label=f'turnover, $s={s}$')
    a.set_xticks(q); a.set_xticklabels([str(v) for v in q]); a.minorticks_off()
    a.set_xlabel('$k$'); a.set_ylabel(r'$P(K\geq k)$'); a.legend(frameon=False, fontsize=6)
    panel(a, '(c)')
    fig.tight_layout(); fig.savefig(f'{F}/fig2.pdf'); plt.close(fig)


# ------------------------------------------------------------------ fig 2: venation
def hld_null():
    cache = f'{F}/cache_hld_null.npz'
    if os.path.exists(cache):
        return dict(np.load(cache))
    sys.path.insert(0, f'{B}/images'); sys.path.insert(0, f'{B}/../beyond')
    import divisions
    from skimage.draw import polygon
    from scipy.spatial import cKDTree
    from planar import label_graph

    def hld(walls, area):
        par = {c: c for c in area}; A = dict(area); sp = []

        def f(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for a, b, w in sorted(walls, key=lambda x: x[2]):
            ra, rb = f(a), f(b)
            if ra == rb:
                continue
            sp.append(min(A[ra], A[rb]) / (A[ra] + A[rb])); par[ra] = rb; A[rb] += A[ra]
        return np.array(sp)
    M = 1536; T = 800; rng = np.random.default_rng(0); out = {}
    srs = []; srn = []; vor = []
    for seed in (2, 3, 4):
        r = divisions.run(T, seed=seed, rule='SR'); Lab = np.zeros((M, M), np.int32)
        for i, p in enumerate(r['leaf_polys']):
            p = np.asarray(p) * (M - 1); rr, cc = polygon(p[:, 1], p[:, 0], (M, M)); Lab[rr, cc] = i + 1
        V, E, vl = label_graph(Lab)
        A = r['A'] * (M - 1); Bq = r['B'] * (M - 1)
        area = {c: float(n) for c, n in zip(*np.unique(Lab[Lab > 0], return_counts=True))}
        births = []
        for u, v, pts, (a, b) in E:
            if a == 0:
                continue
            m = pts[len(pts) // 2][::-1]; d = Bq - A
            t = np.clip(((m - A) * d).sum(1) / (d * d).sum(1), 0, 1)
            births.append((a, b, int(np.argmin(np.hypot(*(A + t[:, None] * d - m).T)))))
        srs += list(hld([(a, b, -np.log(t + 1)) for a, b, t in births], area))
        srn += list(hld([(a, b, -np.log(t + 1) + 1.0 * rng.standard_normal()) for a, b, t in births], area))
        seeds = rng.random((T, 2)) * (M - 1); yy, xx = np.mgrid[0:M, 0:M]
        _, lab = cKDTree(seeds).query(np.c_[xx.ravel(), yy.ravel()])
        Lv = (lab + 1).reshape(M, M).astype(np.int32)
        V, E, vl = label_graph(Lv)
        area = {c: float(n) for c, n in zip(*np.unique(Lv, return_counts=True))}
        vor += list(hld([(a, b, rng.random()) for _, _, _, (a, b) in E if a > 0], area))
    out = dict(sr=np.array(srs), srnoise=np.array(srn), voronoi=np.array(vor))
    np.savez(cache, **out)
    return out


def fig_venation():  # -> figs/fig5.pdf
    files = sorted(glob.glob(f'{B}/data/venation/*.npz'))
    angs = []; thin = [0, 0]
    for f in files:
        d = np.load(f); n1 = d['e_EndNodes_1'].astype(int) - 1; n2 = d['e_EndNodes_2'].astype(int) - 1
        bnd = (d['e_Ai'] == 0) | (d['e_Aj'] == 0); typ = d['e_Type']; W = d['e_Width']

        def u(a):
            a = np.radians(a); return np.c_[-np.sin(a), np.cos(a)]
        ui = u(d['e_Or_ij']); uj = u(d['e_Or_ji']); arms = {}
        for e in range(len(n1)):
            if n1[e] == n2[e]:
                continue
            arms.setdefault(n1[e], []).append((e, ui[e])); arms.setdefault(n2[e], []).append((e, uj[e]))
        for j, a in arms.items():
            if len(a) != 3:
                continue
            es = [x[0] for x in a]
            if bnd[es].any() or not typ[es].all():
                continue
            Vd = np.array([x[1] for x in a]); Cm = Vd @ Vd.T; iu = np.triu_indices(3, 1)
            p = np.argmin(Cm[iu]); angs.append(np.degrees(np.arccos(np.clip(Cm[iu][p], -1, 1))))
            stem = 3 - iu[0][p] - iu[1][p]
            thin[1] += 1; thin[0] += W[es[stem]] <= min(W[es[iu[0][p]]], W[es[iu[1][p]]])
    angs = np.array(angs)
    sp = np.load(f'{R}/venation_pool.npz')['split']
    nul = hld_null()
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.3))
    a = ax[0]
    a.hist(angs, bins=np.arange(90, 181, 5), color=C1, alpha=0.8, density=True)
    a.axvline(180, color='k', lw=0.6); a.axvline(120, color='k', lw=0.6, ls=':')
    a.text(121, a.get_ylim()[1] * 0.9, 'Y', fontsize=7); a.text(172, a.get_ylim()[1] * 0.9, 'T', fontsize=7)
    a.set_xlabel(r'largest arm angle at 3-way node (deg)'); a.set_ylabel('density')
    a.text(0.03, 0.62, f'{len(angs)} junctions\narm opposite largest\nangle thinnest: {thin[0] / thin[1]:.2f}\n(unrelated widths: 0.33)',
           transform=a.transAxes, fontsize=6.5)
    panel(a, '(a)')
    a = ax[1]; bins = np.linspace(0, 0.5, 11); c = 0.5 * (bins[1:] + bins[:-1])
    for arr, lab, col, ls in [(sp, 'leaves (44)', 'k', '-'), (nul['sr'], 'recursive partitioning, birth-ordered widths', C1, '-'),
                              (nul['srnoise'], 'recursive partitioning, noisy widths', C1, '--'),
                              (nul['voronoi'], 'Voronoi, random widths', C2, '-')]:
        h, _ = np.histogram(arr, bins=bins, density=True)
        a.plot(c, h, ls, color=col, marker='o', ms=2.5, label=lab)
    a.axvspan(0.45, 0.5, color=C3, alpha=0.15); a.text(0.41, 0.3, 'Errera-like', fontsize=6, color=C3, rotation=90)
    a.set_xlabel(r'split fraction $\min(A_1,A_2)/(A_1+A_2)$'); a.set_ylabel('density'); a.set_ylim(0, 5)
    a.legend(frameon=False, fontsize=6, loc='upper center'); panel(a, '(b)')
    fig.tight_layout(); fig.savefig(f'{F}/fig5.pdf'); plt.close(fig)
    return len(angs), thin[0] / thin[1]


# ------------------------------------------------------------------ fig 3: meristem
def fig_plant():  # -> figs/fig4.pdf
    rel = []
    for f in sorted(glob.glob(f'{R}/meristem_plant*.json')):
        rel += json.load(open(f))['relax']
    rel = np.array(rel)
    ages = np.arange(0, 40, 4)
    med = []; q1 = []; q3 = []
    for a_ in ages:
        s = rel[:, 0] == a_
        med.append(np.median(rel[s, 1])); q1.append(np.percentile(rel[s, 1], 25)); q3.append(np.percentile(rel[s, 1], 75))
    fit = json.load(open(f'{R}/fit_tension.json'))
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.3))
    a = ax[0]
    a.fill_between(ages + 2, q1, q3, color=C4, alpha=0.25, lw=0)
    a.plot(ages + 2, med, 'ko', ms=3.5, label='shoot meristems (6 plants)')
    for t0, tau, col in [(0.6, 24, C1), (0.75, 48, C2), (0.3, 8, C3)]:
        r = [x for x in fit if x['t0'] == t0 and x['tau_h'] == tau][0]
        a.plot(ages + 2, r['curve'], '-', color=col,
               label=rf"$t_0={t0}$, $\tau={tau}$ h (rms {r['rms']:.1f}$^\circ$)")
    rf = np.array(json.load(open(f'{R}/refahi_stats.json'))['ang'])
    eb = [0, 8, 16, 24, 36, 48, 72, 120]; xc = []; mr = []
    for lo, hi in zip(eb[:-1], eb[1:]):
        s_ = (rf[:, 0] >= lo) & (rf[:, 0] < hi); xc.append(0.5 * (lo + hi)); mr.append(np.median(rf[s_, 1]))
    print('Refahi angle medians', np.round(mr, 1))
    a.plot(xc, mr, 'D', mfc='none', mec='k', ms=4, label='flower meristems (6 films)')
    rm = json.load(open(f'{R}/refahi_model.json'))['relaxed K=20']['angle']
    a.plot(xc, list(rm.values()), ':', color=C1, label=r'model, $T_{cc}=31$ h')
    a.axhline(120, color='k', lw=0.6, ls=':')
    a.set_xlabel('wall age (h)'); a.set_ylabel(r'through angle $\theta$ (deg)'); a.set_ylim(115, 165); a.set_xlim(0, 100)
    a.legend(frameon=False, fontsize=6); panel(a, '(a)')
    a = ax[1]
    mt = json.load(open(f'{R}/meristem_topo.json'))
    n = np.arange(3, 11); wt = np.array([mt[p]['nPn'] for p in mt])
    Pd = np.array([[mt[p]['Pn'][str(k)] for k in n] for p in mt])
    Pm = (Pd * wt[:, None]).sum(0) / wt.sum()
    pr = json.load(open(f'{R}/plant_predictions_N2000.json'))
    w = 0.2
    a.bar(n - 1.5 * w, Pm, w, color='k', label='shoot meristems')
    for i, (k, lab, col) in enumerate([('errera, no relax', 'Errera', C1),
                                       ('plant (errera+relax, size pick)', 'Errera + tension', C2),
                                       ('SR (no relax)', 'random cuts', C4)]):
        a.bar(n + (i - 0.5) * w, [pr[k]['Pn'][str(x)] for x in n], w, color=col, label=lab)
    a.set_xlabel('sides $n$'); a.set_ylabel('$P(n)$'); a.legend(frameon=False, fontsize=6); panel(a, '(b)')
    fig.tight_layout(); fig.savefig(f'{F}/fig4.pdf'); plt.close(fig)
    return dict(ages=(ages + 2).tolist(), med=np.round(med, 1).tolist(), Pn=np.round(Pm, 3).tolist(), nobs=int(wt.sum()))


def fig_actin():  # -> figs/figS1.png
    import shutil
    shutil.copy(f'{R}/overlay_Cell_011.png', f'{F}/figS1.png')


# ------------------------------------------------------------------ schematic
def _clip(poly, p, d):
    """chord of convex polygon poly (m,2, ccw) through p along d: returns (A, B, ia, ib)
    with A, B on edges ia, ib."""
    hits = []
    m = len(poly)
    for i in range(m):
        a, b = poly[i], poly[(i + 1) % m]
        e = b - a
        M = np.array([[d[0], -e[0]], [d[1], -e[1]]])
        if abs(np.linalg.det(M)) < 1e-14:
            continue
        t, u = np.linalg.solve(M, a - p)
        if -1e-12 <= u <= 1 + 1e-12:
            hits.append((t, i, p + t * d))
    hits.sort(key=lambda x: x[0])
    return hits[0][2], hits[-1][2], hits[0][1], hits[-1][1]


def _area(P):
    x, y = P[:, 0], P[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _split(poly, A, B, ia, ib):
    m = len(poly)
    p1 = [A]; i = ia
    while True:
        i = (i + 1) % m; p1.append(poly[i])
        if i == ib:
            break
    p1.append(B)
    p2 = [B]; i = ib
    while True:
        i = (i + 1) % m; p2.append(poly[i])
        if i == ia:
            break
    p2.append(A)
    return np.array(p1), np.array(p2)


def recursive_division(T, seed=3):
    """SR: cell prop. to area, uniform point, uniform angle. Each cell keeps, per edge, the
    id of the segment it lies on (-1 = boundary). Returns segments [(A, B, hostA, hostB)]."""
    rng = np.random.default_rng(seed)
    cells = [(np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float), [-1, -1, -1, -1])]
    segs = []
    for s in range(T):
        ar = np.array([_area(c[0]) for c in cells])
        j = rng.choice(len(cells), p=ar / ar.sum())
        P, lab = cells.pop(j)
        lo, hi = P.min(0), P.max(0)
        while True:   # uniform point (rejection in bbox; convex test)
            q = lo + rng.random(2) * (hi - lo)
            cr = [np.cross(P[(i + 1) % len(P)] - P[i], q - P[i]) for i in range(len(P))]
            if min(cr) > 0:
                break
        th = rng.random() * np.pi
        A, Bp, ia, ib = _clip(P, q, np.array([np.cos(th), np.sin(th)]))
        segs.append((A, Bp, lab[ia], lab[ib]))
        p1, p2 = _split(P, A, Bp, ia, ib)
        m = len(P)
        # edge labels: p1 = A, P[ia+1..ib], B ; edges: A->P[ia+1] (ia), ..., P[ib]->B (ib), B->A (s)
        l1 = []; i = ia
        while True:
            l1.append(lab[i]); i = (i + 1) % m
            if i == (ib + 1) % m:
                break
        l1.append(s)
        l2 = []; i = ib
        while True:
            l2.append(lab[i]); i = (i + 1) % m
            if i == (ia + 1) % m:
                break
        l2.append(s)
        cells += [(p1, l1), (p2, l2)]
    return segs


def fig_schematic():  # -> figs/fig1.pdf
    from matplotlib.patches import Polygon, FancyArrowPatch, Arc
    fig = plt.figure(figsize=(7.0, 2.2))
    gsp = fig.add_gridspec(1, 6, width_ratios=[0.8, 0.8, 0.8, 2.0, 1.5, 2.1], wspace=0.12,
                           left=0.01, right=0.995, top=0.9, bottom=0.03)
    segs = recursive_division(60, seed=3)
    lw = lambda i, n: 3.0 * (1 + i) ** -0.55 + 0.25
    sq = dict(color='k', lw=1.4)

    def draw(ax, n):
        ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], **sq)
        for i, (A, B, _, _) in enumerate(segs[:n]):
            ax.plot([A[0], B[0]], [A[1], B[1]], color=C1, lw=lw(i, n), solid_capstyle='butt')
        ax.set_xlim(-.03, 1.03); ax.set_ylim(-.03, 1.03); ax.set_aspect('equal'); ax.axis('off')
    for c, n in enumerate((1, 3, 10)):
        ax = fig.add_subplot(gsp[:, c])
        draw(ax, n)
        ax.set_title(f'{n} cut' + ('s' if n > 1 else ''), fontsize=7, pad=2)
        if c == 0:
            ax.text(-0.05, 1.1, 'a', transform=ax.transAxes, fontweight='bold', fontsize=9)
    # (b) contact network
    ax = fig.add_subplot(gsp[:, 3])
    draw(ax, len(segs))
    mid = np.array([(A + B) / 2 for A, B, _, _ in segs])
    deg = np.zeros(len(segs), int)
    for i, (_, _, ha, hb) in enumerate(segs):
        for h in (ha, hb):
            if h >= 0:
                ax.plot([mid[i, 0], mid[h, 0]], [mid[i, 1], mid[h, 1]], color=C2, lw=0.45, alpha=0.7, zorder=3)
                deg[i] += 1; deg[h] += 1
    ax.scatter(mid[:, 0], mid[:, 1], s=2 + 4 * deg, color=C2, edgecolor='w', lw=0.4, zorder=4)
    ax.text(-0.02, 1.06, 'b', transform=ax.transAxes, fontweight='bold', fontsize=9)
    # (c) plant cell division
    ax = fig.add_subplot(gsp[:, 4])
    ang = np.radians([10, 70, 125, 185, 240, 300])
    rad = np.array([1.0, 0.9, 1.05, 0.95, 1.0, 0.92])
    P = np.c_[1.25 * rad * np.cos(ang), rad * np.sin(ang)]
    ax.add_patch(Polygon(P, closed=True, fc='#eef3f8', ec='k', lw=1.4))
    # centroid
    x, y = P[:, 0], P[:, 1]; xs, ys = np.roll(x, -1), np.roll(y, -1)
    cr = x * ys - xs * y; Ar = cr.sum() / 2
    cen = np.array([((x + xs) * cr).sum(), ((y + ys) * cr).sum()]) / (6 * Ar)
    best = None
    for th in np.linspace(0, np.pi, 360, endpoint=False):
        A, Bp, ia, ib = _clip(P, cen, np.array([np.cos(th), np.sin(th)]))
        L = np.linalg.norm(Bp - A)
        if best is None or L < best[0]:
            best = (L, A, Bp)
    _, A, Bp = best
    ax.plot([A[0], Bp[0]], [A[1], Bp[1]], color=C3, lw=2.0)
    ax.plot(*cen, 'o', color=C3, ms=3)
    # corner-cutting chord: across vertex 1 (between edges 0 and 1)
    v = P[2]; a = P[2] + 0.3 * (P[1] - P[2]); b = P[2] + 0.55 * (P[3] - P[2])
    ax.plot([a[0], b[0]], [a[1], b[1]], color=C2, lw=1.4, ls='--')
    ax.text(v[0] - 0.1, v[1] + 0.1, 'corner cut', color=C2, fontsize=6.5, ha='center', va='bottom')
    ax.text(cen[0], -1.1, 'shortest wall\nthrough centroid', color=C3, fontsize=6.5, ha='center', va='top')
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.75, 1.3); ax.set_aspect('equal'); ax.axis('off')
    ax.text(-0.02, 1.06, 'c', transform=ax.transAxes, fontweight='bold', fontsize=9)
    # (d) junction relaxation
    ax = fig.add_subplot(gsp[:, 5])
    def junction(x0, th_deg, lab):
        h = np.radians((180 - th_deg) / 2)     # host arms bend downward by h each
        L = 0.62
        for sgn in (-1, 1):
            d = np.array([sgn * np.cos(h), -np.sin(h)])
            ax.plot([x0, x0 + L * d[0]], [0, L * d[1]], color='k', lw=2.4, solid_capstyle='round')
        ax.plot([x0, x0], [0, 0.7], color=C3, lw=1.2, solid_capstyle='round')
        t1 = 180 + np.degrees(h); t2 = 360 - np.degrees(h)
        ax.add_patch(Arc((x0, 0), 0.5, 0.5, theta1=t1, theta2=t2, color=C2, lw=0.9))
        ax.text(x0, -0.42 - 0.25 * np.sin(h) * 0 , rf'$\theta={lab}^\circ$', color=C2, fontsize=7, ha='center', va='top')
    junction(0, 180, 180); junction(2.3, 125, 125)
    ax.add_patch(FancyArrowPatch((0.8, 0.35), (1.5, 0.35), arrowstyle='->', mutation_scale=9, color='k', lw=0.8))
    ax.text(1.15, 0.25, 'new wall\nmatures', fontsize=6.5, ha='center', va='top')
    ax.text(0, 0.78, 'birth', fontsize=6.5, ha='center', va='bottom')
    ax.text(2.3, 0.78, r'$\sim$1 cell cycle', fontsize=6.5, ha='center', va='bottom')
    ax.set_xlim(-0.75, 3.05); ax.set_ylim(-0.95, 1.2); ax.set_aspect('equal'); ax.axis('off')
    ax.text(-0.02, 1.06, 'd', transform=ax.transAxes, fontweight='bold', fontsize=9)
    fig.savefig(f'{F}/fig1.pdf', bbox_inches='tight', pad_inches=0.02); plt.close(fig)
    k = deg
    return dict(T=len(segs), kmax=int(k.max()), kmean=float(k.mean()))


# ------------------------------------------------------------------ meristem wall map
def fig_meristem_map(plant='plant1'):  # -> figs/fig3.pdf
    sys.path.insert(0, f'{B}/images')
    import meristem
    from matplotlib.collections import LineCollection
    out = meristem.run(plant, keep_last=True)
    times = out['times']; tl = times[-1]
    E = out['last_E']; lb = out['last_birth']
    old, new, nb, unk = [], [], [], []
    for u, v, pts, ab in E:
        if ab[0] == 0:
            continue
        b = lb.get(ab, 'unknown')
        seg = np.asarray(pts, float)[:, ::-1]      # (x, y)
        if b is None:
            old.append(seg)
        elif isinstance(b, int):
            new.append(seg); nb.append(b)
        else:
            unk.append(seg)
    lab = out['last_lab']
    fig, ax = plt.subplots(figsize=(3.4, 3.1))
    ax.imshow(lab > 0, cmap='Greys', vmin=0, vmax=6, interpolation='nearest')
    ax.add_collection(LineCollection(old + unk, colors='#9a9a9a', linewidths=0.5))
    lc = LineCollection(new, array=np.array(nb, float), cmap='viridis', linewidths=1.1)
    lc.set_clim(0, tl)
    ax.add_collection(lc)
    cb = fig.colorbar(lc, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label('wall birth time (h)')
    ys, xs = np.nonzero(lab > 0)
    ax.set_xlim(xs.min() - 5, xs.max() + 5); ax.set_ylim(ys.max() + 5, ys.min() - 5)
    ax.set_aspect('equal'); ax.axis('off')
    fig.tight_layout(); fig.savefig(f'{F}/fig3.pdf', dpi=300); plt.close(fig)
    return dict(plant=plant, t_first=times[0], t_last=tl, walls_last=len(old) + len(new) + len(unk),
                walls_dated=len(new), walls_undated=len(old) + len(unk),
                cells_last=int(len(np.unique(lab)) - 1))


if __name__ == '__main__':
    if len(sys.argv) > 1:
        for f in sys.argv[1:]:
            print(f, globals()[f]())
        sys.exit()
    fig_schematic(); fig_theory(); print('fig3', fig_meristem_map()); print('fig4', fig_plant())
    print('fig5', fig_venation()); fig_actin()
