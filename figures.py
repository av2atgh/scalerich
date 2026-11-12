"""Figures for manuscript.tex.  Run: python3 figures.py [1 2 3 4 5]
Simulations are cached in figs/cache_*.npz (git-ignored)."""
import os, sys, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, 'lng'))
sys.path.insert(0, os.path.join(ROOT, 'rg'))
import srcells, surrogate, graphstats as gs
from exp_topology import edges_from_hits
FIG = os.path.join(ROOT, 'figs'); os.makedirs(FIG, exist_ok=True)
RES = os.path.join(ROOT, 'lng', 'results')

C1, C2, C3, C4 = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
INK, MUTED = '#262625', '#8a8980'
plt.rcParams.update({'font.size': 8, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'axes.spines.top': False,
                     'axes.spines.right': False, 'legend.frameon': False,
                     'legend.fontsize': 7, 'lines.linewidth': 1.2, 'lines.markersize': 3.5,
                     'savefig.dpi': 300, 'font.family': 'serif', 'mathtext.fontset': 'cm'})
W1, W2 = 3.4, 7.0   # PRE column / page width (in)


def lab(ax, s):
    ax.text(-0.02, 1.02, s, transform=ax.transAxes, fontweight='bold', va='bottom', ha='right')


def cached(name, fn):
    p = os.path.join(FIG, f'cache_{name}.npz')
    if os.path.exists(p):
        return dict(np.load(p, allow_pickle=True))
    d = fn(); np.savez(p, **d); return d


# ---------------------------------------------------------------- Fig. 1
def fig1():
    o = srcells.run(60, alpha=0.5, lam0=0.0, seed=4)
    T = o['T']; lam = 0.035 * np.arange(1, T + 1) ** -0.5
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.45), gridspec_kw=dict(width_ratios=[1, 1, 1.05]))
    ax = axs[0]
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=INK, lw=1.5)
    segs = [[(e[0], e[1]), (e[2], e[3])] for e in o['ends']]
    ax.add_collection(LineCollection(segs, linewidths=lam / lam[0] * 3.2 + 0.3, colors=INK, capstyle='butt'))
    ax.set_title('SR material, $T=60$', fontsize=8)
    ax = axs[1]
    E = edges_from_hits(o['hit']); A = gs.adjacency(E, T); k = gs.degrees(A)
    ax.add_collection(LineCollection(segs, linewidths=0.5, colors='#c9c8c0'))
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color='#c9c8c0', lw=0.8)
    m = o['mid']
    ax.add_collection(LineCollection([[m[i], m[j]] for i, j in E], linewidths=0.7, colors=C1))
    ax.scatter(m[:, 0], m[:, 1], s=6 + 5 * k, color=C1, edgecolor='white', linewidth=0.6, zorder=3)
    ax.set_title('dual graph: ligament = node', fontsize=8)
    for a in axs[:2]:
        a.set_aspect('equal'); a.set_xlim(-0.02, 1.02); a.set_ylim(-0.02, 1.02); a.axis('off')
    for a, s in zip(axs, 'abc'):
        a.text(0.0, 1.06, f'({s})', transform=a.transAxes, fontweight='bold')
    # (c) the two kinds of split of a quadrilateral cell
    ax = axs[2]; ax.axis('off'); ax.set_xlim(0, 2.3); ax.set_ylim(-0.15, 1.25); ax.set_aspect('equal')
    for x0, chord, name, extra in [(0.0, [(0.0, 0.62), (0.4, 1.0)], 'corner cut\ncloses a triangle', 'tri'),
                                   (1.25, [(0.0, 0.3), (1.0, 0.72)], 'cut across\ncloses a square', 'sq')]:
        sq = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]) + [x0, 0]
        ax.plot(sq[:, 0], sq[:, 1], color=INK, lw=1.4)
        (ax_, ay), (bx, by) = chord
        ax.plot([x0 + ax_, x0 + bx], [ay, by], color=C2, lw=1.8)
        ax.text(x0 + 0.5, -0.12, name, ha='center', va='top', fontsize=7)
        # dual nodes: side midpoints + chord midpoint
        sides = {'b': (x0 + .5, .0), 'r': (x0 + 1.0, .5), 't': (x0 + .5, 1.0), 'l': (x0 + .0, .5)}
        cyc = ['b', 'r', 't', 'l']
        for i in range(4):
            p, q = sides[cyc[i]], sides[cyc[(i + 1) % 4]]
            ax.plot([p[0], q[0]], [p[1], q[1]], color=C1, lw=0.7, ls=':')
        nm = (x0 + (ax_ + bx) / 2, (ay + by) / 2)
        hits = ('l', 't') if extra == 'tri' else ('l', 'r')
        for h in hits:
            ax.plot([nm[0], sides[h][0]], [nm[1], sides[h][1]], color=C2, lw=0.9)
        for s in sides.values():
            ax.scatter(*s, s=16, color=C1, zorder=3, edgecolor='white', linewidth=0.5)
        ax.scatter(*nm, s=16, color=C2, zorder=4, edgecolor='white', linewidth=0.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.04, wspace=0.08)
    fig.savefig(os.path.join(FIG, 'fig1.pdf')); plt.close(fig)


# ---------------------------------------------------------------- Fig. 2 (Cowan / BST)
def fig2():
    import srtree
    def sim():
        r, TA, x, lnh = [], [], [], []
        for s in range(6):
            g = srtree.run(20000, seed=s)
            r.append(g['r'][100:])
            leafA = g['cell_area'][g['leaf']]; TA.append(leafA * len(leafA))
            t = np.arange(1, 20001); h = np.concatenate([g['hostA'], g['hostB']])
            tt = np.concatenate([t, t]); ok = h >= 0
            x.append(np.log(tt[ok] / (h[ok] + 1)))
            lnh.append(tt[ok] > 2000)
        return dict(r=np.concatenate(r), TA=np.concatenate(TA),
                    x=np.concatenate(x), late=np.concatenate(lnh))
    d = cached('bst', sim)
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.1))
    ax = axs[0]
    ax.hist(d['r'], bins=20, range=(0, 1), density=True, color=C1, edgecolor='white', linewidth=0.8)
    ax.axhline(1, color=INK, lw=1, ls='--'); ax.set_ylim(0, 1.3)
    ax.set_xlabel('area fraction $r$ of a daughter cell'); ax.set_ylabel('density')
    ax = axs[1]
    xg = np.arange(0, 11, 0.5); cc = np.array([(d['TA'] > x).mean() for x in xg])
    ax.semilogy(xg[cc > 0], cc[cc > 0], 'o', color=C1, mfc='none', label='SR cells')
    xx = np.linspace(0, 11, 50); ax.semilogy(xx, np.exp(-xx), color=INK, lw=1, label='$e^{-x}$')
    ax.set_xlabel('$x=(T+1)A$'); ax.set_ylabel('$P(TA>x)$'); ax.set_ylim(1e-5, 1.5); ax.legend()
    ax = axs[2]
    xl = d['x'][d['late']]
    ax.hist(xl, bins=40, range=(0, 10), density=True, color=C1, edgecolor='white', linewidth=0.5,
            label='$\\ln(t/h)$, $t>2000$')
    xx = np.linspace(0, 10, 100); ax.plot(xx, 0.5 * np.exp(-xx / 2), color=INK, lw=1, label='$\\frac{1}{2} e^{-x/2}$')
    ax.set_xlabel('$\\ln(t/h)$, $h$ = host birth index'); ax.set_ylabel('density'); ax.legend()
    for a, s in zip(axs, 'abc'):
        lab(a, f'({s})')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig2.pdf')); plt.close(fig)


# ---------------------------------------------------------------- Fig. 3 (degrees, clustering)
def sim_main():
    T = 100000; ks, Ls, Cs, kn = [], [], [], []
    for s in (1, 2, 3):
        o = srcells.run(T, seed=s)
        A = gs.adjacency(edges_from_hits(o['hit']), T); k = gs.degrees(A)
        ks.append(k); Ls.append(o['length']); Cs.append(gs.local_clustering(A, k)[0]); kn.append(gs.knn(A, k))
    return dict(k=np.array(ks), L=np.array(Ls), C=np.array(Cs), knn=np.array(kn))


def fig3():
    d = cached('main', sim_main)
    T = d['k'].shape[1]
    fig, axs = plt.subplots(2, 2, figsize=(W1 * 1.02, 4.6))
    k = d['k'].ravel()
    ax = axs[0, 0]
    K = np.unique(np.logspace(np.log10(2), np.log10(k.max()), 40).astype(int))
    P = np.array([(k >= x).mean() for x in K])
    ax.loglog(K, P, 'o', color=C1, mfc='none', label='SR, $T=10^5$')
    ax.loglog(K, 6 / (K * (K + 1.0)) * (K > 2) + (K == 2), color=INK, lw=1, label='$6/[k(k+1)]$')
    ax.set_xlabel('$k$'); ax.set_ylabel('$P(\\geq k)$'); ax.legend(loc='upper right')
    ax = axs[0, 1]
    i = np.arange(T)
    x = 4 * d['L'] / 2.33 * (np.sqrt(T) - np.sqrt(i + 1)); y = d['k'] - 2
    m = (i > 30) & (i < 0.8 * T)
    xs, ys = x[:, m].ravel(), y[:, m].ravel()
    sub = np.random.default_rng(0).choice(len(xs), 6000, replace=False)
    ok = ys[sub] >= 1
    ax.loglog(xs[sub][ok], ys[sub][ok], '.', color='#9ec5f4', ms=1.5, alpha=0.6)
    xb = np.logspace(0, np.log10(xs.max()), 16)
    mb = np.array([(np.sqrt(xb[j] * xb[j + 1]), ys[(xs >= xb[j]) & (xs < xb[j + 1])].mean())
                   for j in range(15) if ((xs >= xb[j]) & (xs < xb[j + 1])).sum() > 5])
    ax.loglog(mb[:, 0], mb[:, 1], 'o', color=C1, mfc='none', label='binned mean')
    zz = np.logspace(0, 3, 10); ax.loglog(zz, zz, color=INK, lw=1, label='slope 1'); ax.legend(loc='upper left')
    ax.set_xlim(0.8, 1.2e3); ax.set_ylim(0.8, 1.2e3)
    ax.set_xlabel('$4\\ell_i(\\sqrt{T}-\\sqrt{i})/c_L$'); ax.set_ylabel('$k_i-2$')
    ax = axs[1, 0]
    kb = np.unique(np.logspace(np.log10(2), np.log10(k.max()), 18).astype(int))
    C = d['C'].ravel(); kk = np.arange(2, k.max() + 1)
    cb = [(np.exp(np.mean(np.log(kb[j:j + 2]))), C[(k >= kb[j]) & (k < kb[j + 1])].mean())
          for j in range(len(kb) - 1) if ((k >= kb[j]) & (k < kb[j + 1])).sum() > 5]
    cb = np.array(cb)
    ax.loglog(cb[:, 0], cb[:, 1], 'o', color=C1, mfc='none', label='SR')
    ax.loglog(kk, 2 * (0.494 + 0.53 * (kk - 2)) / (kk * (kk - 1.0)), color=INK, lw=1,
              label='$2[p_0+p(k-2)]/[k(k-1)]$')
    ax.set_xlabel('$k$'); ax.set_ylabel('$C(k)$'); ax.legend(loc='lower left')
    ax = axs[1, 1]
    kn = d['knn'].ravel()
    nb = np.array([(np.exp(np.mean(np.log(kb[j:j + 2]))), kn[(k >= kb[j]) & (k < kb[j + 1])].mean())
                   for j in range(len(kb) - 1) if ((k >= kb[j]) & (k < kb[j + 1])).sum() > 5])
    ax.semilogx(nb[:, 0], nb[:, 1], 'o-', color=C1, mfc='none', lw=0.8)
    ax.axhline(np.log(T), color=INK, lw=1, ls='--'); ax.text(60, np.log(T) + 0.4, '$\\ln T$', fontsize=7)
    ax.set_xlabel('$k$'); ax.set_ylabel('$k_{nn}(k)$')
    for a, s in zip(axs.ravel(), 'abcd'):
        lab(a, f'({s})')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig3.pdf')); plt.close(fig)


# ---------------------------------------------------------------- Fig. 4 (what moves gamma)
def fig4():
    import csv
    g = list(csv.DictReader(open(os.path.join(RES, 'geometry.csv'))))
    pol = list(csv.DictReader(open(os.path.join(RES, 'polar.csv'))))
    names = {'SR0': 'SR, zero thickness', 'a1.0_l0.1': '$\\alpha=1,\\lambda_0=0.1$',
             'a0.5_l0.03': '$\\alpha=0.5,\\lambda_0=0.03$', 'a0.25_l0.01': '$\\alpha=0.25$ (jammed)',
             'a0_l0.02': '$\\alpha=0$ (jammed, $T\\approx2000$)', 'set_0_90': 'angles $\\{0,\\pi/2\\}$',
             'set_0_60_120': 'angles $\\{0,\\pi/3,2\\pi/3\\}$', 'vm_k16': 'aligned, $\\kappa=16$',
             'q0': 'cell $\\propto$ area$^0$', 'q2': 'cell $\\propto$ area$^2$',
             'grad_x': 'density gradient', 'field_circ': 'polarized (circular)'}
    rows = []
    for key, lbl in names.items():
        v = [float(r['gmle10']) for r in g if r['name'] == key]
        c = [float(r['Cmean']) for r in g if r['name'] == key]
        if v:
            rows.append((lbl, np.mean(v), np.std(v), np.mean(c)))
    rows.append(('3D plates', 2.495, 0.015, 0.72))
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.6), gridspec_kw=dict(width_ratios=[1.35, 1]))
    ax = axs[0]
    y = np.arange(len(rows))[::-1]
    col = [C2 if ('polar' in r[0] or '3D' in r[0]) else C1 for r in rows]
    ax.errorbar([r[1] for r in rows], y, xerr=[r[2] for r in rows], fmt='none', ecolor=MUTED, lw=0.8)
    ax.scatter([r[1] for r in rows], y, color=col, s=18, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax.axvline(3, color=INK, lw=0.8, ls='--'); ax.axvline(2.5, color=MUTED, lw=0.8, ls=':')
    for yy, r in zip(y, rows):
        ax.text(3.08, yy, f'{r[3]:.2f}', va='center', fontsize=6.5, color=INK)
    ax.text(3.08, len(rows) - 0.2, '$C$', fontsize=7)
    ax.set_xlim(2.3, 3.2); ax.set_xlabel('$\\gamma_{\\rm MLE}$ ($k\\geq10$)')
    ax = axs[1]
    for field, cc, mk in (('circ', C2, 'o'), ('wave', C3, 's')):
        pf = sorted({float(r['pfield']) for r in pol if r['field'] == field and r['T'] == '30000'})
        gv = [np.mean([float(r['gmle10']) for r in pol if r['field'] == field and r['T'] == '30000'
                       and float(r['pfield']) == p]) for p in pf]
        pf = [0.0] + pf; gv = [2.90] + gv
        ax.plot(pf, gv, mk + '-', color=cc, label={'circ': 'circular field', 'wave': 'wave field'}[field])
    ax.axhline(3, color=INK, lw=0.8, ls='--')
    ax.set_xlabel('fraction of ligaments following the field'); ax.set_ylabel('$\\gamma_{\\rm MLE}$')
    ax.legend(loc='lower left')
    for a, s in zip(axs, 'ab'):
        lab(a, f'({s})')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig4.pdf')); plt.close(fig)


# ---------------------------------------------------------------- Fig. 5 (paths)
def fig5():
    import csv
    rows = list(csv.DictReader(open(os.path.join(RES, 'paths2_scaling.csv'))))
    old = list(csv.DictReader(open(os.path.join(RES, 'paths_scaling.csv'))))
    fig, axs = plt.subplots(1, 2, figsize=(W1 * 1.02, 2.1))
    ax = axs[0]
    for model, cc, lbl, src in (('SR', C1, 'SR', rows), ('SR*', C4, 'SR rewired', rows),
                                ('DM', C3, 'DM', rows), ('U-padj', C2, 'face surrogate', old)):
        Ts = sorted({int(r['T']) for r in src if r['model'] == model})
        v = [np.array([float(r['meand']) for r in src if r['model'] == model and int(r['T']) == t]) for t in Ts]
        ax.errorbar(Ts, [x.mean() for x in v], yerr=[x.std(ddof=1) / np.sqrt(len(x)) for x in v],
                    fmt='o-', color=cc, mfc='none', label=lbl, capsize=1.5, lw=1)
    ax.set_xscale('log'); ax.set_xlabel('$T$'); ax.set_ylabel('$\\langle d\\rangle$'); ax.legend(fontsize=6)
    ax = axs[1]
    c = list(csv.DictReader(open(os.path.join(RES, 'paths2_control.csv'))))
    for fam, cc, mk, lbl in (('SR', C1, 'o', 'SR variants'), ('surr', C2, 's', 'surrogate, varied $p_4$')):
        names = [n for n in dict.fromkeys(r['variant'] for r in c) if n.startswith(fam)]
        for n in names:
            v = [r for r in c if r['variant'] == n]
            tri = np.mean([float(r['tri']) for r in v]); d = np.array([float(r['delta']) for r in v])
            ax.errorbar(tri, d.mean(), yerr=d.std(ddof=1) / np.sqrt(len(d)), fmt=mk, color=cc, mfc='none',
                        capsize=1.5, label=lbl if n == names[0] else None)
            if n == 'SR':
                ax.annotate('SR', (tri, d.mean()), xytext=(-12, -9), textcoords='offset points', fontsize=6)
            if n == 'SR_q0':
                ax.annotate('$q=0$', (tri, d.mean()), xytext=(4, 2), textcoords='offset points', fontsize=6)
    ax.set_xlabel('triangle-closing fraction'); ax.set_ylabel('$\\Delta\\langle\\ln\\mu\\rangle$')
    ax.set_ylim(0, 4.2); ax.legend(fontsize=6, loc='upper right', handletextpad=0.2)
    for a, s in zip(axs, 'ab'):
        lab(a, f'({s})')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig5.pdf')); plt.close(fig)


if __name__ == '__main__':
    which = sys.argv[1:] or ['1', '2', '3', '4', '5']
    for w in which:
        globals()['fig' + w](); print('fig', w, 'done', flush=True)
