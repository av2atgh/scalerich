"""Figures for manuscript.tex.  python figs.py [1 2 3]"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HERE, '..')
sys.path.insert(0, B); sys.path.insert(0, os.path.join(B, 'd3'))
plt.rcParams.update({'font.size': 9, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'legend.fontsize': 7.5, 'lines.linewidth': 1.2,
                     'lines.markersize': 4, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})
RULES = ['SR', 'uniform', 'rect', 'centroid', 'balanced', 'regular']
COL = dict(zip(RULES, ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']))
MK = {'SR': '*', 'uniform': 'X', 'rect': 's', 'centroid': 'o', 'balanced': 'v', 'regular': 'p'}
REF = {'square': ('#555555', 's'), 'triangular': ('#555555', '^'), 'honeycomb': ('#555555', 'h'), 'voronoi': ('#555555', 'D')}
LABEL = {'SR': 'SR', 'uniform': 'uniform', 'rect': 'rect', 'centroid': 'centroid',
         'balanced': 'balanced', 'regular': 'regular'}


OPEN = ('centroid', 'regular')                       # open markers, same everywhere
from matplotlib.lines import Line2D


def sty(r, **kw):
    """one style per structure, used in every figure (lines with markers)"""
    d = dict(color=COL[r], marker=MK[r], ms=8.5 if MK[r] == '*' else 5.5, mec=COL[r], mew=0.9,
             mfc='white' if r in OPEN else COL[r], ls='-')
    d.update(kw); return d


def sct(r, scale=1.0):
    """the same symbols for scatter plots"""
    return dict(marker=MK[r], s=(46 if MK[r] == '*' else 20) * scale, edgecolors=COL[r], linewidths=0.8,
                facecolors='white' if r in OPEN else COL[r])


def _top_legend(fig, axes_top, rules, extra=(), ncol=3, lift=0.0):
    """the structure key, repeated at the top of every figure"""
    h = [Line2D([], [], **sty(r)) for r in rules] + [e[0] for e in extra]
    l = [LABEL[r] for r in rules] + [e[1] for e in extra]
    y = axes_top.get_position().y1 + 0.035 + lift
    fig.legend(h, l, loc='lower center', bbox_to_anchor=(0.5, y), ncol=ncol, frameon=False,
               handlelength=2.0, columnspacing=1.0, handletextpad=0.4)


def lab(ax, s):
    ax.text(-0.02, 1.02, s, transform=ax.transAxes, fontweight='bold', va='bottom', ha='right')


def fig1():
    import divisions
    fig, axs = plt.subplots(2, 3, figsize=(3.4, 2.4))
    T = 400
    for ax, r in zip(axs.ravel(), RULES):
        g = divisions.run(T, 3, r)
        lw = 2.2 * np.arange(1, T + 1.) ** -0.5 + 0.15
        segs = np.stack([g['A'], g['B']], 1)
        ax.add_collection(LineCollection(segs, linewidths=lw, colors='#222222', capstyle='butt'))
        ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color='#999999', lw=0.4)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal'); ax.axis('off')
        ax.text(0.5, 1.04, LABEL[r], transform=ax.transAxes, ha='center', va='bottom', fontsize=8)
        ax.plot([0.5 - 0.075 * len(LABEL[r]) / 2 - 0.1], [1.10], transform=ax.transAxes, clip_on=False,
                **sty(r, ls='none', ms=7.5 if MK[r] == '*' else 4.5))
    fig.subplots_adjust(wspace=0.08, hspace=0.24)
    fig.savefig(os.path.join(HERE, 'fig1.pdf'))


def _agg(rows, rule, key, by='rho_t', alpha=0.0):
    v = [x for x in rows if x['rule'] == rule and (alpha is None or x['alpha'] == alpha)]
    xs = sorted({x[by] for x in v})
    return np.array(xs), np.array([np.exp(np.mean(np.log([x[key] for x in v if x[by] == q]))) for q in xs])


def _stack():
    """two panels stacked vertically, full column width (symbols explained in the caption)"""
    fig, (a, b) = plt.subplots(2, 1, figsize=(3.4, 4.4))
    fig.subplots_adjust(hspace=0.42)
    return fig, a, b


def _panel_labels(fig, axes, labels):
    """panel labels left-aligned with the y-axis labels, at the top of each axes"""
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    x0 = min(inv.transform(ax.yaxis.label.get_window_extent())[0, 0] for ax in axes)
    for ax, s in zip(axes, labels):
        y1 = inv.transform(ax.get_window_extent())[1, 1]
        fig.text(x0, y1 + 0.012, s, fontweight='bold', ha='left', va='bottom')


def _legend(fig, a, b):
    h, l = a.get_legend_handles_labels()
    h2, l2 = b.get_legend_handles_labels()
    for hh, ll in zip(h2, l2):
        if ll not in l:
            h.append(hh); l.append(ll)
    fig.legend(h, l, loc='center left', bbox_to_anchor=(0.71, 0.5), frameon=False,
               handlelength=1.6, handletextpad=0.4, labelspacing=0.5)


def _logticks(ax, axis, ticks):
    getattr(ax, f'set_{axis}scale')('log')
    getattr(ax, f'set_{axis}ticks')(ticks)
    getattr(ax, f'set_{axis}ticklabels')([('%g' % t) for t in ticks])
    getattr(ax, f'{axis}axis').set_minor_formatter(matplotlib.ticker.NullFormatter())


def fig3():
    p = json.load(open(os.path.join(B, 'results', 'psweep_T1000.json')))
    vor = json.load(open(os.path.join(B, 'results', 'voronoi.json')))
    bl = json.load(open(os.path.join(B, 'results', 'bloch2_T1000.json')))
    fr = os.path.join(B, 'results', 'bloch2_rect.json')
    if os.path.exists(fr):
        bl = [x for x in bl if x['rule'] != 'rect'] + json.load(open(fr))
    sp_ = json.load(open(os.path.join(B, 'results', 'spanmode_bloch2.json')))
    fig, (a, b, c) = plt.subplots(3, 1, figsize=(3.4, 6.8))
    fig.subplots_adjust(hspace=0.45)
    MAIN = ('SR', 'uniform', 'rect', 'centroid', 'balanced', 'regular')
    # (a) index vs density
    for r in MAIN:
        v = [x for x in p if x['rule'] == r and x['alpha'] == 0.0]
        qs = sorted({x['rho_t'] for x in v})
        rh = np.array([np.mean([x['rho'] for x in v if x['rho_t'] == q]) for q in qs])
        nn = [np.array([x['n_mean'] for x in v if x['rho_t'] == q]) for q in qs]
        gm = np.array([np.exp(np.mean(np.log(z))) for z in nn])
        if r in ('SR', 'balanced'):
            a.fill_between(rh, [z.min() for z in nn], [z.max() for z in nn], color=COL[r], alpha=0.15, lw=0)
        a.plot(rh, gm, **sty(r))
    for name in ('honeycomb', 'triangular', 'square'):
        v = [x for x in p if x['rule'] == name]
        cc, m = REF[name]
        a.plot([x['rho'] for x in v], [x['n_mean'] for x in v], color=cc, marker=m, ms=5, lw=0.7, label=name)
    rv = sorted({x['rho_t'] for x in vor})
    a.plot([np.mean([x['rho'] for x in vor if x['rho_t'] == q]) for q in rv],
           [np.mean([x['n_mean_per'] for x in vor if x['rho_t'] == q]) for q in rv],
           color='#555555', marker='D', ms=4.5, lw=0.7, label='Voronoi')
    a.set_xlabel(r'$\bar\rho$'); a.set_ylabel(r'$\langle n\rangle$')
    a.set_ylim(1, 8); a.set_yticks([1, 2, 3, 4, 5, 6, 7, 8])
    a.set_xticks([0.1, 0.2, 0.3, 0.4])
    # (b) band edge at fixed index: each (realisation, direction) density curve is interpolated in
    # ln(omega_h) at fixed n; mean and standard error over the curves that reach that n (>= 3)
    nall = []
    curves = {}
    for x in bl:
        if x['resolved'] and x['omega_h'] > 0:
            curves.setdefault((x['rule'], x['seed'], x['dir']), []).append((x['n0'], x['omega_h']))
    grid = np.arange(2.0, 8.01, 0.25)
    for r in ('SR', 'rect', 'centroid', 'balanced', 'regular'):
        M = []; E = []; G = []
        for n in grid:
            v = []
            for (rr, _, _), cv in curves.items():
                if rr != r:
                    continue
                cv = sorted(cv); xx = np.array([q[0] for q in cv]); yy = np.log([q[1] for q in cv])
                if xx.min() <= n <= xx.max():
                    v.append(np.interp(n, xx, yy))
            if len(v) >= 3:
                v = np.array(v); G.append(n); M.append(v.mean()); E.append(v.std() / np.sqrt(len(v)))
        G = np.array(G); M = np.array(M); E = np.array(E)
        nall.extend(G)
        if len(G) == 0:
            continue
        b.fill_between(G, np.exp(M - E), np.exp(M + E), color=COL[r], alpha=0.18, lw=0)
        b.plot(G, np.exp(M), markevery=2, **sty(r))
    b.set_xlabel(r'$n$ (propagation direction)'); b.set_ylabel(r'$\omega_h$')
    b.set_yscale('log'); b.set_xlim(min(nall) - 0.15, max(nall) + 0.15)
    b.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
    # (c) band edge vs weakest-span prediction
    for r in ('SR', 'centroid', 'balanced'):
        v = [x for x in sp_ if x['rule'] == r and x['resolved'] and x['omega_h'] > 0]
        c.scatter([x['omega_span'] for x in v], [x['omega_h'] for x in v], **sct(r, 0.8))
    allx = [x['omega_span'] for x in sp_ if x['resolved'] and x['omega_h'] > 0]
    xx = np.array([0.7 * min(allx), 1.4 * max(allx)])
    c.plot(xx, xx, color='#999999', lw=0.6, ls=':')
    c.plot(xx, (np.pi / 4.730) ** 2 * xx, color='#999999', lw=0.6, ls=':')
    c.set_xlim(*xx)
    c.set_xscale('log'); c.set_yscale('log')
    c.set_xlabel(r'$\omega_{\rm span}$ (weakest span)'); c.set_ylabel(r'$\omega_h$')
    _panel_labels(fig, (a, b, c), ('(a)', '(b)', '(c)'))
    lat = [(Line2D([], [], color='#555555', marker=REF[k][1], ms=4.5, lw=0.7), k) for k in ('square', 'triangular', 'honeycomb')]
    lat.append((Line2D([], [], color='#555555', marker='D', ms=4, lw=0.7), 'Voronoi'))
    y = a.get_position().y1 + 0.035
    fig.legend([e[0] for e in lat], [e[1] for e in lat], loc='lower center', bbox_to_anchor=(0.5, y), ncol=4,
               frameon=False, handlelength=1.6, columnspacing=0.8, handletextpad=0.3)
    _top_legend(fig, a, MAIN, lift=0.03)
    fig.savefig(os.path.join(HERE, 'fig3.pdf'))


def fig5():
    pl = json.load(open(os.path.join(B, 'results', 'sweep3d_T300.json')))
    st = json.load(open(os.path.join(B, 'results', 'skel3d_T300.json')))
    rules = [r for r in RULES if r != 'uniform']      # plates thicker than wide (thin-shell model invalid)
    fig, a, b = _stack()
    for r in rules:
        for rows, ls in ((pl, '-'), (st, ':')):
            xs, ys = _agg(rows, r, 'n_mean')
            rhos = [np.mean([x['rho'] for x in rows if x['rule'] == r and x['alpha'] == 0 and x['rho_t'] == q]) for q in xs]
            a.plot(rhos, ys, **sty(r, ls=ls))
            v = [x for x in rows if x['rule'] == r]
            b.scatter([x['n_mean'] for x in v], [x['omega1'] for x in v], **sct(r, 0.8))
    a.set_xlabel(r'$\bar\rho$'); a.set_ylabel(r'$\langle n\rangle$'); a.set_xscale('log')
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xticks([0.01, 0.03, 0.1, 0.2]); a.set_xticklabels(['0.01', '0.03', '0.1', '0.2'])
    a.set_ylim(1, 4); a.set_yticks([1, 2, 3, 4])
    b.set_xlabel(r'$\langle n\rangle$'); b.set_ylabel(r'$\omega_1$')
    b.set_xticks([1.5, 2, 2.5, 3, 3.5, 4]); b.set_yscale('log')
    _panel_labels(fig, (a, b), ('(a)', '(b)'))
    _top_legend(fig, a, rules, extra=[(Line2D([], [], color='#555555', ls='-'), 'plates'),
                                      (Line2D([], [], color='#555555', ls=':'), 'struts')], ncol=4)
    fig.savefig(os.path.join(HERE, 'fig5.pdf'))


def fig_shock_unused():
    """shock routing in a lattice insert (homogenised, cut-angle grading at rho = 0.2)"""
    d = json.load(open(os.path.join(B, 'results', 'shockfig_tau20.json')))
    sU = np.load(os.path.join(B, 'results', 'shockfig_U.npz')); sD = np.load(os.path.join(B, 'results', 'shockfig_D5.npz'))
    fig = plt.figure(figsize=(3.4, 4.9))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.62, 1.0, 1.0], hspace=0.5, wspace=0.08)
    a1 = fig.add_subplot(gs[0, 0]); a2 = fig.add_subplot(gs[0, 1])
    b = fig.add_subplot(gs[1, :]); c = fig.add_subplot(gs[2, :])
    Wd, Hd = 60.0, 30.0
    key = 's100'
    vmax = max(sU[key].max(), sD[key].max())
    for ax, S, ttl in ((a1, sU, 'uniform'), (a2, sD, 'graded')):
        ax.imshow(S[key].T, origin='lower', extent=(-Wd / 2, Wd / 2, 0, Hd), cmap='magma', vmin=0, vmax=vmax, aspect='equal')
        ax.plot([-4, 4], [0.3, 0.3], color='#2a78d6', lw=3)               # component
        ax.annotate('', xy=(0, Hd - 0.5), xytext=(0, Hd + 6), arrowprops=dict(arrowstyle='->', color='k', lw=0.8), annotation_clip=False)
        if ttl == 'graded':
            for xx in (-5, 5):
                ax.axvline(xx, color='white', lw=0.5, ls=':')
        ax.set_title(ttl, fontsize=8, pad=8); ax.set_xticks([]); ax.set_yticks([])
    t = np.array(d['U']['t'])
    for name, col, lab_, ls in (('U', '#555555', 'uniform', '-'), ('D5', COL['centroid'], 'diverging', '-'), ('C5', COL['balanced'], 'converging', '--')):
        b.plot(t, np.array(d[name]['F']), color=col, ls=ls, lw=1.2, label=lab_)
    b.set_xlim(0, 300); b.set_xlabel(r'$t$ ($L/c_0$)'); b.set_ylabel('force on\ncomponent')
    b.legend(frameon=False, loc='lower right', handlelength=2)
    sc = d['scan']
    for dsg, col, mk in (('D', COL['centroid'], 'o'), ('C', COL['balanced'], 'v')):
        v = [x for x in sc if x['design'] == dsg]
        X = [x['X'] for x in v]
        c.plot(X, [x['peak'] for x in v], color=col, marker=mk, ms=5, ls='-', mfc=col)
        c.plot(X, [x['impulse'] for x in v], color=col, marker=mk, ms=5, ls='--', mfc='white')
    c.axhline(1, color='#999999', lw=0.6, ls=':')
    c.set_xlabel(r'grading width $X$ (cells)'); c.set_ylabel('relative to uniform')
    _panel_labels(fig, (a1, b, c), ('(a)', '(b)', '(c)'))
    fig.savefig(os.path.join(HERE, 'fig4.pdf'))


def _style(r):
    return dict(color=COL[r], marker=MK[r], ms=8 if MK[r] == '*' else 5.5)


def fig2():
    """network mappings: split fraction, dual degree law, cell elongation"""
    A = json.load(open(os.path.join(B, 'results', 'net_asym_T20000.json')))
    F = json.load(open(os.path.join(B, 'results', 'net_fig.json')))['dist']
    fig, (a, b, c) = plt.subplots(3, 1, figsize=(3.4, 6.8))
    fig.subplots_adjust(hspace=0.45)
    edges = np.linspace(0, 0.5, 21); mid = 0.5 * (edges[1:] + edges[:-1])
    for r in ('SR', 'centroid', 'balanced', 'regular'):
        h = np.sum([x['frac_hist'] for x in A if x['rule'] == r], axis=0).astype(float)
        h /= h.sum() * (edges[1] - edges[0])
        a.plot(mid, h, **sty(r))
    a.axvline(4 / 9, color='#999999', lw=0.6, ls=':')
    a.set_xlabel(r'smaller daughter area fraction $f$'); a.set_ylabel(r'$p(f)$'); a.set_yscale('log')
    a.set_xlim(0, 0.5); a.set_ylim(0.5, 200)
    Ks = np.array([3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128])
    for r in RULES:
        cc = np.mean([[x['ccdf'][str(K)] for K in Ks] for x in A if x['rule'] == r], axis=0)
        m = cc * 4 * 20000 >= 20
        b.plot(Ks[m], (Ks ** 2 * cc)[m], **sty(r))
    KK = np.geomspace(3, 128, 50)
    b.axhline(6, color='#555555', lw=0.8, ls=':')
    b.set_xscale('log'); b.set_xlabel(r'degree $K$'); b.set_ylabel(r'$K^2P(k\geq K)$'); b.set_ylim(0, 9)
    _logticks(b, 'x', [3, 10, 30, 100])
    for r in RULES:
        d = F[r]; g = np.array(d['elong_grid']); y = np.array(d['elong_ccdf']); m = y > 0
        c.plot(g[m], y[m], markevery=6, **sty(r))
    c.set_xscale('log'); c.set_yscale('log'); c.set_xlim(1, 1e4); c.set_ylim(1e-5, 1.5)
    c.set_xlabel(r'cell elongation $e$'); c.set_ylabel(r'$P(e\geq x)$')
    _panel_labels(fig, (a, b, c), ('(a)', '(b)', '(c)'))
    _top_legend(fig, a, RULES)
    fig.savefig(os.path.join(HERE, 'fig2.pdf'))


def fig4():
    """from network to physics: span tail, index collapse, cut-angle family"""
    F = json.load(open(os.path.join(B, 'results', 'net_fig.json')))
    Cl = json.load(open(os.path.join(B, 'results', 'collapse.json')))
    fig, (a, b, c) = plt.subplots(3, 1, figsize=(3.4, 6.8))
    fig.subplots_adjust(hspace=0.45)
    for r in RULES:
        d = F['dist'][r]; g = np.array(d['span_grid']); y = np.array(d['span_ccdf']); m = (y > 0) & (g >= 0.3)
        a.plot(g[m], y[m], markevery=8, **sty(r))
    a.set_xscale('log'); a.set_yscale('log'); a.set_ylim(1e-5, 1.5); a.set_xlim(0.3, 300)
    a.set_xlabel(r'span length $\ell\sqrt{T}$'); a.set_ylabel(r'$P(\ell\sqrt{T}\geq x)$')
    kcol = plt.get_cmap('Greys')
    for p in Cl['pts']:
        r = p['rule']
        if r == 'uniform':
            continue
        if r.startswith('k'):
            kv = float(r[1:]); b.scatter(p['x'], p['n'] ** -2, s=9, marker='.', color=kcol(0.95 - 0.6 * kv), linewidths=0)
        else:
            b.scatter(p['x'], p['n'] ** -2, **sct(r, 0.8))
    xx = np.geomspace(min(p['x'] for p in Cl['pts']) * 0.8, max(p['x'] for p in Cl['pts']) * 1.2, 100)
    b.plot(xx, Cl['s'] + Cl['c'] * xx ** Cl['q'], color='#333333', lw=0.8)
    b.axhline(Cl['s'], color='#999999', lw=0.6, ls=':')
    b.set_xscale('log'); b.set_yscale('log')
    b.set_xlabel(r'$\lambda/\ell_*$'); b.set_ylabel(r'$\langle n\rangle^{-2}$')
    K = F['kappa']; ks = sorted({x['kappa'] for x in K})
    mean = lambda f: np.array([np.mean([x[f] for x in K if x['kappa'] == k]) for k in ks])
    c.plot(ks, mean('n') / mean('n')[0], color='#222222', marker='^', ms=5, label=r'$\langle n\rangle$')
    c.plot(ks, mean('l_star') / mean('l_star')[0], color='#222222', marker='^', ms=5, ls='--', mfc='white', label=r'$\ell_*$')
    c.plot(ks, mean('eta'), color='#777777', marker='>', ms=5, label=r'$\eta_n$')
    c.plot(ks, mean('eta_aff'), color='#777777', marker='>', ms=5, ls='--', mfc='white', label=r'$\eta_{\rm aff}$')
    c.plot(ks, mean('fabric4'), color='#aaaaaa', marker='d', ms=5, ls=':', label=r'$a_4$')
    c.set_xlabel(r'cut-angle spread $\kappa$'); c.set_ylabel('value')
    c.legend(loc='upper left', frameon=False, ncol=3, handlelength=2.0, borderaxespad=0.2, columnspacing=0.8)
    c.set_ylim(0, 2.3)
    _panel_labels(fig, (a, b, c), ('(a)', '(b)', '(c)'))
    _top_legend(fig, a, RULES)
    fig.savefig(os.path.join(HERE, 'fig4.pdf'))


if __name__ == '__main__':
    which = sys.argv[1:] or ['1', '2', '3', '4', '5']
    for w in which:
        globals()['fig' + w]()
