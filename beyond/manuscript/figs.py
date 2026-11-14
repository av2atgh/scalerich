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
        ax.set_title(LABEL[r], fontsize=8, pad=2)
    fig.subplots_adjust(wspace=0.08, hspace=0.18)
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


def fig2():
    p = json.load(open(os.path.join(B, 'results', 'psweep_T1000.json')))
    vor = json.load(open(os.path.join(B, 'results', 'voronoi.json')))
    bl = json.load(open(os.path.join(B, 'results', 'bloch2_T1000.json')))
    sp_ = json.load(open(os.path.join(B, 'results', 'spanmode_bloch2.json')))
    fig, (a, b, c) = plt.subplots(3, 1, figsize=(3.4, 6.8))
    fig.subplots_adjust(hspace=0.45)
    MAIN = ('SR', 'centroid', 'balanced', 'regular')
    # (a) index vs density
    for r in MAIN:
        v = [x for x in p if x['rule'] == r and x['alpha'] == 0.0]
        qs = sorted({x['rho_t'] for x in v})
        rh = np.array([np.mean([x['rho'] for x in v if x['rho_t'] == q]) for q in qs])
        nn = [np.array([x['n_mean'] for x in v if x['rho_t'] == q]) for q in qs]
        gm = np.array([np.exp(np.mean(np.log(z))) for z in nn])
        dashed = r in ('centroid', 'regular')       # coincide with SR / balanced: same cut angles
        if not dashed:
            a.fill_between(rh, [z.min() for z in nn], [z.max() for z in nn], color=COL[r], alpha=0.15, lw=0)
        a.plot(rh, gm, color=COL[r], marker=MK[r], ms=8 if MK[r] == '*' else 5.5, ls='--' if dashed else '-',
               mfc='white' if dashed else COL[r], label=LABEL[r])
    for name in ('honeycomb', 'triangular', 'square'):
        v = [x for x in p if x['rule'] == name]
        cc, m = REF[name]
        a.plot([x['rho'] for x in v], [x['n_mean'] for x in v], color=cc, marker=m, ms=5, lw=0.7, label=name)
    rv = sorted({x['rho_t'] for x in vor})
    a.plot([np.mean([x['rho'] for x in vor if x['rho_t'] == q]) for q in rv],
           [np.mean([x['n_mean_per'] for x in vor if x['rho_t'] == q]) for q in rv],
           color='#555555', marker='D', ms=4.5, lw=0.7, label='Voronoi')
    a.set_xlabel(r'$\bar\rho$'); a.set_ylabel(r'$\langle n\rangle$')
    a.set_ylim(1, 6.5); a.set_yticks([1, 2, 3, 4, 5, 6])
    a.set_xticks([0.1, 0.2, 0.3, 0.4])
    # (b) band edge at fixed index: each (realisation, direction) density curve is interpolated in
    # ln(omega_h) at fixed n; mean and standard error over the curves that reach that n (>= 3)
    nall = []
    curves = {}
    for x in bl:
        if x['resolved'] and x['omega_h'] > 0:
            curves.setdefault((x['rule'], x['seed'], x['dir']), []).append((x['n0'], x['omega_h']))
    grid = np.arange(2.0, 8.01, 0.25)
    for r in ('SR', 'centroid', 'balanced'):
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
        dashed = r == 'centroid'
        b.fill_between(G, np.exp(M - E), np.exp(M + E), color=COL[r], alpha=0.18, lw=0)
        b.plot(G, np.exp(M), color=COL[r], marker=MK[r], ms=8 if MK[r] == '*' else 5.5, ls='--' if dashed else '-',
               mfc='white' if dashed else COL[r], label=LABEL[r], markevery=2)
    b.set_xlabel(r'$n$ (propagation direction)'); b.set_ylabel(r'$\omega_h$')
    b.legend(loc='upper right', frameon=False, handlelength=2.2, borderaxespad=0.2)
    b.set_yscale('log'); b.set_xlim(min(nall) - 0.15, max(nall) + 0.15)
    b.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
    # (c) band edge vs weakest-span prediction
    for r in ('SR', 'centroid', 'balanced'):
        v = [x for x in sp_ if x['rule'] == r and x['resolved'] and x['omega_h'] > 0]
        c.scatter([x['omega_span'] for x in v], [x['omega_h'] for x in v], s=40 if MK[r] == '*' else 18, marker=MK[r],
                  facecolors='none' if r == 'centroid' else COL[r], edgecolors=COL[r], linewidths=0.5)
    allx = [x['omega_span'] for x in sp_ if x['resolved'] and x['omega_h'] > 0]
    xx = np.array([0.7 * min(allx), 1.4 * max(allx)])
    c.plot(xx, xx, color='#999999', lw=0.6, ls=':')
    c.plot(xx, (np.pi / 4.730) ** 2 * xx, color='#999999', lw=0.6, ls=':')
    c.set_xlim(*xx)
    c.set_xscale('log'); c.set_yscale('log')
    c.set_xlabel(r'$\omega_{\rm span}$ (weakest span)'); c.set_ylabel(r'$\omega_h$')
    _panel_labels(fig, (a, b, c), ('(a)', '(b)', '(c)'))
    fig.savefig(os.path.join(HERE, 'fig2.pdf'))


def fig3():
    pl = json.load(open(os.path.join(B, 'results', 'sweep3d_T300.json')))
    st = json.load(open(os.path.join(B, 'results', 'skel3d_T300.json')))
    rules = [r for r in RULES if r != 'uniform']      # plates thicker than wide (thin-shell model invalid)
    fig, a, b = _stack()
    for r in rules:
        for rows, ls, filled in ((pl, '-', True), (st, '--', False)):
            xs, ys = _agg(rows, r, 'n_mean')
            rhos = [np.mean([x['rho'] for x in rows if x['rule'] == r and x['alpha'] == 0 and x['rho_t'] == q]) for q in xs]
            a.plot(rhos, ys, color=COL[r], marker=MK[r], ms=8 if MK[r] == '*' else 5.5, ls=ls,
                   mfc=COL[r] if filled else 'white', label=LABEL[r] if filled else None)
            v = [x for x in rows if x['rule'] == r]
            b.scatter([x['n_mean'] for x in v], [x['omega1'] for x in v], s=36 if MK[r] == '*' else 16, marker=MK[r],
                      facecolors=COL[r] if filled else 'none', edgecolors=COL[r], linewidths=0.6)
    a.plot([], [], color='#555555', ls='-', label='plates'); a.plot([], [], color='#555555', ls='--', mfc='white', label='struts')
    a.set_xlabel(r'$\bar\rho$'); a.set_ylabel(r'$\langle n\rangle$'); a.set_xscale('log')
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xticks([0.01, 0.03, 0.1, 0.2]); a.set_xticklabels(['0.01', '0.03', '0.1', '0.2'])
    a.set_ylim(1, 4); a.set_yticks([1, 2, 3, 4])
    b.set_xlabel(r'$\langle n\rangle$'); b.set_ylabel(r'$\omega_1$')
    b.set_xticks([1.5, 2, 2.5, 3, 3.5, 4]); b.set_yscale('log')
    _panel_labels(fig, (a, b), ('(a)', '(b)'))
    fig.savefig(os.path.join(HERE, 'fig3.pdf'))


if __name__ == '__main__':
    which = sys.argv[1:] or ['1', '2', '3']
    for w in which:
        globals()['fig' + w]()
