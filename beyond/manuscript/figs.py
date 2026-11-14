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
MK = dict(zip(RULES, ['o', 'v', 's', '^', 'D', 'P']))
REF = {'square': ('#555555', 'x'), 'triangular': ('#555555', '+'), 'honeycomb': ('#555555', '*')}
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
    """two panels stacked vertically, legend outside on the right"""
    fig, (a, b) = plt.subplots(2, 1, figsize=(3.4, 3.7))
    fig.subplots_adjust(hspace=0.42, right=0.70)
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
    f2 = os.path.join(B, 'results', 'bloch2_T1000.json')
    bl = json.load(open(f2)) if os.path.exists(f2) else json.load(open(os.path.join(B, 'results', 'bloch_T1000.json')))
    fig, a, b = _stack()
    for r in RULES:
        v = [x for x in p if x['rule'] == r and x['alpha'] == 0.0]
        a.scatter([x['rho'] for x in v], [x['n_x'] for x in v], s=5, marker=MK[r], color=COL[r], alpha=0.3, linewidths=0)
        xs, ys = _agg(p, r, 'n_x')
        rhos = [np.mean([x['rho'] for x in v if x['rho_t'] == q]) for q in xs]
        a.plot(rhos, ys, color=COL[r], marker=MK[r], ms=3.5, label=LABEL[r])
    for name in ('square', 'triangular', 'honeycomb'):
        v = [x for x in p if x['rule'] == name]
        c, m = REF[name]
        a.plot([x['rho'] for x in v], [x['n_x'] for x in v], color=c, marker=m, ms=4, lw=0.8, label=name)
    rv = sorted({x['rho_t'] for x in vor})
    a.plot([np.mean([x['rho'] for x in vor if x['rho_t'] == q]) for q in rv],
           [np.mean([x['n_x_per'] for x in vor if x['rho_t'] == q]) for q in rv],
           color='#555555', marker='o', mfc='white', ms=3.5, lw=0.8, label='Voronoi')
    a.set_xlabel(r'$\bar\rho$'); a.set_ylabel(r'$n_x$')
    _logticks(a, 'y', [1, 2, 5, 10])
    for r in RULES:
        v = [x for x in bl if x['rule'] == r and x.get('alpha', 0.0) == 0.0 and x['omega_h'] > 0
             and x.get('resolved', x['k_h'] > 0.3)]
        if not v:
            continue
        iso = np.array([x['eta_n'] >= 0.7 for x in v])
        n0 = np.array([x['n0'] for x in v]); wh = np.array([x['omega_h'] for x in v])
        b.scatter(n0[iso], wh[iso], s=10, marker=MK[r], color=COL[r], linewidths=0)
        b.scatter(n0[~iso], wh[~iso], s=10, marker=MK[r], facecolors='none', edgecolors=COL[r], linewidths=0.6)
    b.set_xlabel(r'$n$ (propagation direction)'); b.set_ylabel(r'$\omega_h$')
    _logticks(b, 'x', [2, 3, 4, 6, 8]); b.set_yscale('log')
    b.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    _legend(fig, a, b)
    _panel_labels(fig, (a, b), ('(a)', '(b)'))
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
            a.plot(rhos, ys, color=COL[r], marker=MK[r], ms=3.5, ls=ls,
                   mfc=COL[r] if filled else 'white', label=LABEL[r] if filled else None)
            v = [x for x in rows if x['rule'] == r]
            b.scatter([x['n_mean'] for x in v], [x['omega1'] for x in v], s=10, marker=MK[r],
                      facecolors=COL[r] if filled else 'none', edgecolors=COL[r], linewidths=0.6)
    a.plot([], [], color='#555555', ls='-', label='plates'); a.plot([], [], color='#555555', ls='--', mfc='white', label='struts')
    a.set_xlabel(r'$\bar\rho$'); a.set_ylabel(r'$\langle n\rangle$'); a.set_xscale('log')
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xticks([0.01, 0.03, 0.1, 0.2]); a.set_xticklabels(['0.01', '0.03', '0.1', '0.2'])
    _logticks(a, 'y', [1, 2, 3, 4])
    b.set_xlabel(r'$\langle n\rangle$'); b.set_ylabel(r'$\omega_1$')
    _logticks(b, 'x', [1.5, 2, 3, 4]); b.set_yscale('log')
    _legend(fig, a, b)
    _panel_labels(fig, (a, b), ('(a)', '(b)'))
    fig.savefig(os.path.join(HERE, 'fig3.pdf'))


if __name__ == '__main__':
    which = sys.argv[1:] or ['1', '2', '3']
    for w in which:
        globals()['fig' + w]()
