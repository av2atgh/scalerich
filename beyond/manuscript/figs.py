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


def panels(rows, rules, refs, fname, xlab_rho=r'$\bar\rho$', aff=False):
    fig, (a, b) = plt.subplots(1, 2, figsize=(3.4, 1.75))
    for r in rules:
        v = [x for x in rows if x['rule'] == r]
        a.scatter([x['rho'] for x in v], [x['n_x'] for x in v], s=6, marker=MK[r], color=COL[r],
                  alpha=0.35, linewidths=0)
        xs, ys = _agg(rows, r, 'n_x')
        a.plot(xs, ys, color=COL[r], marker=MK[r], ms=3.5, label=LABEL[r])
        iso = np.array([x['eta_n'] >= 0.7 for x in v])
        nm = np.array([x['n_mean'] for x in v]); w1 = np.array([x['omega1'] for x in v])
        b.scatter(nm[iso], w1[iso], s=10, marker=MK[r], color=COL[r], linewidths=0)
        b.scatter(nm[~iso], w1[~iso], s=10, marker=MK[r], facecolors='none', edgecolors=COL[r], linewidths=0.6)
    for name in refs:
        v = [x for x in rows if x['rule'] == name]
        if not v:
            continue
        c, m = REF[name]
        a.plot([x['rho'] for x in v], [x['n_x'] for x in v], color=c, marker=m, ms=4, lw=0.8, label=name)
        b.scatter([x['n_mean'] for x in v], [x['omega1'] for x in v], s=14, marker=m, color=c, linewidths=0.8)
    if aff:
        for r in rules:
            xs, ys = _agg(rows, r, 'n_aff_mean')
            a.plot(xs, ys, color=COL[r], lw=0.6, ls=':')
    a.set_xlabel(xlab_rho); a.set_ylabel(r'$n_x$'); a.set_yscale('log')
    b.set_xlabel(r'$\langle n\rangle$'); b.set_ylabel(r'$\omega_1$', labelpad=1); b.set_yscale('log'); b.set_xscale('log')
    for ax in (a, b):
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_yticks([1, 2, 5, 10, 20]); a.set_yticklabels(['1', '2', '5', '10', '20'])
    b.set_xticks([2, 3, 5, 10]); b.set_xticklabels(['2', '3', '5', '10'])
    b.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    lab(a, '(a)'); lab(b, '(b)')
    h, l = a.get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.93),
               handlelength=1.4, columnspacing=0.8, handletextpad=0.3)
    fig.subplots_adjust(wspace=0.55)
    fig.savefig(os.path.join(HERE, fname))


def fig2():
    rows = json.load(open(os.path.join(B, 'results', 'sweep_T1000.json')))
    panels(rows, RULES, ['square', 'triangular', 'honeycomb'], 'fig2.pdf')


def fig3():
    rows = json.load(open(os.path.join(B, 'results', 'sweep3d_T300.json')))
    # 'uniform' excluded: plates thicker than wide at matched density (thin-shell model invalid)
    panels(rows, [r for r in RULES if r != 'uniform'], [], 'fig3.pdf', aff=True)


if __name__ == '__main__':
    which = sys.argv[1:] or ['1', '2', '3']
    for w in which:
        globals()['fig' + w]()
