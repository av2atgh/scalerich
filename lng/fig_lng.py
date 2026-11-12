"""Figures for the LNG section (results/fig_lng.pdf/.png).
(a) CCDF P(k>=K): SR, uniform-face surrogate, DM vs 6/(K(K+1)) [p_k = 12/(k(k+1)(k+2))]
(b) C(k): SR vs 2[p0+p(k-2)]/(k(k-1))
(c) k_nn(k): SR, surrogate, DM; ln T line
(d) <mu>(d): SR, SR* (rewired), DM, T = 1e5 (fresh runs, 40 sources)
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, surrogate, graphstats as gs
from exp_topology import edges_from_hits

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
C1, C2, C3, C4 = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
T = 100000
mods = {'SR': srcells.run(T, seed=21)['hit'],
        'face surrogate': surrogate.run(T, face='uniform', side='padj', padj=0.29, seed=21)['hit'],
        'DM / BB(n,1)': surrogate.run_dm(T, seed=21)['hit']}
cols = {'SR': C1, 'face surrogate': C2, 'DM / BB(n,1)': C3}
G = {m: gs.adjacency(edges_from_hits(h), T) for m, h in mods.items()}
fig, ax = plt.subplots(2, 2, figsize=(9, 7.5))
K = np.unique(np.logspace(np.log10(2), np.log10(1500), 40).astype(int))
for m, A in G.items():
    k = gs.degrees(A)
    ax[0, 0].loglog(K, [(k >= v).mean() for v in K], 'o', ms=4, color=cols[m], label=m, mfc='none')
ax[0, 0].loglog(K, 6 / (K * (K + 1.0)), '-', color='0.3', lw=1.5, label=r'$6/[K(K+1)]$ (BA $m=2$)')
ax[0, 0].set_xlabel('K'); ax[0, 0].set_ylabel(r'$P(k\geq K)$'); ax[0, 0].legend(frameon=False, fontsize=8)
# C(k)
bins = np.unique(np.logspace(np.log10(2), np.log10(800), 18).astype(int))
for m, A in G.items():
    k = gs.degrees(A); C, _ = gs.local_clustering(A, k)
    xs, ys = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        s = (k >= lo) & (k < hi)
        if s.sum() > 3: xs.append(k[s].mean()); ys.append(C[s].mean())
    ax[0, 1].loglog(xs, ys, 'o', ms=4, color=cols[m], mfc='none', label=m)
kk = np.logspace(np.log10(2), np.log10(800), 50)
ax[0, 1].loglog(kk, 2 * (0.494 + 0.53 * (kk - 2)) / (kk * (kk - 1)), '-', color='0.3', lw=1.5,
                label=r'$2[p_0+p(k-2)]/[k(k-1)]$, $p_0{=}0.49,p{=}0.53$')
ax[0, 1].set_xlabel('k'); ax[0, 1].set_ylabel('C(k)'); ax[0, 1].legend(frameon=False, fontsize=8)
for m, A in G.items():
    k = gs.degrees(A); kn = gs.knn(A, k)
    xs, ys = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        s = (k >= lo) & (k < hi)
        if s.sum() > 3: xs.append(k[s].mean()); ys.append(kn[s].mean())
    ax[1, 0].semilogx(xs, ys, 'o-', ms=4, lw=1, color=cols[m], mfc='none', label=m)
ax[1, 0].axhline(np.log(T), color='0.3', lw=1, ls='--', label=r'$\ln T$')
ax[1, 0].set_xlabel('k'); ax[1, 0].set_ylabel(r'$k_{nn}(k)$'); ax[1, 0].legend(frameon=False, fontsize=8)
A = G['SR']; k = gs.degrees(A)
E2, _ = gs.configuration_rewire(A, k, seed=3, nswap_per_edge=5)
for m, AA, c in [('SR', A, C1), ('SR rewired', gs.adjacency(E2, T), C4), ('DM / BB(n,1)', G['DM / BB(n,1)'], C3)]:
    ps = gs.path_stats(AA, nsrc=40, seed=5)
    d = [p[0] for p in ps['per'] if p[1] > 200]; mu = [p[2] for p in ps['per'] if p[1] > 200]
    ax[1, 1].semilogy(d, mu, 'o-', ms=4, lw=1, color=c, mfc='none', label=f"{m}: <d>={ps['meand']:.1f}, <mu>={ps['mu']:.1f}")
kk2 = gs.degrees(A); c2 = (kk2 * (kk2 - 1.0)).mean() / kk2.mean()
ax[1, 1].axhline(c2 / np.e, color='0.3', lw=1, ls='--', label=r'$c_2/e$ (random baseline)')
ax[1, 1].set_xlabel('distance d'); ax[1, 1].set_ylabel(r'$\langle\mu\rangle(d)$'); ax[1, 1].legend(frameon=False, fontsize=7)
for a, l in zip(ax.ravel(), 'abcd'):
    a.text(-0.13, 1.02, f'({l})', transform=a.transAxes, fontsize=10, va='top')
    for s in ('top', 'right'): a.spines[s].set_visible(False)
fig.suptitle(f'SR dual vs face-local growth rules, T = {T:,}', fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_lng.png'), dpi=150)
fig.savefig(os.path.join(OUT, 'fig_lng.pdf'))
print('ok')
