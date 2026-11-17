"""Fig. S3 of supplement.tex: spin waves on magnetic wire networks (classes/spinwave.py).
(a) participation fraction of the modes, geometric mean over eight frequency windows
(omega = 4..8, 60 modes each, median per window), against the number of cuts T; random rules:
mean of three seeds, bars = range. (b) median participation per window at T = 16000."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9})
d = {}
for s in 'AB': d.update(json.load(open('../results/spinwave_avg%s.json' % s)))
BLUE, RED, GREEN, PURPLE, GREY = '#2c6e9e', '#c0392b', '#1e8449', '#7d3c98', '#8c8c8c'
S = [('periodic', 'periodic', dict(color=BLUE, marker='s')),
     ('golden-aligned', 'golden, aligned', dict(color=RED, marker='o')),
     ('r0.3-aligned', r'$r_0=0.3$, aligned', dict(color=RED, marker='^', ls='--')),
     ('golden-anti', 'golden, anti-aligned', dict(color=GREEN, marker='o', mfc='none')),
     ('uniform', 'uniform split', dict(color=PURPLE, marker='D')),
     ('SR', 'SR', dict(color=GREY, marker='v'))]
fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.8))
for key, lab, st in S:
    rows = [v for v in d.values() if v['rule'] == key]
    Ts = sorted({v['T'] for v in rows})
    G = [[v['gm'] for v in rows if v['T'] == T] for T in Ts]
    m = np.array([np.exp(np.mean(np.log(g))) for g in G])
    lo = m - np.array([min(g) for g in G]); hi = np.array([max(g) for g in G]) - m
    a.errorbar(Ts, m, yerr=[lo, hi], ms=4, lw=1.1, capsize=0, label=lab, **st)
    r16 = [v for v in rows if v['T'] == 16000]
    med = np.exp(np.mean([np.log(v['med']) for v in r16], axis=0))
    b.plot(r16[0]['oms'], med, ms=4, lw=1.1, **st)
x = np.array([1000, 16000]); a.plot(x, 0.0022 * (x / 1000.0) ** -1, 'k:', lw=0.8)
a.text(1150, 0.0008, r'$\propto T^{-1}$', fontsize=7)
a.set_xscale('log'); a.set_yscale('log')
a.set_xlabel('number of cuts $T$'); a.set_ylabel('participation fraction')
fig.legend(*a.get_legend_handles_labels(), frameon=False, fontsize=7, loc='lower center', ncol=6, bbox_to_anchor=(0.5, -0.06))
b.set_yscale('log'); b.set_xlabel(r'frequency $\omega$'); b.set_ylabel('participation fraction')
b.set_title('$T=16000$', fontsize=8)
a.set_title('(a)', loc='left', fontsize=9); b.set_title('(b)', loc='left', fontsize=9)
fig.tight_layout(pad=0.4, rect=(0, 0.06, 1, 1))
fig.savefig('figSW.pdf', bbox_inches='tight'); fig.savefig('figSW.png', dpi=150, bbox_inches='tight')
