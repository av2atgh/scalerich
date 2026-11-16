"""Figures for classes/manuscript/manuscript.tex. Run from classes/manuscript."""
import sys, json, numpy as np
sys.path[:0] = ['..', '../../beyond', '../../rg']
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import qdiv, divisions
from diffraction import RULES
plt.rcParams.update({'font.size': 8, 'axes.linewidth': 0.6, 'lines.linewidth': 1.0,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'legend.frameon': False})
R = '../results/'
COL = {'golden-avoid': '#d35400', '(5,1)': '#16a085', '(5,1)-avoid': '#48c9b0', '(3,2)-avoid': '#bb8fce', 'uniform': '#7f8c8d', 'plastic51-avoid': '#48c9b0', 'regular': 'k', 'golden': '#c0392b', 'supergolden': '#e67e22', 'theta4': '#d4a017', 'plastic': '#8e44ad',
       'r0.3': '#2980b9', 'golden-side': '#27ae60', 'largest-uniform': '#7f8c8d', 'balanced': '#16a085', 'SR': '#34495e'}
LAB = {'golden-avoid': r'$r=\tau^{-2}$, anti', 'regular': r'$r=1/2$', 'golden': r'$r=\tau^{-2}$', 'plastic': r'$r=\rho^{-3}$', 'r0.3': r'$r=0.3$',
       'golden-side': r'$r=\tau^{-2}$, random side', 'largest-uniform': r'$r\sim U(0,1)$', 'balanced': 'balanced', 'SR': 'SR'}
# short legend names
LEG = {'golden-avoid': 'golden, anti-aligned', '(5,1)': '(5,1)', '(5,1)-avoid': '(5,1), anti-aligned', '(3,2)': 'plastic (3,2)',
       '(3,2)-avoid': 'plastic, anti-aligned', 'uniform': 'uniform split', 'golden': 'golden', 'golden-side': 'golden, random side',
       'plastic51-avoid': '(5,1), anti-aligned', 'r=0.3': r'$r=0.3$', '(2,1)': 'golden (2,1)', '(3,1)': '(3,1)', '(4,1)': '(4,1)', 'regular': 'periodic', 'golden': 'golden', 'plastic': 'plastic', 'r0.3': r'$r=0.3$',
       'golden-side': 'golden, rnd side', 'largest-uniform': 'uniform split', 'balanced': 'balanced', 'SR': 'SR',
       'supergolden': 'Pisot (3,1)', 'theta4': 'Pisot (4,1)'}


from matplotlib.lines import Line2D
TAU = qdiv.TAU
def panels3():
    fig = plt.figure(figsize=(7.0, 2.35))
    ax = [fig.add_axes([0.08 + 0.255 * i, 0.18, 0.175, 0.74]) for i in range(3)]
    for a, l in zip(ax, 'abc'):
        a.text(-0.4, 1.0, '(%s)' % l, transform=a.transAxes, fontsize=11, fontweight='bold')
    return fig, ax

# ---- Fig 1: six members, 2 x 3, one title per panel
T = 233
RULE1 = [('(a) $r=1/2$', ('largest', 0.5, 'perp'), None),
         (r'(b) $\tau^{-2}$ aligned', ('largest', TAU ** -2, 'perp'), None),
         (r'(c) $\tau^{-2}$ anti', ('largest', ('avoid', TAU ** -2, 'low'), 'perp'), None),
         (r'(d) $\tau^{-2}$ random', ('largest', ('side', TAU ** -2), 'perp'), None),
         (r'(e) $r\sim U(0,1)$', ('largest', 'uniform', 'perp'), None),
         ('(f) SR', None, 'SR')]
fig, ax = plt.subplots(2, 3, figsize=(3.4, 2.6))
plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.92, wspace=0.08, hspace=0.24)
for a, (tit, rule, nm) in zip(ax.ravel(), RULE1):
    g = qdiv.run(T, rule, seed=2) if rule else divisions.run(T, 3, nm)
    a.add_collection(LineCollection(np.stack([g['A'], g['B']], 1), colors='k', linewidths=0.35))
    a.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], 'k', lw=0.7)
    a.set_xlim(-0.01, 1.01); a.set_ylim(-0.01, 1.01); a.set_aspect('equal'); a.axis('off')
    a.set_title(tit, fontsize=8.5, pad=3, loc='left')
plt.savefig('fig1.pdf'); plt.close()

# ---- Fig 2: order. (a) 1D split sequences, (b) 2D junction sets for three side rules
fig = plt.figure(figsize=(7.0, 2.5))
ax = [fig.add_axes([0.095, 0.18, 0.29, 0.74]), fig.add_axes([0.485, 0.18, 0.29, 0.74])]
for a, l in zip(ax, 'ab'):
    a.text(-0.25, 1.0, '(%s)' % l, transform=a.transAxes, fontsize=11, fontweight='bold')
od = json.load(open(R + 'oned.json'))
C1 = {'(2,1)': '#c0392b', '(3,1)': '#e67e22', '(4,1)': '#d4a017', '(5,1)': '#16a085', '(3,2)': '#8e44ad'}
for nm, c in C1.items():
    d = [o for o in od if o['case'] == nm]
    ax[0].plot([o['N'] for o in d], [o['Imax'] for o in d], 'o-', ms=3, color=c, label=nm)
for nm in ('(4,3)', '(5,3)'):
    d = [o for o in od if o['case'] == nm]
    ax[0].plot([o['N'] for o in d], [o['Imax'] for o in d], 's:', ms=2.5, color='0.55')
d = [o for o in od if o['case'] == 'r=0.3']
ax[0].plot([o['N'] for o in d], [o['Imax'] for o in d], 'v--', ms=2.5, color='#2980b9')
Nn = np.array([1e3, 3e6]); ax[0].plot(Nn, 2.0 / Nn, 'k-', lw=0.5)
ax[0].text(3e4, 3e-6, r'$\propto N^{-1}$', fontsize=9)
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlabel('$N$'); ax[0].set_ylabel(r'$I_{\max}$ (1D)')
b2 = json.load(open(R + 'bragg2.json'))
STY = {'golden': ('#c0392b', 'o-'), 'golden-side': ('#27ae60', 'o-'), 'golden-avoid': ('#d35400', 'D-'),
       '(5,1)': ('#16a085', 'o-'), '(5,1)-avoid': ('#16a085', 'D--'), '(3,2)': ('#8e44ad', 'o-'), '(3,2)-avoid': ('#8e44ad', 'D--'),
       'r0.3': ('#2980b9', 'v--'), 'uniform': ('#7f8c8d', 's:')}
for nm, (c, st) in STY.items():
    d = [o for o in b2 if o['rule'] == nm]
    if d: ax[1].plot([o['T'] for o in d], [o['Imax'] for o in d], st, ms=3, color=c, mfc=c if 'avoid' not in nm else 'w')
Tt = np.array([3e3, 1.5e6]); ax[1].plot(Tt, 6.0 / Tt, 'k-', lw=0.5); ax[1].text(1.5e4, 1.5e-5, r'$\propto T^{-1}$', fontsize=9)
ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].set_xlabel('$T$'); ax[1].set_ylabel(r'$I_{\max}$ (2D)')
h = [Line2D([], [], color=c, marker='o', ms=3) for c in C1.values()] + \
    [Line2D([], [], color='0.55', ls=':', marker='s', ms=2.5), Line2D([], [], color='#2980b9', ls='--', marker='v', ms=2.5),
     Line2D([], [], color='#7f8c8d', ls=':', marker='s', ms=2.5), Line2D([], [], color='#27ae60', marker='o', ms=3),
     Line2D([], [], color='k', ls='-', marker='o', ms=3), Line2D([], [], color='k', ls='--', marker='D', ms=3, mfc='w')]
l = ['golden (2,1)', '(3,1)', '(4,1)', '(5,1)', 'plastic (3,2)', 'not Pisot', r'$r=0.3$', 'uniform split', 'golden, random side',
     'aligned', 'anti-aligned']
fig.legend(h, l, loc='center left', bbox_to_anchor=(0.80, 0.55), fontsize=7, handlelength=1.8)
plt.savefig('fig2.pdf', bbox_inches='tight', pad_inches=0.02); plt.close()

# ---- Fig 3: elasticity. (a) spanning lines, (b) exact pinned ceiling, (c) index vs density
fig, ax = panels3()
ce = json.load(open(R + 'ceiling.json'))
SC = {'golden': ('#c0392b', 'o-', 'golden, aligned'), 'r0.3': ('#2980b9', 'v--', r'$r=0.3$, aligned'),
      'golden-side': ('#27ae60', 'o-', 'golden, random side'), 'golden-avoid': ('#d35400', 'D-', 'golden, anti-aligned'),
      'plastic51-avoid': ('#16a085', 'D--', '(5,1), anti-aligned'), 'uniform': ('#7f8c8d', 's:', 'uniform split')}
for nm, (c, st, lab) in SC.items():
    Ts = sorted(set(o['T'] for o in ce if o['rule'] == nm))
    m = [np.mean([o['m'] for o in ce if o['rule'] == nm and o['T'] == T]) for T in Ts]
    nx = [np.mean([o['nmax_x'] for o in ce if o['rule'] == nm and o['T'] == T]) for T in Ts]
    mf = c if 'avoid' not in nm else 'w'
    ax[0].plot(Ts, m, st, ms=3, color=c, mfc=mf, label=lab); ax[1].plot(Ts, nx, st, ms=3, color=c, mfc=mf)
Tg = np.array([1e3, 1e6]); ax[0].plot(Tg, 1.6 * Tg ** 0.5, 'k-', lw=0.5); ax[0].text(1.5e4, 500, r'$T^{1/2}$', fontsize=9)
ax[1].plot(Tg, 1.45 * Tg ** 0.25, 'k-', lw=0.5); ax[1].text(2e3, 22, r'$T^{1/4}$', fontsize=9)
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlabel('$T$'); ax[0].set_ylabel('$m$')
ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].set_xlabel('$T$'); ax[1].set_ylabel(r'$n_{\max,x}$')
m3 = json.load(open(R + 'mech3.json'))
M3 = {'golden': ('#c0392b', 'o'), 'golden-avoid': ('#d35400', 'D'), '(5,1)-avoid': ('#16a085', 'D'),
      '(3,2)-avoid': ('#8e44ad', 'D'), 'uniform': ('#7f8c8d', 's')}
def m3get(nm, rho, key):
    Ts = sorted(set(o['T'] for o in m3 if o['rule'] == nm))
    return Ts, [np.mean([o[key] for o in m3 if o['rule'] == nm and o['T'] == T and o['rho_t'] == rho and o['model'] == 'EB']) for T in Ts]
for nm in ('golden', 'golden-avoid', 'uniform'):
    c, mk = M3[nm]; Ts, nx = m3get(nm, 0.1, 'n_x')
    ax[1].plot(Ts, nx, ':', marker='s', ms=3.5, color=c, mfc=c if 'avoid' not in nm else 'w')
for nm, (c, mk) in M3.items():
    for rho in (0.05, 0.1, 0.2):
        Ts, s = m3get(nm, rho, 'pin_ratio'); s = np.array(s); x = rho ** 2 * np.sqrt(Ts)
        ax[2].plot(x, s / (1 - s), mk, ms=3, color=c, mfc=c if 'avoid' not in nm else 'w', ls='none')
xg = np.array([0.02, 6]); ax[2].plot(xg, 0.17 / xg, 'k-', lw=0.5); ax[2].text(0.03, 0.06, r'$\propto x^{-1}$', fontsize=9)
ax[2].set_xscale('log'); ax[2].set_yscale('log'); ax[2].set_xlabel(r'$x=\bar\rho^{2}T^{1/2}$'); ax[2].set_ylabel(r'$s/(1-s)$', labelpad=1)
h, l = ax[0].get_legend_handles_labels()
h += [Line2D([], [], color='#8e44ad', ls='none', marker='D', ms=3, mfc='w'),
      Line2D([], [], color='k', ls=':', marker='s', ms=3.5)]
l += ['plastic, anti-aligned', r'rigid, $\bar\rho=0.1$']
fig.legend(h, l, loc='center left', bbox_to_anchor=(0.77, 0.55), fontsize=7, handlelength=1.8)
plt.savefig('fig3.pdf', bbox_inches='tight', pad_inches=0.02); plt.close()
print('ok')
