"""Fig. S1 of supplement.tex: recursive cell division step by step. (a) aligned (top) vs
anti-aligned (bottom) golden split; red: newest cut; dashed: the other candidate position.
(b) scale-rich rule; shaded: the cell picked (probability proportional to area); dot: the uniform
random point the cut passes through; red: newest cut, at a uniform random angle."""
import sys; sys.path[:0] = ['..', '../../beyond']
import numpy as np, qdiv, divisions
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
plt.rcParams.update({'font.size': 9})
r0 = qdiv.TAU ** -2
ROWS = [('aligned', ('largest', r0, 'perp')), ('anti-aligned', ('largest', ('avoid', r0, 'low'), 'perp'))]
STEPS = (1, 2, 3, 4, 6, 9, 15, 40)
fig, ax = plt.subplots(3, len(STEPS), figsize=(7.0, 3.45))
plt.subplots_adjust(left=0.1, right=0.995, bottom=0.02, top=0.92, wspace=0.12, hspace=0.32)
for i, (name, rule) in enumerate(ROWS):
    ax[i, 0].text(-0.2, 0.5, name, transform=ax[i, 0].transAxes, rotation=90, va='center', ha='center')
    for j, T in enumerate(STEPS):
        a = ax[i, j]; g = qdiv.run(T, rule)
        a.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], 'k', lw=0.8)
        a.add_collection(LineCollection(np.stack([g['A'][:-1], g['B'][:-1]], 1), colors='k', linewidths=0.6))
        a.plot(*np.stack([g['A'][-1], g['B'][-1]], 1), color='#c0392b', lw=1.4)
        x0, y0, w, h = g['cellrect'][-1]; f = 1 - g['frac'][-1]
        if g['orient'][-1]:
            a.plot([x0 + f * w] * 2, [y0, y0 + h], color='#c0392b', lw=0.8, ls=(0, (2, 1.5)))
        else:
            a.plot([x0, x0 + w], [y0 + f * h] * 2, color='#c0392b', lw=0.8, ls=(0, (2, 1.5)))
        a.set_xlim(-0.02, 1.02); a.set_ylim(-0.02, 1.02); a.set_aspect('equal'); a.axis('off')
        if i == 0: a.set_title('$T=%d$' % T, fontsize=9, pad=3)
# annotate r0 on the first panel
a = ax[0, 0]
a.annotate('', xy=(0, -0.08), xytext=(r0, -0.08), arrowprops=dict(arrowstyle='<->', lw=0.6), annotation_clip=False)
a.text(r0 / 2, -0.2, '$r_0$', ha='center', va='top', fontsize=8)
a.annotate('', xy=(r0, -0.08), xytext=(1, -0.08), arrowprops=dict(arrowstyle='<->', lw=0.6), annotation_clip=False)
a.text((1 + r0) / 2, -0.2, '$1-r_0$', ha='center', va='top', fontsize=8)
# (b) scale-rich rule
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 24
ax[2, 0].text(-0.2, 0.5, 'scale-rich', transform=ax[2, 0].transAxes, rotation=90, va='center', ha='center')
for j, T in enumerate(STEPS):
    a = ax[2, j]; g = divisions.run(T, SEED, 'SR')
    a.fill(*np.asarray(g['cell'][-1]).T, color='#eadcef', lw=0)
    a.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], 'k', lw=0.8)
    a.add_collection(LineCollection(np.stack([g['A'][:-1], g['B'][:-1]], 1), colors='k', linewidths=0.6))
    a.plot(*np.stack([g['A'][-1], g['B'][-1]], 1), color='#c0392b', lw=1.4)
    a.plot(*g['point'][-1], 'o', color='#c0392b', ms=2.5)
    a.set_xlim(-0.02, 1.02); a.set_ylim(-0.02, 1.02); a.set_aspect('equal'); a.axis('off')
ax[0, 0].text(-0.32, 1.25, '(a)', transform=ax[0, 0].transAxes, fontsize=9)
ax[2, 0].text(-0.32, 1.12, '(b)', transform=ax[2, 0].transAxes, fontsize=9)
plt.savefig('figS1.pdf'); plt.savefig('figS1.png', dpi=150); print('ok')
