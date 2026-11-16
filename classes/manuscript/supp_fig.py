"""Fig. S1 of supplement.tex: recursive cell division step by step, aligned (top) vs anti-aligned
(bottom) golden split. Red: newest cut; dashed: the other candidate position."""
import sys; sys.path[:0] = ['..']
import numpy as np, qdiv
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
plt.rcParams.update({'font.size': 9})
r0 = qdiv.TAU ** -2
ROWS = [('aligned', ('largest', r0, 'perp')), ('anti-aligned', ('largest', ('avoid', r0, 'low'), 'perp'))]
STEPS = (1, 2, 3, 4, 6, 9, 15, 40)
fig, ax = plt.subplots(2, len(STEPS), figsize=(7.0, 2.25))
plt.subplots_adjust(left=0.1, right=0.995, bottom=0.02, top=0.88, wspace=0.12, hspace=0.12)
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
plt.savefig('figS1.pdf'); print('ok')
