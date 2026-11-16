"""Fig. S2 of supplement.tex: four property circles (ordered, disordered; stretch dominated,
bending dominated) with the divisions placed in their overlaps."""
import sys; sys.path[:0] = ['..', '../../beyond']
import numpy as np, qdiv, divisions
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.collections import LineCollection
plt.rcParams.update({'font.size': 9})
TAU = qdiv.TAU
fig = plt.figure(figsize=(7.0, 4.45))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 13.6); ax.set_ylim(0, 8.6); ax.set_aspect('equal'); ax.axis('off')
# order ovals (filled), elastic ovals (outlined)
ax.add_patch(Ellipse((3.3, 4.0), 5.0, 7.7, fc='#d6e4f0', ec='#2c6e9e', lw=1.2, alpha=0.9))
ax.add_patch(Ellipse((8.4, 4.0), 3.8, 7.7, fc='#eadcef', ec='#7d3c98', lw=1.2, alpha=0.9))
ax.add_patch(Ellipse((2.1, 4.0), 1.8, 6.2, fc='none', ec='#2c6e9e', lw=0.9, ls=(0, (3, 2))))
ax.add_patch(Ellipse((6.0, 6.15), 11.6, 2.5, fc='none', ec='#c0392b', lw=1.6))
ax.add_patch(Ellipse((6.0, 1.85), 11.6, 2.5, fc='none', ec='#1e8449', lw=1.6))
ax.text(3.3, 7.95, 'ordered (Bragg peaks)', ha='center', va='bottom', color='#2c6e9e')
ax.text(8.4, 7.95, 'disordered (no Bragg peaks)', ha='center', va='bottom', color='#7d3c98')
ax.text(2.1, 4.0, 'periodic', ha='center', va='center', rotation=90, color='#2c6e9e', fontsize=8)
ax.text(11.9, 6.15, 'aligned\n$\\mu=1/2$\n(stretch\ndominated)', ha='left', va='center', color='#c0392b', fontsize=8.5)
ax.text(11.9, 1.85, 'non-aligned\n$\\mu=0$\n(bending\ndominated)', ha='left', va='center', color='#1e8449', fontsize=8.5)

def thumb(cx, cy, g, label, size=1.05):
    a = fig.add_axes([(cx - size / 2) / 13.6, (cy - size / 2 + 0.12) / 8.6, size / 13.6, size / 8.6])
    a.add_collection(LineCollection(np.stack([g['A'], g['B']], 1), colors='k', linewidths=0.3))
    a.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], 'k', lw=0.5)
    a.set_xlim(-0.02, 1.02); a.set_ylim(-0.02, 1.02); a.set_aspect('equal'); a.axis('off')
    ax.text(cx, cy - size / 2 + 0.05, label, ha='center', va='top', fontsize=7.5)
T = 89
thumb(2.1, 6.15, qdiv.run(63, ('largest', 0.5, 'perp')), '$r_0=1/2$', size=0.85)
thumb(3.9, 6.15, qdiv.run(T, ('largest', TAU ** -2, 'perp')), 'aligned Pisot')
thumb(3.9, 1.85, qdiv.run(T, ('largest', ('avoid', TAU ** -2, 'low'), 'perp')), 'anti-aligned Pisot')
thumb(3.9, 4.0, qdiv.run(T, ('largest', ('side', TAU ** -2), 'perp'), seed=2), 'random side\n' + r'partially aligned, $\mu\approx0.2$', size=0.9)
ax.text(2.1, 1.85, r'$\emptyset$', ha='center', va='center', fontsize=12, color='#2c6e9e')
thumb(8.4, 6.15, qdiv.run(T, ('largest', 0.3, 'perp')), r'aligned $r_0=0.3$')
thumb(7.75, 1.85, qdiv.run(T, ('largest', 'uniform', 'perp'), seed=2), 'uniform split', size=0.95)
thumb(9.05, 1.85, divisions.run(T, 3, 'SR'), 'SR', size=0.95)
plt.savefig('figS2.pdf'); print('ok')
