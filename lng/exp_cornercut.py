"""Triangle closure from integral geometry.

For a convex cell K, a chord through a uniform point at a uniform angle has line
density sigma(h, theta) dh dtheta / (pi A) (sigma = chord length).  The chord hits
two ADJACENT sides iff it separates one vertex v from the rest; for fixed theta
this happens for h between the extreme projection and the next one, and the
integral of sigma over that range is the area of the cut-off corner.  Hence

   P_corner(K) = (1/(pi A)) int_0^pi [A_lo(theta) + A_hi(theta)] dtheta,

A_lo/hi = area of K below/above the second-smallest/second-largest vertex
projection.  Square: P = (4/pi) int_0^{pi/4} tan = 2 ln 2 / pi = 0.4413.
Triangle: 1.  The SR step uses area-weighted cells, so
   P_corner(SR) = sum_K (A_K/A) P_corner(K).
Check against the measured next-step corner-cut frequency, and tabulate by n.
Output: results/cornercut.txt
"""
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')


def pcorner(P, nth=180):
    A = srcells.poly_area(P)
    if A <= 0:
        return np.nan
    L = [0] * len(P)
    tot = 0.0
    ths = (np.arange(nth) + 0.5) * math.pi / nth
    for th in ths:
        nx, ny = -math.sin(th), math.cos(th)
        h = sorted(nx * x + ny * y for x, y in P)
        # area with projection <= h[1]
        Q, _ = srcells.clip(P, L, -nx, -ny, -h[1], 0)
        alo = srcells.poly_area(Q) if Q else 0.0
        Q, _ = srcells.clip(P, L, nx, ny, h[-2], 0)
        ahi = srcells.poly_area(Q) if Q else 0.0
        tot += alo + ahi
    return tot * (math.pi / nth) / (math.pi * A)


if __name__ == '__main__':
    lines = []
    sq = pcorner([(0, 0), (1, 0), (1, 1), (0, 1)], nth=2000)
    lines.append(f"square: formula {sq:.5f}  exact 2ln2/pi = {2*math.log(2)/math.pi:.5f}")
    lines.append(f"rectangle 1x4: {pcorner([(0,0),(4,0),(4,1),(0,1)], 2000):.5f}")
    lines.append(f"regular pentagon: {pcorner([(math.cos(2*math.pi*j/5), math.sin(2*math.pi*j/5)) for j in range(5)], 2000):.5f}")
    lines.append(f"regular hexagon: {pcorner([(math.cos(2*math.pi*j/6), math.sin(2*math.pi*j/6)) for j in range(6)], 2000):.5f}")
    lines.append(f"equilateral triangle: {pcorner([(0,0),(1,0),(0.5,0.8660254)], 500):.5f}")
    for T in (3000, 20000):
        for seed in (1, 2):
            o = srcells.run(T, seed=seed, keep_cells=True)
            P = o['cellP']; A = o['cellA']
            n = np.array([len(p) for p in P])
            pc = np.array([pcorner(p, nth=90) for p in P])
            w = A / A.sum()
            pred = np.nansum(w * pc)
            # measured: next-step corner-cut frequency from a continuation is the same
            # quantity in expectation; use the last 20% of the run as an estimate
            h = o['hit']; nb = (h >= 0).all(1)
            meas = o['adj'][int(0.8 * T):].mean()
            lines.append(f"T={T} seed={seed}: area-weighted P_corner (formula) = {pred:.4f}; measured (last 20% of steps, all splits) = {meas:.4f}")
            lines.append(f"   number fractions n=3..7: " + ' '.join(f"{(n==m).mean():.3f}" for m in range(3, 8)) + f"  <n>={n.mean():.4f}")
            lines.append(f"   area fractions   n=3..7: " + ' '.join(f"{w[n==m].sum():.3f}" for m in range(3, 8)))
            lines.append(f"   mean P_corner|n  n=3..7: " + ' '.join(f"{np.nanmean(pc[n==m]):.3f}" for m in range(3, 8)))
    txt = '\n'.join(lines); print(txt)
    open(os.path.join(OUT, 'cornercut.txt'), 'w').write(txt + '\n')
