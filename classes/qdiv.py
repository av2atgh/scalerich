"""Recursive divisions of the unit square with a prescribed split law.

Same rule set as beyond/divisions.py, (pick, split, angle), with the cut position
given by the area fraction r of one daughter instead of a point:
  pick   'largest' (deterministic, ties by birth order) | 'area' (prob. ∝ area)
  split  float r            daughter on the low side (left/bottom) gets fraction r
         ('side', r)        fraction r on a random side (fair coin)
         ('choice', r1, r2) r1 or r2 with prob. 1/2, low side
         ('avoid', r)       fraction r on the side (r or 1 - r) that creates fewer X-crossings
                            (a cut end coinciding with a junction already on the host from the
                            other side); ties by a fair coin
         ('avoid', r, 'low') same, ties to the low side (fully deterministic)
         'uniform'          r ~ U(0, 1)
  angle  'perp'  cut perpendicular to the long side; squares cut vertically
On rectangles 'perp' keeps every cell a rectangle, so the division is exact
(no angle noise, unlike divisions.py). Output has the keys of divisions.run.
"""
import heapq
import numpy as np

TAU = (1 + 5 ** 0.5) / 2
PLASTIC = 1.3247179572447460  # rho^3 = rho + 1, so rho^-2 + rho^-3 = 1


def run(T, rule=('largest', 1 / TAU, 'perp'), seed=0):
    pick, split, angle = rule
    assert angle == 'perp'
    rng = np.random.default_rng(seed)
    # cell = (x0, y0, w, h, labs) with labs = (bottom, right, top, left) ligament ids
    cells = {0: (0.0, 0.0, 1.0, 1.0, (-1, -2, -3, -4))}
    on = {}                      # ligament id -> set of rounded positions of junctions on it
    depth = {0: 0}
    born = {0: -1}                # cell id -> cut that created it (-1: the square)
    heap = [(-1.0, 0)]
    nxt = 1
    A = np.zeros((T, 2)); B = np.zeros((T, 2))
    hostA = np.zeros(T, int); hostB = np.zeros(T, int)
    frac = np.zeros(T); cut_depth = np.zeros(T, int); orient = np.zeros(T, int)
    cellrect = np.zeros((T, 4)); parent = np.full(T, -1); cand = np.zeros((T, 2), int)
    for t in range(T):
        if pick == 'largest':
            _, c = heapq.heappop(heap)
        else:  # area-proportional: point location by linear scan of leaves (small T only)
            keys = list(cells); a = np.array([cells[k][2] * cells[k][3] for k in keys])
            c = keys[rng.choice(len(keys), p=a / a.sum())]
        x0, y0, w, h, (lb, lr, lt, ll) = cells.pop(c)
        if split == 'uniform':
            r = rng.random()
        elif isinstance(split, tuple) and split[0] == 'side':
            r = split[1] if rng.random() < 0.5 else 1 - split[1]
        elif isinstance(split, tuple) and split[0] == 'avoid':
            r = None                     # set below, once the cut direction is known
        elif isinstance(split, tuple) and split[0] == 'choice':
            r = split[1] if rng.random() < 0.5 else split[2]
        else:
            r = float(split)
        vert = w >= h * (1 - 1e-12)
        if isinstance(split, tuple) and split[0] == 'avoid':
            def ncross(rr):
                if vert:
                    k = round(x0 + rr * w, 9); hs = (lb, lt)
                else:
                    k = round(y0 + rr * h, 9); hs = (ll, lr)
                return sum(1 for hh in hs if hh >= 0 and k in on.get(hh, ()))
            c0, c1 = ncross(split[1]), ncross(1 - split[1])
            tie = (rng.random() < 0.5) if len(split) < 3 else (split[2] == 'low')
            r = split[1] if (c0 < c1 or (c0 == c1 and tie)) else 1 - split[1]
            cand[t] = (c0, c1)
        if vert:
            xc = x0 + r * w
            A[t] = (xc, y0); B[t] = (xc, y0 + h); hostA[t] = lb; hostB[t] = lt
            ch = [(x0, y0, r * w, h, (lb, t, lt, ll)), (xc, y0, (1 - r) * w, h, (lb, lr, lt, t))]
        else:
            yc = y0 + r * h
            A[t] = (x0, yc); B[t] = (x0 + w, yc); hostA[t] = ll; hostB[t] = lr
            ch = [(x0, y0, w, r * h, (lb, lr, t, ll)), (x0, yc, w, (1 - r) * h, (t, lr, lt, ll))]
        frac[t] = r; cut_depth[t] = depth[c]; orient[t] = int(vert)
        cellrect[t] = (x0, y0, w, h); parent[t] = born[c]
        k = round(xc if vert else yc, 9)
        for hh in ((lb, lt) if vert else (ll, lr)):
            on.setdefault(hh, set()).add(k)
        for cc in ch:
            cells[nxt] = cc; depth[nxt] = depth[c] + 1; born[nxt] = t
            if pick == 'largest':
                heapq.heappush(heap, (-cc[2] * cc[3], nxt))
            nxt += 1
    keys = sorted(cells)
    polys = [np.array([[x, y], [x + w, y], [x + w, y + h], [x, y + h]])
             for (x, y, w, h, _) in (cells[k] for k in keys)]
    return dict(T=T, A=A, B=B, hostA=hostA, hostB=hostB,
                leaf_area=np.array([cells[k][2] * cells[k][3] for k in keys]),
                leaf_wh=np.array([[cells[k][2], cells[k][3]] for k in keys]),
                leaf_xy=np.array([[cells[k][0], cells[k][1]] for k in keys]),
                leaf_polys=polys, leaf_depth=np.array([depth[k] for k in keys]),
                frac=frac, cut_depth=cut_depth, orient=orient, cellrect=cellrect, parent=parent, cand=cand)
