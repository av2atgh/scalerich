"""Face-based (cell-division) simulator of the scale-rich (SR) model.

Zero thickness: the SR step "random point + random angle, grow both ways to the
first hit" is exactly: pick a cell of the current tessellation with probability
proportional to its area, draw the chord of that (convex) cell through a uniform
point at a uniform angle.  Every cell is a convex polygon whose edges carry the
label of the ligament (dual node) they lie on (-1 = outer square).

Thickness: ligament t has width lam_t = lam0 * t**(-alpha) (t = 1, 2, ...).  The
free region of a cell is the convex polygon bounded by the faces of the
neighbouring rectangles; a new centreline grows to the free-polygon boundary
(= the surface of the ligament it hits) and each child is the parent clipped by
the half-plane offset lam_t/2 from the centreline.  Nucleation only in free area.

Knobs (all generalisations of the paper's Fig. 6 rules):
  q        cell weight = area**q  (q = 1 is the SR model; q = 0 uniform cell)
  angle    'uniform' | 'set' (angles drawn from `angles`) | 'vonmises' (kappa
           around angle0) | 'field' (theta = field(x,y) + vonmises noise, with prob. pfield; else uniform)
  dens     None or callable g(x, y) >= 0: nucleation density; cell weight is
           area**q * g(centroid) and the point is drawn by rejection inside the
           cell.

Returns a dict of numpy arrays (see `run`).
"""
import math
import numpy as np


class Fenwick:
    def __init__(self, n):
        self.n = n
        self.t = [0.0] * (n + 1)
        self.w = [0.0] * n
        self.tot = 0.0

    def set(self, i, w):
        d = w - self.w[i]
        self.w[i] = w
        self.tot += d
        i += 1
        t = self.t
        n = self.n
        while i <= n:
            t[i] += d
            i += i & (-i)

    def find(self, u):
        """smallest i with prefix(i) > u"""
        pos = 0
        t = self.t
        step = 1 << (self.n.bit_length())
        while step:
            nx = pos + step
            if nx <= self.n and t[nx] <= u:
                pos = nx
                u -= t[nx]
            step >>= 1
        return min(pos, self.n - 1)


def poly_area(P):
    a = 0.0
    n = len(P)
    for i in range(n):
        x0, y0 = P[i]
        x1, y1 = P[(i + 1) % n]
        a += x0 * y1 - x1 * y0
    return 0.5 * a


def poly_centroid(P, A):
    cx = cy = 0.0
    n = len(P)
    for i in range(n):
        x0, y0 = P[i]
        x1, y1 = P[(i + 1) % n]
        c = x0 * y1 - x1 * y0
        cx += (x0 + x1) * c
        cy += (y0 + y1) * c
    return cx / (6 * A), cy / (6 * A)


def sample_in_poly(P, rng):
    # fan triangulation from vertex 0
    x0, y0 = P[0]
    areas = []
    for i in range(1, len(P) - 1):
        x1, y1 = P[i]
        x2, y2 = P[i + 1]
        areas.append(0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)))
    s = sum(areas)
    u = rng.random() * s
    j = 0
    while j < len(areas) - 1 and u > areas[j]:
        u -= areas[j]
        j += 1
    x1, y1 = P[j + 1]
    x2, y2 = P[j + 2]
    r1, r2 = rng.random(), rng.random()
    if r1 + r2 > 1:
        r1, r2 = 1 - r1, 1 - r2
    return (x0 + r1 * (x1 - x0) + r2 * (x2 - x0), y0 + r1 * (y1 - y0) + r2 * (y2 - y0))


def clip(P, L, nx, ny, c, newlab):
    """keep the part of convex polygon P with nx*x+ny*y >= c; edge i of P
    (P[i]->P[i+1]) has label L[i]; the cut edge gets label newlab."""
    n = len(P)
    f = [nx * x + ny * y - c for (x, y) in P]
    outP = []
    outL = []
    for i in range(n):
        j = (i + 1) % n
        fi, fj = f[i], f[j]
        if fi >= 0:
            outP.append(P[i])
            if fj >= 0:
                outL.append(L[i])
            else:  # leaving: edge i up to crossing, then cut edge
                s = fi / (fi - fj)
                xi, yi = P[i]
                xj, yj = P[j]
                outL.append(L[i])
                outP.append((xi + s * (xj - xi), yi + s * (yj - yi)))
                outL.append(newlab)
        else:
            if fj >= 0:  # entering: crossing point then rest of edge i
                s = fi / (fi - fj)
                xi, yi = P[i]
                xj, yj = P[j]
                outP.append((xi + s * (xj - xi), yi + s * (yj - yi)))
                outL.append(L[i])
    if len(outP) < 3:
        return None, None
    return outP, outL


def run(T, alpha=0.5, lam0=0.0, q=1.0, angle='uniform', angles=None, kappa=0.0,
        angle0=0.0, field=None, pfield=1.0, dens=None, seed=0, keep_cells=False, maxtry=200):
    rng = np.random.default_rng(seed)
    pyrng = rng  # numpy Generator: .random()
    fw = Fenwick(T + 2)
    cellP = [[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]]
    cellL = [[-1, -1, -1, -1]]
    cellA = [1.0]
    cellB = [-1]          # birth time of cell
    if dens is None:
        fw.set(0, 1.0)
    else:
        cx, cy = 0.5, 0.5
        fw.set(0, dens(cx, cy))
    # per ligament records
    born_cell_n = np.zeros(T, np.int32)      # corners of the split cell
    hit = np.full((T, 2), -2, np.int64)      # labels hit (-1 boundary)
    adjflag = np.zeros(T, np.int8)           # the two hit sides were adjacent in the cell
    length = np.zeros(T)
    lam = np.zeros(T)
    theta = np.zeros(T)
    mid = np.zeros((T, 2))
    ends = np.zeros((T, 4))
    cellarea = np.zeros(T)                   # area of the split cell
    child_n = np.zeros((T, 2), np.int32)
    freeA = np.zeros(T)                      # free area after step
    tot = 1.0
    Tdone = T
    for t in range(T):
        if fw.tot <= 1e-14:
            Tdone = t
            break
        lt = lam0 * (t + 1) ** (-alpha) if lam0 > 0 else 0.0
        # choose cell
        for _ in range(maxtry):
            c = fw.find(rng.random() * fw.tot)
            if fw.w[c] > 0:
                break
        P = cellP[c]
        L = cellL[c]
        n = len(P)
        if dens is None:
            px, py = sample_in_poly(P, rng)
        else:
            # rejection on g inside the cell (g bounded by gmax on the cell vertices / centroid)
            gm = max(dens(x, y) for (x, y) in P)
            gm = max(gm, dens(*poly_centroid(P, cellA[c])))
            for _ in range(1000):
                px, py = sample_in_poly(P, rng)
                if rng.random() * gm <= dens(px, py):
                    break
        # angle
        if angle == 'uniform':
            th = rng.random() * math.pi
        elif angle == 'set':
            th = angles[rng.integers(len(angles))]
        elif angle == 'vonmises':
            th = (angle0 + rng.vonmises(0.0, kappa) / 2.0) % math.pi  # axial
        elif angle == 'field':
            if pfield < 1.0 and rng.random() >= pfield:
                th = rng.random() * math.pi
            else:
              th = (field(px, py) + (rng.vonmises(0.0, kappa) / 2.0 if kappa > 0 else 0.0)) % math.pi
        else:
            raise ValueError(angle)
        dx, dy = math.cos(th), math.sin(th)
        # chord: intersect with edges
        smax = math.inf
        smin = -math.inf
        ep = em = -1
        for i in range(n):
            ax, ay = P[i]
            bx, by = P[(i + 1) % n]
            ex, ey = bx - ax, by - ay
            den = dx * ey - dy * ex
            if den == 0.0:
                continue
            wx, wy = ax - px, ay - py
            s = (wx * ey - wy * ex) / den
            u = (wx * dy - wy * dx) / den
            if -1e-12 <= u <= 1 + 1e-12:
                if s > 0 and s < smax:
                    smax = s
                    ep = i
                elif s < 0 and s > smin:
                    smin = s
                    em = i
        if ep < 0 or em < 0:
            # numerical failure; skip this step (count it)
            hit[t] = (-3, -3)
            continue
        lp, lm = L[ep], L[em]
        hit[t] = (lp, lm)
        d = (ep - em) % n
        adjflag[t] = 1 if (d == 1 or d == n - 1) else 0
        born_cell_n[t] = n
        length[t] = smax - smin
        lam[t] = lt
        theta[t] = th
        mid[t] = (px + 0.5 * (smax + smin) * dx, py + 0.5 * (smax + smin) * dy)
        ends[t] = (px + smin * dx, py + smin * dy, px + smax * dx, py + smax * dy)
        cellarea[t] = cellA[c]
        # children: left normal (-dy, dx)
        nx, ny = -dy, dx
        c0 = nx * px + ny * py
        P1, L1 = clip(P, L, nx, ny, c0 + 0.5 * lt, t)
        P2, L2 = clip(P, L, -nx, -ny, -c0 + 0.5 * lt, t)
        A1 = poly_area(P1) if P1 is not None else 0.0
        A2 = poly_area(P2) if P2 is not None else 0.0
        tot += A1 + A2 - cellA[c]
        child_n[t] = (len(P1) if P1 else 0, len(P2) if P2 else 0)
        # replace c by child 1, append child 2
        cellP[c] = P1 if P1 is not None else [P[0], P[0], P[0]]
        cellL[c] = L1 if L1 is not None else [t, t, t]
        cellA[c] = A1
        cellB[c] = t
        cellP.append(P2 if P2 is not None else [P[0], P[0], P[0]])
        cellL.append(L2 if L2 is not None else [t, t, t])
        cellA.append(A2)
        cellB.append(t)
        k2 = len(cellP) - 1
        for cc, Ac in ((c, A1), (k2, A2)):
            if Ac <= 0:
                fw.set(cc, 0.0)
            else:
                w = Ac ** q if q != 1.0 else Ac
                if dens is not None:
                    w *= dens(*poly_centroid(cellP[cc], Ac))
                fw.set(cc, w)
        freeA[t] = tot
    out = dict(T=Tdone, hit=hit[:Tdone], adj=adjflag[:Tdone], ncell=born_cell_n[:Tdone],
               length=length[:Tdone], lam=lam[:Tdone], theta=theta[:Tdone], mid=mid[:Tdone],
               ends=ends[:Tdone], cellarea=cellarea[:Tdone], child_n=child_n[:Tdone],
               freeA=freeA[:Tdone])
    if keep_cells:
        out['cellP'] = cellP
        out['cellL'] = cellL
        out['cellA'] = np.array(cellA)
    return out


def dual_edges(out):
    """edge list (i, j) of the dual graph (ligament = node, junction = edge)."""
    h = out['hit']
    T = len(h)
    e = []
    for col in (0, 1):
        m = h[:, col] >= 0
        e.append(np.stack([np.arange(T)[m], h[m, col]], 1))
    E = np.concatenate(e)
    # a ligament could hit the same ligament at both ends only if non-convex: drop dups
    E = np.unique(np.sort(E, 1), axis=0)
    return E
