"""Parallel iterated division: SR with finite growth speed (Gilbert <-> Cowan/SR).

Unit square, growth speed v = 1. Segments nucleate as a Poisson process of rate
R per unit area per unit time at a uniform point with angle theta (uniform, or
von Mises about 0 with concentration kappa on the doubled angle). Each segment
grows both ways at speed 1; a tip stops where it first meets a segment that is
already there (static or still growing) or the boundary. At a crossing of two
growing segments the later arrival stops (Gilbert race). Exact event-driven
resolution inside each time step; the step only bounds the candidate search.

Pi = R L^3 / v (L = 1) = nucleations per domain-crossing time.
Pi -> 0 : one segment at a time = SR (Cowan 2010, cell picked prop. to area).
Pi >> N^2: everything nucleates before anything stops = Gilbert tessellation.
Local crossover expected at n* ~ Pi^2 segments (cell size n^-1/2 crossed in
time n^-1/2 while R n^-1/2 = Pi n^-1/2 nucleations occur).

Turnover (delta > 0): each segment is removed at rate delta; tips that had
stopped on it resume growth from their stop point (the obstacle is gone).
Then the network reaches a steady state with mean number ~ R/delta-limited.

run(N, Pi, seed, kappa=0, delta=0, Tmax=None) -> dict with
  P (N,2), D (N,2), t0, tend (removal time or inf), sp, sm (final extents +/-),
  hitp, hitm (segment hit at + / - tip; -1 boundary, -2 still growing/free),
  E (dual edges between segments present at the end).
"""
import math
import numpy as np
from numba import njit

BIG = 1e18


@njit(cache=True)
def _cross(px, py, dx, dy, qx, qy, ex, ey):
    """Signed distances (s along d from p, u along e from q) of line crossing."""
    den = dx * ey - dy * ex
    if abs(den) < 1e-14:
        return BIG, BIG
    wx = qx - px; wy = qy - py
    s = (wx * ey - wy * ex) / den
    u = (wx * dy - wy * dx) / den
    return s, u


@njit(cache=True)
def _bdist(px, py, dx, dy):
    """distance along +d to the unit-square boundary."""
    best = BIG
    if dx > 1e-15:
        best = min(best, (1 - px) / dx)
    elif dx < -1e-15:
        best = min(best, -px / dx)
    if dy > 1e-15:
        best = min(best, (1 - py) / dy)
    elif dy < -1e-15:
        best = min(best, -py / dy)
    return best


@njit(cache=True)
def _cells_on_ray(px, py, dx, dy, s0, s1, G, out):
    """grid cells (index) met by points p + s d, s in [s0, s1]; returns count."""
    n = 0
    x0 = px + s0 * dx; y0 = py + s0 * dy
    x1 = px + s1 * dx; y1 = py + s1 * dy
    i = min(G - 1, max(0, int(x0 * G))); j = min(G - 1, max(0, int(y0 * G)))
    i1 = min(G - 1, max(0, int(x1 * G))); j1 = min(G - 1, max(0, int(y1 * G)))
    out[n] = i * G + j; n += 1
    if i == i1 and j == j1:
        return n
    stepx = 1 if dx > 0 else -1
    stepy = 1 if dy > 0 else -1
    h = 1.0 / G
    if abs(dx) > 1e-15:
        nx = (i + (1 if dx > 0 else 0)) * h
        tmx = (nx - x0) / dx; tdx = h / abs(dx)
    else:
        tmx = BIG; tdx = BIG
    if abs(dy) > 1e-15:
        ny = (j + (1 if dy > 0 else 0)) * h
        tmy = (ny - y0) / dy; tdy = h / abs(dy)
    else:
        tmy = BIG; tdy = BIG
    L = s1 - s0
    while n < out.shape[0]:
        if tmx < tmy:
            if tmx > L:
                break
            i += stepx; tmx += tdx
        else:
            if tmy > L:
                break
            j += stepy; tmy += tdy
        if i < 0 or i >= G or j < 0 or j >= G:
            break
        out[n] = i * G + j; n += 1
        if i == i1 and j == j1:
            break
    return n


@njit(cache=True)
def _remove(j, n, present, alive, hit, ext, bnd, lastcell, kdeg):
    """remove segment j; tips stopped on j resume. Returns change in alive tips."""
    dn = 0
    present[j] = False
    for sd in range(2):
        if alive[j, sd]:
            alive[j, sd] = False; dn -= 1
    for sd in range(2):
        if hit[j, sd] >= 0:
            kdeg[hit[j, sd]] -= 1
    kdeg[j] = 0
    for i in range(n):
        if present[i]:
            for sd in range(2):
                if hit[i, sd] == j:
                    hit[i, sd] = -2; kdeg[i] -= 1
                    if ext[i, sd] < bnd[i, sd] - 1e-12:
                        alive[i, sd] = True; dn += 1
                        lastcell[i, sd] = -1
    return dn


@njit(cache=True)
def _sim(N, R, seed, kappa, delta, Tmax, G, maxent, cap, stab):
    np.random.seed(seed)
    P = np.zeros((N, 2)); D = np.zeros((N, 2)); t0 = np.full(N, BIG)
    tend = np.full(N, BIG)
    ext = np.zeros((N, 2))          # current extent of side 0 (+d) and 1 (-d)
    alive = np.zeros((N, 2), np.bool_)
    hit = np.full((N, 2), -2, np.int64)
    bnd = np.zeros((N, 2))
    present = np.zeros(N, np.bool_)
    kdeg = np.zeros(N, np.int64)
    # grid linked lists
    head = np.full(G * G, -1, np.int64)
    nxt = np.empty(maxent, np.int64); eseg = np.empty(maxent, np.int64)
    ne = 0
    lastcell = np.full((N, 2), -1, np.int64)
    cellbuf = np.empty(4 * G + 8, np.int64)
    # candidate events
    ev_t = np.empty(cap); ev_i = np.empty(cap, np.int64); ev_side = np.empty(cap, np.int64)
    ev_s = np.empty(cap); ev_j = np.empty(cap, np.int64); ev_jside = np.empty(cap, np.int64)
    ev_u = np.empty(cap)
    # dependents: tips stopped on segment j (for turnover resume) -- found by scan
    n = 0; t = 0.0; nalive = 0; npresent = 0
    tnext = np.random.exponential(1.0 / R)
    nrem_total = 0
    step = 0
    while True:
        if n >= N and nalive == 0:
            break
        if t >= Tmax:
            break
        # time step: bounded by typical cell size
        dt = 0.25 / math.sqrt(npresent + 4.0)
        if nalive == 0 and n < N and tnext > t:
            if delta > 0 and npresent > 0:
                trem = t + np.random.exponential(1.0 / (delta * npresent))
                if trem < tnext and trem < Tmax:
                    # one removal, then continue (memoryless)
                    t = trem
                    j = np.random.randint(0, n)
                    if present[j] and np.random.random() < (2.0 / max(kdeg[j], 2)) ** stab:
                        tend[j] = t; npresent -= 1; nrem_total += 1
                        nalive += _remove(j, n, present, alive, hit, ext, bnd, lastcell, kdeg)
                    continue
            t = tnext - 1e-12 if tnext > t else t
        t1 = t + dt
        if t1 > Tmax:
            t1 = Tmax
        step += 1
        # turnover: removals in (t, t1]
        if delta > 0 and npresent > 0:
            pr = 1.0 - math.exp(-delta * (t1 - t))
            for j in range(n):
                if present[j] and np.random.random() < pr * (2.0 / max(kdeg[j], 2)) ** stab:
                    tend[j] = t; npresent -= 1; nrem_total += 1
                    nalive += _remove(j, n, present, alive, hit, ext, bnd, lastcell, kdeg)
        # nucleations in (t, t1]
        while n < N and tnext <= t1:
            th = 0.0
            if kappa > 0:
                th = 0.5 * np.random.vonmises(0.0, kappa)
            else:
                th = np.random.random() * math.pi
            P[n, 0] = np.random.random(); P[n, 1] = np.random.random()
            D[n, 0] = math.cos(th); D[n, 1] = math.sin(th)
            t0[n] = tnext
            bnd[n, 0] = _bdist(P[n, 0], P[n, 1], D[n, 0], D[n, 1])
            bnd[n, 1] = _bdist(P[n, 0], P[n, 1], -D[n, 0], -D[n, 1])
            alive[n, 0] = True; alive[n, 1] = True; nalive += 2
            present[n] = True; npresent += 1
            n += 1
            tnext += np.random.exponential(1.0 / R)
        # register sweeps (new cells along [ext, ext_end]) for alive tips
        for i in range(n):
            if not present[i]:
                continue
            for sd in range(2):
                if not alive[i, sd]:
                    continue
                sgn = 1.0 if sd == 0 else -1.0
                tb = max(t, t0[i])
                send = min(ext[i, sd] + (t1 - tb), bnd[i, sd])
                if send < ext[i, sd]:
                    send = ext[i, sd]
                m = _cells_on_ray(P[i, 0], P[i, 1], sgn * D[i, 0], sgn * D[i, 1],
                                  ext[i, sd], send, G, cellbuf)
                for q in range(m):
                    c = cellbuf[q]
                    if c == lastcell[i, sd]:
                        continue
                    # avoid immediate duplicates of the last few cells
                    if ne >= maxent:
                        return P, D, t0, tend, ext, hit, present, n, -1
                    eseg[ne] = i; nxt[ne] = head[c]; head[c] = ne; ne += 1
                    lastcell[i, sd] = c
        # candidate events
        nev = 0
        for i in range(n):
            if not present[i]:
                continue
            for sd in range(2):
                if not alive[i, sd]:
                    continue
                sgn = 1.0 if sd == 0 else -1.0
                dx = sgn * D[i, 0]; dy = sgn * D[i, 1]
                tb = max(t, t0[i])
                sa = ext[i, sd]
                sb = sa + (t1 - tb)
                # boundary
                if bnd[i, sd] <= sb:
                    if nev >= cap:
                        return P, D, t0, tend, ext, hit, present, n, -2
                    ev_t[nev] = tb + (bnd[i, sd] - sa)
                    ev_i[nev] = i; ev_side[nev] = sd; ev_s[nev] = bnd[i, sd]
                    ev_j[nev] = -1; ev_jside[nev] = -1; ev_u[nev] = 0.0
                    nev += 1
                sbb = min(sb, bnd[i, sd])
                m = _cells_on_ray(P[i, 0], P[i, 1], dx, dy, sa, sbb, G, cellbuf)
                for q in range(m):
                    e = head[cellbuf[q]]
                    while e >= 0:
                        j = eseg[e]; e = nxt[e]
                        if j == i or not present[j]:
                            continue
                        s, u = _cross(P[i, 0], P[i, 1], dx, dy, P[j, 0], P[j, 1], D[j, 0], D[j, 1])
                        if s <= sa or s > sbb:
                            continue
                        jsd = 0 if u >= 0 else 1
                        au = abs(u)
                        if au > bnd[j, jsd]:
                            continue
                        tau = tb + (s - sa)
                        # does j cover the crossing by tau?
                        if au <= ext[j, jsd]:
                            pass  # already there at step start
                        elif alive[j, jsd]:
                            tbj = max(t, t0[j])
                            sig = tbj + (au - ext[j, jsd])
                            if sig > tau:
                                continue
                        else:
                            continue
                        if nev >= cap:
                            return P, D, t0, tend, ext, hit, present, n, -2
                        ev_t[nev] = tau; ev_i[nev] = i; ev_side[nev] = sd; ev_s[nev] = s
                        ev_j[nev] = j; ev_jside[nev] = jsd; ev_u[nev] = au
                        nev += 1
        # process events in time order
        order = np.argsort(ev_t[:nev])
        for q in range(nev):
            k = order[q]
            i = ev_i[k]; sd = ev_side[k]
            if not alive[i, sd]:
                continue
            j = ev_j[k]
            if j >= 0:
                jsd = ev_jside[k]; au = ev_u[k]
                # j must reach au by tau: either still growing (it passes au at
                # sig <= tau; a stop of j short of au would be an earlier event,
                # already processed) or stopped at a point >= au.
                if au > ext[j, jsd] and not alive[j, jsd]:
                    continue
            alive[i, sd] = False; nalive -= 1
            ext[i, sd] = ev_s[k]; hit[i, sd] = j
            if j >= 0:
                kdeg[i] += 1; kdeg[j] += 1
        # advance surviving tips
        for i in range(n):
            if not present[i]:
                continue
            for sd in range(2):
                if alive[i, sd]:
                    tb = max(t, t0[i])
                    ext[i, sd] = min(ext[i, sd] + (t1 - tb), bnd[i, sd])
        t = t1
    return P, D, t0, tend, ext, hit, present, n, step


def run(N, Pi, seed=0, kappa=0.0, delta=0.0, Tmax=None, G=None, stab=0.0):
    if G is None:
        G = max(8, int(math.sqrt(N) / 2))
    if Tmax is None:
        Tmax = BIG
    maxent = int(60 * N + 200 * G * G) + 100000
    cap = 1 << 16
    while True:
        P, D, t0, tend, ext, hit, present, n, st = _sim(N, float(Pi), seed, float(kappa),
                                                       float(delta), float(Tmax), G, maxent, cap, float(stab))
        if st == -2:
            cap *= 4; continue
        if st == -1:
            maxent *= 2; continue
        break
    P, D, t0, tend, ext, hit, present = P[:n], D[:n], t0[:n], tend[:n], ext[:n], hit[:n], present[:n]
    E = []
    for i in range(n):
        if not present[i]:
            continue
        for sd in range(2):
            j = hit[i, sd]
            if j >= 0 and present[j]:
                E.append((i, j))
    return dict(P=P, D=D, t0=t0, tend=tend, sp=ext[:, 0], sm=ext[:, 1], hitp=hit[:, 0],
                hitm=hit[:, 1], present=present, E=np.array(E, dtype=np.int64).reshape(-1, 2),
                steps=st)
