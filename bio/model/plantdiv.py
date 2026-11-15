"""Plant-like iterated division in a vertex model.

Tissue = planar polygon partition of the unit square (shared vertices, straight walls).
Division (rule):
  'errera'   wall through the centroid, shortest among K directions, subject to
             no corner cut (the two cut walls are not adjacent around the cell) and
             no four-way junction (end points >= dmin x wall length from wall ends).
  'errera_cc' same, corner cuts allowed (only the four-way constraint).
  'sr'       random point, random angle (SR / Cowan), no constraints.
Cell choice (pick): 'area' (prop. to area) or 'size' (prop. to area^beta, beta = 4:
size-threshold-like, plants).
Mechanics (relax=True): quasi-static vertex model,
  E = sum_e T_e l_e + sum_c (K/2) (A_c - A0_c)^2 / A0_c^{3/2},
  A0_c = area at birth; new-wall tension matures T(a) = 1 - (1 - t0) exp(-a / tau),
  a = wall age in cell cycles (time advances 1/N per division). Gradient descent,
  boundary vertices slide on the square, corners fixed.
Bookkeeping: each wall carries a line id (the division that made it); splitting a
wall keeps its id; the dual edge (new line -> host line) is recorded at each end.
"""
import math
import numpy as np


class Tissue:
    def __init__(self, seed=0, rule='errera', pick='size', beta=4.0, relax=True,
                 t0=0.6, tau=1.0, K=20.0, nsteps=30, K_dir=36, dmin=0.15):
        self.rng = np.random.default_rng(seed)
        self.rule, self.pick, self.beta, self.relax = rule, pick, beta, relax
        self.t0, self.tau, self.K, self.nsteps, self.K_dir, self.dmin = t0, tau, K, nsteps, K_dir, dmin
        self.pos = [np.array(p, float) for p in [(0, 0), (1, 0), (1, 1), (0, 1)]]
        self.cells = [[0, 1, 2, 3]]
        self.A0 = [1.0]
        self.vcells = {0: {0}, 1: {0}, 2: {0}, 3: {0}}
        self.eline = {}                       # frozenset(u,v) -> line id (-1 boundary)
        for i in range(4):
            self.eline[frozenset((i, (i + 1) % 4))] = -1
        self.birth = []                       # line id -> time
        self.hosts = []                       # line id -> [host, host]
        self.ends = []                        # line id -> [vertex a, vertex b]
        self.t = 0.0
        self.corner_cuts = 0
        self.fails = 0

    # ---------------------------------------------------------------- geometry
    def P(self, c):
        return np.array([self.pos[v] for v in self.cells[c]])

    @staticmethod
    def area_centroid(P):
        x, y = P[:, 0], P[:, 1]; x1, y1 = np.roll(x, -1), np.roll(y, -1)
        cr = x * y1 - x1 * y; A = 0.5 * cr.sum()
        return A, np.array([((x + x1) * cr).sum(), ((y + y1) * cr).sum()]) / (6 * A)

    def crossings(self, P, p, d):
        """intersections of line p + s d with polygon edges: list of (s, edge index, u)."""
        Q = np.roll(P, -1, 0); e = Q - P; w = P - p
        den = d[0] * e[:, 1] - d[1] * e[:, 0]
        with np.errstate(divide='ignore', invalid='ignore'):
            s = (w[:, 0] * e[:, 1] - w[:, 1] * e[:, 0]) / den
            u = (w[:, 0] * d[1] - w[:, 1] * d[0]) / den
        ok = np.isfinite(s) & (u >= 0) & (u < 1)
        return s, u, ok

    def plan(self, c):
        P = self.P(c); n = len(P)
        A, C = self.area_centroid(P)
        if self.rule == 'sr':
            for _ in range(50):
                lo, hi = P.min(0), P.max(0)
                p = lo + self.rng.random(2) * (hi - lo)
                if self._inside(P, p):
                    break
            th = self.rng.random() * math.pi
            cands = [(p, th)]
        else:
            off = self.rng.random() * math.pi / self.K_dir
            cands = [(C, off + k * math.pi / self.K_dir) for k in range(self.K_dir)]
        best = None
        for p, th in cands:
            d = np.array([math.cos(th), math.sin(th)])
            s, u, ok = self.crossings(P, p, d)
            sp = np.where(ok & (s > 0), s, np.inf); sn = np.where(ok & (s < 0), s, -np.inf)
            i = int(np.argmin(sp)); j = int(np.argmax(sn))
            if not np.isfinite(sp[i]) or not np.isfinite(sn[j]) or i == j:
                continue
            L = sp[i] - sn[j]
            if self.rule != 'sr':
                if min(u[i], 1 - u[i]) < self.dmin or min(u[j], 1 - u[j]) < self.dmin:
                    continue
                adjacent = (j == (i + 1) % n) or (i == (j + 1) % n)
                if self.rule == 'errera' and adjacent:
                    continue
            if best is None or L < best[0]:
                best = (L, i, u[i], j, u[j])
        return best

    @staticmethod
    def _inside(P, p):
        x, y = P[:, 0], P[:, 1]; x1, y1 = np.roll(x, -1), np.roll(y, -1)
        c = ((y > p[1]) != (y1 > p[1])) & (p[0] < (x1 - x) * (p[1] - y) / (y1 - y + 1e-300) + x)
        return c.sum() % 2 == 1

    # ---------------------------------------------------------------- topology
    def _split_edge(self, c, k, u):
        """insert a vertex on edge k of cell c at fraction u; returns vertex id."""
        vs = self.cells[c]; a, b = vs[k], vs[(k + 1) % len(vs)]
        x = (1 - u) * self.pos[a] + u * self.pos[b]
        v = len(self.pos); self.pos.append(x)
        key = frozenset((a, b)); ln = self.eline.pop(key)
        self.eline[frozenset((a, v))] = ln; self.eline[frozenset((v, b))] = ln
        self.vcells[v] = set()
        for cc in self.vcells[a] & self.vcells[b]:
            L = self.cells[cc]; m = len(L)
            for q in range(m):
                if {L[q], L[(q + 1) % m]} == {a, b}:
                    L.insert(q + 1, v); break
            self.vcells[v].add(cc)
        return v, ln

    def divide(self, c):
        pl = self.plan(c)
        if pl is None:
            self.fails += 1
            return False
        _, i, ui, j, uj = pl
        vs = self.cells[c]; n = len(vs)
        if (j == (i + 1) % n) or (i == (j + 1) % n):
            self.corner_cuts += 1
        ei = (vs[i], vs[(i + 1) % n]); ej = (vs[j], vs[(j + 1) % n])
        va, hA = self._split_edge(c, i, ui)
        # locate edge j again (indices shifted)
        vs = self.cells[c]
        k = [q for q in range(len(vs)) if {vs[q], vs[(q + 1) % len(vs)]} == set(ej)][0]
        vb, hB = self._split_edge(c, k, uj)
        vs = self.cells[c]; m = len(vs)
        ia, ib = vs.index(va), vs.index(vb)
        c1 = [vs[(ia + q) % m] for q in range((ib - ia) % m + 1)]
        c2 = [vs[(ib + q) % m] for q in range((ia - ib) % m + 1)]
        A1 = abs(self.area_centroid(np.array([self.pos[v] for v in c1]))[0])
        A2 = abs(self.area_centroid(np.array([self.pos[v] for v in c2]))[0])
        a0 = self.A0[c]
        self.cells[c] = c1; self.A0[c] = a0 * A1 / (A1 + A2)
        new = len(self.cells); self.cells.append(c2); self.A0.append(a0 * A2 / (A1 + A2))
        for v in c2:
            self.vcells[v].add(new)
            if v not in (va, vb):
                self.vcells[v].discard(c)
        lid = len(self.birth)
        self.eline[frozenset((va, vb))] = lid
        self.birth.append(self.t); self.hosts.append([hA, hB]); self.ends.append([va, vb])
        return True

    # ---------------------------------------------------------------- mechanics
    def _arrays(self):
        E = list(self.eline.items())
        eu = np.array([tuple(k)[0] for k, _ in E]); ev = np.array([tuple(k)[1] for k, _ in E])
        el = np.array([l for _, l in E])
        bt = np.array([self.birth[l] if l >= 0 else -1e9 for l in el])
        cv = np.concatenate([np.array(c) for c in self.cells])
        cn = np.array([len(c) for c in self.cells]); cs = np.r_[0, np.cumsum(cn)[:-1]]
        cid = np.repeat(np.arange(len(self.cells)), cn)
        nxt = cv[(np.arange(len(cv)) - cs[cid] + 1) % cn[cid] + cs[cid]]
        prv = cv[(np.arange(len(cv)) - cs[cid] - 1) % cn[cid] + cs[cid]]
        return eu, ev, el, bt, cv, cid, nxt, prv

    def relax_steps(self, n):
        X = np.array(self.pos)
        eu, ev, el, bt, cv, cid, nxt, prv = self._arrays()
        age = self.t - bt
        T = np.where(el < 0, 1.0, 1 - (1 - self.t0) * np.exp(-age / self.tau))
        A0 = np.array(self.A0)
        # boundary vertices: on x = 0/1 or y = 0/1
        onx = (np.abs(X[:, 0]) < 1e-9) | (np.abs(X[:, 0] - 1) < 1e-9)
        ony = (np.abs(X[:, 1]) < 1e-9) | (np.abs(X[:, 1] - 1) < 1e-9)
        lmean = math.sqrt(1.0 / len(self.cells))
        eta = 0.02 * lmean
        for _ in range(n):
            F = np.zeros_like(X)
            d = X[ev] - X[eu]; L = np.linalg.norm(d, axis=1) + 1e-15
            f = (T / L)[:, None] * d
            np.add.at(F, eu, f); np.add.at(F, ev, -f)
            x, y = X[cv, 0], X[cv, 1]
            A = 0.5 * np.bincount(cid, weights=x * X[nxt, 1] - X[nxt, 0] * y)
            pr = -self.K * (A - A0) / A0 ** 1.5            # -dE/dA
            gx = 0.5 * (X[nxt, 1] - X[prv, 1]); gy = 0.5 * (X[prv, 0] - X[nxt, 0])
            np.add.at(F[:, 0], cv, pr[cid] * gx); np.add.at(F[:, 1], cv, pr[cid] * gy)
            F[onx, 0] = 0; F[ony, 1] = 0
            step = eta * F
            sm = np.linalg.norm(step, axis=1); cap = 0.1 * lmean
            step *= np.minimum(1, cap / (sm + 1e-15))[:, None]
            X += step
        for i in range(len(self.pos)):
            self.pos[i] = X[i]

    # ---------------------------------------------------------------- driver
    def choose(self):
        A = np.array([abs(self.area_centroid(self.P(c))[0]) for c in range(len(self.cells))])
        w = A if self.pick == 'area' else A ** self.beta
        return int(self.rng.choice(len(A), p=w / w.sum()))

    def grow(self, N, record_every=None):
        snaps = []
        while len(self.cells) < N:
            c = self.choose()
            if self.divide(c):
                self.t += 1.0 / len(self.cells)
                if self.relax:
                    self.relax_steps(self.nsteps)
        return self

    # ---------------------------------------------------------------- observables
    def sides(self, interior=True):
        """topological sides (= number of vertices, all trivalent) of interior cells."""
        out = []
        for c, vs in enumerate(self.cells):
            if interior and any(self.eline.get(frozenset((vs[q], vs[(q + 1) % len(vs)]))) == -1
                                for q in range(len(vs))):
                continue
            out.append(len(vs))
        return np.array(out)

    def junction_angles(self):
        """for each end of each line: (age, through-angle between the two host pieces)."""
        nb = {}
        for k, l in self.eline.items():
            u, v = tuple(k)
            nb.setdefault(u, []).append((v, l)); nb.setdefault(v, []).append((u, l))
        out = []
        for lid, (ends, hs) in enumerate(zip(self.ends, self.hosts)):
            for v, h in zip(ends, hs):
                if h < 0:
                    continue
                arms = [w for w, l in nb.get(v, []) if l == h]
                if len(arms) != 2:
                    continue
                a = self.pos[arms[0]] - self.pos[v]; b = self.pos[arms[1]] - self.pos[v]
                ang = math.degrees(math.acos(np.clip(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1)))
                out.append((self.t - self.birth[lid], ang, self.birth[h] if h >= 0 else None))
        return out

    def dual(self):
        n = len(self.birth); k = np.zeros(n, int); adj = [set() for _ in range(n)]
        for lid, hs in enumerate(self.hosts):
            for h in hs:
                if h >= 0:
                    k[lid] += 1; k[h] += 1; adj[lid].add(h); adj[h].add(lid)
        C = []
        for i in range(n):
            m = list(adj[i])
            if len(m) >= 2:
                t = sum(1 for x in range(len(m)) for y in range(x + 1, len(m)) if m[y] in adj[m[x]])
                C.append((len(m), 2 * t / (len(m) * (len(m) - 1))))
        return k, C
