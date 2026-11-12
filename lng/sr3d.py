"""3D SR model with plates (zero thickness): iterated division of convex polyhedra.

Step: pick a cell with probability proportional to its volume, a uniform point in
it and a uniform (isotropic) plane normal; the plate is the section of the cell by
that plane; it ends on every face it crosses (T-junctions along lines).  Dual
graph: plate = node, contact = edge.  Faces carry the plate label (-1 = cube).

Prediction (face-incidence linear preferential attachment): if a newcomer hits
m plates (sides of its section polygon), incidences grow by 2 + m per step and
degree by 2m, so beta = <m>/(2 + <m>) and gamma = 1 + 1/beta = 2 + 2/<m>.
In 2D m = 2 exactly -> gamma = 3.
Output: results/sr3d.txt
"""
import sys, os, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srcells import Fenwick
import graphstats as gs

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')


def cube():
    V = [np.array(v, float) for v in [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]]
    F = [[0,3,2,1],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,4,7,3]]
    return [[V[i] for i in f] for f in F], [-1] * 6


def poly_vol(faces):
    c = np.mean([v for f in faces for v in f], axis=0)
    vol = 0.0
    tets = []
    for f in faces:
        a = f[0]
        for j in range(1, len(f) - 1):
            v = abs(np.dot(a - c, np.cross(f[j] - c, f[j + 1] - c))) / 6.0
            vol += v
            tets.append((v, c, a, f[j], f[j + 1]))
    return vol, tets


def sample_tet(tets, vol, rng):
    u = rng.random() * vol
    for v, a, b, c, d in tets:
        if u <= v:
            break
        u -= v
    s, t, w = rng.random(3)
    if s + t > 1: s, t = 1 - s, 1 - t
    if t + w > 1: t, w = 1 - w, 1 - s - t
    elif s + t + w > 1: s, w = 1 - t - w, s + t + w - 1
    return a + s * (b - a) + t * (c - a) + w * (d - a)


def clip_poly(faces, labs, n, c, newlab):
    """keep part with n.x >= c. faces: list of vertex lists (ordered)."""
    outF, outL, cut = [], [], []
    hitlabs = set()
    for f, L in zip(faces, labs):
        d = [np.dot(n, v) - c for v in f]
        nf = []
        m = len(f)
        crossed = False
        for i in range(m):
            j = (i + 1) % m
            if d[i] >= 0:
                nf.append(f[i])
            if (d[i] >= 0) != (d[j] >= 0):
                s = d[i] / (d[i] - d[j])
                p = f[i] + s * (f[j] - f[i])
                nf.append(p); cut.append(p); crossed = True
        if crossed:
            hitlabs.add(L)
        if len(nf) >= 3:
            outF.append(nf); outL.append(L)
    if len(cut) < 3:
        return None, None, hitlabs
    # order cut points around their centroid in the plane
    P = np.array(cut)
    cen = P.mean(0)
    e1 = P[0] - cen; e1 /= np.linalg.norm(e1) + 1e-300
    e2 = np.cross(n, e1)
    ang = np.arctan2((P - cen) @ e2, (P - cen) @ e1)
    order = np.argsort(ang)
    uniq = []
    for i in order:
        if not uniq or np.linalg.norm(P[i] - uniq[-1]) > 1e-12:
            uniq.append(P[i])
    if len(uniq) >= 2 and np.linalg.norm(uniq[0] - uniq[-1]) < 1e-12:
        uniq.pop()
    if len(uniq) >= 3:
        outF.append(uniq); outL.append(newlab)
    return outF, outL, hitlabs


def run(T, seed=0):
    rng = np.random.default_rng(seed)
    F0, L0 = cube()
    cells = [(F0, L0)]
    vol0, tets0 = poly_vol(F0)
    tets = [tets0]
    vols = [vol0]
    fw = Fenwick(T + 2); fw.set(0, vol0)
    edges = []
    mhit = np.zeros(T, int)
    for t in range(T):
        c = fw.find(rng.random() * fw.tot)
        F, L = cells[c]
        p = sample_tet(tets[c], vols[c], rng)
        n = rng.normal(size=3); n /= np.linalg.norm(n)
        cc = float(np.dot(n, p))
        F1, L1, hl = clip_poly(F, L, n, cc, t)
        F2, L2, _ = clip_poly(F, L, -n, -cc, t)
        if F1 is None or F2 is None:
            continue
        plates = [h for h in hl if h >= 0]
        mhit[t] = len(hl)
        for h in plates:
            edges.append((t, h))
        for (Fi, Li), idx in (((F1, L1), c), ((F2, L2), len(cells))):
            v, tt = poly_vol(Fi)
            if idx == len(cells):
                cells.append((Fi, Li)); tets.append(tt); vols.append(v)
            else:
                cells[idx] = (Fi, Li); tets[idx] = tt; vols[idx] = v
            fw.set(idx, v)
    return np.array(edges), mhit


if __name__ == '__main__':
    lines = []
    for T in (5000, 20000):
        for seed in (1, 2):
            t0 = time.time()
            E, mhit = run(T, seed)
            A = gs.adjacency(E, T)
            k = gs.degrees(A)
            # interior plates (all hits on plates) : m = sides of section polygon
            mm = mhit[T // 10:].mean()
            ii = np.unique(np.logspace(np.log10(20), np.log10(T / 10), 12).astype(int))
            mk = [k[a:int(a * 1.3) + 1].mean() for a in ii]
            beta = -np.polyfit(np.log(ii), np.log(mk), 1)[0]
            C, _ = gs.local_clustering(A, k)
            line = (f"T={T} seed={seed}: <k>={k.mean():.3f} <m> (sides of new plate, late 90%)={mm:.3f} "
                    f"pred beta=m/(2+m)={mm/(2+mm):.3f} gamma=2+2/m={2+2/mm:.3f} | measured beta={beta:.3f} "
                    f"gamma_MLE(k>=10)={gs.gamma_mle(k,10)[0]:.3f} (k>=20)={gs.gamma_mle(k,20)[0]:.3f} kmax={k.max()} C={C[k>=2].mean():.3f} [{time.time()-t0:.0f}s]")
            print(line, flush=True); lines.append(line)
    open(os.path.join(OUT, 'sr3d.txt'), 'w').write('\n'.join(lines) + '\n')
