"""Purely topological surrogates of the SR dual (no geometry).

State: a list of faces, each a cyclic list of side labels (dual nodes; -1 = outer
boundary).  Step t: pick a face f (rule `face`), pick two distinct side positions
i, j of f (rule `side`), add node t joined to the nodes on sides i and j, and
split the cycle:  child1 = f[i..j] + [t],  child2 = f[j..i] + [t]  (sides i, j
appear in both children: the hit ligament is cut into two portions).
Corner number bookkeeping: n1 + n2 = n + 4, exactly as in the T-tessellation.

face rules:  'uniform'  each face equally likely
             'nsides'   weight = number of sides
             'area'     each face carries an area; children split it with a
                        fraction u ~ Beta(a, a) (a = `abeta`); weight = area
side rules:  'uniform'  uniform pair of distinct positions
             'padj'     adjacent pair w.p. `padj` (always if n == 3), otherwise a
                        uniform non-adjacent pair
             'length'   each side-portion carries a length; pick two positions
                        with prob. proportional to length (without replacement);
                        portions are cut at a uniform point (needs face='area'
                        style bookkeeping of lengths; lengths are a mean-field
                        stand-in for geometry)
"""
import numpy as np
from srcells import Fenwick


def run(T, face='uniform', side='uniform', padj=0.45, abeta=1.0, seed=0):
    rng = np.random.default_rng(seed)
    faces = [[-1, -1, -1, -1]]
    flen = [[1.0, 1.0, 1.0, 1.0]]
    area = [1.0]
    fw = Fenwick(T + 2) if face != 'uniform' else None
    if fw is not None:
        fw.set(0, 4.0 if face == 'nsides' else 1.0)
    hit = np.full((T, 2), -2, np.int64)
    adj = np.zeros(T, np.int8)
    ncell = np.zeros(T, np.int32)
    for t in range(T):
        if fw is None:
            c = int(rng.integers(len(faces)))
        else:
            c = fw.find(rng.random() * fw.tot)
        F = faces[c]
        n = len(F)
        if side == 'uniform':
            i, j = rng.choice(n, 2, replace=False)
        elif side == 'padj':
            if n == 3 or rng.random() < padj:
                i = int(rng.integers(n)); j = (i + 1) % n
            else:
                while True:
                    i, j = rng.choice(n, 2, replace=False)
                    if (j - i) % n not in (1, n - 1):
                        break
        elif side == 'length':
            w = np.array(flen[c])
            i = rng.choice(n, p=w / w.sum())
            w2 = w.copy(); w2[i] = 0
            j = rng.choice(n, p=w2 / w2.sum())
        i, j = int(min(i, j)), int(max(i, j))
        d = j - i
        adj[t] = 1 if d in (1, n - 1) else 0
        ncell[t] = n
        hit[t] = (F[i], F[j])
        c1 = F[i:j + 1] + [t]
        c2 = F[j:] + F[:i + 1] + [t]
        Lf = flen[c]
        # lengths: cut sides i and j at uniform points; new side length ~ mean of the two
        ui, uj = rng.random(), rng.random()
        li1, li2 = Lf[i] * ui, Lf[i] * (1 - ui)
        lj1, lj2 = Lf[j] * uj, Lf[j] * (1 - uj)
        lnew = 0.5 * (sum(Lf) / n) * 2 if side == 'length' else 1.0
        l1 = [li2] + Lf[i + 1:j] + [lj1] + [lnew]
        l2 = [lj2] + Lf[j + 1:] + Lf[:i] + [li1] + [lnew]
        faces[c] = c1
        flen[c] = l1
        faces.append(c2)
        flen.append(l2)
        if fw is not None:
            if face == 'nsides':
                fw.set(c, float(len(c1)))
                fw.set(len(faces) - 1, float(len(c2)))
            else:
                u = rng.beta(abeta, abeta)
                a = area[c]
                area[c] = a * u
                area.append(a * (1 - u))
                fw.set(c, area[c])
                fw.set(len(faces) - 1, area[-1])
        else:
            area.append(0.0)
    return dict(T=T, hit=hit, adj=adj, ncell=ncell)


def run_dm(T, seed=0):
    """Dorogovtsev-Mendes / BB(n,1): new node joins both ends of a uniform random edge.
    Started from a triangle (nodes 0,1,2); returns the same dict layout (hit = the two ends)."""
    rng = np.random.default_rng(seed)
    eu = [0, 1, 0]
    ev = [1, 2, 2]
    hit = np.full((T, 2), -2, np.int64)
    hit[1] = (0, -1); hit[2] = (0, 1)
    for t in range(3, T):
        e = int(rng.integers(len(eu)))
        a, b = eu[e], ev[e]
        hit[t] = (a, b)
        eu += [t, t]; ev += [a, b]
    return dict(T=T, hit=hit, adj=np.ones(T, np.int8), ncell=np.full(T, 3))
