"""Graph statistics for the SR dual and its surrogates (sparse, single-threaded)."""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components


def adjacency(E, N):
    E = np.asarray(E)
    A = sp.coo_matrix((np.ones(len(E)), (E[:, 0], E[:, 1])), shape=(N, N))
    A = ((A + A.T) > 0).astype(np.float64).tocsr()
    A.setdiag(0)
    A.eliminate_zeros()
    return A


def degrees(A):
    return np.asarray(A.sum(1)).ravel().astype(int)


def triangles(A):
    A2 = A @ A
    t = np.asarray(A2.multiply(A).sum(1)).ravel() / 2.0
    return t


def local_clustering(A, k=None):
    if k is None:
        k = degrees(A)
    t = triangles(A)
    C = np.zeros(len(k))
    m = k >= 2
    C[m] = 2 * t[m] / (k[m] * (k[m] - 1))
    return C, t


def knn(A, k=None):
    if k is None:
        k = degrees(A)
    s = A @ k.astype(float)
    out = np.zeros(len(k))
    m = k > 0
    out[m] = s[m] / k[m]
    return out


def binned(x, y, edges):
    res = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi)
        if m.sum():
            res.append((lo, hi, int(m.sum()), float(np.mean(y[m]))))
    return res


def gamma_mle(k, kmin):
    """discrete power-law MLE (Clauset et al. approx) for k >= kmin."""
    x = k[k >= kmin].astype(float)
    return 1 + len(x) / np.sum(np.log(x / (kmin - 0.5))), len(x)


def ccdf_slope(k, kmin, kmax):
    ks = np.arange(kmin, kmax + 1)
    cc = np.array([(k >= v).mean() for v in ks])
    m = cc > 0
    s = np.polyfit(np.log(ks[m]), np.log(cc[m]), 1)[0]
    return 1 - s


def assortativity(A, k):
    Ac = sp.triu(A).tocoo()
    x = k[Ac.row].astype(float)
    y = k[Ac.col].astype(float)
    xx = np.concatenate([x, y])
    yy = np.concatenate([y, x])
    return np.corrcoef(xx, yy)[0, 1]


def bfs_counts(A, src):
    """distances and number of shortest paths from src (level-synchronous)."""
    N = A.shape[0]
    dist = np.full(N, -1, np.int32)
    sig = np.zeros(N)
    dist[src] = 0
    sig[src] = 1.0
    front = np.array([src])
    d = 0
    AT = A  # symmetric
    while len(front):
        v = np.zeros(N)
        v[front] = sig[front]
        nxt = AT @ v
        new = (nxt > 0) & (dist < 0)
        idx = np.nonzero(new)[0]
        d += 1
        dist[idx] = d
        sig[idx] = nxt[idx]
        front = idx
    return dist, sig


def path_stats(A, nsrc=100, seed=0, maxd=200):
    """mean distance, <mu>, <ln mu>, and per-distance <mu>, <ln mu>, counts."""
    rng = np.random.default_rng(seed)
    N = A.shape[0]
    ncomp, lab = connected_components(A, directed=False)
    big = np.bincount(lab).argmax()
    nodes = np.nonzero(lab == big)[0]
    srcs = rng.choice(nodes, size=min(nsrc, len(nodes)), replace=False)
    D = np.zeros(maxd)
    M = np.zeros(maxd)
    LM = np.zeros(maxd)
    totd = totm = totlm = cnt = 0.0
    for s in srcs:
        dist, sig = bfs_counts(A, s)
        m = dist > 0
        dd = dist[m]
        mu = sig[m]
        np.add.at(D, np.minimum(dd, maxd - 1), 1)
        np.add.at(M, np.minimum(dd, maxd - 1), mu)
        np.add.at(LM, np.minimum(dd, maxd - 1), np.log(mu))
        totd += dd.sum()
        totm += mu.sum()
        totlm += np.log(mu).sum()
        cnt += m.sum()
    per = [(d, int(D[d]), M[d] / D[d], LM[d] / D[d]) for d in range(maxd) if D[d] > 0]
    return dict(meand=totd / cnt, mu=totm / cnt, lnmu=totlm / cnt, per=per,
                giant=len(nodes) / N)


def configuration_rewire(A, k, seed=0, nswap_per_edge=10):
    """degree-preserving double-edge swaps (Maslov-Sneppen), simple graph kept."""
    rng = np.random.default_rng(seed)
    Ac = sp.triu(A).tocoo()
    u = Ac.row.copy()
    v = Ac.col.copy()
    m = len(u)
    S = set(zip(u.tolist(), v.tolist()))
    S |= set((b, a) for a, b in S)
    nsw = nswap_per_edge * m
    i1 = rng.integers(m, size=nsw)
    i2 = rng.integers(m, size=nsw)
    fl = rng.random(nsw) < 0.5
    acc = 0
    for a, b, f in zip(i1.tolist(), i2.tolist(), fl.tolist()):
        if a == b:
            continue
        x, y = u[a], v[a]
        w, z = (u[b], v[b]) if f else (v[b], u[b])
        if x == z or y == w or x == w or y == z:
            continue
        if (x, z) in S or (w, y) in S:
            continue
        S.discard((x, y)); S.discard((y, x)); S.discard((w, z)); S.discard((z, w))
        S.add((x, z)); S.add((z, x)); S.add((w, y)); S.add((y, w))
        u[a], v[a] = x, z
        u[b], v[b] = w, y
        acc += 1
    E = np.stack([u, v], 1)
    return E, acc / nsw


def summary(A, nsrc=60, seed=0, paths=True):
    k = degrees(A)
    C, t = local_clustering(A, k)
    kn = knn(A, k)
    m2 = k >= 2
    res = dict(N=A.shape[0], meank=k.mean(), kmax=int(k.max()),
               c2=(k * (k - 1.0)).mean() / k.mean(),
               Cmean=C[m2].mean(), transitivity=3 * t.sum() / max(1, (k * (k - 1) / 2).sum()),
               r=assortativity(A, k),
               gmle6=gamma_mle(k, 6)[0], gmle10=gamma_mle(k, 10)[0])
    res['Ck'] = binned(k, C, [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 128, 256, 100000])
    res['knnk'] = binned(k, kn, [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 128, 256, 100000])
    if paths:
        res.update(path_stats(A, nsrc=nsrc, seed=seed))
    return res, k, C, kn
