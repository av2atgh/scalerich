"""Goal 4: triangles, clustering spectrum C(k), k_nn(k) for 3D plates, and the
same statistics for the 3D random Apollonian network (RAN: pick a tetrahedron
uniformly, join a new node to its 4 vertices, replace it by 4 tetrahedra).

Triangle bookkeeping: a newcomer's section polygon is a cycle of faces; two
consecutive faces meet along an edge of the cell, so their plates are in
contact.  A newcomer with m plate sides closes one triangle per consecutive
plate pair (m if no cube face), plus any non-consecutive hit pairs already
linked.  When plate i is hit, its two polygon neighbours are linked to i, so i
gains ~2 triangles per received link:  e_i = e0_i + q recv_i,
C(k) = 2[e0 + q(k - o)]/[k(k-1)] -> 2q/k.

usage: python3 exp_cluster.py run [run ...] -> results/cluster_<run>.txt ; 'ran' for the RAN.
"""
import sys, os
import numpy as np
import scipy.sparse as sp
from common import load, plate_arrays, gamma_mle, cohort_beta, RES


def adjacency(src, dst, N):
    A = sp.coo_matrix((np.ones(len(src)), (src, dst)), shape=(N, N))
    A = ((A + A.T) > 0).astype(np.float64).tocsr()
    A.setdiag(0); A.eliminate_zeros()
    return A


def graph_stats(A, out, k=None):
    k = np.asarray(A.sum(1)).ravel().astype(int) if k is None else k
    A2 = A @ A
    tri = np.asarray(A2.multiply(A).sum(1)).ravel() / 2
    C = np.zeros(len(k)); mk = k >= 2
    C[mk] = 2 * tri[mk] / (k[mk] * (k[mk] - 1))
    knn = np.zeros(len(k)); s = A @ k.astype(float); knn[k > 0] = s[k > 0] / k[k > 0]
    Ac = sp.triu(A).tocoo()
    x = k[Ac.row].astype(float); y = k[Ac.col].astype(float)
    r = np.corrcoef(np.r_[x, y], np.r_[y, x])[0, 1]
    out(f'  <k>={k[k > 0].mean():.4f}  <C>={C[mk].mean():.4f} (k>=2)  transitivity={tri.sum() / np.sum(k * (k - 1) / 2):.4f}  r={r:.4f}  kmax={k.max()}')
    edges = [3, 4, 5, 6, 8, 11, 16, 23, 32, 45, 64, 90, 128, 181, 256, 362, 512, 1024, 4096]
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        b = (k >= lo) & (k < hi)
        if b.sum() >= 5:
            rows.append((lo, hi, b.sum(), C[b].mean(), (C[b] * k[b]).mean(), knn[b].mean(), tri[b].mean()))
    out('  k-bin      n      C(k)    k*C(k)   knn(k)   triangles')
    for lo, hi, n, c, kc, kn, tr in rows:
        out(f'  [{lo},{hi})  {n:7d}  {c:.4f}  {kc:.3f}  {kn:8.2f}  {tr:9.1f}')
    return k, C, tri, knn


def plates(name):
    lines = []

    def out(s):
        print(s, flush=True); lines.append(s)

    S, E, Cc = load(name)
    P = plate_arrays(S, E)
    T = P['T']
    pe = E['h'] >= 0
    A = adjacency(E['t'][pe], E['h'][pe], T + 1)
    out(f'== {name}: T={T}')
    # newcomer closures (labels are written in polygon order)
    Et = E['t']; Eh = E['h']
    starts = np.r_[0, np.where(np.diff(Et) != 0)[0] + 1]; ends = np.r_[starts[1:], len(Et)]
    linked = set()
    cons = cons_link = noncons_link = noncons = 0
    closed = np.zeros(T + 1)
    for a, b in zip(starts, ends):
        t = int(Et[a]); L = Eh[a:b]; m = len(L)
        c_here = 0
        for i in range(m):
            for j in range(i + 1, m):
                u, v = L[i], L[j]
                if u < 0 or v < 0:
                    continue
                adj = (j == i + 1) or (i == 0 and j == m - 1)
                key = (u, v) if u < v else (v, u)
                isl = key in linked
                if adj:
                    cons += 1; cons_link += isl
                else:
                    noncons += 1; noncons_link += isl
                c_here += isl
        closed[t] = c_here
        for u in L:
            if u >= 0:
                linked.add((u, t))
    nt = len(starts)
    out(f'[closures] per newcomer: consecutive plate pairs {cons / nt:.4f} (linked fraction {cons_link / max(cons, 1):.5f}); '
        f'non-consecutive pairs {noncons / nt:.4f} (linked {noncons_link / max(noncons, 1):.4f}, i.e. {noncons_link / nt:.4f} per newcomer); '
        f'triangles closed per newcomer {closed.sum() / nt:.4f}')
    k, C, tri, knn = graph_stats(A, out, P['k'])
    # e_i = e0 + q recv
    recv = P['recv']; o = P['out']
    sel = P['interior'] & (recv >= 20)
    q = np.polyfit(recv[sel], tri[sel], 1)
    e0 = closed[P['interior']].mean()
    out(f'[e_i] linear fit triangles = {q[1]:.2f} + {q[0]:.4f} recv (recv>=20, interior); <e0> (closed at birth, interior) = {e0:.4f}; '
        f'C(k) -> 2q/k = {2 * q[0]:.3f}/k')
    mo = o[P['interior']].mean()
    for K in (8, 16, 32, 64, 128, 256):
        b = (k >= K / 1.2) & (k < K * 1.2) & P['interior']
        if b.sum() > 5:
            pred = 2 * (e0 + q[0] * (k[b] - o[b])) / (k[b] * (k[b] - 1))
            out(f'  C(k~{K}) measured {C[b].mean():.4f}  predicted {pred.mean():.4f}')
    out(f'[knn] hubs: <knn | k>=100> = {knn[k >= 100].mean():.2f}; lnT = {np.log(T):.2f}; T^(1/3)={T**(1/3):.1f}')
    open(os.path.join(RES, f'cluster_{name}.txt'), 'w').write('\n'.join(lines) + '\n')


def ran(N=200_000, seed=7):
    lines = []

    def out(s):
        print(s, flush=True); lines.append(s)

    rng = np.random.default_rng(seed)
    tets = np.zeros((3 * N + 1, 4), np.int64)
    tets[0] = [0, 1, 2, 3]; nt = 1
    src, dst = [0, 0, 0, 1, 1, 2], [1, 2, 3, 2, 3, 3]
    v = 4
    r = rng.integers(0, 1 << 62, size=N)
    for step in range(N - 4):
        j = int(r[step] % nt)
        a, b, c, d = tets[j]
        src += [v, v, v, v]; dst += [a, b, c, d]
        tets[j] = (v, b, c, d); tets[nt] = (a, v, c, d); tets[nt + 1] = (a, b, v, d); tets[nt + 2] = (a, b, c, v)
        nt += 3; v += 1
    A = adjacency(np.array(src), np.array(dst), N)
    out(f'== 3D random Apollonian network, N={N}')
    k, C, tri, knn = graph_stats(A, out)
    s = ' '.join(f'k>={q}:{gamma_mle(k, q)[0]:.3f}' for q in (20, 40, 80, 160, 320))
    out(f'[MLE] {s}; cohort beta [1e2,N/1e2] = {cohort_beta(k.astype(float), np.ones(N, bool), 100, N / 100):.4f} (2/3)')
    q = np.polyfit(k[k >= 20], tri[k >= 20], 1)
    out(f'[e_i] triangles = {q[1]:.2f} + {q[0]:.4f} k  (exact: 3 per new link -> C(k) ~ 6/k)')
    out(f'[knn] hubs <knn | k>=100> = {knn[k >= 100].mean():.2f}')
    open(os.path.join(RES, 'cluster_ran3d.txt'), 'w').write('\n'.join(lines) + '\n')


def knn_scaling(names):
    """k_nn(k) only (no A^2), to compare T: expected k_nn ~ A T^(1/3) k^(-1/2) + B."""
    lines = []
    for name in names:
        S, E, Cc = load(name)
        P = plate_arrays(S, E); T = P['T']; k = P['k']
        pe = E['h'] >= 0
        A = adjacency(E['t'][pe], E['h'][pe], T + 1)
        knn = np.zeros(T + 1); s = A @ k.astype(float); knn[k > 0] = s[k > 0] / k[k > 0]
        row = []
        for K in (4, 8, 16, 32, 64, 128, 256):
            b = (k >= K / 1.2) & (k < K * 1.2)
            row.append(knn[b].mean())
        Ac = sp.triu(A).tocoo(); x = k[Ac.row].astype(float); y = k[Ac.col].astype(float)
        r = np.corrcoef(np.r_[x, y], np.r_[y, x])[0, 1]
        l = (f'[knn] {name} T={T}: knn at k~4,8,16,32,64,128,256 = ' + ' '.join(f'{v:.1f}' for v in row)
             + f'; knn/T^(1/3) at k~4: {row[0] / T**(1/3):.3f}; r={r:.4f}')
        print(l, flush=True); lines.append(l)
    open(os.path.join(RES, 'knn_scaling.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'knn':
        knn_scaling(sys.argv[2:])
    for name in sys.argv[1:] if sys.argv[1] != 'knn' else []:
        ran() if name == 'ran' else plates(name)
