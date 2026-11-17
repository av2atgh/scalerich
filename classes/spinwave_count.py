"""IDOS by Sylvester inertia: N(omega) = number of modes below omega (counts positive pivots of
[[0,G],[G^T,0]] - blockdiag(Ku,Kw)/omega). Mode count in [w1, w2] relative to the Weyl count."""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import spinwave as sw, spinwave_scale as ss
a, b = 1.0, 4.0
def kof(om): return np.sqrt(max(0, (-(a + b) + np.sqrt((a - b) ** 2 + 4 * om ** 2)) / 2))
def counter(Ku, Kw, G):
    A = sp.block_diag([Ku, Kw]).tocsc(); Bb = sp.bmat([[None, G], [G.T, None]]).tocsc()
    def N(w):
        lu = spl.splu((Bb - A / w).tocsc(), permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0, options=dict(SymmetricMode=True))
        return int((lu.U.diagonal() > 0).sum())
    return N
if __name__ == '__main__':
    w1, w2 = float(sys.argv[3]), float(sys.argv[4])
    for nm in sys.argv[1].split(','):
        for T in [int(x) for x in sys.argv[2].split(',')]:
            X, S = sw.network(ss.build(nm, T, 0), T); Ku, Kw, G, info = sw.assemble(X, S); N = counter(Ku, Kw, G)
            n1, n2 = N(w1), N(w2); weyl = info['Ltot'] * (kof(w2) - kof(w1)) / np.pi
            print('%-16s T %6d  modes in [%.2f,%.2f]: %5d  Weyl %7.1f  ratio %.3f  N/(T+1) %.3f' % (nm, T, w1, w2, n2 - n1, weyl, (n2 - n1) / weyl, n1 / (T + 1)), flush=True)
