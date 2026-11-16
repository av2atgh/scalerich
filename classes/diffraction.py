"""Diffraction of the junction set: S(q) = |sum_j exp(i q.x_j)|^2 / N over distinct vertices
(ligament endpoints + square corners). Bragg test: the brightest peak I_max = max S(q)/N in the
annulus 2*pi*sqrt(T)*[0.25, 2] stays O(1) for pure-point order and falls as ~ln N / N otherwise.
FFT on a fine grid finds candidates; each is refined by direct summation."""
import sys, json, numpy as np
from scipy.optimize import minimize
sys.path[:0] = ['../beyond']
import qdiv, divisions

def vertices(g):
    P = np.vstack([g['A'], g['B'], [[0, 0], [1, 0], [1, 1], [0, 1]]])
    return np.unique(np.round(P, 10), axis=0)

def S_direct(P, q):
    return abs(np.exp(1j * (P @ q)).sum()) ** 2 / len(P)

def bragg(P, T, ncand=30):
    N = len(P); M = int(2 ** np.ceil(np.log2(24 * np.sqrt(T))))
    H = np.zeros((M, M)); ix = np.minimum((P * M).astype(int), M - 1)
    np.add.at(H, (ix[:, 0], ix[:, 1]), 1.0)
    F = np.abs(np.fft.fft2(H)) ** 2 / N
    k = np.fft.fftfreq(M, 1.0 / M)                       # q = 2 pi k
    KX, KY = np.meshgrid(k, k, indexing='ij'); K = np.hypot(KX, KY)
    ann = (K > 0.25 * np.sqrt(T)) & (K < 2 * np.sqrt(T))
    F[~ann] = 0
    idx = np.argsort(F.ravel())[-ncand:]
    best = 0.0; qbest = None
    for i in idx:
        q0 = 2 * np.pi * np.array([KX.ravel()[i], KY.ravel()[i]])
        res = minimize(lambda q: -S_direct(P, q), q0, method='Nelder-Mead',
                       options=dict(xatol=1e-3, fatol=1e-6, initial_simplex=q0 + 2 * np.pi * np.array([[0, 0], [1, 0], [0, 1]])))
        if -res.fun > best:
            best = -res.fun; qbest = res.x
    return best / N, (qbest / (2 * np.pi * np.sqrt(T))).tolist()

RULES = {
    'regular': ('largest', 0.5, 'perp'),
    'golden': ('largest', 1 / qdiv.TAU ** 2, 'perp'),
    'plastic': ('largest', qdiv.PLASTIC ** -3, 'perp'),
    'r0.3': ('largest', 0.3, 'perp'),
    'golden-side': ('largest', ('side', 1 / qdiv.TAU ** 2), 'perp'),
    'largest-uniform': ('largest', 'uniform', 'perp'),
}

if __name__ == '__main__':
    out = []
    for name, rule in list(RULES.items()) + [('SR', None), ('balanced', None)]:
        for T in (1023, 4095, 16383):
            g = qdiv.run(T, rule, seed=1) if rule else divisions.run(T, seed=1, rule=name)
            P = vertices(g); I, q = bragg(P, T)
            out.append(dict(rule=name, T=T, N=len(P), Imax=I, q=q))
            print(name, T, len(P), round(I, 4), np.round(q, 3), flush=True)
    json.dump(out, open('results/bragg.json', 'w'), indent=1)
