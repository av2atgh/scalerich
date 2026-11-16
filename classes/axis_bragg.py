"""Axis diffraction S(0, q) of the vertex set: 1D, binned on a fine grid (M = 2^22) and FFT.
I_max = max S/N over q in 2 pi [0.25, 2] sqrt(T). Larger T than the 2D search."""
import sys, json, numpy as np
sys.path[:0] = ['../beyond']
import qdiv
from diffraction import vertices, RULES

def axis_I(P, T, M=2 ** 22):
    y = P[:, 1]; N = len(y)
    H = np.bincount(np.minimum((y * M).astype(int), M - 1), minlength=M).astype(float)
    F = np.abs(np.fft.rfft(H)) ** 2 / N
    k = np.arange(len(F)); s = np.sqrt(T)
    m = (k > 0.25 * s) & (k < 2 * s)
    i = np.argmax(np.where(m, F, 0))
    return F[i] / N, int(k[i])

if __name__ == '__main__':
    out = []
    for name in ('golden', 'plastic', 'r0.3', 'golden-side', 'largest-uniform', 'regular'):
        for T in (2 ** 12 - 1, 2 ** 14 - 1, 2 ** 16 - 1, 2 ** 18 - 1):
            P = vertices(qdiv.run(T, RULES[name], seed=1))
            I, k = axis_I(P, T)
            out.append(dict(rule=name, T=T, N=len(P), Imax=I, k=k))
            print(name, T, len(P), round(I, 4), k, flush=True)
    json.dump(out, open('results/axis_bragg.json', 'w'), indent=1)
