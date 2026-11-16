"""Commensurate splits r = theta^-p, 1 - r = theta^-q (theta root of x^p - x^(p-q) - 1):
finite shape set for every (p, q); Bragg peaks along the axis only for Pisot theta.
-> results/pisot.json"""
import json, numpy as np, qdiv
from diffraction import vertices
from math import gcd

def theta(p, q):
    c = np.zeros(p + 1); c[0] = 1; c[q] -= 1; c[p] -= 1
    r = np.roots(c); th = max(r.real[abs(r.imag) < 1e-9])
    return th, float(max(sorted(abs(r))[:-1]))

def axis_top(P, T, M=2 ** 20, pad=8):
    y = P[:, 1]; N = len(y)
    H = np.zeros(M * pad); np.add.at(H, np.minimum((y * M).astype(int), M - 1), 1.0)
    F = np.abs(np.fft.rfft(H)) ** 2 / N; k = np.arange(len(F)) / pad; s = np.sqrt(T)
    m = (k > 0.25 * s) & (k < 2 * s); i = np.argmax(np.where(m, F, 0)); return F[i] / N, k[i]

if __name__ == '__main__':
    out = []
    cases = [(1, 1), (2, 1), (3, 1), (3, 2), (4, 1), (4, 3), (5, 2), (5, 3), (7, 2)]
    for p, q in cases:
        th, mo = (2.0, 0.0) if p == q else theta(p, q); r = th ** -p
        for T in (2 ** 14 - 1, 2 ** 16 - 1, 2 ** 18 - 1, 2 ** 20 - 1):
            g = qdiv.run(T, ('largest', r, 'perp'))
            wh = np.sort(g['leaf_wh'], axis=1) * np.sqrt(T + 1)
            ns = len(set(map(tuple, np.round(wh, 6))))
            I, k = axis_top(vertices(g), T)
            out.append(dict(p=p, q=q, theta=th, r=r, maxconj=mo, T=T, shapes=ns, Imax=I, k=k))
            print(p, q, round(th, 4), round(r, 4), round(mo, 3), T, ns, round(I, 4), k, flush=True)
    json.dump(out, open('results/pisot.json', 'w'), indent=1)
