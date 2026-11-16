"""1D r-splitting (low side gets r): refine every interval longer than eps; endpoints set X.
Axis diffraction of X at N up to ~1e7: I_max = max |sum exp(2 pi i k x)|^2 / N^2 over
k in [0.25, 2] N (window follows the mean spacing). -> results/oned.json"""
import json, numpy as np
from pisot import theta

def split_points(r, eps):
    pts = [0.0, 1.0]; stack = [(0.0, 1.0)]
    while stack:
        a, L = stack.pop()
        if L <= eps: continue
        c = a + r * L; pts.append(c); stack += [(a, r * L), (c, (1 - r) * L)]
    return np.sort(np.array(pts))

def imax(x, pad=2):
    N = len(x); M = 1 << int(np.ceil(np.log2(16 * N)))
    H = np.zeros(M * pad); np.add.at(H, np.minimum((x * M).astype(int), M - 1), 1.0)
    F = np.abs(np.fft.rfft(H)) ** 2 / N ** 2; k = np.arange(len(F)) / pad
    msk = (k > 0.25 * N) & (k < 2 * N); i = np.argmax(np.where(msk, F, 0)); return float(F[i]), float(k[i])

if __name__ == '__main__':
    out = []
    cases = {'(2,1)': theta(2, 1)[0] ** -2, '(3,1)': theta(3, 1)[0] ** -3, '(4,1)': theta(4, 1)[0] ** -4,
             '(3,2)': theta(3, 2)[0] ** -3, '(5,1)': theta(3, 2)[0] ** -5, '(4,3)': theta(4, 3)[0] ** -4,
             '(5,3)': theta(5, 3)[0] ** -5, 'r=0.3': 0.3}
    for nm, r in cases.items():
        for eps in (1e-3, 1e-4, 1e-5, 1e-6):
            x = split_points(r, eps); I, k = imax(x)
            out.append(dict(case=nm, r=r, eps=eps, N=len(x), Imax=I, k=k))
            print(nm, round(r, 4), eps, len(x), round(I, 4), k, flush=True)
    json.dump(out, open('results/oned.json', 'w'), indent=1)
