"""Loaders and small helpers for the d3 runs (raw data in data/)."""
import os
import numpy as np

D = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(D, 'data')
RES = os.path.join(D, 'results')
os.makedirs(RES, exist_ok=True)

CELL_DT = np.dtype([('nf', 'i4'), ('nv', 'i4'), ('v', 'f8'), ('tc', 'i4'), ('depth', 'i4')])
EDGE_DT = np.dtype([('t', 'i4'), ('h', 'i4'), ('fa', 'f4'), ('fp', 'f4'), ('vol', 'f4'), ('nf', 'i4')])
STEP_COLS = ['t', 'depth', 'nfp', 'm', 'ncube', 'area', 'perim', 'vparent', 'u', 'free']


def load(name):
    S = np.fromfile(os.path.join(DATA, name + '.steps')).reshape(-1, 10)
    E = np.fromfile(os.path.join(DATA, name + '.edges'), dtype=EDGE_DT)
    C = np.fromfile(os.path.join(DATA, name + '.cells'), dtype=CELL_DT)
    return S, E, C


def plate_arrays(S, E):
    """per-plate arrays indexed by plate label t (index 0 unused).
    Returns dict with area, perim, out (plate hits at birth), ncube, recv, k, interior, ok."""
    T = len(S)
    ok = np.zeros(T + 1, bool); ok[S[S[:, 0] > 0, 0].astype(int)] = True
    area = np.zeros(T + 1); perim = np.zeros(T + 1); m = np.zeros(T + 1, int); nc = np.zeros(T + 1, int)
    idx = S[:, 0].astype(int); g = idx > 0
    area[idx[g]] = S[g, 5]; perim[idx[g]] = S[g, 6]; m[idx[g]] = S[g, 3]; nc[idx[g]] = S[g, 4]
    P = E['h'] >= 0
    recv = np.bincount(E['h'][P], minlength=T + 1)
    out = m - nc
    k = out + recv
    return dict(area=area, perim=perim, m=m, ncube=nc, out=out, recv=recv, k=k,
                interior=ok & (nc == 0), ok=ok, T=T)


def gamma_mle(k, kmin):
    x = k[k >= kmin].astype(float)
    return 1 + len(x) / np.sum(np.log(x / (kmin - 0.5))), len(x)


def cohort_beta(k, sel, lo, hi, nb=14):
    """slope of log <k>_cohort vs log i over log-spaced cohorts in [lo, hi]."""
    edges = np.unique(np.logspace(np.log10(lo), np.log10(hi), nb + 1).astype(int))
    xs, ys = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        s = sel[a:b]
        if s.sum() > 5:
            xs.append(np.log(np.sqrt(a * b))); ys.append(np.log(k[a:b][s].mean()))
    return -np.polyfit(xs, ys, 1)[0]
