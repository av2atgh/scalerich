"""Does area-weighted cell selection route hits to high-degree ligaments?
At snapshot time t, for ligaments born in [t/40, t/4): current degree k_i(t),
birth length l_i, hits received in [t, 2t). Poisson regression
   E[hits] = exp(c0 + a ln l_i + b ln k_i(t) + d ln i).
Length-only hit rate (no degree link) => b = 0, a = 1.
Also: catchment area per unit length vs degree at the end of the run
(cells kept): area of adjacent cells weighted by the share of their perimeter
carried by ligament i, divided by l_i.
Usage: python exp_area_degree.py T seed q
Output: results/area_degree.txt (appended)"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells

T = int(sys.argv[1]); seed = int(sys.argv[2]); q = float(sys.argv[3])
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'area_degree.txt')
f = open(OUT, 'a')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); f.write(s + '\n'); f.flush()

def poisson_fit(X, y, it=50):
    beta = np.zeros(X.shape[1]); beta[0] = np.log(y.mean() + 1e-9)
    for _ in range(it):
        mu = np.exp(X @ beta)
        W = mu
        z = X @ beta + (y - mu) / mu
        XtW = X.T * W
        cov = np.linalg.inv(XtW @ X)
        nb = cov @ (XtW @ z)
        if np.max(np.abs(nb - beta)) < 1e-9:
            beta = nb; break
        beta = nb
    mu = np.exp(X @ beta)
    phi = np.sum((y - mu) ** 2 / mu) / (len(y) - X.shape[1])
    se = np.sqrt(np.diag(cov) * phi)
    return beta, se

o = srcells.run(T, seed=seed, q=q, keep_cells=True)
T = o['T']; h = o['hit']; L = o['length']
tt = np.repeat(np.arange(T), 2); tg = h.ravel(); m = tg >= 0
tt, tg = tt[m], tg[m]
P(f"== T={T} seed={seed} q={q}")
for t in (T // 8, T // 4, T // 2):
    pre = tt < t
    k = np.bincount(tg[pre], minlength=T)[:t] + 2
    win = (tt >= t) & (tt < 2 * t) if 2 * t <= T else (tt >= t)
    fut = np.bincount(tg[win], minlength=T)[:t].astype(float)
    i = np.arange(t)
    s = (i >= t // 40) & (i < t // 4) & (L[:t] > 0)
    X = np.column_stack([np.ones(s.sum()), np.log(L[:t][s]), np.log(k[s]), np.log(i[s] + 1)])
    b, se = poisson_fit(X, fut[s])
    # partial correlation of residual future hits with k at fixed (l, i)
    P(f"t={t} n={s.sum()}  a(ln l)={b[1]:.3f}±{se[1]:.3f}  b(ln k)={b[2]:.3f}±{se[2]:.3f}  d(ln i)={b[3]:.3f}±{se[3]:.3f}")

# catchment at end: each cell's area shared among its sides by side length
cellP, cellL, cellA = o['cellP'], o['cellL'], o['cellA']
catch = np.zeros(T)
for P_, L_, A_ in zip(cellP, cellL, cellA):
    if A_ <= 0: continue
    Pa = np.asarray(P_); n = len(Pa)
    seg = np.linalg.norm(np.roll(Pa, -1, axis=0) - Pa, axis=1)
    per = seg.sum()
    for lab, sl in zip(L_, seg):
        if lab >= 0:
            catch[lab] += A_ * sl / per
k = np.bincount(tg, minlength=T) + 2
i = np.arange(T)
s = (i >= 30) & (i < T // 4) & (L > 0)
r = np.log(catch[s] / L[s])
# within age-and-length bins: correlation of catchment per unit length with degree
ab = np.digitize(np.log(i[s] + 1), np.linspace(np.log(30), np.log(T // 4), 6))
lb = np.zeros(s.sum(), int)
for a_ in np.unique(ab):
    w = ab == a_
    lb[w] = np.digitize(np.log(L[s][w]), np.quantile(np.log(L[s][w]), [0.2, 0.4, 0.6, 0.8]))
cs = []
for a_ in np.unique(ab):
    for l_ in range(5):
        w = (ab == a_) & (lb == l_)
        if w.sum() > 30:
            cs.append((np.corrcoef(r[w], np.log(k[s][w]))[0, 1], w.sum()))
cs = np.array(cs)
P(f"catchment/length vs k at fixed (age,length): mean r={np.average(cs[:,0], weights=cs[:,1]):.3f} "
  f"(range {cs[:,0].min():.2f}..{cs[:,0].max():.2f}, {len(cs)} bins); "
  f"overall corr(ln catch/l, ln k)={np.corrcoef(r, np.log(k[s]))[0,1]:.3f}")
