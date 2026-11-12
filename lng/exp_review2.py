"""Second-round referee items (2D).
(a) dispersion around each ligament's own prediction, D = <(k-2-mu)^2>/<mu> by mu bin, the mean ratio,
    and the binning artefact: Fano within a mu bin vs D. Reinforcement for the oldest ligaments.
(b) tail amplitude at T=1e6: measured vs Poisson(mu) vs Poisson(r(mu) mu) with the binned mean ratio r.
(c) birth-length scaling under orientation fields.
Output: results/review2.txt"""
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
f = open(os.path.join(OUT, 'review2.txt'), 'w')
def P(s): print(s, flush=True); f.write(s + '\n'); f.flush()
def prep(o, T):
    h = o['hit']; L = o['length']; i = np.arange(T)
    t = np.repeat(np.arange(T), 2); tg = h.ravel(); m = tg >= 0
    k = np.bincount(tg[m], minlength=T) + 2
    cL = L.sum() / np.sqrt(T)
    mu = 4 * L / cL * (np.sqrt(T) - np.sqrt(i + 1))
    return t[m], tg[m], k, mu, L, cL, i
bins = (2, 4, 8, 16, 32, 64)
T = 100000
for s in (1, 2, 3):
    o = srcells.run(T, seed=s); t, tg, k, mu, L, cL, i = prep(o, T); y = k - 2.0
    out = []
    for lo in bins:
        sel = (mu >= lo) & (mu < 2 * lo) & (i > 30)
        if sel.sum() > 50:
            out.append(f"mu[{lo},{2*lo}) n={sel.sum()}: ratio={y[sel].mean()/mu[sel].mean():.3f} D={np.mean((y[sel]-mu[sel])**2)/mu[sel].mean():.2f} Fano-in-bin={y[sel].var()/y[sel].mean():.2f}")
    P(f"seed {s} (c_L={cL:.3f})"); [P("   " + x) for x in out]
    early = np.zeros(T); late = np.zeros(T); e = t < 2 * tg
    np.add.at(early, tg[e], 1); np.add.at(late, tg[~e], 1)
    ee = 4 * L / cL * (np.sqrt(2 * (i + 1)) - np.sqrt(i + 1)); el = 4 * L / cL * (np.sqrt(T) - np.sqrt(2 * (i + 1)))
    for name, sel in (("mu>16", (mu > 16) & (i > 10)), ("mu in [4,16]", (mu >= 4) & (mu <= 16) & (i > 10))):
        re, rl = early[sel] - ee[sel], late[sel] - el[sel]
        P(f"   reinforcement {name}: n={sel.sum()} corr = {np.corrcoef(re, rl)[0,1]:.3f}")
T = 1000000; Ks = (32, 64, 128); rng = np.random.default_rng(3)
res = []
for s in (11, 12, 13):
    o = srcells.run(T, seed=s); t, tg, k, mu, L, cL, i = prep(o, T); y = k - 2.0
    edges = np.array([0, 2, 4, 8, 16, 32, 64, 128, 256, 1e9]); r = np.ones(len(edges) - 1)
    for j in range(len(edges) - 1):
        sel = (mu >= edges[j]) & (mu < edges[j + 1]) & (i > 0)
        if sel.sum() > 20: r[j] = y[sel].sum() / mu[sel].sum()
    rr = r[np.clip(np.searchsorted(edges, mu, 'right') - 1, 0, len(r) - 1)]
    k1 = 2 + rng.poisson(mu); k2 = 2 + rng.poisson(rr * mu)
    res.append([[K * K * (x >= K).mean() for K in Ks] for x in (k, k1, k2)])
    P(f"T=1e6 seed {s}: ratio by mu bin " + " ".join(f"{a:.3f}" for a in r))
res = np.array(res).mean(0)
for j, K in enumerate(Ks):
    P(f"K={K}: measured {res[0][j]:.3f}  Poisson(mu) {res[1][j]:.3f}  Poisson(r(mu) mu) {res[2][j]:.3f}  BA {6*K/(K+1):.3f}")
circ = lambda x, y: math.atan2(y - 0.5, x - 0.5) + math.pi / 2
wave = lambda x, y: math.pi * (x + y) / 0.125
for name, kw in (('uniform', {}), ('circ', dict(angle='field', field=circ)), ('wave_w0.125', dict(angle='field', field=wave))):
    o = srcells.run(100000, seed=1, **kw); L = o['length']
    ii = np.unique(np.logspace(1, 4.5, 15).astype(int))
    ml = [L[a:int(a * 1.3) + 1].mean() for a in ii]
    P(f"{name}: <l_i> ~ i^{np.polyfit(np.log(ii), np.log(ml), 1)[0]:.3f} (i=10..3e4)")
f.close()
