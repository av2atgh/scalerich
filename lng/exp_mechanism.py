"""Referee test (2D): reinforcement (preferential attachment) vs length x age.
(a) conditional Fano factor of k-2 at fixed predicted mean mu_i = 4 l_i (sqrt T - sqrt i)/c_L;
(b) reinforcement: do hits in [i,2i] predict hits in [2i,T] after controlling for l_i?
Positive control: topological surrogate (incidence PA, no lengths).
(c) tail amplitude K^2 P(k>=K) at T = 1e6 vs 6 (BA) and 4<x^2> (length picture).
Output: results/mechanism.txt"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, surrogate
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
f = open(os.path.join(OUT, 'mechanism.txt'), 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); f.write(s + '\n'); f.flush()

def hit_times(h, T):
    """list of (time, target) for every hit on an existing ligament"""
    t = np.repeat(np.arange(T), 2); tg = h.ravel()
    m = tg >= 0
    return t[m], tg[m]

def analyse(name, h, T, L=None):
    t, tg = hit_times(h, T)
    k = np.bincount(tg, minlength=T) + 2
    i = np.arange(T)
    # reinforcement: early = hits in [i, 2i), late = hits in [2i, T)
    early = np.zeros(T); late = np.zeros(T)
    e = t < 2 * tg
    np.add.at(early, tg[e], 1); np.add.at(late, tg[~e], 1)
    sel = (i >= 100) & (i < T // 20)
    res = {}
    if L is not None:
        mu = 4 * L / 2.33 * (np.sqrt(T) - np.sqrt(i + 1))
        y = k - 2.0
        # Fano factor in bins of mu
        out = []
        for lo, hi in [(2, 4), (4, 8), (8, 16), (16, 32), (32, 64), (64, 128)]:
            s = (mu >= lo) & (mu < hi) & (i > 30)
            if s.sum() > 50:
                out.append(f"mu[{lo},{hi}) n={s.sum()} mean(k-2)/mu={y[s].mean()/mu[s].mean():.3f} Fano={y[s].var()/y[s].mean():.2f}")
        P(name, 'conditional Fano of k-2 at fixed predicted mean:'); [P('   ', o) for o in out]
        # residual Fano: var(y - mu)/mean(mu) overall for mid-age ligaments
        s = sel
        P(name, f'  overall var(k-2-mu)/<mu> (i in [100,T/20)) = {np.var(y[s]-mu[s])/mu[s].mean():.2f}')
        # reinforcement controlling for length: regress late on [log-age expectation pieces]
        # expected early and late under length picture: rate 2 l_i /(c_L sqrt t)
        ee = 4 * L / 2.33 * (np.sqrt(2 * (i + 1)) - np.sqrt(i + 1))
        el = 4 * L / 2.33 * (np.sqrt(T) - np.sqrt(2 * (i + 1)))
        re, rl = early[s] - ee[s], late[s] - el[s]
        P(name, f'  corr(early,late) raw = {np.corrcoef(early[s], late[s])[0,1]:.3f}; '
                f'after removing length*age expectation = {np.corrcoef(re, rl)[0,1]:.3f}')
        # within narrow length bins and narrow age bins
        cs = []
        for a0 in (100, 300, 1000, 3000):
            for q in range(5):
                ss = (i >= a0) & (i < 1.3 * a0)
                if ss.sum() < 50: continue
                Lq = np.quantile(L[ss], [q / 5, (q + 1) / 5])
                s2 = ss & (L >= Lq[0]) & (L <= Lq[1])
                if s2.sum() > 30 and early[s2].std() > 0 and late[s2].std() > 0:
                    cs.append(np.corrcoef(early[s2], late[s2])[0, 1])
        P(name, f'  corr(early,late) within (age, length-quintile) cells: mean {np.mean(cs):.3f} sd {np.std(cs):.3f} ({len(cs)} cells)')
    else:
        s = sel
        P(name, f'  corr(early,late) raw = {np.corrcoef(early[s], late[s])[0,1]:.3f}')
        cs = []
        for a0 in (100, 300, 1000, 3000):
            ss = (i >= a0) & (i < 1.3 * a0)
            if ss.sum() > 30: cs.append(np.corrcoef(early[ss], late[ss])[0, 1])
        P(name, f'  corr(early,late) within age cells: mean {np.mean(cs):.3f}')
        # Fano at fixed age
        for a0 in (100, 1000):
            ss = (i >= a0) & (i < 1.1 * a0); y = k[ss] - 2.0
            P(name, f'  Fano of k-2 at age~{a0}: {y.var()/y.mean():.2f} (mean {y.mean():.1f})')
    return k

T = 100000
for s in (1, 2, 3):
    o = srcells.run(T, seed=s)
    k = analyse(f'SR seed{s}', o['hit'], T, o['length'])
    i = np.arange(T)
    for a0 in (100, 1000):
        ss = (i >= a0) & (i < 1.1 * a0); y = k[ss] - 2.0
        P(f'SR seed{s}', f'  Fano of k-2 at age~{a0} (not conditioning on length): {y.var()/y.mean():.2f} (mean {y.mean():.1f})')
for s in (1, 2):
    su = surrogate.run(T, face='uniform', side='padj', padj=0.29, seed=s)
    analyse(f'surrogate seed{s}', su['hit'], T)

# tail amplitude at T=1e6
T = 1000000
o = srcells.run(T, seed=11)
t, tg = hit_times(o['hit'], T); k = np.bincount(tg, minlength=T) + 2
L = o['length']; i = np.arange(T); x = L * np.sqrt(i + 1) / 1.163
m = i > 1000
P(f'T=1e6: <x^2> (birth-length scatter, i>1000) = {np.mean(x[m]**2):.3f} -> 4<x^2> = {4*np.mean(x[m]**2):.2f}; BA = 6')
for K in (8, 16, 32, 64, 128, 256):
    P(f'  K={K}: K^2 P(k>=K) = {K*K*(k>=K).mean():.2f}   BA exact 6K/(K+1) = {6*K/(K+1):.2f}')
f.close()
