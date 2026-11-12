"""Elastic experiments on SR beam frames (Goal 2c/3).

E1 'lam0':  fixed T, vary lam0 -> E(rho) at fixed architecture; local
            m_lam0 = dlnE/dlnrho (1 = stretching, 3 = bending).
E2 'flow':  fixed (alpha, lam0), truncations tau of one T realisation ->
            E(tau), rho(tau), elastic affine bound aff(tau) = sum lam l cos^4,
            window increments psi_el = dE/daff and the newest-chord slenderness
            s_tau = l/lam ~ (c_l/lam0) tau^(alpha-1/2).  RG prediction:
            psi_el collapses on s (marginal alpha-dependence only).
            Also m_T = dlnE/dlnrho along tau (the paper's density route if rho is
            varied through T) vs eq. (3) with p = 0: (1-alpha/2)/(1/2-alpha).
E3 'lig':   per-ligament energies at T=3000, rho=0.4: fit U_i ~ i^-x and
            U ~ k^a lam^b (paper: x = 3/4+alpha/2, a = 3/2, b = 1/2).
"""
import json, sys, time
import numpy as np
import srtree, frame


def rho_of(g, lam):
    return float(np.sum(lam * np.linalg.norm(g['B'] - g['A'], axis=1)))


def trunc(g, n):
    h = dict(T=n)
    for k in ('A', 'B', 'hostA', 'hostB'):
        h[k] = g[k][:n]
    return h


def Eboth(g, lam):
    nodes, seg, sl, bn = srtree.network(g)
    pack = frame.assemble(nodes, seg, lam[sl])
    Ex, Uax, Ubx = frame.stiffness(nodes, seg, lam[sl], bn, 'x', Kpack=pack)
    Ey, Uay, Uby = frame.stiffness(nodes, seg, lam[sl], bn, 'y', Kpack=pack)
    fa = (Uax.sum() + Uay.sum()) / (Uax.sum() + Ubx.sum() + Uay.sum() + Uby.sum())
    return Ex, Ey, fa, sl, (Uax + Ubx), (Uay + Uby)


def e1(Ts=(1000, 3000), NS=3, alphas=(0.0, 0.2, 0.4, 0.6, 0.8)):
    rt = np.logspace(np.log10(0.01), np.log10(0.6), 10)
    out = {}
    for T in Ts:
        for a in alphas:
            rows = []
            for s in range(NS):
                g = srtree.run(T, 200 + s)
                S = rho_of(g, srtree.thickness(T, 1.0, a))
                for r in rt:
                    lam = srtree.thickness(T, r / S, a)
                    try:
                        Ex, Ey, fa, *_ = Eboth(g, lam)
                    except Exception as e:
                        print('FAILED e1', T, a, s, r, e); continue
                    rows.append((s, r, r / S, Ex, Ey, fa))
            rows = np.array(rows)
            Em = np.array([np.exp(np.mean(np.log(0.5 * (rows[rows[:, 1] == r, 3] + rows[rows[:, 1] == r, 4])))) for r in rt])
            fa = np.array([rows[rows[:, 1] == r, 5].mean() for r in rt])
            m = np.gradient(np.log(Em), np.log(rt))
            out[f'{T}_{a}'] = dict(rho=rt.tolist(), E=Em.tolist(), m=m.round(3).tolist(), axial_frac=fa.round(3).tolist(),
                                   aniso=float(np.mean(np.minimum(rows[:, 3], rows[:, 4]) / np.maximum(rows[:, 3], rows[:, 4]))))
            print('E1 T', T, 'alpha', a, 'm', m.round(2), 'axial', fa.round(2), flush=True)
            json.dump(out, open('results/elastic_e1.json', 'w'), indent=1)


def e2(T=10000, NS=4, alphas=(0.0, 0.1, 0.2, 0.3, 0.4, 0.6), rho_end=(0.1, 0.4)):
    taus = np.unique(np.round(np.logspace(np.log10(30), np.log10(T), 14)).astype(int))
    out = {'taus': taus.tolist()}
    for a in alphas:
        for re in rho_end:
            recs = []
            for s in range(NS):
                g = srtree.run(T, 300 + s)
                S = rho_of(g, srtree.thickness(T, 1.0, a))
                lam0 = re / S
                lam = srtree.thickness(T, lam0, a)
                d = g['B'] - g['A']; l = np.linalg.norm(d, axis=1)
                c2 = (d[:, 0] / l) ** 2
                affx = np.cumsum(lam * l * c2 ** 2); affy = np.cumsum(lam * l * (1 - c2) ** 2)
                rho = np.cumsum(lam * l)
                for tau in taus:
                    try:
                        Ex, Ey, fa, *_ = Eboth(trunc(g, int(tau)), lam[:tau])
                    except Exception as e:
                        print('FAILED e2', a, re, s, tau, e); continue
                    recs.append(dict(seed=s, tau=int(tau), E=float(0.5 * (Ex + Ey)), aff=float(0.5 * (affx[tau - 1] + affy[tau - 1])),
                                     rho=float(rho[tau - 1]), axial=float(fa), lam0=float(lam0),
                                     s_new=float(1.165 * tau ** -0.5 / lam[tau - 1])))
            out[f'{a}_{re}'] = recs
            # summary
            E = np.array([[r['E'] for r in recs if r['seed'] == s] for s in range(NS)])
            A = np.array([[r['aff'] for r in recs if r['seed'] == s] for s in range(NS)])
            R = np.array([[r['rho'] for r in recs if r['seed'] == s] for s in range(NS)])
            psi = (np.diff(E, axis=1).sum(0) / np.diff(A, axis=1).sum(0))
            mT = np.diff(np.log(E.mean(0))) / np.diff(np.log(R.mean(0)))
            print('E2 alpha', a, 'rho_end', re, 'psi_el', psi.round(3), 'm_T', mT.round(2), flush=True)
            json.dump(out, open('results/elastic_e2.json', 'w'))


def e3(T=3000, NS=3, alphas=(0.0, 0.3, 0.6), rho=0.4):
    out = {}
    for a in alphas:
        X = []
        for s in range(NS):
            g = srtree.run(T, 400 + s)
            S = rho_of(g, srtree.thickness(T, 1.0, a))
            lam = srtree.thickness(T, rho / S, a)
            Ex, Ey, fa, sl, Ux, Uy = Eboth(g, lam)
            Ui = np.bincount(sl, Ux + Uy, minlength=T)
            k = np.zeros(T)
            for h in np.r_[g['hostA'], g['hostB']]:
                if h >= 0: k[h] += 1
            k += (g['hostA'] >= 0) + (g['hostB'] >= 0)
            l = np.linalg.norm(g['B'] - g['A'], axis=1)
            X.append(np.c_[np.arange(1, T + 1), k, lam, l, Ui])
        X = np.concatenate(X)
        m = (X[:, 0] > 30) & (X[:, 0] < 0.8 * T) & (X[:, 4] > 0) & (X[:, 1] > 2)
        x_i = -np.polyfit(np.log(X[m, 0]), np.log(X[m, 4]), 1)[0]
        M = np.c_[np.log(X[m, 1]), np.log(X[m, 2]), np.ones(m.sum())]
        coef = np.linalg.lstsq(M, np.log(X[m, 4]), rcond=None)[0]
        M2 = np.c_[np.log(X[m, 1]), np.log(X[m, 2]), np.log(X[m, 3]), np.ones(m.sum())]
        coef2 = np.linalg.lstsq(M2, np.log(X[m, 4]), rcond=None)[0]
        # share of total energy by birth decade
        tot = X[:, 4].sum()
        share = {f'{lo}-{3*lo}': float(X[(X[:, 0] >= lo) & (X[:, 0] < 3 * lo), 4].sum() / tot) for lo in (1, 3, 10, 30, 100, 300, 1000)}
        out[str(a)] = dict(x_i=float(x_i), paper_x=0.75 + a / 2, U_k_lam=coef[:2].tolist(),
                           U_k_lam_l=coef2[:3].tolist(), share=share)
        print('E3 alpha', a, out[str(a)], flush=True)
    json.dump(out, open('results/elastic_e3.json', 'w'), indent=1)


if __name__ == '__main__':
    which = sys.argv[1]
    t0 = time.time()
    {'e1': e1, 'e2': e2, 'e3': e3}[which]()
    print('time', time.time() - t0)
