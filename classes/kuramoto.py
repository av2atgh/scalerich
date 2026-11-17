"""Kuramoto oscillators on the ligament graph (ligaments = oscillators, linked where they meet;
frame excluded): d theta_i/dt = omega_i + K sum_j A_ij sin(theta_j - theta_i), omega ~ N(0, 1).
Heterogeneous mean field: K_c = 2/(pi g(0)) <k>/<k^2> = 2 sqrt(2 pi)/pi <k>/<k^2>.
Adiabatic sweep in kappa = K / K_c^MF (state carried over); RK4, dt = min(0.05, 2/(K lambda_max(L))); transient t_tr,
average of r = |<exp(i theta)>| over t_av; also r among the top-5% degree ligaments and the
largest fraction with equal effective frequency (tolerance 0.02). argv[5] = 'abs': common grid K = 0.02..2. -> results/kuramoto.json (appends per case/T/seed)."""
import sys, json, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import supercond as sc, ligaments as lg

def adjacency(g):
    X, S = sc.wires(g)
    # rebuild ligament ids and edges as in ligaments.ligament_graph, returning the sparse matrix
    d = X[S[:, 1]] - X[S[:, 0]]; u = d / np.linalg.norm(d, axis=1)[:, None]
    par = np.arange(len(S))
    def find(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    inc = [[] for _ in range(len(X))]
    for e, (a, b) in enumerate(S): inc[a].append(e); inc[b].append(e)
    for v in range(len(X)):
        es = inc[v]
        for i in range(len(es)):
            for j in range(i + 1, len(es)):
                e, f = es[i], es[j]
                se = u[e] if S[e, 0] == v else -u[e]; sf = u[f] if S[f, 0] == v else -u[f]
                if se @ sf < -1 + 1e-9: par[find(e)] = find(f)
    lig = np.array([find(e) for e in range(len(S))])
    onf = np.zeros(len(S), bool)
    for ax in (0, 1):
        for c in (0.0, 1.0):
            onf |= (np.abs(X[S[:, 0], ax] - c) < 1e-9) & (np.abs(X[S[:, 1], ax] - c) < 1e-9)
    frame = set(lig[onf])
    ids = [l for l in np.unique(lig) if l not in frame]; ren = {l: i for i, l in enumerate(ids)}
    edges = set()
    for v in range(len(X)):
        ls = sorted({ren[l] for l in {lig[e] for e in inc[v]} if l not in frame})
        for i in range(len(ls)):
            for j in range(i + 1, len(ls)): edges.add((ls[i], ls[j]))
    E = np.array(sorted(edges)); n = len(ids)
    A = sp.coo_matrix((np.ones(2 * len(E)), (np.r_[E[:, 0], E[:, 1]], np.r_[E[:, 1], E[:, 0]])), shape=(n, n)).tocsr()
    return A

def sweep(A, kappas, seed=0, dt=0.05, t_tr=200.0, t_av=200.0, absK=False):
    n = A.shape[0]; k = np.asarray(A.sum(1)).ravel()
    Kc = 2 * np.sqrt(2 * np.pi) / np.pi * k.mean() / (k ** 2).mean()
    rng = np.random.default_rng(seed); om = rng.standard_normal(n); th = rng.uniform(0, 2 * np.pi, n)
    def F(th, K):
        z = np.exp(1j * th); return om + K * np.imag(np.conj(z) * (A @ z))
    L = sp.diags(k) - A                               # RK4 stable for dt K lambda_max(L) < 2.78
    lmax = float(spl.eigsh(L.astype(float), k=1, which='LA', return_eigenvectors=False)[0])
    hub = k >= np.quantile(k, 0.95)
    rs = []
    for kap in kappas:
        K = kap * Kc if not absK else kap; acc = 0.0; acch = 0.0; m = 0
        h = min(dt, 2.0 / (K * lmax))
        nsteps = int((t_tr + t_av) / h); ntr = int(t_tr / h); th0 = None
        for step in range(nsteps):
            a = F(th, K); b = F(th + 0.5 * h * a, K); c = F(th + 0.5 * h * b, K); d = F(th + h * c, K)
            th = th + h / 6 * (a + 2 * b + 2 * c + d)
            if step == ntr: th0 = th.copy()
            if step >= ntr and step % max(1, int(dt / h)) == 0:
                z = np.exp(1j * th); acc += abs(z.mean()); acch += abs(z[hub].mean()); m += 1
        weff = np.sort((th - th0) / ((nsteps - 1 - ntr) * h))
        tol = 0.02; j = np.searchsorted(weff, weff + tol, side='right')
        locked = float((j - np.arange(n)).max() / n)        # largest group of equal effective frequency
        rs.append((acc / m, acch / m, locked))
    return Kc, rs, k

if __name__ == '__main__':
    names = sys.argv[1].split(',')
    Ts = [int(x) for x in sys.argv[2].split(',')]
    seeds = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [0]
    fout = sys.argv[4] if len(sys.argv) > 4 else 'results/kuramoto.json'
    absK = len(sys.argv) > 5 and sys.argv[5] == 'abs'
    kappas = np.round(np.geomspace(0.02, 2, 16), 4) if absK else np.round(np.geomspace(0.25, 8, 16), 4)
    for nm in names:
        for T in Ts:
            for sd in (seeds if nm in ('uniform', 'SR') else [0]):
                try: out = json.load(open(fout))
                except Exception: out = {}
                key = '%s_T%d_s%d' % (nm, T, sd)
                if key in out: continue
                t0 = time.time()
                try:
                    if nm == 'uniform': g = sc.qdiv.run(T, ('largest', 'uniform', 'perp'), seed=sd)
                    elif nm == 'SR': g = sc.divisions.run(T, sd, 'SR')
                    else: g = sc.CASES[nm](T)
                    A = adjacency(g)
                    Kc, rs, k = sweep(A, kappas, seed=sd, absK=absK)
                except Exception as err:          # skip and keep going
                    print(key, 'FAILED', repr(err), flush=True); continue
                out[key] = dict(rule=nm, T=T, seed=sd, n=A.shape[0], k1=float(k.mean()),
                                k2k1=float((k ** 2).mean() / k.mean()), KcMF=float(Kc),
                                kappa=kappas.tolist(), absK=absK, r=[x[0] for x in rs],
                                r_hub=[x[1] for x in rs], locked=[x[2] for x in rs])
                json.dump(out, open(fout, 'w'))
                print(key, 'n', A.shape[0], 'KcMF %.4f' % Kc, 'r:', ' '.join('%.2f' % x[0] for x in rs), 'locked:', ' '.join('%.2f' % x[2] for x in rs),
                      '%.0fs' % (time.time() - t0), flush=True)
