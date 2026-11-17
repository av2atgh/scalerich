"""Single-mode quantum-wire networks on the divisions at B = 0: Anderson localisation test.
Quantum graph (Schrodinger on each wire, Kirchhoff/Neumann junctions), open sample with frame,
linear FEM along wires (element <= hfac <A>^{1/2}). Eigenstates in the window around
k^2 <A> = E0 (wavelength ~ 4 cell sides) by shift-invert. Per state: relative participation
p = 1 / (N sum |psi_i|^4) (|psi|^2 normalised on FEM nodes). Per window: level-spacing ratio
<r> (GOE 0.531, Poisson 0.386). Junction scattering is blind to angles, so alignment enters only
through topology. -> results/localisation.json (appends per case/T/seed)."""
import sys, json, time, numpy as np
import scipy.sparse as sp, scipy.sparse.linalg as spl
import supercond as sc

def run(g, T, nev=200, E0=2.0, hfac=1 / 6):
    X, S = sc.wires(g)
    Y, E = sc.discretise(X, S, hfac / np.sqrt(T + 1))
    K, M = sc.matrices(Y, E, 0.0)
    K = K.real.tocsc(); M = M.real.tocsc()
    sigma = E0 * (T + 1)
    w, v = spl.eigsh(K, k=nev, M=M, sigma=sigma, which='LM')
    o = np.argsort(w); w = w[o]; v = v[:, o]
    P = v ** 2; P /= P.sum(0)
    pr = 1.0 / (len(Y) * (P ** 2).sum(0))
    s = np.diff(w); s = s[s > 1e-9 * sigma]
    r = np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
    return dict(N=len(Y), k2A=(w / (T + 1)).tolist(), pr=pr.tolist(),
                r=float(r.mean()), pr_med=float(np.median(pr)), ndeg=int(nev - 1 - len(s)))

if __name__ == '__main__':
    names = sys.argv[1].split(',')
    Ts = [int(x) for x in sys.argv[2].split(',')]
    seeds = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [0]
    fout = sys.argv[4] if len(sys.argv) > 4 else 'results/localisation.json'
    for nm in names:
        for T in Ts:
            for sd in (seeds if nm in ('uniform', 'SR') else [0]):
                try: out = json.load(open(fout))
                except Exception: out = {}
                key = '%s_T%d_s%d' % (nm, T, sd)
                if key in out: continue
                t0 = time.time()
                try:
                    if nm == 'uniform':
                        g = sc.qdiv.run(T, ('largest', 'uniform', 'perp'), seed=sd)
                    elif nm == 'SR':
                        g = sc.divisions.run(T, sd, 'SR')
                    else:
                        g = sc.CASES[nm](T)
                    res = run(g, T)
                except Exception as err:          # skip and keep going
                    print(key, 'FAILED', repr(err), flush=True); continue
                res.update(rule=nm, T=T, seed=sd)
                out[key] = res
                json.dump(out, open(fout, 'w'))
                print(key, 'N', res['N'], '<r> %.3f' % res['r'], 'p_med %.4f' % res['pr_med'],
                      'deg', res['ndeg'], '%.0fs' % (time.time() - t0), flush=True)
