"""Linear spin waves on magnetic wire networks drawn on the divisions (frame excluded, free ends).
Each wire is magnetised along its tangent t; the deviation is m - t = u n + w z, n = z x t.
Energy  E = 1/2 sum_wires int [ (u'^2 + w'^2) + a u^2 + b w^2 ] ds  (exchange + shape anisotropy,
in-plane a, out-of-plane b > a), lengths in units of the mean cell side l = <A>^{1/2}.
Junctions: the deviation vector is continuous. Out of plane, w takes one value at the junction
(scalar Kirchhoff, blind to angles). In plane, the junction carries a vector V and each wire
takes u = V . n_wire, so collinear wires (straight continuation) are joined and perpendicular ones
are not: a T-junction is a free end for the stem and transparent for the bar.
Precession: K_u u = omega G w,  K_w w = omega G^T u  ->  K_u u = omega^2 (G K_w^-1 G^T) u.
FEM: linear elements of length <= h along every wire."""
import sys, json, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, scipy.linalg as sl
import supercond as sc

def network(g, T):
    X, S = sc.wires(g)
    onf = np.zeros(len(S), bool)
    for ax in (0, 1):
        for c in (0.0, 1.0):
            onf |= (np.abs(X[S[:, 0], ax] - c) < 1e-9) & (np.abs(X[S[:, 1], ax] - c) < 1e-9)
    S = S[~onf]
    used = np.unique(S); ren = -np.ones(len(X), int); ren[used] = np.arange(len(used))
    return X[used] * np.sqrt(T + 1), ren[S]

def assemble(X, S, h=0.25, a=1.0, b=4.0):
    nJ = len(X)
    d = X[S[:, 1]] - X[S[:, 0]]; L = np.linalg.norm(d, axis=1); t = d / L[:, None]
    nrm = np.c_[-t[:, 1], t[:, 0]]
    # junction in-plane basis: 1 dof if all incident normals are collinear, else 2
    inc = [[] for _ in range(nJ)]
    for s, (p, q) in enumerate(S): inc[p].append(s); inc[q].append(s)
    ub = []; nu = 0; ubase = np.zeros(nJ, int); bas = []
    for j in range(nJ):
        N = nrm[inc[j]]
        if np.all(np.abs(N @ np.array([-N[0, 1], N[0, 0]])) < 1e-9): B = N[:1]       # rows = basis vectors
        else: B = np.eye(2)
        bas.append(B); ubase[j] = nu; nu += len(B)
    nw = nJ
    # elements
    Eu = []; Ew = []; El = []; Eseg = []
    for s, (p, q) in enumerate(S):
        k = max(1, int(np.ceil(L[s] / h))); l = L[s] / k
        # u coefficients of node p seen by this wire: V . n = sum_r c_r V_r, V = sum_r c_r B_r
        cp = bas[p] @ nrm[s]; cq = bas[q] @ nrm[s]
        up = [(ubase[p] + r, cp[r]) for r in range(len(cp))]
        uq = [(ubase[q] + r, cq[r]) for r in range(len(cq))]
        ui = [[(nu + i, 1.0)] for i in range(k - 1)]; nu += k - 1
        wi = list(range(nw, nw + k - 1)); nw += k - 1
        un = [up] + ui + [uq]; wn = [p] + wi + [q]
        for e in range(k):
            Eu.append((un[e], un[e + 1])); Ew.append((wn[e], wn[e + 1])); El.append(l); Eseg.append(s)
    ke = np.array([[1, -1], [-1, 1]]); me = np.array([[2, 1], [1, 2]]) / 6
    Ru, Cu, Vku, Vmu = [], [], [], []; Rw, Cw, Vkw, Vmw = [], [], [], []; Rg, Cg, Vg = [], [], []
    for (eu, ew, l) in zip(Eu, Ew, El):
        for al in (0, 1):
            for be in (0, 1):
                for (i, ci) in eu[al]:
                    for (j, cj) in eu[be]:
                        Ru.append(i); Cu.append(j); Vku.append(ci * cj * ke[al, be] / l); Vmu.append(ci * cj * me[al, be] * l)
                    Rg.append(i); Cg.append(ew[be]); Vg.append(ci * me[al, be] * l)
                Rw.append(ew[al]); Cw.append(ew[be]); Vkw.append(ke[al, be] / l); Vmw.append(me[al, be] * l)
    Ku = sp.coo_matrix((np.array(Vku) + a * np.array(Vmu), (Ru, Cu)), shape=(nu, nu)).tocsc()
    Kw = sp.coo_matrix((np.array(Vkw) + b * np.array(Vmw), (Rw, Cw)), shape=(nw, nw)).tocsc()
    G = sp.coo_matrix((Vg, (Rg, Cg)), shape=(nu, nw)).tocsc()
    info = dict(Eu=Eu, Ew=Ew, El=np.array(El), Eseg=np.array(Eseg), Ltot=float(L.sum()), nJ=nJ,
                deg=np.array([len(x) for x in inc]))
    return Ku, Kw, G, info

def spectrum_dense(Ku, Kw, G):
    """all omega; Ku u = omega^2 (G Kw^-1 G^T) u solved as Bm u = mu Ku u, mu = omega^-2"""
    W = spl.splu(Kw.tocsc()).solve(G.T.toarray())
    Bm = G @ W; Bm = 0.5 * (Bm + Bm.T)
    mu, U = sl.eigh(Bm, Ku.toarray())
    o = mu > 1e-12
    om = 1 / np.sqrt(mu[o]); U = U[:, o]
    s = np.argsort(om); return om[s], U[:, s], W

def spectrum_window(Ku, Kw, G, om0, nev=200):
    """nev modes nearest om0: blockdiag(Ku, Kw) x = omega [[0, G], [G^T, 0]] x, shift-invert in 1/omega"""
    nu = Ku.shape[0]
    A = sp.block_diag([Ku, Kw]).tocsc(); Bb = sp.bmat([[None, G], [G.T, None]]).tocsc()
    mu, V = spl.eigsh(Bb, k=nev, M=A, sigma=1.0 / om0, which='LM')
    o = mu > 0; om = 1 / mu[o]; U = V[:nu, o]; s = np.argsort(om)
    U = U[:, s]; om = om[s]
    return om, U, V[nu:, o][:, s]

def density(U, W, om, info, Wv=None):
    """energy-like weight u^2 + w^2 per element; returns participation fraction per mode"""
    Eu, Ew, El = info['Eu'], info['Ew'], info['El']
    if Wv is None: Wv = (W @ U) * om[None, :]       # w = omega Kw^-1 G^T u
    def pad(end):
        I = np.array([[t[0] for t in x] + [0] * (2 - len(x)) for x in end])
        C = np.array([[t[1] for t in x] + [0.0] * (2 - len(x)) for x in end]); return I, C
    I0, C0 = pad([e[0] for e in Eu]); I1, C1 = pad([e[1] for e in Eu])
    u0 = C0[:, :1] * U[I0[:, 0]] + C0[:, 1:] * U[I0[:, 1]]
    u1 = C1[:, :1] * U[I1[:, 0]] + C1[:, 1:] * U[I1[:, 1]]
    Jw = np.array(Ew); w0 = Wv[Jw[:, 0]]; w1 = Wv[Jw[:, 1]]
    rho = El[:, None] * (u0 ** 2 + u0 * u1 + u1 ** 2 + w0 ** 2 + w0 * w1 + w1 ** 2) / 3
    rho /= rho.sum(0)
    pr = 1.0 / (El.sum() * (rho ** 2 / El[:, None]).sum(0))
    return pr

CASES = sc.CASES

if __name__ == '__main__':
    names = sys.argv[1].split(','); Ts = [int(x) for x in sys.argv[2].split(',')]
    fout = sys.argv[3] if len(sys.argv) > 3 else 'results/spinwave.json'
    for nm in names:
        for T in Ts:
            try: out = json.load(open(fout))
            except Exception: out = {}
            key = '%s_T%d' % (nm, T)
            if key in out: continue
            t0 = time.time()
            try:
                g = CASES[nm](T); X, S = network(g, T)
                Ku, Kw, G, info = assemble(X, S)
                om, U, W = spectrum_dense(Ku, Kw, G)
                pr = density(U, W, om, info)
            except Exception as err:
                print(key, 'FAILED', repr(err), flush=True); continue
            out[key] = dict(rule=nm, T=T, nu=Ku.shape[0], Ltot=info['Ltot'], nJ=info['nJ'],
                            om=om.tolist(), pr=pr.tolist())
            json.dump(out, open(fout, 'w'))
            print(key, 'nu', Ku.shape[0], 'L', '%.0f' % info['Ltot'], 'modes', len(om),
                  'om %.2f..%.2f' % (om[0], om[-1]), '%.0fs' % (time.time() - t0), flush=True)
