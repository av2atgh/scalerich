"""Referee 3a: isotropy without mirror symmetry. Single cell with affine (kinematic) boundary
displacements, full C including C16, C26; cut-angle family and the rules at rho = 0.2, seeds 100-103.
-> results/kinematic_eta.json"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [B, os.path.join(B, '..', 'rg'), HERE]
import divisions, waves, netstats
out = []
for name in ['kappa%g' % k for k in (0, 0.2, 0.4, 0.6, 0.8, 1)] + ['SR', 'centroid', 'balanced', 'regular']:
    for s in range(4):
        g = divisions.run(1000, 100 + s, netstats.rule_of(name))
        L = np.linalg.norm(g['B'] - g['A'], axis=1); lam = np.full(1000, 0.2 / L.sum())
        r = waves.analyse(g, lam, resonance=False)
        nodes, seg, lam_s, Ls, K, bnd, rho = waves.system(g, lam)
        C = waves.stiffness_tensor(nodes, K, bnd)[0]
        off = float(np.hypot(C[0, 2], C[1, 2]) / np.sqrt(C[0, 0] * C[1, 1]))
        out.append(dict(rule=name, seed=s, eta_kin=r['n_min'] / max(r['n_max'], 1e-9) if 'n_max' in r else None,
                        n_kin=r['n_mean'], off=off))
        print(name, s, 'eta_kin %.3f off-axis %.3f' % (out[-1]['eta_kin'], off), flush=True)
json.dump(out, open(os.path.join(B, 'results', 'kinematic_eta.json'), 'w'), indent=1)
