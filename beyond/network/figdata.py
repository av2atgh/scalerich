"""Distributions for the network figures (T = 2e4, seeds 1-2): cell elongation (diameter / width)
and span lengths in units of T^(-1/2). Also the cut-angle family summary at T = 1e3 (seeds 100-103)
joined with the physics of results/anglegrade.json.  -> results/net_fig.json"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [B, os.path.join(B, '..', 'rg'), HERE]
import divisions, srtree, netstats

TH = np.linspace(0, np.pi, 180, endpoint=False); DIR = np.stack([np.cos(TH), np.sin(TH)], 1)
def elong(p):
    pr = p @ DIR.T; w = pr.max(0) - pr.min(0); return w.max() / max(w.min(), 1e-300)

def ccdf(x, grid):
    x = np.sort(x); return (len(x) - np.searchsorted(x, grid, side='left')) / len(x)

def main():
  out = {"dist": {}, "kappa": []}
  EG = np.geomspace(1, 1e4, 81); SG = np.geomspace(1e-3, 300, 111)
  for name in ['SR', 'uniform', 'rect', 'centroid', 'balanced', 'regular', 'kappa0.4']:
      e = []; s = []
      for sd in (1, 2):
          T = 20000
          g = divisions.run(T, sd, netstats.rule_of(name))
          e += [elong(p) for p in g['leaf_polys']]
          nodes, seg, sl, bn = srtree.network(g)
          ls = np.linalg.norm(nodes[seg[:, 1]] - nodes[seg[:, 0]], axis=1); ls = ls[ls > 1e-14]
          s += list(ls * np.sqrt(T))
      out['dist'][name] = dict(elong_grid=EG.tolist(), elong_ccdf=ccdf(np.array(e), EG).tolist(),
                               span_grid=SG.tolist(), span_ccdf=ccdf(np.array(s), SG).tolist(),
                               elong_q=[float(np.quantile(e, q)) for q in (0.5, 0.95, 0.999)])
      print(name, out['dist'][name]['elong_q'], flush=True)
  ag = json.load(open(os.path.join(B, 'results', 'anglegrade.json')))
  for k in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
      for sd in range(4):
          g = divisions.run(1000, 100 + sd, netstats.rule_of('kappa%g' % k))
          nodes, seg, sl, bn = srtree.network(g, frame=True)
          d = nodes[seg[:, 1]] - nodes[seg[:, 0]]; l = np.linalg.norm(d, axis=1); m = l > 1e-12; d = d[m]; l = l[m]
          phi = np.arctan2(d[:, 1], d[:, 0]); th = np.linspace(0, np.pi, 181)
          Caff = (l[None, :] * np.cos(phi[None, :] - th[:, None]) ** 4).sum(1)
          e = np.array([elong(p) for p in g['leaf_polys']])
          row = dict(kappa=k, seed=sd, l_star=float(((l ** 11).sum() / l.sum()) ** 0.1) * np.sqrt(1000),
                     l_max=float(l.max()) * np.sqrt(1000), eta_aff=float(np.sqrt(Caff.min() / Caff.max())),
                     elong_q95=float(np.quantile(e, 0.95)),
                     fabric4=float((l * np.cos(4 * phi)).sum() / l.sum()))
          for x in ag:
              if abs(x['kappa'] - k) < 1e-9 and x['seed'] == sd and x['rho_t'] == 0.2:
                  row.update(n=x['n_mean'], eta=x['eta_n'], omega1=x['omega1'])
          out['kappa'].append(row)
  json.dump(out, open(os.path.join(B, 'results', 'net_fig.json'), 'w'))


if __name__ == '__main__':
    main()
