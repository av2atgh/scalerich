import json, sys, numpy as np, collections
rows = json.load(open(sys.argv[1]))
g = collections.defaultdict(list)
for r in rows: g[(r['rule'], r['alpha'], r['rho_t'])].append(r)
print(f"{'rule':10s} {'a':>4s} {'rho':>4s} {'n_x':>6s} {'n_mean':>6s} {'eta_n':>5s} {'w1':>6s}")
for k in sorted(g, key=lambda k: (k[0], k[2], k[1])):
    v = g[k]
    f = lambda q: np.exp(np.mean(np.log([x[q] for x in v])))
    print(f"{k[0]:10s} {k[1]:4.2f} {k[2]:4.2f} {f('n_x'):6.2f} {f('n_mean'):6.2f} {np.mean([x['eta_n'] for x in v]):5.2f} {f('omega1'):6.3f}")
print()
print('ranges over (alpha, rho), per-sample; isotropic = eta_n >= 0.8')
for rule in dict.fromkeys(r['rule'] for r in rows):
    v = [r for r in rows if r['rule'] == rule]
    nx = [r['n_x'] for r in v]; iso = [r['n_mean'] for r in v if r['eta_n'] >= 0.8]
    print(f"{rule:10s} n_x {min(nx):5.2f}-{max(nx):5.2f}  ratio {max(nx)/min(nx):5.2f} | iso n {('%5.2f-%5.2f'%(min(iso),max(iso))) if iso else '   none    '} ({len(iso)}/{len(v)}) | w1 {min(r['omega1'] for r in v):.3f}-{max(r['omega1'] for r in v):.3f}")
