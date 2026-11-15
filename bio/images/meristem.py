"""Shoot apical meristem L1 (Willis et al. 2016, 4 h time-lapse with lineage) as an
iterated-division partition. Run: ../.venv/bin/python -I meristem.py <plant> [...]

Per time point: L1 projection (last tissue voxel along z; label 1 = background),
walls and triple junctions (planar.label_graph). Wall identity through time from the
lineage: piece (a|b) at t+1 belongs to the line of (mom a | mom b) at t, or is a new
division wall when mom a == mom b. For each division:
  hosts  = lines of the mother's walls on which the new wall ends,
  corner = both hosts are sides of the mother that meet at a vertex of the mother
           (SR 0.50, uniform clocks 0.63, centroid 0.41, balanced/Errera 0.07),
  T-angle at each new end: angle between the two host pieces (180 = straight T),
          followed through later time points by line identity (relaxation vs age).
Dual degree of a line = new walls ending on it + its own ends on hosts.
"""
import sys, glob, gzip, re, ast, json, numpy as np
sys.path.insert(0, '/Users/avazquez/av2atg/scalerich/bio/images')
import tifffile
from scipy import ndimage as ndi
from planar import label_graph, arm_dir

ROOT = '/Users/avazquez/av2atg/scalerich/bio/data/raw/plant_willis2016/PNAS'


def load_proj(f):
    a = tifffile.imread(gzip.open(f))
    fg = a != 1
    zr = fg[::-1].argmax(0); has = fg.any(0)
    ztop = a.shape[0] - 1 - zr
    lab = np.take_along_axis(a, ztop[None], 0)[0].astype(np.int32)
    lab[~has] = 0; lab[lab == 1] = 0
    h = np.where(has, ztop, np.nan)
    return lab, h


def tracking(plant, t0, t1):
    f = f'{ROOT}/{plant}/tracking_data/{t0}hrs_to_{t1}hrs.txt'
    mom = {}; kids = {}
    for line in open(f):
        if ':' not in line:
            continue
        m, ds = line.split(':', 1)
        m = int(m); ds = [int(x) for x in ds.replace(',', ' ').split()]
        kids[m] = ds
        for d in ds:
            mom[d] = m
    return mom, kids


def frame(lab):
    V, E, vl = label_graph(lab)
    walls = {}
    for i, (u, v, pts, ab) in enumerate(E):
        walls.setdefault(ab, []).append(i)
    return V, E, vl, walls


def arms_at(V, E, vertex, r1=2.0, r2=8.0):
    """{(a,b): direction} of walls leaving vertex."""
    out = {}
    for u, v, pts, ab in E:
        if u == vertex:
            d = arm_dir(np.asarray(pts, float), r1, r2)
        elif v == vertex:
            d = arm_dir(np.asarray(pts, float)[::-1], r1, r2)
        else:
            continue
        if d is not None:
            out[ab] = d
    return out


def run(plant, slope_max=0.7):
    files = glob.glob(f'{ROOT}/{plant}/segmentation_tiffs/*hrs_*.tif.gz')
    times = sorted({int(re.match(r'(\d+)hrs', f.split('/')[-1]).group(1)) for f in files})
    res = [l for l in open(f'{ROOT}/{plant}/{plant}_tiff_resolutions.txt') if not l.startswith('#')]
    resl = ast.literal_eval(''.join(res).strip())
    resd = {int(r[0][:-3]): r for r in resl}
    frames = {}
    for t in times:
        f = [x for x in files if x.split('/')[-1].startswith(f'{t}hrs_')][0]
        lab, h = load_proj(f)
        r = resd[t]; xr = r[1][0]; zr = r[1][2] / r[2]
        gy, gx = np.gradient(ndi.gaussian_filter(np.nan_to_num(h, nan=0.0), 6))
        slope = np.hypot(gx, gy) * zr / xr
        V, E, vl, walls = frame(lab)
        frames[t] = dict(lab=lab, slope=slope, V=V, E=E, vl=vl, walls=walls)
        print(plant, t, 'cells', len(np.unique(lab)) - 1, 'walls', len(walls), flush=True)
    # line identity
    nid = [0]
    line = {}          # (t, (a,b)) -> line id
    birth = {}         # line id -> birth time (None = present at start)
    hosts = {}         # line id -> list of host line ids
    t0 = times[0]
    for ab in frames[t0]['walls']:
        if ab[0] == 0:
            continue
        line[(t0, ab)] = nid[0]; birth[nid[0]] = None; nid[0] += 1
    divs = []
    for ta, tb in zip(times[:-1], times[1:]):
        if tb - ta != 4:
            continue
        try:
            mom, kids = tracking(plant, ta, tb)
        except FileNotFoundError:
            continue
        Fa, Fb = frames[ta], frames[tb]
        for ab in Fb['walls']:
            if ab[0] == 0:
                continue
            a, b = ab
            ma, mb = mom.get(a), mom.get(b)
            if ma is None or mb is None:
                lid = None
            elif ma == mb:
                lid = nid[0]; birth[lid] = tb; nid[0] += 1
            else:
                lid = line.get((ta, (min(ma, mb), max(ma, mb))))
            if lid is None:
                lid = nid[0]; birth[lid] = 'unknown'; nid[0] += 1
            line[(tb, ab)] = lid
        # divisions
        for M, ds in kids.items():
            if len(ds) != 2:
                continue
            d1, d2 = sorted(ds); ab = (d1, d2)
            if ab not in Fb['walls']:
                continue
            new = line[(tb, ab)]
            # mother geometry at ta
            mpix = Fa['lab'] == M
            if mpix.sum() == 0:
                continue
            cy, cx = np.argwhere(mpix).mean(0)
            if Fa['slope'][int(cy), int(cx)] > slope_max:
                continue
            ends = []
            for v, labs in Fb['vl'].items():
                if d1 in labs and d2 in labs:
                    X = [x for x in labs if x not in (d1, d2)]
                    ends.append((v, X))
            hl = []; ok = True; tang = []
            for v, X in ends:
                if len(X) != 1 or X[0] == 0:
                    ok = False; continue
                x = X[0]; mx = mom.get(x)
                if mx is None:
                    ok = False; continue
                h = line.get((ta, (min(M, mx), max(M, mx))))
                if h is None:
                    ok = False; continue
                hl.append((h, mx))
                arms = arms_at(Fb['V'], Fb['E'], v)
                p = (min(d1, x), max(d1, x)); q = (min(d2, x), max(d2, x))
                if p in arms and q in arms:
                    tang.append(float(np.degrees(np.arccos(np.clip(arms[p] @ arms[q], -1, 1)))))
            corner = None
            if ok and len(hl) == 2:
                x1, x2 = hl[0][1], hl[1][1]
                corner = any({M, x1, x2} <= set(l) for l in Fa['vl'].values())
                hosts[new] = [hl[0][0], hl[1][0]]
            divs.append(dict(t=tb, M=int(M), new=int(new), area=int(mpix.sum()), nends=len(ends),
                             hosts=[int(h) for h, _ in hl], corner=corner, Tangle=tang))
    # T-angle relaxation: at every frame, vertices where a dated new line ends on a host
    relax = []
    for t in times:
        F = frames[t]
        for v, labs in F['vl'].items():
            labs = [l for l in labs if l != 0]
            if len(labs) != 3:
                continue
            arms = arms_at(F['V'], F['E'], v)
            ids = {ab: line.get((t, ab)) for ab in arms}
            for ab, lid in ids.items():
                if lid is None or not isinstance(birth.get(lid), int) or lid not in hosts:
                    continue
                others = [(ab2, l2) for ab2, l2 in ids.items() if ab2 != ab]
                if len(others) == 2 and others[0][1] == others[1][1] and others[0][1] in hosts[lid]:
                    ang = float(np.degrees(np.arccos(np.clip(arms[others[0][0]] @ arms[others[1][0]], -1, 1))))
                    relax.append((t - birth[lid], ang))
    # dual degree over the observation
    deg = {}
    for new, hs in hosts.items():
        for h in hs:
            deg[h] = deg.get(h, 0) + 1
            deg[new] = deg.get(new, 0) + 1
    return dict(plant=plant, times=times, divs=divs, relax=relax,
                deg={int(k): v for k, v in deg.items()},
                birth={int(k): (v if v is None or isinstance(v, int) else str(v)) for k, v in birth.items()})


if __name__ == '__main__':
    for p in sys.argv[1:]:
        out = run(p)
        json.dump(out, open(f'/Users/avazquez/av2atg/scalerich/bio/results/meristem_{p}.json', 'w'))
        dv = [d for d in out['divs'] if d['corner'] is not None]
        cf = np.mean([d['corner'] for d in dv]) if dv else float('nan')
        ta = [a for d in out['divs'] for a in d['Tangle']]
        print(p, 'divisions', len(out['divs']), 'with two hosts', len(dv), f'corner-cut {cf:.3f}',
              f'new-T angle median {np.median(ta):.1f}' if ta else '', flush=True)
        rl = np.array(out['relax'])
        if len(rl):
            for age in sorted(set(rl[:, 0].astype(int)))[:10]:
                s = rl[:, 0] == age
                print('   age', age, 'h  n', s.sum(), ' through-angle median', np.median(rl[s, 1]).round(1))
