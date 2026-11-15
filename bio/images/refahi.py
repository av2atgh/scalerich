"""Flower meristem L1 (Refahi et al. 2021, wt FM1-FM6) as an iterated-division partition.
Like meristem.py, with: variable intervals, mothers with >= 2 daughters (all walls between
sisters are new), mother 1 (background) dropped, stack orientation detected per file.
Outputs per film: divisions (new line, hosts, corner cut for 2-daughter divisions, T-angles),
host relations for every new wall, line births, through-angle vs age, division fraction per interval.
Run: ../.venv/bin/python -I refahi.py FM1 [FM5 ...]"""
import sys, glob, gzip, re, json, numpy as np
sys.path.insert(0, '/Users/avazquez/av2atg/scalerich/bio/images')
import tifffile
from scipy import ndimage as ndi
from planar import label_graph
from meristem import arms_at

ROOT = '/Users/avazquez/av2atg/scalerich/bio/data/raw/plant_refahi2021/wt'


def load_proj(f):
    a = tifffile.imread(gzip.open(f))
    fg = a != 1
    top_high = (a[-1] == 1).mean() >= (a[0] == 1).mean()     # apex toward the emptier end
    if top_high:
        zr = fg[::-1].argmax(0); z = a.shape[0] - 1 - zr
    else:
        z = fg.argmax(0)
    has = fg.any(0)
    lab = np.take_along_axis(a, z[None], 0)[0].astype(np.int32)
    lab[~has] = 0; lab[lab == 1] = 0
    return lab, np.where(has, z, np.nan)


def tracking(fm, ta, tb):
    mom, kids = {}, {}
    for line in open(f'{ROOT}/{fm}/tracking_data/{ta}hrs_to_{tb}hrs.txt'):
        if ':' not in line:
            continue
        m_, ds = line.split(':', 1); m_ = int(m_)
        if m_ == 1:
            continue
        ds = [int(x) for x in ds.replace(',', ' ').split()]
        kids[m_] = ds
        for d in ds:
            mom[d] = m_
    return mom, kids


def run(fm, slope_max=0.7):
    files = glob.glob(f'{ROOT}/{fm}/segmentation_tiffs/*h_segmented.tif.gz')
    times = sorted(int(re.match(r'(\d+)h_', f.split('/')[-1]).group(1)) for f in files)
    F = {}
    for t in times:
        lab, h = load_proj(f'{ROOT}/{fm}/segmentation_tiffs/{t}h_segmented.tif.gz')
        gy, gx = np.gradient(ndi.gaussian_filter(np.nan_to_num(h, nan=0.0), 6))
        V, E, vl = label_graph(lab)
        walls = {}
        for i, e in enumerate(E):
            walls.setdefault(e[3], []).append(i)
        F[t] = dict(lab=lab, slope=np.hypot(gx, gy), V=V, E=E, vl=vl, walls=walls)
        print(fm, t, 'cells', len(np.unique(lab)) - 1, flush=True)
    nid = [0]; line = {}; birth = {}; hostrel = {}
    for ab in F[times[0]]['walls']:
        if ab[0]:
            line[(times[0], ab)] = nid[0]; birth[nid[0]] = None; nid[0] += 1
    divs = []; fracs = []
    for ta, tb in zip(times[:-1], times[1:]):
        try:
            mom, kids = tracking(fm, ta, tb)
        except FileNotFoundError:
            continue
        fracs.append((ta, tb, float(np.mean([len(v) >= 2 for v in kids.values()])),
                      float(np.mean([np.log2(len(v)) for v in kids.values()]))))
        Fa, Fb = F[ta], F[tb]
        for ab in Fb['walls']:
            if not ab[0]:
                continue
            ma, mb = mom.get(ab[0]), mom.get(ab[1])
            if ma is None or mb is None:
                lid = None
            elif ma == mb:
                lid = nid[0]; birth[lid] = tb; nid[0] += 1
            else:
                lid = line.get((ta, (min(ma, mb), max(ma, mb))))
            if lid is None:
                lid = nid[0]; birth[lid] = 'unknown'; nid[0] += 1
            line[(tb, ab)] = lid
        # host relations of every new wall
        for ab in Fb['walls']:
            if not ab[0]:
                continue
            a, b = ab; M = mom.get(a)
            if M is None or mom.get(b) != M:
                continue
            new = line[(tb, ab)]; hs = []
            for v, labs in Fb['vl'].items():
                if a in labs and b in labs:
                    X = [x for x in labs if x not in (a, b)]
                    if len(X) != 1 or X[0] == 0:
                        hs.append(None); continue
                    x = X[0]; mx = mom.get(x)
                    if mx is None:
                        hs.append(None)
                    elif mx != M:
                        hs.append(line.get((ta, (min(M, mx), max(M, mx)))))
                    else:                              # ends on a sister wall born in the same interval
                        hs.append(line.get((tb, (min(a, x), max(a, x)))))
            hostrel[new] = hs
        # 2-daughter divisions: corner cut, T-angles
        for M, ds in kids.items():
            if len(ds) != 2:
                continue
            d1, d2 = sorted(ds)
            if (d1, d2) not in Fb['walls']:
                continue
            mp = Fa['lab'] == M
            if not mp.any():
                continue
            cy, cx = np.argwhere(mp).mean(0)
            if Fa['slope'][int(cy), int(cx)] > slope_max:
                continue
            ends = [(v, [x for x in labs if x not in (d1, d2)]) for v, labs in Fb['vl'].items()
                    if d1 in labs and d2 in labs]
            hl = []; ok = True; tang = []
            for v, X in ends:
                if len(X) != 1 or X[0] == 0 or mom.get(X[0]) in (None, M):
                    ok = False; continue
                x = X[0]; mx = mom[x]
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
                corner = any({M, hl[0][1], hl[1][1]} <= set(l) for l in Fa['vl'].values())
            divs.append(dict(ta=ta, tb=tb, M=int(M), new=int(line[(tb, (d1, d2))]),
                             hosts=[int(h) for h, _ in hl], corner=corner, Tangle=tang))
    relax = []
    for t in times:
        Fr = F[t]
        for v, labs in Fr['vl'].items():
            labs = [l for l in labs if l]
            if len(labs) != 3:
                continue
            arms = arms_at(Fr['V'], Fr['E'], v)
            ids = {ab: line.get((t, ab)) for ab in arms}
            for ab, lid in ids.items():
                if lid is None or not isinstance(birth.get(lid), int) or lid not in hostrel:
                    continue
                o = [(ab2, l2) for ab2, l2 in ids.items() if ab2 != ab]
                if len(o) == 2 and o[0][1] == o[1][1] and o[0][1] in hostrel[lid]:
                    relax.append((t - birth[lid], float(np.degrees(np.arccos(np.clip(arms[o[0][0]] @ arms[o[1][0]], -1, 1))))))
    return dict(fm=fm, times=times, divs=divs, relax=relax, fracs=fracs,
                hostrel={int(k): [None if h is None else int(h) for h in v] for k, v in hostrel.items()},
                birth={int(k): (v if v is None or isinstance(v, int) else str(v)) for k, v in birth.items()})


if __name__ == '__main__':
    for fm in sys.argv[1:]:
        out = run(fm)
        json.dump(out, open(f'/Users/avazquez/av2atg/scalerich/bio/results/refahi_{fm}.json', 'w'))
        dv = [d for d in out['divs'] if d['corner'] is not None]
        print(fm, 'new walls', len(out['hostrel']), '2-daughter divisions', len(out['divs']),
              'corner-decidable', len(dv), 'corner cuts', sum(d['corner'] for d in dv), flush=True)
