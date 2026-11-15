# bio/ — do scale-rich (SR) partitions occur in living systems?

Started 2026-11-15. Question: the SR metamaterial = iterated cell division
(Cowan 2010): a segment nucleates inside a cell, grows both ways to the first
hit (T-junctions only), older segments thicker. Does this architecture arise in
cells and tissues, e.g. the actin network of flat (spread) human cells as a 2D example?

Relaxation of the exact model: divisions need not be sequential. Many
compartments split at the same time; segments grow at a finite speed v from
nucleation events of rate r per unit area per time, so growing segments can hit
each other. One dimensionless control:
    Π = r ℓ³ / v  ~  number of nucleations while a segment crosses a cell,
    sequential SR (Π → 0, v → ∞)  ...  Gilbert tessellation (all nucleated at once).
Plus turnover (segments removed at rate δ) since actin is dynamic.

## Workstreams
- `lit/`    literature: cytoskeletal (actin stress fibres, transverse arcs, cortex,
            keratin, vimentin) architecture as a planar partition; junction types;
            segment-length/thickness statistics; other candidates (plant cell
            division / Errera rule, leaf venation, crack-like skin patterns);
            prior art linking biology to Gilbert/STIT tessellations.
- `data/`   public F-actin images of flat human cells at network-resolving resolution.
- `model/`  parallel-division model: dual degree law p_k, gamma, C(k), thickness vs
            age/length, junction statistics vs Π and turnover; discriminators vs
            Poisson lines, Gilbert, Voronoi, Arp2/3 branching (70°), random fibres (X-crossings).
- `images/` skeleton → segments/junctions → dual graph pipeline on real cells; test
            the model's signatures.

Signatures of the SR class to test (from ../lng): T-junction dominance
(no X crossings), dual p_k ~ k^-3 in 2D, C(k) ~ 2p/k, thickness decreasing with
birth order, degree ∝ length × age.

## Status

### 2026-11-15 results (first pass)
Model `model/pargrow.py` (exact event-driven Gilbert race, numba; checks in model/check.py:
no X crossings, all tips land on their targets, Pi -> 0 reproduces SR).

1. **Parallel nucleation keeps the SR class** (`model/scan.py`, N = 2e4, 3 seeds,
   results/scan_N20000.log). Pi = 1e-3 ... 1e3: gamma10 2.82-2.87, C 0.376, C(k) ~ 1/k,
   triangle closure 0.49, corr(k-2, L x age) 0.94 — indistinguishable from SR. Only the
   oldest cohorts change. Local criterion: nucleations per cell crossing = Pi n^{-3/2},
   so the top of the hierarchy (first n* segments) is Gilbert-like:
   n* ~ 0.04 Pi^{2/3} (measured 0.024-0.055 for Pi = 1e3-1e6; `model/crossover.py`).
   Pi = 1e4: k_max 65 (SR 280); Pi = 1e6 (n* ~ N): Gilbert, gamma > 4, k_max 16.
   => SR statistics hold at scales below l* = (v/R)^{1/3} whatever happens above.
2. **Per-compartment clocks don't matter for gamma** (`model/percell.py`, T = 2e4):
   SR, uniform-pick (every compartment divides at the same rate), centroid, balanced
   (Errera-like) all give gamma10 = 2.9-3.1. Discriminators are clustering
   (SR 0.38, uniform 0.45, centroid 0.28, balanced 0.055) and length spread
   (sd ln L: 0.86, 3.3, 0.69, 0.59).
3. **Turnover truncates the tail unless loaded segments are stabilised**
   (`model/turnover.py`, n_ss = 1500). Uniform removal: steady state P(k>=16) 0.005-0.009,
   k_max 23-34 (SR at n = 1500: 0.021, 82) — a characteristic degree, almost
   n_ss-independent because a lifetime always spans n_ss nucleations. Removal rate
   prop. to (2/k)^s, s >= 1 (load-bearing bundles protected): P(k>=16) 0.018-0.020,
   k_max 70-90 — SR tail restored. Prediction: scale-rich actin needs tension-
   (degree-) dependent stability; corr(k, age) ~ 0.7 in all turnover runs.
4. **Pipeline** (`images/pipeline.py`, validated in images/validate.py on rendered networks,
   1024 px): T share 0.87-0.97 for SR/Gilbert; X/junction 0.04-0.06 (SR) vs 0.16-0.32
   (Mikado rods). Degree tail undercounted (P(k>=4) 0.16 vs true 0.29; k_max ~ 1/3 of truth):
   pooling over cells needed; gamma from images is biased up.
5. **BioSR F-actin (COS-7, Lifeact, TIRF-SIM, 51 cells; data/DATA.md)** — `images/run_biosr.py`,
   results/biosr_summary.json. 11389 junctions: T 4192, Y 2652, X 772, >=5-arm clusters
   ~1900; T/(T+Y) = 0.61; skeleton fragmented (62791 fibre pieces, P(k>=2) 0.07);
   corr(k, width | L) = 0.05. Overlay (results/overlay_Cell_011.png): the ventral
   Lifeact network is a mesh of long curved bundles crossing at shallow angles, not a
   T-partition. First pass says: not SR in COS-7 ventral actin at this resolution.
   Segmentation untuned (fragmentation), so not a final negative.

### 2026-11-15 literature (bio/lit/LIT.md)
- No junction census or dual degree law published for flat-cell actin; no power-law degree
  reported for any image-derived cytoskeletal network; Gilbert/STIT/Cowan never applied.
- Caution: SOAX tracing *creates* T-junctions (tips stop on hit, Xu 2015) — biased toward SR.
- Closest prior art: Sadegh 2017 PRX 7, 011031 — cortical actin compartments (HEK, STORM)
  log-normal, explained as compartments split by growing filaments; never tested vs STIT.
- Actin: Arp2/3 = one-sided Y (70°); dorsal fibres grow unidirectionally from adhesions;
  arcs thicken by lateral fusion (older = thicker, as in SR); stress-fibre actin t1/2 ~ 1 min
  => turnover regime, where our model says the tail is cut unless loaded bundles are stabilised.
- Keratin mean vertex degree 3.5-4.2 (Lück 2010): not T-dominated — negative control.
- Plant cytokinesis is literally "grow both ways, stop at first wall", 4-way junctions avoided;
  angle = Errera/Gibbs rule (our 'balanced': gamma ~ 3, C ~ 0.05). Leaf venation: T-junctions +
  radius hierarchy (Bohn 2002, Couder 2002); segmented data Ronellenfitsch 2015, Blonder 2020.
- Rock-fracture dual graph already built (Andresen 2013): broad degree, no T/X split.
- Cheap first test: mean sides per face = 4 for successive division without reorganisation
  (vs 6 Voronoi-like; Bohn, Douady & Couder PRL 94, 054503).
- Ranked targets: leaf minor venation > meristem time-lapse > stress fibres > crocodile
  scales > cortical actin > keratin (control). Order of tests: T/Y/X share, sides/face,
  width vs age, then dual p_k and C(k).

## Next
1. Leaf venation (Blonder 2020 / Ronellenfitsch 2015 data): T share, sides/face, dual p_k,
   C(k), width vs degree at fixed length — the best biological SR candidate.
2. Meristem lineage time-lapse: divisions as Errera-rule iterated partition; test clustering ~0.05.
3. Actin: retune segmentation on nonlinear SIM; stress-fibre/arc regions; request U2OS/RPE1
   images (SFEX, Vignaud 2021). Cortical-compartment sizes (Sadegh 2017) vs model cell-area law.

### 2026-11-15 plants (data: bio/data/VENATION.md, bio/data/PLANT_WALLS.md)
Tools: images/planar.py (label image / vein graph -> junction classes, lines, dual), validated on
rasterised model tissues (images/validate_tissue.py: dual p_k and C recovered, corners 3.97-3.99).

**Leaf venation** (LeafVeinCNN, Xu et al. 2021; 44 SLF leaf segments, 25 species, 283k vein edges;
images/venation*.py, results/venation_*):
- 3-way nodes: max arm angle median ~147 deg (T > 155 deg: 31%); 34% of nodes are free ends
  (free-ending veinlets: not in SR).
- At geometric T's the stem is the thinnest arm 0.78 (random 0.33); at all loop 3-way nodes,
  the arm opposite the widest gap is thinnest 0.71. Geometry and width hierarchy agree.
- HLD (width-ordered loop fusion) split fraction, 20k fusions: mean 0.28, quartiles .16/.30/.41,
  histogram rising toward 0.5 with no small-split spike. SR exact: flat, .12/.25/.37; SR with width
  noise and Voronoi with random widths: spike at 0 (images/hld_null.py); centroid/balanced: all > 0.45.
  => nested sequential subdivision, mildly biased to even splits.
- Dual via width continuation, loop network: P(k>=4,8,16) = .31/.066/.008, C(k) .43,.24,.15,.09,.05,.03,
  pc(log k, log W | log L) = 0.40. BUT the same construction on Voronoi with random widths gives
  P(k>=4)=.30 and C(k) .34->.04, pc 0.32 (images/width_lines_null.py): degree law and C(k) are NOT
  evidence here; only pc slightly above null (SR exact 0.57).

**Shoot apical meristem L1** (Willis et al. 2016, 6 plants, 4 h steps, 76 h; images/meristem.py,
results/meristem_*.json): 1117 divisions, 672 with both hosts identified.
- Corner cuts (new wall on two adjacent sides of the mother): 0/672 (95% upper bound 0.45%).
  SR 0.50, uniform clocks 0.63, centroid 0.41, balanced/Errera 0.07 (model refs). Detector checked:
  44% of neighbour pairs of a cell are adjacent around it (~2/(n-1)).
- New junctions are not T at birth: through-angle median 143 deg (quartiles 135/143/151; >165: 3.6%),
  relaxing to 125-134 deg by 32-36 h in all six plants (toward Y).
- Dual degree: wall identity breaks where lineage is missing (only central L1 tracked); <= 3 cell
  generations in 76 h — too short for a degree law (max k = 6).

**Reading.** Plant cell walls are an iterated division but not SR: the rule forbids corner cuts
(dual clustering -> 0) and junctions relax from ~145 to ~128 deg, so straight lines dissolve.
The class prediction gamma ~ 3 is rule-independent but untestable on these films.
Leaf veins carry the SR hierarchy (thinner veins end on thicker ones; ordered nested loops) with a
balanced bias, but their dual degree law cannot be separated from the line construction.

## Next (plants)
- Longer/whole-organ data: UCBG whole leaves (Dryad 1g1jwsv36, needs browser download) for the full
  vein hierarchy (HLD on whole leaves, oldest veins not cut by the sample boundary).
- Refahi 2021 flower meristems (173 h, all layers) for more generations and a degree law.
- Model: iterated division with zero corner cuts + junction relaxation, to give plant-specific
  predictions (C ~ 0, gamma ~ 3, sides distribution).
