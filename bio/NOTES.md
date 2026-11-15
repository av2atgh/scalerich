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
