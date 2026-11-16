# classes/ — universality classes of recursive cell division (periodic, quasiperiodic, random)

Started 2026-11-15. Question: do quasiperiodic tessellations sit in the SR family with no change
of model? Answer: yes, through the split law (area fraction r of the low-side daughter).

## Code
- `qdiv.py` — exact rectangle divisions, rule (pick, split, 'perp'); pick 'largest' | 'area';
  split float r | ('side', r) | ('choice', r1, r2) | 'uniform'. Output keys = beyond/divisions.run.
- `shapes.py` (N_s vs T), `diffraction.py` (2D Bragg search), `axis_bragg.py`, `pisot.py`
  ((p,q) commensurate splits), `netq.py` (netstats features + crossing fraction chi),
  `mech.py` (framed periodic index, merged crossings, spanning lines m). Results in results/.
- `manuscript/` — manuscript.tex (revtex prl, no sections, pdflatex + bibtex), figs.py -> fig1-3.pdf.
  4 pp (2026-11-15 PRL format; abstract ~600 chars; refs abbreviated with DOIs verified via Crossref).

## Results
- largest-first + fixed r: areas r^a (1-r)^b in window [r_min A_max, A_max]. Commensurate
  r = theta^-p, 1-r = theta^-q (theta root of x^p - x^(p-q) - 1) -> finite shapes N_s
  (2 periodic, 4-5 golden, 8 plastic, 12-32 others); generic r = 0.3: N_s, areas ∝ ln T.
- Pisot (irrational) for p <= 7 only (2,1) tau, (3,1) 1.466, (3,2) plastic 1.325, (4,1) 1.380.
- Axis Bragg I_max (window q/2pi in [0.25, 2] sqrt T), T 4e3..2e6: golden 0.78 flat (peaks at
  Fibonacci q), (3,1) 0.41-0.47, (4,1) 0.19-0.27; plastic 0.11 -> 0.037 (undecided); non-Pisot
  and r = 0.3 decay to 0.01-0.03; uniform split ~1/N. Random-side golden keeps 0.63-0.71.
  (Fixed-q tracking is ill-posed: the division refines inside a fixed box; dropped.)
- Dual graph identical across orders: K^2 P(k>=64) 3.4-6.4, h/t ~ U^2 moments, C = 0 for all
  rectangle rules. What differs: crossing fraction chi (0.98 periodic, 0.79 golden, 0.85 plastic,
  0.75 r=0.3, 0.54 random side, 0 random splits) -> spanning lines m ≈ 1.4-1.9 sqrt T.
- Index (T = 1e3, alpha 0): deterministic splits <n> 1.6-1.9 at every density = pinned value
  (stretch dominated), T-independent; random side 2.8 (rho .1), ceiling 3.1/3.6/4.5 at
  T = 300/1e3/3e3 (m ~ T^0.22 -> n_max ~ T^0.14, measured 0.16); uniform split 5.5, ceiling ∝ T^1/4.
- Manuscript reframed 2026-11-15 as "Universality classes of recursive cell division": one dual-graph
  class; four order classes (periodic / quasiperiodic (Pisot) / multiscale / random); elastic classes by
  mu in m ~ T^mu, n_max ~ T^((1/2-mu)/2): mu fitted .51/.46/.50/.51 (det.), .23 (random side), 0;
  pinned-n exponents -.01 (golden), .16 (side, pred .14), .249 (uniform, pred .25).
- Message: order and high index exchanged through m. Rectangle members fourfold only (eta .6-.7).

## Open
- plastic Bragg weight at larger T; Bloch edges on Fibonacci approximants; cut angles from a
  finite set beyond the axes (pinwheel-type rotational order).
- New refs: DOI, title, journal, volume, pages verified via Crossref 2026-11-15 (not the cited statements).
