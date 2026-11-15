# Plant cell-wall partitions: public segmented data

Purpose: test whether plant cell walls form an iterated-division partition (each new wall
inserted as a T-junction into its mother cell; Errera-like). Best data = surface (L1) or 2D
segmentations with cell lineage from time-lapse, so each wall's birth interval is known.

Survey date 2026-11-15. All downloads are untrusted data: nothing in `raw/` is executed;
inspection only with `python3 -I` / `bio/.venv/bin/python -I` from scripts in `bio/data/`.

## Candidate sources

| # | Source | URL | Species / organ | 2D/3D | Lineage | Time points | Cells | Format | Size | License | Suitability |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **Willis et al. 2016 PNAS** 113:E8238 (data: Jönsson, Willis, Refahi) | https://doi.org/10.17863/CAM.7793 (= repository.cam.ac.uk/handle/1810/262530) | *A. thaliana* shoot apical meristem (NPA-treated, pin-like, no organs), 6 plants | 3D label stacks (several layers); lineage tracked for L1 | **yes**, every 4 h | 19–22 per plant, 0–76/84 h (plant18 lacks 40 h) | ~750 labels/stack (all layers); ~200–400 tracked L1 mothers per step; ~3100 multi-daughter events in total | uint16 TIFF labels + txt lineage + txt resolutions | PNAS.zip 4.89 GB (seg+tracking members only ≈ 215 MB compressed) | CC BY 4.0 | **Best.** Densest time sampling (4 h < cell cycle), so nearly every wall's birth step is known; flat-ish dome surface; 6 replicates. Needs L1 surface extraction from 3D. |
| 2 | **Refahi et al. 2021 Dev Cell** 56:540 | https://doi.org/10.17863/CAM.61991 (= handle/1810/318119) | *A. thaliana* early flower meristem FM1–FM6 (+3 *lfy-12*) | 3D label stacks; lineage for all cells (all layers) | **yes**, 4–13 h steps | FM1 18 (0–132 h), FM2 6, FM3 7, FM4 7, FM5 15 (to 173 h), FM6 5 | 446 labels at FM1 0 h → >1000 tracked daughters late; 1000–2900 per late stack | uint16 TIFF labels + txt lineage | wt.zip 1.73 GB (seg+tracking ≈ 111 MB compressed); leafy.zip 250 MB | CC BY 4.0 | **Very good.** Long series, curved organ (sepal primordia, boundaries); intervals sometimes ≥ cell cycle (3–6 daughters per mother), so some wall birth orders are ambiguous. |
| 3 | **Jackson et al. 2019 Cell Systems** 8:53 | https://osf.io/mxgrh/ (WT 1.59 GB zip, *ktn1* 2.68 GB zip); https://osf.io/ucd76/ (empty) | *A. thaliana* SAM WT Col (4 reps) and *katanin1* | 3D (4 cell layers); MorphoGraphX | partial: parents for ~140 central cells per step | 3 (0, 11, 22 h) | ~133 cells in topology table (~45 in L1); network ~13k directed rows | CSV (adjacency with wall areas; per-cell topology; parents) + `.mgxs` stacks + `.mgxm` meshes (75–116 MB compressed each) + raw OME-TIFF | 1.59 GB WT | none stated on OSF node | Good ready-made adjacency graphs and layer labels; lineage sparse, 11 h steps. Meshes would give junction coordinates. |
| 4 | Lin et al. (Kierzkowski lab) 2026, moss phyllids | https://doi.org/10.5061/dryad.kkwh70sjp (+ companion Zenodo, models/meshes) | *Physcomitrium patens* (and other mosses) leaf-like phyllids — single-layer lamina | 2.5D MorphoGraphX surface | yes (`PARENT` CSVs, Fig S1 A–C) | 3 series/genotype | n/a | MorphoGraphX CSV exports + XLSX; no segmented meshes in Dryad | 3.8 MB | Dryad (CC0 by Dryad policy) | Lineage and per-cell scalars only, **no wall geometry**; flat lamina is ideal topologically if the meshes were obtainable (ask authors). |
| 5 | PlantSeg (Wolny et al. 2020 eLife) | https://osf.io/uzq3w (children: LateralRootPrimordia 2rszy, Ovules w38uf, InflorescenceMeristem 295su, CLSMLeaf kfx3d); BioImage Archive S-BIAD1367 | *A. thaliana* LRP (light-sheet time-lapse), ovules, SAM, leaf | 3D | no (LRP is time-lapse but untracked); CLSMLeaf has T1/T2 pairs (lineage not verified) | LRP every 30 min | thousands | HDF5 (`raw` uint8, `label` uint16); leaf: TIFF + `.mgxm` | LRP/ovules several GB; CLSMLeaf ≈ 1.8 GB | CC BY 4.0 (paper) | Benchmark-quality 3D labels; no lineage → only static topology tests. LRP is 3D and highly anisotropic. |
| 6 | Galvan-Ampudia et al. 2020 eLife | https://zenodo.org/records/3737795 | *A. thaliana* SAM | 3D raw CZI | no | time series | — | raw CZI only | large | CC BY | Raw only; would need own segmentation. Not useful now. |
| 7 | 3DCoordX (BioStudies S-BSST734) | https://www.ebi.ac.uk/biostudies/studies/S-BSST734 | *Marchantia polymorpha*, *Utricularia gibba* organs | 3D meshes | no | single | — | meshes | — | per BioStudies | Static Marchantia meshes; possible flat-tissue check without lineage. |
| 8 | SSBD repos-000290 | https://ssbd.riken.jp/repository/290/ | Marchantia gemmalings | 2D bright-field | no | time-lapse | — | zipped TIFF | 32.2 GB | per SSBD | No segmentation, low cellular contrast. Not suitable. |
| 9 | Fernandez et al. 2010 Nat Methods (MARS-ALT) | — | *A. thaliana* flower | 3D | yes in paper | — | — | — | — | — | Public deposit not located; superseded by #2 (same group). |
| 10 | Hong et al. 2016 Dev Cell; Hervieux et al. 2016 Curr Biol; Tsugawa et al. 2017 (sepals) | — | *A. thaliana* sepal abaxial epidermis | 2.5D MorphoGraphX | yes in papers | 6–12 h steps | — | — | — | — | **No public deposit found** (data statements point to authors). Would be ideal (flat 2.5D, lineage). Request from Roeder / Hamant labs. |
| 11 | Besson & Dumais 2011 PNAS 108:6294 | PMC3076879 | many species (Chara, Coleochaete, Arabidopsis etc.), single dividing cells | 2D | n/a | — | — | figures only | — | — | No data deposit; figures only. |
| 12 | Dupuy, Mackenzie & Haseloff 2010 PNAS 107:2711 | PMC2823878 | *Coleochaete scutata/orbicularis* disc thallus | 2D | yes (tracked 80 h in figure) | — | — | SI movies / figures only | — | — | Ideal geometry (2D disc, straight walls) but no segmented deposit found; lab page cert broken. Digitising SI figure is the only route. |
| 13 | Harline & Roeder 2023 Plant Methods (leaf pipeline); Kierzkowski et al. 2019 Cell (leaf) | — | *A. thaliana* / *C. hirsuta* leaf epidermis | 2.5D MorphoGraphX | yes in papers | — | — | — | — | — | No public mesh/lineage deposit located. |

## Downloads (pilot, total stored 344 MB; 4.9 + 1.7 + 1.6 GB archives were *not* fetched whole)

Method: `bio/data/remote_zip_gz.py` reads only the selected members of each remote zip via
HTTP range requests (no full archive download), writes `.txt/.csv` as-is and every other member
gzipped (`<member>.gz`). Each `raw/plant_*/MANIFEST.tsv` lists
`member  uncompressed_bytes  sha256(uncompressed member)  stored_path` (tab-separated);
`SHA256SUMS.stored` holds `shasum -a 256` of every stored file.

| Dir | Source archive URL | Members fetched | Stored | Uncompressed | Index sha256 |
|---|---|---|---|---|---|
| `raw/plant_willis2016/` | https://www.repository.cam.ac.uk/bitstreams/95babfae-ae16-4dae-82f5-8656ceace4d6/download (PNAS.zip, 4.89 GB); README also at …/bitstreams/569a9ec3-7d44-4abe-8f45-3f5004c1f200/download | `PNAS/README.txt`, `PNAS/plantX/plantX_tiff_resolutions.txt` (6), `PNAS/plantX/segmentation_tiffs/*.tif` (125), `PNAS/plantX/tracking_data/*.txt` (121); X ∈ {1,2,4,13,15,18} | 208 MB, 253 files | 20.85 GB | MANIFEST.tsv d036bfcd…7da0d7f34; SHA256SUMS.stored f48cc244…e5dd7bc |
| `raw/plant_refahi2021/` | https://www.repository.cam.ac.uk/bitstreams/59bea6a3-06b4-40db-9eaa-2a4981dfb00d/download (wt.zip, 1.73 GB); README …/bitstreams/7402f9cd-ff36-4375-ba96-cdfd3b1fa66b/download | `wt/FM{1..6}/segmentation_tiffs/*` (59 tif + `FM4/.../L1Width.pkl`), `wt/FM{1..6}/tracking_data/*.txt` (51) | 118 MB, 111 files | 6.27 GB | MANIFEST.tsv e0b84d95…c2eae; SHA256SUMS.stored ee18b23a…cd738 |
| `raw/plant_jackson2019/` | https://osf.io/download/h587k/ → https://files.osf.io/v1/resources/mxgrh/providers/osfstorage/5b4ca1dbd75a110010d29405 (WT Arabidopsis.zip, 1 589 927 599 B, OSF sha256 cbf8d2f9e46f4de4988eae8e0c9a444bee7571525397f80d385788bbcff552de) | `WT Arabidopsis/README.txt`, all `*.csv` (28), all `*.mgxs` (12); not `.mgxm` meshes nor raw OME-TIFF | 9.6 MB, 41 files | 14 MB | MANIFEST.tsv 937282db…f0b1ef4; SHA256SUMS.stored 5c598f1d…5bb31e7 |

Not downloaded (available if needed): Willis raw membrane `processed_tiffs` and nuclear
segmentations; Refahi raw `tiffs/`, `FM1/geneExpression/`, `leafy.zip` (250 MB); Jackson
`.mgxm` meshes (12 × 75–116 MB compressed) and *ktn1* zip (2.68 GB).

**Do not unpickle** `raw/plant_refahi2021/wt/FM4/segmentation_tiffs/L1Width.pkl.gz`
(Python pickle = arbitrary code execution); it is an auxiliary L1-width estimate and not needed.

## Formats

Reading a stored label stack: `tifffile.TiffFile(io.BytesIO(gzip.open(p).read())).asarray()`
(see `bio/data/inspect_plant_seg.py`, which prints shape, dtype, label count, tags, corner labels).

### Willis 2016 (`raw/plant_willis2016/PNAS/`)
- Plant folder names are the original ids; article numbering: plant1→1, plant2→2, plant4→3,
  plant13→4, plant15→5, plant18→6.
- `segmentation_tiffs/<T>hrs_plant<X>_trim-acylYFP_hmin_2_asf_<a>_s_<s>_clean_3.tif.gz`:
  multipage uint16 TIFF, array order (z, y, x), e.g. plant18 0 h = (134, 512, 512).
  **Label 1 = background** (all corners are 1, largest label); cells = unique integers 2…N
  (plant18 0 h: 756 cells, max label 757). Watershed segmentation (h-min 2) of acyl-YFP plasma
  membrane; walls are not a separate label (labels tile the tissue, membrane voxels assigned
  to cells). Several cell layers (L1, L2, L3…) are present; L1 = cells touching label 1 on the
  dome side. TIFF resolution tags are meaningless (written by pylibtiff) — use the
  resolutions file.
- `plantX_tiff_resolutions.txt`: Python-literal list `[["0hrs", (x, y, z), zfactor], …]`
  (parse with `ast.literal_eval` after dropping `#` lines, never `eval`). x = y ≈ 0.2196 or
  0.2396 µm/px, z = 0.26 µm nominal; true z spacing = z / zfactor (zfactor 0.97–2.8,
  Table S10 of the paper), i.e. stacks are not isotropic.
- `tracking_data/<T>hrs_to_<T+4>hrs.txt` (one 36→44 h step in plant18): one line per tracked
  mother, `mother_label: d1[, d2, …]`, labels referring to the segmentation at T and T+4.
  One daughter = same cell; two = division during the interval. A few lines have 3 daughters
  (two divisions in 4 h or segmentation over-split; filter). Only tracked cells are listed
  (≈200–400 per step, L1 cells in the central region per the paper).
  Lines per plant (all steps): 7163/4815/3782/6164/5472/6079; multi-daughter lines
  964/487/340/478/390/460 for plants 1/2/4/13/15/18.

### Refahi 2021 (`raw/plant_refahi2021/wt/FMk/`)
- `segmentation_tiffs/<T>h_segmented.tif.gz`: multipage uint16 TIFF, (z, y, x), e.g. FM1 0 h
  = (280, 400, 400), written by tifffile (`ImageDescription {"shape": …}`).
  **Label 1 = background**; FM1 0 h has 446 cells (labels 2…447). X/Y resolution tag
  3154279/761859 px per unit → 0.2415 µm/px; stacks were resampled (MARS fusion) and appear
  isotropic, but the authors applied a later voxel-size correction (code
  `estimateVoxelSize` in https://gitlab.com/slcu/teamHJ/publications/refahi_etal_2020); for
  topology/angles this scale is irrelevant. All layers segmented and tracked.
  Time points: FM1 0,10,18,24,32,40,48,57,64,72,81,88,96,104,112,120,128,132 h; FM2
  0,10,23,33,45,55; FM3 10,23,33,45,50,55,65; FM4 0,10,23,33,45,50,55; FM5 16…173 h (15);
  FM6 6 stacks to 69 h.
- `tracking_data/<T1>hrs_to_<T2>hrs.txt`: per line `mother: d1, d2, …` (comma+space separated),
  mother label in image T1, daughters in image T2. **Drop the line for mother 1**: it maps
  background to newly appearing (untracked/entering) cells. Remaining lines in FM1 0→10 h:
  313 with 1 daughter, 61 with 2, 13 with 3, single lines with 4–6 (multiple rounds or
  errors).
- The code repo's `data/FMk/tv/*_tvformat_neighbors.txt` and `*_volume_position.txt` hold
  precomputed neighbour lists and cell volumes/positions for some time points (not
  downloaded; text, small).

### Jackson 2019 (`raw/plant_jackson2019/WT Arabidopsis/Rep r/T{0,1,2}/`)
- T0 = 0 h, T1 = 11 h, T2 = 22 h; files prefixed with image ids 02/06/08/09 (Reps 1–4).
- `<id>_T<t>_network.csv`: header `label, neighbour, wall area`; one row per directed pair
  (both orientations listed), wall area in µm² from the MorphoGraphX 3D mesh; 4 cell layers.
- `<id>_T<t>_topology.csv` (T0, T1 only): `Label, Log_Betweenness_Centrality,
  Random_Walk_Centrality, Degree, Division, Layer, Area, Volume, Anisotropy`; ~130–150
  central cells in layers 1–3 (~45 per layer, so only ~45 L1 cells); `Division` (0/1) matches
  the 2-child parents in the next parents file (Rep 1 T0: 5 divisions), i.e. divides before t+1.
- `<id>_T<t>-T<t+1>_parents.csv` (CRLF line ends): `Label, Parent Label, Layer` = label at
  t+1, its parent label at t, layer. Rep 1 T0→T1: 128 parents with 1 child, 5 with 2.
- `*.mgxs.gz` → MorphoGraphX stack `MGXS 1.3`: binary header (contents not decoded here;
  presumably voxel steps, origin, dimensions), then a zlib stream (starts `78 5e` at byte 0x3f) holding uint16 labels.
  Parse manually if needed (MorphoGraphX source documents the header); do not run MGX
  scripts from the archive.

## Assessment

- **Best source: Willis 2016** — six SAMs, 4 h steps over 76–84 h, CC BY, and lineage for L1
  cells: each L1 wall born during the series can be assigned a birth interval, which is the
  datum the iterated-division test needs. Work needed: extract the L1 surface (cells adjacent
  to background label 1 on the dome side), build the L1 adjacency/junction graph from shared
  boundaries (3-cell junctions = vertices), and map new walls through `mother: d1, d2`.
- **Second: Refahi 2021** — longer and larger series; whole-organ lineage including L2/L3,
  curved geometry, but coarser steps (multi-round divisions) and mother-1 artefact lines.
- **Jackson 2019** gives ready adjacency graphs with wall areas and layer labels (useful for a
  quick topological check, e.g. 3-valent T-junction statistics in L1), but only 3 time points.
- Truly flat 2D tissues with lineage (sepal, leaf, moss phyllid, Coleochaete) are not public
  as geometry; those would require asking the Roeder, Hamant, Kierzkowski or Haseloff labs.

## Helper scripts (in `bio/data/`, never in `raw/`)
- `remote_zip_gz.py` — range-request member extraction, gzip storage, manifest with sha256.
- `inspect_plant_seg.py` — print shape, dtype, labels, TIFF tags, corner labels of a label stack.
