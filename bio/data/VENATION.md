# Public leaf-venation network data for the scale-rich partition test

Purpose: test whether leaf minor venation is an iterated-division partition (Cowan/STIT-like):
veins inserted as T-junctions into existing areoles, older veins thicker. Needs per leaf either
(a) an extracted spatial graph (node coordinates, edge list, per-edge widths) or (b) cleared-leaf
images / binary segmentations with resolvable vein widths.

Surveyed 2026-11-15. "Verified" = record opened (API, file listing or file contents) that day.

## Survey

| # | Source | URL | Species / leaves | Contents | Resolution | Size | License | Suitability |
|---|--------|-----|------------------|----------|-----------|------|---------|-------------|
| 1 | **LeafVeinCNN results** (Xu, Blonder, Jodra, Malhi & Fricker 2021 New Phytol 229:631, doi:10.1111/nph.16923; used by Blonder et al. 2020 New Phytol 228:1796, doi:10.1111/nph.16830) | https://doi.org/10.5281/zenodo.4008361 | 728 leaf segments (724 analysed), ~295 species, 50 families, 8 Bornean plots (Sabah) | **(a) extracted graphs**: one `.xlsx` per leaf with sheets Vein Edges / Vein Nodes (pixel coordinates, edge list, **per-edge width**, length, tortuosity, orientations, loop/tree class), Areole + Polygon (faces) tables, dual graph, **Hierarchical Loop Decomposition** tree with fusion width thresholds; `results.xlsx` = per-leaf summary (80+ metrics) | 3.3614 µm/px (images 1.68 µm/px, downsampled x2) | 10 files; per-site `Results_<SITE>.zip` 149 MB (SLF) - 1.16 GB (ESA), ≈ 6.05 GB all sites; plus ~4.2 GB of PDFs | CC BY 4.0 | **best**. Ready graphs with widths over ~1 x 1.2 cm segments (median useful area 115 mm²), ~6500 edges per leaf; minor veins plus 2°/3° veins crossing the segment. Verified (pilot downloaded, below). |
| 2 | **LeafVeinCNN images** (same paper) | https://doi.org/10.5281/zenodo.4008614 | same 728 leaves | **(b)** per leaf: grey image, CNN vein-probability map, binary masks (useful tissue, ground-truth ROI, hand-traced veins in ROI, big veins) | 3.3614 µm/px | per-site zips 0.60-1.64 GB, ≈ 11.8 GB total | CC BY 4.0 | **best companion**: lets us re-threshold / re-skeletonize the CNN map and check the graph. Verified (pilot downloaded). |
| 3 | Blonder et al. 2019 Ecology 100:e02844, "Leaf venation networks of Bornean trees: images and hand-traced segmentations" | https://ora.ox.ac.uk/objects/uuid:de65fc07-4b8f-4277-a6c4-82836afbdeb3 | 726 leaves, 295 spp., 50 families | full-resolution images `<CODE>_img.png`, `_roi.png`, `_seg.png` (hand-traced veins at apparent width, ROI ≈ 700x700 px ≈ 1.2 x 1.2 mm), `_slc.png` (useful tissue), `_big.png` (veins > ~0.5 mm), `Metadata.csv` (CODE, Family, Species) | 595 px/mm = 1.68 µm/px | one zip, 24.3 GB (24846.9 MB), 3216 files | ODC-BY | good but heavy: same leaves at 2x resolution; hand tracings only cover ~1.4 mm² per leaf (too small for a hierarchy test). **ORA does not honour HTTP Range** (tested: 200, full stream), so no partial extraction; whole 24 GB needed. Verified (metadata, data-paper SI). |
| 4 | Fricker/Blonder et al. 2024, UC Berkeley Botanical Garden, 122 spp. (Dryad; Nature Plants 2025 paper) | https://doi.org/10.5061/dryad.1g1jwsv36 | 120 cleared leaves (ferns, basal angiosperms, monocots, eudicots) | `UCBG_venation_form_data.zip` (1.38 GB): per species 8 CSVs from LeafVeinCNN (veins, nodes, areole edges, areoles, polygon edges, polygons, summary, HLD) with widths (README says mm); `UCBG_venation_form_images.zip` (29.4 GB) cleared images + CNN outputs; traits CSVs | whole leaves, 100 mm macro lens on full-frame camera (µm/px not read; README gives a mm/px `Calibration` column) | 44.4 GB total | CC0 (Dryad) / ODC-BY stated in abstract | **very good second source**: whole leaves (all vein orders, phylogenetically broad), graphs with widths as CSV. Not downloaded: Dryad file downloads now need a logged-in session / API token (API returned 401); download manually in a browser. Verified (Dryad API metadata + file list + README via landing page). |
| 5 | Ronellenfitsch, Lasser, Daly & Katifori 2015 PLoS Comput Biol 11:e1004680 | https://doi.org/10.1371/journal.pcbi.1004680 | 186 leaves, 137 spp. (mostly Burseraceae) | Data statement: "within the paper and SI": SI holds only nesting-fingerprint tables (S1, S2); no images or graphs. Software on GitHub (network_extraction, nesting) | scans 6400 dpi (≈ 4 µm/px) | - | CC BY (paper) | n/a: no network or image deposit; author request only. Verified. |
| 6 | Katifori & Magnasco 2012 PLoS One 7:e37994 | https://doi.org/10.1371/journal.pone.0037994 | few leaves | HLD method paper; no data deposit found | - | - | - | n/a (method; HLD is what #1/#4 already compute). Not verified beyond search. |
| 7 | Green, Little, Price, Wing, Smith, Kotrc & Doria 2014 APPS "Reading the leaves" (Dryad) | https://doi.org/10.5061/dryad.8h022 | 120 leaves, 118 genera, 80 families | one tarball: images, LEAF GUI areole data, scripts | not read | 1.95 GB | CC0 | medium-poor: images + areole measurements; graphs/widths not stated; Dryad download needs browser. Verified (API metadata). |
| 8 | Blonder et al. 2018 (Andes-Amazon reticulation, 136 spp., Dryad) | https://doi.org/10.5061/dryad.33bf108 | 136 spp. | only `Data S1.csv` (species-level trait means) | - | 261 kB | CC0 | n/a: no networks or images. Verified. |
| 9 | Wilf et al. 2021 PhytoKeys 187:93 (Smithsonian NCLC + Axelrod cleared leaves), Figshare+ v2.0 | https://doi.org/10.25452/figshare.plus.14980698 | 30,252 cleared/x-rayed leaf images, family-vetted | images only (whole leaves; image format and pixel size not checked) | whole-leaf photos; minor-vein widths probably not resolvable (not checked) | 65 GB extant leaves (3 zips 19-23 GB) | CC BY 4.0 | poor for widths; good only for coarse whole-leaf topology. Verified (Figshare API). |
| 10 | LEAF GUI (Price et al. 2011 Plant Physiol), NEFI (Dirnberger et al. 2015) | http://www.leafgui.org ; https://nefi.mpi-inf.mpg.de | example images only | software + sample images | - | small | - | n/a. Not verified beyond prior knowledge. |

### Ranking

1. **#1 + #2 (Zenodo LeafVeinCNN)**: the only public set with extracted graphs *and* per-edge widths
   for hundreds of species, plus the CNN maps to redo the segmentation. Downloaded as the pilot.
2. #4 (Dryad UCBG 122 spp.): same pipeline on **whole leaves**, CSV format; best for the age (vein
   order) axis; needs a manual browser download of `UCBG_venation_form_data.zip` (1.38 GB).
3. #3 (ORA full-resolution images, 24 GB): only if 1.68 µm/px widths are needed for the thinnest veins.

## Downloads (pilot: site SLF of #1 and #2)

Directory: `raw/venation_leafveincnn/` (804 MB). All md5 values match the Zenodo API checksums.
Nothing from the archives was executed; inspection used `python3 -I venation_inspect.py`
(helper in this directory: `sheets`, `tocsv`, `pnginfo`, `sha256`).

| File | URL | Bytes | sha256 | md5 (Zenodo) |
|------|-----|------:|--------|--------------|
| Results_SLF.zip | https://zenodo.org/api/records/4008361/files/Results_SLF.zip/content | 148815070 | e48d0a0247f284944f071224902f1847b3071ee4d590b0f2f9212526d701889f | fbf2a4e0c934859008c92f6276f991af |
| results.xlsx | https://zenodo.org/api/records/4008361/files/results.xlsx/content | 1594457 | 723da503f6cec0536723619214e0236aab78fa4e79c8bb19163c62f9cc84c42a | 98983ed434778c61dc089999641efae6 |
| SLF-downsampled_images.zip | https://zenodo.org/api/records/4008614/files/SLF-downsampled_images.zip/content | 600973689 | ea4600c61593066b97205901e5fa9077fa335e8a37f900c53efa1b4c0bbd4099 | be88abf118ef4ba54a98b79500ac6c04 |
| LeafVein_Manual_v2.pdf | https://zenodo.org/api/records/4589053/files/LeafVein_Manual%20v2.pdf/content | 72622848 | 537e87c07cdb3b249632d6469a8aa1f9d13e3aaa493ae395de626d2643cb7a5b | e1555991b48af2381913268b4ab41284 |
| Both_tree_functional_traits.xlsx (species per tree; Zenodo 3247631, CC BY 4.0) | https://zenodo.org/api/records/3247631/files/Both_tree_functional_traits.xlsx/content | 662602 | 429efcf163d5cd4c7cbf8a362f68ae510fd55b6ef9220049f4e7fa1d93ff53de | 7789eeac6e0fc36ff589fcebe5e1c21d |

Pilot content: 44 leaf segments from plot SLF (SAFE LF, logged forest), 25 species; 283,380 vein
edges in total (1340-18218 per leaf, median 6554). Remaining sites (same format):
`Results_{BEL,DAS1,DAF2,BSO,BNT,SER,ESA}.zip` (0.46-1.16 GB each) and
`{BEL,DAS1,DAF2,BSO,BNTa,BNTb,SER,ESAa,ESAb}-downsampled_images.zip` (0.93-1.64 GB each) at the
same Zenodo records. Zenodo honours HTTP Range, so `remote_zip.py list|get` can pull single leaves.

## File formats

### Leaf code
`<PLOT>-T<tree>-B<n><S|SH>[ (R)][BIS]-CLAHEm`, e.g. `SLF-T114-B2S-CLAHEm`: plot, tree tag,
branch number, sunlit (S) / shaded (SH). `CLAHEm` = contrast-equalised image. Species: match
(`PLOT`, `tree`) to (`plot_name_trait_campaign`, `tree_id`) in sheet `Tree_functional_traits` of
`Both_tree_functional_traits.xlsx` (header row is the row whose first cell is `field_name`;
rows above it are SAFE metadata; species written `Genus.species`). All 44 SLF leaves match.

### `Results_<SITE>.zip` -> `<CODE>_results.xlsx` (one per leaf, 8 sheets)

**Units in this deposit: lengths and widths in µm, areas µm², volumes µm³** (pipeline v1.0.7,
`Calibration` = 3.3614 µm/px). The v2 manual's tables say "mm"; that applies to newer runs only.
Checked: Euclidean distance between end nodes x 3.3614 = `Length/Tortuosity`; summed polygon
area = useful leaf area; vein density ≈ 7 mm/mm².

**Vein Edges** (one row per vein segment between junctions/free ends; 30 columns):
- `EndNodes_1`, `EndNodes_2`: `node_ID`s (1-based) of the two ends (undirected).
- `node_Idx_1`, `node_Idx_2`: MATLAB linear pixel index of the ends (column-major, 1-based).
- `Name`: edge id. `Type`: `EL` = edge in a loop, `ET` = edge in a tree part (free-ending
  veinlets and bridges); `EB` (boundary) and `F` (petiole) do not occur in the SLF pilot.
- `Calibration`: µm per pixel (3.3614).
- `Width_initial`, `Length_initial`, `Orientation_initial` (rad): whole segment including the
  part overlapping the thicker vein it joins.
- `N_pix`: skeleton pixels; `M_idx`: linear index of the middle pixel.
- `Intensity`, `Intensity_cv`: image grey level; `Probability`, `Probability_cv`: CNN output,
  0-255 scale.
- `Length` (µm), **`Width` (µm)**, `Width_cv`: excluding the overlap region. Width = mean over
  the skeleton pixels of 2 x Euclidean distance to the edge of the thresholded (full-width) CNN
  vein mask; at each junction the overlap of the thinner vein with the thicker one is excluded
  (manual §2.12). `Weight` = `Width_initial` (internal).
- `Or_ij`, `Or_ji`: orientation (deg) from each end node to the midpoint.
- `Tortuosity` = Length / end-to-end distance. `Area` = π(W/2)² (cross-section, µm²);
  `SurfaceArea` = πWL; `Volume` = Area x L. `R2` = L/r² = 4L/W²; `R4` = L/r⁴ = 16L/W⁴.
- `Ai`, `Aj`: `polygon_ID` of the faces on either side; `Ai == Aj` for tree-like edges inside one
  face, 0 = outside the analysed region.
- Missing widths: `Width` = 0 on 109 / 283,380 pilot edges, with `R2 = R4 = 65535` (sentinel).

**Vein Nodes** (20 columns): `node_ID`, `node_Idx` (linear index), `node_Type` (all `E`),
**`node_X_pix` (column), `node_Y_pix` (row)**, 1-based, on the x2-downsampled image (verified:
`node_Idx = (X-1)*H + Y` with H = image height in px, e.g. 3818), `Calibration`,
`node_Degree`, `node_D0 >= node_D1 >= node_D2` (widths of incident veins, µm),
`node_Strength` (sum of widths), `node_Average`, `node_OD0/1/2` (orientation of the incident
veins, deg, same order), `node_AreaRatio` = (D1+D2)/D0, `node_Symmetry` = D2/D1,
`theta_D1_D0`, `theta_D2_D0`, `theta_D2_D1` (branch angles in deg, 3-way junctions). Pilot degrees:
3 → 160,179; 1 → 83,538 (free-ending veinlets *and* cuts at the segment boundary); 2 → 610;
4 → 365; 5 → 1. For a T-junction test: a T has D0 ≈ D1 (through-vein), D2 smaller, and
`theta_D1_D0` ≈ 180°, the other two ≈ 90°.

**Areole Nodes** / **Polygon Nodes** (18 columns): one row per face. Areole = vein-free region
(complement of the full-width mask); polygon = region up to the skeleton. `node_X_pix`,
`node_Y_pix` (centroid, px), `Area` (µm²), `MajorAxisLength`, `MinorAxisLength`, `Perimeter`,
`MaxDistance`, `MeanDistance` (µm), `Eccentricity`, `Orientation` (deg), `ConvexArea`, `Solidity`,
`EquivDiameter`, `Circularity`, `Elongation`, `Roughness`, `areole_ID` / `polygon_ID`
(= `Ai`/`Aj` of Vein Edges).

**Areole Edges** / **Polygon Edges** (5 columns): dual graph, `EndNodes_1`, `EndNodes_2`
(face ids), `Width` (µm, of the vein separating them), `Name` (= vein edge `Name`), `Calibration`.

**HLD Nodes** (21 columns) and **HLD Edges** (`EndNodes_1` = child, `EndNodes_2` = parent):
Hierarchical Loop Decomposition (Katifori & Magnasco 2012). Nodes 1..P are the P polygons
(`width_threshold` = 0); nodes P+1..2P-1 are fusion events in removal order, with
`width_threshold` (µm) = width of the vein removed, non-decreasing; every fusion node has exactly
2 children. Other columns: `node_Area` (µm², fused polygon), `node_Degree` (number of original
polygons in the subtree), `degree_Asymmetry`, `subtree_degree_Asymmetry`, `area_Asymmetry`,
`subtree_area_Asymmetry`, `node_HS` (Horton-Strahler), shape metrics of the fused polygon,
`VTotLen`, `VTotVol` (remaining network after removal), `MSTRatio`, `Boundary` (1 = touches the
segment edge). This tree is directly the "nested partition ordered by vein width" needed for the
iterated-division test.

### `results.xlsx` (Zenodo 4008361, one sheet, 728 rows)
`File` (leaf code), `TimeStamp`, `MicronPerPixel` (3.3614), `DownSample` (2), `Threshold`
(CNN threshold, 0.25), then summary metrics with three-letter codes: `V`(veins) / `N`(nodes) /
`A`(areoles) / `P`(polygons) + `Tot`/`Loop`/`Tree` + `av`/`md`/`sd`/... + `Len`/`Wid`/`Vol`/...
(e.g. `VNE` edges, `VNN` nodes, `VFE` free ends, `Valpha` meshedness, `VTotL` total length µm,
`TotA` analysed area µm², 13.6-177 mm², median 115 mm²). 4 rows are empty (not analysable).

### `SLF-downsampled_images.zip` -> `<CODE>_ds/` (7 PNGs per leaf, 44 leaves)
All PNGs of one leaf share one size (e.g. 3524 x 3818 px = 11.8 x 12.8 mm at 3.3614 µm/px) and
the pixel grid of the node coordinates.
- `<CODE>_ds_img.png`: 8-bit grey, contrast-equalised cleared-leaf image (veins dark, safranin).
- `<CODE>_ds_cnn.png`: 8-bit, ensemble CNN vein probability x 255. Threshold 0.25 (64) gives the
  full-width mask used for the graph (results.xlsx `Threshold`).
- `<CODE>_ds_mask.png`: 1-bit, useful leaf tissue (= `_slc` of #3); analysis restricted to it.
- `<CODE>_ds_cnn_mask.png`: 1-bit, region where CNN prediction is valid (excludes image borders).
- `<CODE>_ds_roi.png`: 1-bit, hand-traced ground-truth ROI (~560 x 600 px).
- `<CODE>_ds_seg.png`: 1-bit, hand-traced veins at apparent width (only inside the ROI).
- `<CODE>_ds_big_mask.png`: 1-bit, hand-delineated large veins (> ~0.5 mm); often empty.

## Caveats for the analysis

- Segments are ~1 cm² cut midway between midrib and margin: minor veins (areoles, FEVs) plus the
  2°/3° veins that cross the segment, not whole leaves. Vein-order range is therefore limited;
  whole leaves are in #4.
- Widths at 3.36 µm/px: median minor-vein width 27.5 µm ≈ 8 px, 10% quantile 15 µm ≈ 4.6 px;
  thinnest veins are quantised and threshold-dependent (re-threshold `_cnn.png` to test).
- Degree-1 nodes mix true free-ending veinlets and cuts by the `mask` boundary; drop nodes on or
  next to the mask boundary (or edges with `Ai`/`Aj` = 0) before counting junction classes.
- Degree-2 nodes (0.25%) are artefacts of splitting duplicate edges; merge them.
