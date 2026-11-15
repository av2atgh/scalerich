# Public F-actin image data for the scale-rich partition test

Purpose: skeletonize the actin network of flat spread cells, classify junctions (T / X / Y),
build the segment graph, and test for an iterated-division partition (T-dominated, broad
segment-length distribution, thicker = older).

Surveyed 2026-11-15. "Verified" means the record was opened (API or file listing) on that date.
"Suitability" is for this purpose only (bundles and their junctions traceable over a large part
of one spread cell).

## Survey

| # | Source | URL | Cell type | Modality | Pixel (nm) | Cells / ROIs | Format | Download size | License | Suitability |
|---|--------|-----|-----------|----------|-----------:|-------------|--------|--------------|---------|-------------|
| 1 | **BioSR, F-actin (linear SIM)** (Qiao et al. 2021 Nat Methods 18:194) | https://doi.org/10.6084/m9.figshare.13264793 (file `F-actin.zip`, id 25944599) | COS-7 (African green monkey kidney, fibroblast-like), live, LifeAct-mEmerald | TIRF-SIM, 1.49 NA; reconstructed SIM ground truth `SIM_gt.mrc` + 9-frame raw SIM at 12 signal levels | 31.3 (SIM, header); 62.6 raw | 51 ROIs (+ `GT_all_a/b.mrc` stacks) | MRC, uint16, 1004x1004 (Cell_024: 1516x1516) | 2.23 GB zip; SIM_gt only = 97 MB (range-extracted) | CC BY 4.0 on figshare (BioSR+ page adds "non-commercial use, cite both papers") | **good** for junction statistics: thin bundles and meshwork resolved at ~100 nm, ventral (TIRF) plane only. FOV 31 µm, so most ROIs show part of a cell; ~5 show a near-whole cell (e.g. Cell_011, 015, 024). Few classical ventral stress fibres; mostly cortical meshwork + some straight bundles (Cell_009/010/040-042). Not human. |
| 2 | **BioSR, F-actin_Nonlinear** (same record) | same, file `F-actin_Nonlinear.zip` (id 25714772) | COS-7, LifeAct-mEmerald | nonlinear SIM (~60 nm); `SIM_gt_a/b.mrc` | 20.1 (header) | 51 ROIs | MRC uint16, 1506x1506 | 4.48 GB zip; SIM_gt_a only ≈ 3.7 MB/cell | CC BY 4.0 | **good**: highest-resolution light data found; FOV 30 µm; denser, noisier skeletons expected; useful to check junction classes are not resolution artefacts (compare with #1). |
| 3 | BioSR+ F-actin and Myosin-IIA | https://zenodo.org/records/7115540 | COS-7 (presumed, same lab/system) | TIRF-SIM raw only (9 phase/angle frames, `RawSIMData_gt.mrc`), no reconstructed SIM | 62.6 raw (31 after reconstruction) | F-actin 28 cells x up to 12 levels; Myosin-IIA ~40 cells x 8 levels | MRC, 512x512x9 | 2.1 GB (F-actin), 1.8 GB (Myosin-IIA) | CC BY 4.0, "non-commercial, cite" | medium: Myosin-IIA marks contractile bundles (age / maturity proxy) but needs a SIM reconstruction step (e.g. fairSIM / HiFi-SIM) first. |
| 4 | BioTISR 2D F-actin (Qiao et al. Nat Biotechnol 2025) | https://zenodo.org/records/13843670 | COS-7, LifeAct | TIRF-SIM and nonlinear SIM, time-lapse, ≥50 ROIs | ~31 / ~20 | ≥50 ROIs, time series | MRC in zip | 15.5 GB (linear), 13.4 GB (nonlinear); range extraction possible | CC BY 4.0 | medium-good: the only SIM set with time; could test "thicker = older" directly by tracking bundle birth. Large. |
| 5 | **EMDB EMD-11870** (Fäßler et al. 2020 Nat Commun 11:6437) | https://www.ebi.ac.uk/emdb/EMD-11870 | mouse fibroblast lamellipodium (extracted, fixed) | cryo-ET tomogram | 1.71 (17.096 Å voxel) | 1 tomogram, 512x720x150 vox (0.88 x 1.23 x 0.26 µm) | MRC (.map.gz) | 35 MB | EMDB: CC0 | bonus/good for single-filament Y-junction (Arp2/3 branch) statistics in a lamellipodium; not stress fibres; filaments must be traced in 3-D. |
| 6 | EMDB EMD-15137 / 15139 / 15141 (Fäßler et al. 2023, B16-F1 WT / ArpC5 KO / ArpC5L KO) | https://www.ebi.ac.uk/emdb/EMD-15137 | mouse B16-F1 melanoma lamellipodium | cryo-ET tomogram | 1.35 | 1 tomogram each, 512x720x200 | MRC .map.gz | 273 MB each | CC0 | medium (bonus): branched network vs Arp2/3 isoform; one tomogram each. |
| 7 | EMDB EMD-1919 (Vinzenz/Small lab) | https://www.ebi.ac.uk/emdb/EMD-1919 | NIH 3T3 lamellipodium | negative-stain ET | 0.75 | 1 tomogram 1866x1792x81 | MRC | not checked (~0.5-1 GB) | CC0 | medium (bonus): classic branched-lamellipodium tomogram. |
| 8 | EMPIAR-12292 | https://www.ebi.ac.uk/empiar/EMPIAR-12292/ | optogenetically induced lamellipodia (cell type not checked) | cryo-ET tilt series + "Final models" | ~Å | not checked | MRC/TIFF | 57.8 GB total | EMPIAR: CC0 | medium: raw tilt series; check whether `data/Final_model` holds filament models (small) before downloading. |
| 9 | EMPIAR-10989 | https://www.ebi.ac.uk/empiar/EMPIAR-10989/ | human RPE1 | cryo-ET of FIB lamellae, with actin and MT annotations | 1.35 (bin 4) | 3 tomograms | MRC | 32.9 GB total (tomograms subset smaller) | CC0 | poor: lamellae cut through cell body, not the spread plane. |
| 10 | Human Protein Atlas, subcellular | e.g. https://www.proteinatlas.org/ENSG00000100345-MYH9/subcellular ; image https://images.proteinatlas.org/1644/62_G1_1_green.jpg | U2OS, U-251MG, A-431 (+ others) | confocal 63x/1.4; channels = nucleus, microtubules, ER, antibody. **No phalloidin channel**; stress fibres only via antibodies to actin-associated proteins (MYH9, TPM, ACTN, PDLIM, MYL) | 83 (from XML `scale`) | many fields, several cells each | JPEG 8-bit (old images 1728 px; TIFF download refused for this image) | ~0.3 MB/field | CC BY-SA 3.0 | medium-poor: human, whole spread cells, but lossy JPEG, crowded fields and the stain is a binding partner, not F-actin. |
| 11 | Zenodo 4011248 (Mazloom-Farsibaf et al. 2021, Lifeact vs phalloidin) | https://zenodo.org/records/4011248 | HeLa, RBL | dSTORM / reversible-binding SMLM, raw camera frames | ~100 camera; ~20 after localization | 6+6 HeLa, 3+3 RBL | HDF5 raw 128x128x(2000-3000) frame sequences | 2-4.7 GB per file, 60 GB total | CC0 | medium-poor: true super-resolution of human cells but needs SMLM fitting and each ROI is ~13 µm. |
| 12 | ShareLoc / NanTex SMLM I | https://shareloc.xyz ; https://zenodo.org/records/17120603 | mixed (U2OS/COS-7 actin dSTORM from ShareLoc) | rendered SMLM crops | ~10-20 | not checked | RAR of images | 106 MB | CC BY 4.0 | medium: crops, not whole cells; ShareLoc itself has full localization tables for actin (render at chosen pixel). |
| 13 | FilamentSensor 2.0 data (Hauke et al. 2023) | https://zenodo.org/records/7314056 | human MSC (Fig2: 32x widefield; Pos12 time series lower mag) | widefield; phalloidin-Atto (Fig2) | ~200 (32x, Axiocam) | 156 fields Fig2 + time series + ground-truth sets | TIF/PNG in one zip | 8.8 GB zip (range extraction works) | CC BY 4.0 | poor: bundles visible but junctions not resolved; contains filament-tracing ground truth that could benchmark a skeletonizer. |
| 14 | Okimura et al. 2024 (stress-fibre contraction rotates cell body) | https://zenodo.org/records/14554179 | keratocytes (Yamaguchi Univ. group; species not verified) | time-lapse fluorescence | not checked | few | TIF in zip | 158 MB | CC BY 4.0 | medium-poor: keratocyte stress fibres, small set, imaging details unverified. |
| 15 | Allen Cell Collection (hiPSC single-cell image dataset) | https://www.allencell.org/data-downloading.html | hiPSC (WTC-11) ACTB-mEGFP, ACTN1-mEGFP lines | spinning-disk 100x/1.25 W, 3-D | 108 | >200k cells (all lines) | OME-TIFF / OME-Zarr (S3, no login) | multi-TB; per-cell crops a few MB | Allen Institute terms of use (non-commercial) | poor: compact epithelial colonies, cells not spread; ventral stress fibres only in edge cells. |
| 16 | JUMP Cell Painting (cpg0016) | https://github.com/jump-cellpainting/datasets ; s3://cellpainting-gallery | U2OS | widefield/confocal 20x; phalloidin mixed with WGA + ConA in AGP channel | ~600-650 | >100k compounds x wells | TIFF on S3 | per image ~2-5 MB | CC0 | poor: 20x, phalloidin not separable from Golgi/membrane stain; bundles not resolved. |
| 17 | IDR | https://idr.openmicroscopy.org (idr0008 Rohn actinome; several YAP/TAZ and Cell Painting screens with phalloidin channel) | HeLa, Drosophila BG2, various | high-content widefield screens | ~150-650 | 10^4-10^6 images | OME-TIFF via OMERO | per image MB | CC BY 4.0 | poor: screening magnification; no high-resolution actin study in the IDR study list (checked `studies.tsv` and search-engine "phalloidin"). |
| 18 | BBBC (Broad) | https://bbbc.broadinstitute.org | MCF-7 (BBBC021 F-actin channel), U2OS (Cell Painting sets) | 20x widefield | ~650 | thousands | TIFF | GB | CC BY 3.0 / CC0 | poor: no high-resolution actin set. |
| 19 | Cell Image Library | e.g. CIL:35283 (NIH 3T3 phalloidin, CC0), CIL:36067 (human skin fibroblast, actin IF) | fibroblasts, various | mostly widefield/confocal, varied | varied | single images | TIFF/JPEG | MB | mostly CC0 / CC BY | medium-poor: individual high-quality images exist but heterogeneous, no series; site blocks automated access (Anubis), manual download. |
| 20 | Confocal/STED live F-actin (Lavoie-Cardinal lab) | https://zenodo.org/records/7908914 | rat hippocampal neurons, SiR-actin | STED | 20 | 884 image pairs | TIFF in zip | 742 MB | CC BY 4.0 | poor: axonal/dendritic actin rings, not spread-cell bundles. |
| 21 | SFEX (Zhang et al. 2017 BMC Bioinf 18:268), Vignaud et al. 2021 Nat Mater, FAST (Aljapur et al. 2025 bioRxiv) | papers only | U2OS TIRF on micropatterns; RPE1 on micropatterns; HeLa confocal | TIRF / confocal | - | - | - | **no public image deposit found** | - | n/a: would need a request to the authors (best human whole-cell stress-fibre images exist here). |

### Ranking

1. BioSR F-actin linear TIRF-SIM (#1): ready-to-skeletonize reconstructed SIM, 51 ROIs, 31 nm pixels.
2. BioSR F-actin nonlinear SIM (#2): same cells' system at 20 nm pixels; a resolution control for junction classes.
3. BioTISR 2D F-actin (#4) for the age question (time-lapse); EMD-11870 (#5) as the filament-level lamellipodium bonus.

Gap: no public, downloadable, high-resolution set of *whole spread human* cells (U2OS/RPE1/fibroblasts)
with phalloidin-resolved ventral stress fibres was found. The closest are HPA (lossy, actin-binding proteins)
and the unpublished image sets behind SFEX / Vignaud 2021 / FAST.

## Pilot download (2026-11-15) — `raw/`, 185.7 MB total (185,724,699 bytes)

All files are data only (MRC, gzip-MRC, XML); nothing was executed from them.

### `raw/biosr/F-actin/Cell_001 … Cell_051/SIM_gt.mrc` (51 files, 105,450,336 bytes)
- Source archive: `F-actin.zip`, https://ndownloader.figshare.com/files/25944599 (2,233,338,886 bytes,
  figshare md5 `80fafab68f26fb71de12be0141face74` — refers to the whole zip, which was not downloaded).
- Members extracted remotely with HTTP range requests by `bio/data/remote_zip.py`
  (`python3 -I remote_zip.py get <url> raw/biosr <members>`); zipfile CRC-32 is checked on extraction.
- Format: MRC header 1024 B, mode 6 (uint16), 1004x1004x1 (Cell_024: 1516x1516), header pixel 31.30 nm.

### `raw/biosr_nonlinear/F-actin_Nonlinear/Cell_001 … Cell_010/SIM_gt_a.mrc` (10 files, 45,370,960 bytes)
- Source archive: `F-actin_Nonlinear.zip`, https://ndownloader.figshare.com/files/25714772
  (4,480,187,050 bytes, figshare md5 `b2884114e5fb034b9261855b3c8dbf77` for the whole zip). Same method.
- Format: MRC uint16, 1506x1506x1, header pixel 20.13 nm.

### `raw/emdb_11870/`
- `emd_11870.map.gz` (34,890,610 bytes), https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-11870/map/emd_11870.map.gz ; `gzip -t` OK.
- `emd-11870.xml` (12,793 bytes), https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-11870/header/emd-11870.xml
- Tomogram 512x720x150, 17.096 Å voxel, mouse fibroblast lamellipodium (extracted, fixed).

### Helper scripts (outside raw/)
- `bio/data/remote_zip.py`: list / extract members of a remote zip by HTTP range requests.
- `bio/data/preview_mrc.py`: montage PNG of 2-D MRC files (numpy + matplotlib, reads bytes only).
  Previews were written to the session scratchpad, not to this directory.

Visual check (montages): linear SIM ROIs show a dense ventral meshwork with straight bundles crossing it;
cell edges visible in about half the ROIs; bright punctate clusters in some (Cell_003, 005, 009, 023, 028)
will need masking before skeletonization.

### SHA-256 of every downloaded file (paths relative to `raw/`)
```
cee0a6afba5f0446a97ce013a229a841d4944ed63384da79f0972cdbed0b084f  biosr_nonlinear/F-actin_Nonlinear/Cell_001/SIM_gt_a.mrc
60f1eb8ff353c461ec105f3da8590d716739dfc60458fe562c5b294cc64cba42  biosr_nonlinear/F-actin_Nonlinear/Cell_002/SIM_gt_a.mrc
5d72f822a4fff9d361f1c7576436af845184ef8b684426c2dda08ed4c78db6e6  biosr_nonlinear/F-actin_Nonlinear/Cell_003/SIM_gt_a.mrc
6d8849dda73c6abcb413bad05536ec910ef6d7d2cea26b91f3487a4fefc34c2f  biosr_nonlinear/F-actin_Nonlinear/Cell_004/SIM_gt_a.mrc
ac2233ae265e80fc374531d2f94aeb82010ea16a29a68bc2ee997db7cda32d59  biosr_nonlinear/F-actin_Nonlinear/Cell_005/SIM_gt_a.mrc
b253db7b77bc5656a37a1059f1f2cc70e3559b74a46d384d8219237e8d1b3b66  biosr_nonlinear/F-actin_Nonlinear/Cell_006/SIM_gt_a.mrc
2076a606b91848ad92b64fd94208971970d97463ee2998a574d75220074077ff  biosr_nonlinear/F-actin_Nonlinear/Cell_007/SIM_gt_a.mrc
a309ec517453bd8320dbd457ae7db1bdb4725c8b6b5b6ee22efafa854def3db6  biosr_nonlinear/F-actin_Nonlinear/Cell_008/SIM_gt_a.mrc
b03183fc87f0e038e82e2553e5cc727ea1fdbdb1529d5c259bb06a82ea37b662  biosr_nonlinear/F-actin_Nonlinear/Cell_009/SIM_gt_a.mrc
504250c9792b8145c1db72ae07d67ef45b9247ca00ae040f8c5d9be9d1a3a13d  biosr_nonlinear/F-actin_Nonlinear/Cell_010/SIM_gt_a.mrc
6af36ad11cdc5e88941395839542e5e491c9058676f68d650c1f60bdecf5b053  biosr/F-actin/Cell_001/SIM_gt.mrc
2542571a03214b857aff48dfd18b37851f9227fd9784b624de2649cebf0f0477  biosr/F-actin/Cell_002/SIM_gt.mrc
bc67307faf2e073960895f91a2dc6575208cf338b71223d257e1190a590ff5f6  biosr/F-actin/Cell_003/SIM_gt.mrc
571c421440eefe7c06f8671415ea115677f14c54c2be52ef41a4b39c92d3684f  biosr/F-actin/Cell_004/SIM_gt.mrc
1be44d82cf835ab5031bc443ee27f6f1b798df99cae0c8e69e884499e35ab5d4  biosr/F-actin/Cell_005/SIM_gt.mrc
eac78ff958e3396e85f8548fa3c23bde8939d5aa2898697cc2e01bb017f3e64e  biosr/F-actin/Cell_006/SIM_gt.mrc
689f58c60cf8bb1409aa1b085bbaed74c9f2f9227b5d5cbfbd0b9d2d16f51725  biosr/F-actin/Cell_007/SIM_gt.mrc
37309c689079ec50c4e714fbf3c2793d4aedfd60d81692c6bd2c8e1277c5cc8c  biosr/F-actin/Cell_008/SIM_gt.mrc
4f39a26da095363fe9a823d6bd9f2397972c7e4f1e851a6990551d1b262564b5  biosr/F-actin/Cell_009/SIM_gt.mrc
e3d451339594ecc80a1610226c9985c192bcb7b154a8c09b027f5ebdcaf7ebe7  biosr/F-actin/Cell_010/SIM_gt.mrc
52b50768fa777edaf238018ea8efde82aab7d92bee94103391e131f50fbb577b  biosr/F-actin/Cell_011/SIM_gt.mrc
5915ed71e76fcfa576f1b321471437d6866a1c355dd96a92d608895918e3c6df  biosr/F-actin/Cell_012/SIM_gt.mrc
3b768e3572f38af533e83d2615195585e7a8dd7ece786382722c349b1b8ec7ba  biosr/F-actin/Cell_013/SIM_gt.mrc
f2cd7eb1c4cab78d24b2eea883ca1a83e7fd6a2bb068c2231a78116b021c1971  biosr/F-actin/Cell_014/SIM_gt.mrc
e7f9e32e6ebf83220d286a921506ece9890f33cd915fbd26fc9fa1170085e3d1  biosr/F-actin/Cell_015/SIM_gt.mrc
4218fc6ade79bdfefe3e7be32f824ee62b460a5bb77376b918bb7d263201c291  biosr/F-actin/Cell_016/SIM_gt.mrc
6cf14297d5b6d461b8247904de7aabce5c3fae172a350afe35197d099ba456f0  biosr/F-actin/Cell_017/SIM_gt.mrc
41db322bf7bfa36fa2de2d54da0d86e8e0a84ba1e2d2e72be64e0c1b515241e2  biosr/F-actin/Cell_018/SIM_gt.mrc
1b738bbf22d7ecbbba35842d60639149a295010c289e95c35e231e06b932001e  biosr/F-actin/Cell_019/SIM_gt.mrc
6e6d1dc16fe096427d7e81f8e754c7a6d4f27e32439190fc6f75472522a43d3c  biosr/F-actin/Cell_020/SIM_gt.mrc
819d3a82d393fa31e2af271d08b8cf53d1bdfd7fe07d8c240cac88f5b1d99f67  biosr/F-actin/Cell_021/SIM_gt.mrc
8dd6543d298d3c37cc1bba79a61e11bc4c4954b857d0720cf185047652992235  biosr/F-actin/Cell_022/SIM_gt.mrc
12e3b79bcaaf1112ce693806715fe3e3d78b7c958c39ebfcbf9b1d8c775c20a9  biosr/F-actin/Cell_023/SIM_gt.mrc
3c342f5c21c2bb1ae0f8d69a7425eab45cefbb5b8753bbb946fbb9712bb58e1d  biosr/F-actin/Cell_024/SIM_gt.mrc
80f8da2ce8cd133653d6d9af384288e780ddaf20b1487001dee33f79a9382f99  biosr/F-actin/Cell_025/SIM_gt.mrc
cf2208d4aca48520845fa0fe18e3ef05227d7307ca2267a390930219da4fc811  biosr/F-actin/Cell_026/SIM_gt.mrc
b64305c8ef132be300b0dc155774696b37b8a486ffe139bfe7caa049d2ff89d8  biosr/F-actin/Cell_027/SIM_gt.mrc
adf2afedf1084adde3d7fa9bb1a6db67a03c3d04e5b70649af157f545f075d63  biosr/F-actin/Cell_028/SIM_gt.mrc
9ebeaa576836c9aee336f61b68c05bbcb335418e1545dd19caf82ba21592379a  biosr/F-actin/Cell_029/SIM_gt.mrc
0df718ac0e1f43273a65e0aead8aab49888f8cf0b1e3496d9f54889827957ec4  biosr/F-actin/Cell_030/SIM_gt.mrc
c504f3c01ff6a4852fabfe6bbfd52fdef48ec54f33852f56e250b71232d452e1  biosr/F-actin/Cell_031/SIM_gt.mrc
a4965c16411c1ddecca014ced0b2bbd9413d51508b43a3a43991364c1e7ddc33  biosr/F-actin/Cell_032/SIM_gt.mrc
0a808c6b7d04d3c7713be6c2daea63d1a54dfbe000412759f44deeb5f71da140  biosr/F-actin/Cell_033/SIM_gt.mrc
5c15bf473015331559f866ed0774702e29413a0836823c1bcab868d678c65db7  biosr/F-actin/Cell_034/SIM_gt.mrc
17e5e9fd713d0bfbb389995469ad451bc8be0738cba8d9152eacf6f3d732f922  biosr/F-actin/Cell_035/SIM_gt.mrc
3ca835bd84f84ab56c02ac9026c75209cc8e912465f79afdcdee4dbdde11a429  biosr/F-actin/Cell_036/SIM_gt.mrc
33f2bf5d787624eb3d576b02efa5fddf90b2297fb1e1b833f27423c9cdd1212a  biosr/F-actin/Cell_037/SIM_gt.mrc
7d7442c0ee3c3b64963715dfbb6d4eff01b7b65955c9357f0fa6253ddd91f540  biosr/F-actin/Cell_038/SIM_gt.mrc
7558437e361969924bc2e445c19042260c2466397e00c410c5d6467a1e36342a  biosr/F-actin/Cell_039/SIM_gt.mrc
1b6741529f0d33feabba971aabafe50631b2090eb44d7fc46520a9f52d4e7fec  biosr/F-actin/Cell_040/SIM_gt.mrc
9a7871ff5fc5d468b92fb074ffaf22f14dd3c1ff7afc56fed13888756fd816fa  biosr/F-actin/Cell_041/SIM_gt.mrc
c39b77664620176d92ce61e1423e37a9b3c30a1742479f7e737bb44f176869f5  biosr/F-actin/Cell_042/SIM_gt.mrc
aa6c5767e130d4d1b9057c0f5b7e007a56d0baf436b51e1b235588cfd81822ad  biosr/F-actin/Cell_043/SIM_gt.mrc
531784a8701f23c2927dd2a79e53c2ac4d5d2b387bce201812c65f48cf9b9fe2  biosr/F-actin/Cell_044/SIM_gt.mrc
de21df825a9d18cbd4a2c74a81bd9549944d0439af68e5b92a573c7c750e8bb2  biosr/F-actin/Cell_045/SIM_gt.mrc
c2bfba9a7fbf0a3dd66a691a16b9609b07d6e4cc78cefe94b6b968c43ff966ee  biosr/F-actin/Cell_046/SIM_gt.mrc
7822659acc982e7e64bac135eee50bd268a622536a68b2a6656f3f67c29325d4  biosr/F-actin/Cell_047/SIM_gt.mrc
b287e3fd324d6eab0615bedc34eb4dadf619ab1d2191749a21b43fdba922ec45  biosr/F-actin/Cell_048/SIM_gt.mrc
ce3eedbb68aa8f8050fe0dbdd1c29bf2f6d42f9935c5f12695704bc6b4379cc8  biosr/F-actin/Cell_049/SIM_gt.mrc
800472c435901d085fd034fe931e91bbf4a3a21af174f7f27434e0b60a28047f  biosr/F-actin/Cell_050/SIM_gt.mrc
31124a0b8b8004bebd5a2e0f0a854361432ef9ef2a554d2e1dc6743ae78bc23c  biosr/F-actin/Cell_051/SIM_gt.mrc
6a1821456a6c5f384d740d07d1d294343ce602dd5a6f65ebfaeb1ba2e6dc941f  emdb_11870/emd_11870.map.gz
d478261ab54831273817fff1ec15068cb9ac57c37ebd69631fdd93528797180c  emdb_11870/emd-11870.xml
```
