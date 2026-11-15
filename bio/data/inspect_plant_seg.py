"""Inspect a (gzipped) 3D label TIFF: shape, dtype, labels, tags, L1 estimate.

Usage: python -I inspect_plant_seg.py FILE.tif[.gz] [FILE2 ...]
Reads bytes only (tifffile); nothing from the data is executed.
"""
import gzip
import io
import sys

import numpy as np
import tifffile


def load(path):
    raw = gzip.open(path).read() if path.endswith(".gz") else open(path, "rb").read()
    return tifffile.TiffFile(io.BytesIO(raw))


for p in sys.argv[1:]:
    with load(p) as tf:
        a = tf.asarray()
        pg = tf.pages[0]
        tags = {t.name: t.value for t in pg.tags.values()
                if t.name in ("XResolution", "YResolution", "ResolutionUnit", "ImageDescription",
                              "BitsPerSample", "SampleFormat", "Software")}
        print("==", p)
        print("  pages", len(tf.pages), "series shape", a.shape, a.dtype,
              "imagej" if tf.is_imagej else "", tf.imagej_metadata if tf.is_imagej else "")
        for k, v in tags.items():
            s = str(v)
            print("  tag", k, s[:300])
        lab, cnt = np.unique(a, return_counts=True)
        print("  n labels", len(lab), "min", lab.min(), "max", lab.max())
        order = np.argsort(cnt)[::-1][:3]
        print("  largest labels (label, voxels):", [(int(lab[i]), int(cnt[i])) for i in order])
        # corner voxels tell which label is background
        corners = [a[z, y, x] for z in (0, -1) for y in (0, -1) for x in (0, -1)] if a.ndim == 3 else []
        print("  corner labels", sorted(set(int(c) for c in corners)))
