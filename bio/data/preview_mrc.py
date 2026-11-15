"""Montage preview of 2-D uint16 MRC files (reads raw bytes only). Usage: python -I preview_mrc.py out.png f1.mrc ..."""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
out, files = sys.argv[1], sys.argv[2:]
n = len(files); c = int(np.ceil(np.sqrt(n))); r = int(np.ceil(n / c))
fig, ax = plt.subplots(r, c, figsize=(3*c, 3*r))
for a in np.ravel(ax): a.axis("off")
for a, f in zip(np.ravel(ax), files):
    b = open(f, "rb").read()
    h = np.frombuffer(b[:1024], dtype="<i4")
    nx, ny, nz, mode = h[:4]; ext = h[23]
    dt = {0: "i1", 1: "<i2", 2: "<f4", 6: "<u2"}[int(mode)]
    im = np.frombuffer(b[1024 + ext:], dtype=dt)[: nx * ny].reshape(ny, nx).astype(float)
    lo, hi = np.percentile(im, [1, 99.7])
    a.imshow(np.clip((im - lo) / (hi - lo), 0, 1), cmap="gray")
    a.set_title(f.split("/")[-2] + f" {nx}x{ny}", fontsize=8)
plt.tight_layout(); plt.savefig(out, dpi=60)
