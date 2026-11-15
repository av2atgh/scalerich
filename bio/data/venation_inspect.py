"""Inspect the LeafVeinCNN venation downloads (Zenodo 4008361 / 4008614) without executing them.

Usage (always with python3 -I, run from outside raw/):
  python3 -I venation_inspect.py sheets  RESULTS_ZIP MEMBER.xlsx      # sheet names, sizes, headers, 2 rows
  python3 -I venation_inspect.py tocsv   RESULTS_ZIP MEMBER.xlsx OUTDIR  # one CSV per sheet
  python3 -I venation_inspect.py pnginfo IMAGES_ZIP PREFIX              # size/mode/unique values of PNGs
  python3 -I venation_inspect.py sha256  FILE [FILE ...]

Reads bytes only (zipfile + openpyxl + PIL); nothing in the archives is executed.
"""
import csv
import hashlib
import io
import os
import sys
import zipfile


def _wb(zpath, member):
    import openpyxl
    data = zipfile.ZipFile(zpath).read(member)
    return openpyxl.load_workbook(io.BytesIO(data), read_only=True)


def sheets(zpath, member):
    wb = _wb(zpath, member)
    for ws in wb.worksheets:
        rows = list(ws.iter_rows(values_only=True))
        print(f"== {ws.title}: {len(rows) - 1} rows x {len(rows[0]) if rows else 0} cols")
        for r in rows[:3]:
            print("  ", r)


def tocsv(zpath, member, outdir):
    os.makedirs(outdir, exist_ok=True)
    stem = os.path.basename(member).replace("_results.xlsx", "")
    wb = _wb(zpath, member)
    for ws in wb.worksheets:
        fn = os.path.join(outdir, f"{stem}__{ws.title.replace(' ', '_')}.csv")
        with open(fn, "w", newline="") as f:
            w = csv.writer(f)
            for r in ws.iter_rows(values_only=True):
                w.writerow(r)
        print("wrote", fn)


def pnginfo(zpath, prefix):
    from PIL import Image
    import numpy as np
    Image.MAX_IMAGE_PIXELS = None
    z = zipfile.ZipFile(zpath)
    for n in z.namelist():
        if n.startswith(prefix) and n.endswith(".png"):
            im = Image.open(io.BytesIO(z.read(n)))
            a = np.asarray(im)
            u = np.unique(a)
            print(n, im.size, im.mode, a.dtype, a.shape,
                  "unique:", u[:6].tolist(), "... n_unique =", len(u))


def sha256(paths):
    for p in paths:
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk)
        print(h.hexdigest(), os.path.getsize(p), p)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "sheets":
        sheets(sys.argv[2], sys.argv[3])
    elif cmd == "tocsv":
        tocsv(sys.argv[2], sys.argv[3], sys.argv[4])
    elif cmd == "pnginfo":
        pnginfo(sys.argv[2], sys.argv[3])
    elif cmd == "sha256":
        sha256(sys.argv[2:])
