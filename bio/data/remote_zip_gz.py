"""Fetch selected members of a remote zip via HTTP range requests, store each gzipped.

Usage:
  python -I remote_zip_gz.py URL OUTDIR REGEX [--nogz REGEX2]

Members whose name matches REGEX (re.search) are written to OUTDIR/<member>.gz
(or uncompressed if they also match REGEX2). A manifest line per member is
appended to OUTDIR/MANIFEST.tsv: member, uncompressed bytes, sha256 of the
uncompressed member, stored path.  Only reads bytes; never executes anything.
"""
import gzip
import hashlib
import io
import os
import re
import subprocess
import sys
import zipfile


def ranged(url, start, end):
    return subprocess.run(["curl", "-sfL", "--retry", "5", "-r", f"{start}-{end}", url],
                          check=True, capture_output=True).stdout


def total_size(url):
    hdr = subprocess.run(["curl", "-sL", "-r", "0-0", "-D", "-", "-o", "/dev/null", url],
                         check=True, capture_output=True, text=True).stdout
    for line in hdr.splitlines():
        if line.lower().startswith("content-range"):
            return int(line.split("/")[-1])
    raise RuntimeError("no content-range")


class RemoteFile(io.RawIOBase):
    def __init__(self, url):
        self.url, self.size, self.pos = url, total_size(url), 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def readinto(self, b):
        n = len(b)
        if n == 0 or self.pos >= self.size:
            return 0
        end = min(self.pos + n, self.size) - 1
        d = ranged(self.url, self.pos, end)
        b[:len(d)] = d
        self.pos += len(d)
        return len(d)


def main():
    url, outdir, pat = sys.argv[1], sys.argv[2], re.compile(sys.argv[3])
    nogz = re.compile(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[4] == "--nogz" else None
    zf = zipfile.ZipFile(io.BufferedReader(RemoteFile(url), buffer_size=8 << 20))
    man = open(os.path.join(outdir, "MANIFEST.tsv"), "a")
    root = os.path.realpath(outdir)
    for info in zf.infolist():
        m = info.filename
        if m.endswith("/") or "__MACOSX" in m or not pat.search(m):
            continue
        plain = nogz is not None and nogz.search(m)
        dest = os.path.join(outdir, m if plain else m + ".gz")
        if not os.path.realpath(dest).startswith(root + os.sep):
            raise RuntimeError("unsafe path " + m)
        if os.path.exists(dest):
            continue
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        h = hashlib.sha256()
        tmp = dest + ".part"
        opener = (lambda p: open(p, "wb")) if plain else (lambda p: gzip.open(p, "wb", compresslevel=6))
        with zf.open(info) as src, opener(tmp) as dst:
            while True:
                chunk = src.read(1 << 22)
                if not chunk:
                    break
                h.update(chunk)
                dst.write(chunk)
        os.replace(tmp, dest)
        man.write(f"{m}\t{info.file_size}\t{h.hexdigest()}\t{os.path.relpath(dest, outdir)}\n")
        man.flush()
        print("wrote", dest, info.file_size, flush=True)


if __name__ == "__main__":
    main()
