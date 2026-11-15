"""List / extract selected members of a remote zip via HTTP range requests.

Usage:
  python -I remote_zip.py list URL
  python -I remote_zip.py get URL OUTDIR member1 [member2 ...]

Only reads bytes; never executes anything from the archive.
"""
import io
import os
import subprocess
import sys
import zipfile


def ranged(url, start, end):
    out = subprocess.run(
        ["curl", "-sfL", "--retry", "3", "-r", f"{start}-{end}", url],
        check=True, capture_output=True).stdout
    return out


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
        self.cache = {}

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        key = (self.pos, end)
        if key not in self.cache:
            self.cache = {key: ranged(self.url, self.pos, end)}
        data = self.cache[key]
        self.pos += len(data)
        return data

    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)


def main():
    cmd, url = sys.argv[1], sys.argv[2]
    rf = io.BufferedReader(RemoteFile(url), buffer_size=1 << 20)
    zf = zipfile.ZipFile(rf)
    if cmd == "list":
        for i in zf.infolist():
            print(i.file_size, i.compress_size, i.filename)
    elif cmd == "get":
        outdir = sys.argv[3]
        for m in sys.argv[4:]:
            info = zf.getinfo(m)
            dest = os.path.join(outdir, m)
            if not os.path.realpath(dest).startswith(os.path.realpath(outdir)):
                raise RuntimeError("unsafe path " + m)
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with zf.open(info) as src, open(dest, "wb") as dst:
                while True:
                    chunk = src.read(1 << 22)
                    if not chunk:
                        break
                    dst.write(chunk)
            print("wrote", dest, info.file_size)


if __name__ == "__main__":
    main()
