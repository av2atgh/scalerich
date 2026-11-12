"""Compile sr3d.c and run the 3D plate configurations, one at a time.

usage: python3 run.py [name ...]   (default: all). Raw output in data/ (not committed).
"""
import os, subprocess, sys, time

D = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(D, 'build', 'sr3d')

RUNS = {
    # name: (T, seed, orient, lam0, alpha)
    'iso_s1': (2_000_000, 1, 'iso', 0, 0),
    'iso_s2': (2_000_000, 2, 'iso', 0, 0),
    'iso_s3': (2_000_000, 3, 'iso', 0, 0),
    'iso_s4': (200_000, 4, 'iso', 0, 0),
    'axes_s1': (1_000_000, 1, 'axes', 0, 0),
    'axes_s2': (1_000_000, 2, 'axes', 0, 0),
    # thickness: tunable (alpha > 1/3) and jammed (alpha < 1/3)
    'lam_a0.5_l0.02': (300_000, 1, 'iso', 0.02, 0.5),
    'lam_a0.8_l0.05': (300_000, 1, 'iso', 0.05, 0.8),
    'lam_a0.4_l0.02': (300_000, 1, 'iso', 0.02, 0.4),
    'lam_a0.25_l0.02': (300_000, 1, 'iso', 0.02, 0.25),
    'lam_a0.0_l0.01': (300_000, 1, 'iso', 0.01, 0.0),
}


def build():
    os.makedirs(os.path.dirname(BIN), exist_ok=True)
    subprocess.check_call(['cc', '-O2', '-o', BIN, os.path.join(D, 'sr3d.c'), '-lm'])


def run(name):
    T, seed, orient, lam0, alpha = RUNS[name]
    os.makedirs(os.path.join(D, 'data'), exist_ok=True)
    t0 = time.time()
    subprocess.check_call(['nice', '-n', '15', BIN, str(T), str(seed),
                           os.path.join(D, 'data', name), orient, str(lam0), str(alpha)])
    print(f'{name}: {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    build()
    for name in (sys.argv[1:] or list(RUNS)):
        run(name)
