"""First steps of three recursive-division rules of beyond/ Fig. 1:
SR (area-weighted cell, uniform point, uniform angle), uniform (cell chosen uniformly),
rect (area-weighted cell, axis-parallel cut). Uses beyond/divisions.run, same seed.
Each step: the picked cell is shaded, the point appears, the chord grows to the cell
boundary. Width ~ t^-1/2.  Usage: python make_gifs_beyond.py"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path
OUT = Path(__file__).parent
sys.path.insert(0, str(OUT.parent / "beyond"))
import divisions

T, PICK, NGROW, HOLD, SEED = 40, 3, 6, 2, 3
TITLE = {"SR": "SR: cell ∝ area, uniform angle",
         "uniform": "uniform: cell chosen uniformly",
         "rect": r"rect: cell ∝ area, angle in $\{0,\pi/2\}$"}

def frames_for():
    per = PICK + NGROW + HOLD
    return [(t, f) for t in range(T) for f in range(per)] + [(T, 0)] * 14

def width(t):
    return 1.0 + 7.0 / np.sqrt(t + 1)

def draw(ax, rule, g, t, f, fs=10):
    A, B, P, C = g["A"], g["B"], g["point"], g["cell"]
    ax.clear()
    ax.set_xlim(-.02, 1.02); ax.set_ylim(-.02, 1.02); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, lw=1.5, ec="0.3"))
    if t < T:
        ax.add_patch(plt.Polygon(C[t], closed=True, fc="#a8c8f0", ec="none", alpha=.6))
    for i in range(min(t, T)):
        ax.plot([A[i, 0], B[i, 0]], [A[i, 1], B[i, 1]], c="k", lw=width(i), solid_capstyle="butt")
    if t == T:
        ax.set_title(f"{TITLE[rule]}\n{T} ligaments", fontsize=fs); return
    ax.set_title(f"{TITLE[rule]}\nstep t = {t + 1}", fontsize=fs)
    if f < PICK: return                                   # cell picked
    p = P[t]; sa = np.linalg.norm(A[t] - p); sb = np.linalg.norm(B[t] - p)
    r = min((f - PICK) / NGROW, 1.0) * max(sa, sb)
    ua = p + (A[t] - p) * min(1, r / sa if sa else 1); ub = p + (B[t] - p) * min(1, r / sb if sb else 1)
    ax.plot([ua[0], ub[0]], [ua[1], ub[1]], c="#e8710a", lw=width(t), solid_capstyle="butt")
    ax.plot(*p, "o", c="#e8710a", ms=5, mec="k", mew=.6)

def animate(rule):
    g = divisions.run(T, SEED, rule)
    fig, ax = plt.subplots(figsize=(4.2, 4.6), dpi=110)
    FuncAnimation(fig, lambda fr: draw(ax, rule, g, *fr), frames=frames_for(), interval=90).save(
        OUT / f"div_{rule}.gif", writer=PillowWriter(fps=11))
    plt.close(fig)

def animate_all(rules=("SR", "uniform", "rect")):
    gs = [divisions.run(T, SEED, r) for r in rules]
    fig, axs = plt.subplots(1, len(rules), figsize=(4.0 * len(rules), 4.4), dpi=100)
    fig.subplots_adjust(left=.01, right=.99, bottom=.01, top=.88, wspace=.05)
    def upd(fr):
        for ax, r, g in zip(axs, rules, gs): draw(ax, r, g, *fr, fs=11)
    FuncAnimation(fig, upd, frames=frames_for(), interval=90).save(
        OUT / "div_all.gif", writer=PillowWriter(fps=11))
    plt.close(fig)

if __name__ == "__main__":
    for r in ("SR", "uniform", "rect"): animate(r)
    animate_all()
