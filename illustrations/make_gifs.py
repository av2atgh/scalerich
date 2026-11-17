"""First growth steps of the SR construction under the three orientation rules of
Both et al. (Nature 659, 345 (2026), Fig. 6a): random, aligned, polarized.
Each step: uniform nucleation point, angle from the rule, ligament grows both ways
until it hits the structure (T-junction). Width ~ lambda_t = lambda0 t^-1/2.
Usage: python make_gifs.py  -> sr_random.gif, sr_aligned.gif, sr_polarized.gif"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

T, NGROW, HOLD = 40, 6, 2           # ligaments, growth frames per step, still frames
OUT = Path(__file__).parent

def angle(rule, p, rng):
    if rule == "random":
        return rng.random() * np.pi
    if rule == "aligned":                       # prescribed angle set {0, pi/2}
        return rng.choice([0.0, np.pi / 2])
    if rule == "polarized":                     # circular field around the centre
        return np.arctan2(p[1] - .5, p[0] - .5) + np.pi / 2

def grow(rule, seed):
    rng = np.random.default_rng(seed); A, B, P = [], [], []
    for t in range(T):
        p = rng.random(2); th = angle(rule, p, rng); d = np.array([np.cos(th), np.sin(th)])
        with np.errstate(divide="ignore"):
            ss = np.concatenate([(np.array([0, 1]) - p[0]) / d[0], (np.array([0, 1]) - p[1]) / d[1]])
        smax = ss[ss > 0].min(); smin = ss[ss < 0].max()
        if t:
            a = np.array(A); e = np.array(B) - a; w = a - p
            den = d[0] * e[:, 1] - d[1] * e[:, 0]
            with np.errstate(divide="ignore", invalid="ignore"):
                s = (w[:, 0] * e[:, 1] - w[:, 1] * e[:, 0]) / den
                u = (w[:, 0] * d[1] - w[:, 1] * d[0]) / den
            ok = (u >= 0) & (u <= 1) & np.isfinite(s)
            if (m := ok & (s > 0) & (s < smax)).any(): smax = s[m].min()
            if (m := ok & (s < 0) & (s > smin)).any(): smin = s[m].max()
        A.append(p + smin * d); B.append(p + smax * d); P.append((p, d, smin, smax))
    return P

def animate(rule, seed, title):
    P = grow(rule, seed)
    fig, ax = plt.subplots(figsize=(4.2, 4.6), dpi=110)
    frames = [(t, f) for t in range(T) for f in range(NGROW + HOLD)] + [(T, 0)] * 14
    width = lambda t: 1.0 + 7.0 / np.sqrt(t + 1)
    def draw(fr):
        t, f = fr; ax.clear()
        ax.set_xlim(-.02, 1.02); ax.set_ylim(-.02, 1.02); ax.set_aspect("equal"); ax.axis("off")
        box = ax.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, lw=1.5, ec="0.3"))
        if rule == "polarized":
            for r in (.15, .3, .45, .6):
                c = ax.add_patch(plt.Circle((.5, .5), r, fill=False, ls=":", lw=.8, ec="0.7"))
                c.set_clip_path(box)
        elif rule == "aligned":
            for x in np.linspace(.1, .9, 5):
                ax.plot([x, x], [0, 1], ":", c="0.85", lw=.6); ax.plot([0, 1], [x, x], ":", c="0.85", lw=.6)
        for i in range(min(t, T)):
            p, d, a, b = P[i]; q = np.array([p + a * d, p + b * d])
            ax.plot(q[:, 0], q[:, 1], c="k", lw=width(i), solid_capstyle="butt")
        if t == T:
            ax.set_title(f"{title}\n{T} ligaments", fontsize=10); return
        p, d, a, b = P[t]; g = min(f / NGROW, 1.0); r = g * max(-a, b)
        q = np.array([p + max(a, -r) * d, p + min(b, r) * d])
        ax.plot(q[:, 0], q[:, 1], c="#e8710a", lw=width(t), solid_capstyle="butt")
        if f < NGROW + HOLD:
            ax.plot(*p, "o", c="#e8710a", ms=5, mec="k", mew=.6)
            if f == 0: ax.arrow(*(p - .06 * d), *(.12 * d), color="#e8710a", width=.004, head_width=0)
        ax.set_title(f"{title}\nstep t = {t + 1}", fontsize=10)
    anim = FuncAnimation(fig, draw, frames=frames, interval=90)
    anim.save(OUT / f"sr_{rule}.gif", writer=PillowWriter(fps=11)); plt.close(fig)

if __name__ == "__main__":
    animate("random", 3, "Random: uniform angle")
    animate("aligned", 3, r"Aligned: angle in $\{0,\pi/2\}$")
    animate("polarized", 3, "Polarized: circular orientation field")
