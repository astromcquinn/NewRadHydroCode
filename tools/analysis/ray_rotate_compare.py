#!/usr/bin/env python3
"""Figure for the B16 fix (ray_rotate, docs/port_notes.md): x_HII slices
through the source of tests/problems/uniform_1src (N = 32, source cell 16)
for rt_driver runs with different ray_rotate, plus the axis deficit of
axis_deficit.py against the step.

    ray_rotate_compare.py OUT.png LABEL=RUNDIR ... [--step 20]
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 32; s0 = 15  # source cell (0-based)
i = np.arange(N)
X, Y, Z = np.meshgrid(i, i, i, indexing="ij")
R = np.sqrt((X - s0) ** 2 + (Y - s0) ** 2 + (Z - s0) ** 2)
AXES = {"x": (Y == s0) & (Z == s0), "y": (X == s0) & (Z == s0), "z": (X == s0) & (Y == s0)}


def xhii(run, step):
    p = os.path.join(run, "step_%04d" % step, "fHI.bin")
    return 1 - np.fromfile(p, "<f4").reshape(N, N, N, order="F")


def axis_deficit(x):
    """most negative (axis cell - shell mean) over the x and y axes"""
    d = []
    for name in ("x", "y"):
        m = AXES[name]
        for rr in range(2, 14):
            on = m & (np.round(R) == rr)
            sh = (np.abs(R - rr) < 0.5) & ~m
            if on.any():
                d.append(x[on].mean() - x[sh].mean())
    return min(d)


def shell_rms(x):
    tot = n = 0
    for rr in range(1, 14):
        v = x[np.abs(R - rr) < 0.5]
        tot += ((v - v.mean()) ** 2).sum(); n += v.size
    return np.sqrt(tot / n)


ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("runs", nargs="+")
ap.add_argument("--step", type=int, default=20)
a = ap.parse_args()
runs = [r.split("=", 1) for r in a.runs]

fig = plt.figure(figsize=(3.3 * len(runs), 6.6))
gs = fig.add_gridspec(2, len(runs), height_ratios=[1, 0.8])
for k, (lab, run) in enumerate(runs):
    ax = fig.add_subplot(gs[0, k])
    x = xhii(run, a.step)
    im = ax.imshow(x[:, :, s0].T, origin="lower", cmap="cividis", vmin=0, vmax=1,
                   extent=(-0.5, N - 0.5, -0.5, N - 0.5))
    ax.set_xlim(6.5, 23.5); ax.set_ylim(6.5, 23.5)
    ax.set_title("%s\nstep %d, slice z = source" % (lab, a.step), fontsize=9)
    ax.set_xlabel("x")
    if k == 0: ax.set_ylabel("y")
fig.colorbar(im, ax=fig.axes, label="x_HII", shrink=0.45, anchor=(0, 1))
steps = sorted(int(d[5:]) for d in os.listdir(runs[0][1]) if d.startswith("step_"))
ax1 = fig.add_subplot(gs[1, : (len(runs) + 1) // 2])
ax2 = fig.add_subplot(gs[1, (len(runs) + 1) // 2:])
for lab, run in runs:
    xs = [xhii(run, s) for s in steps]
    ax1.plot(steps, [axis_deficit(x) for x in xs], label=lab)
    ax2.plot(steps, [shell_rms(x) for x in xs], label=lab)
ax1.set_xlabel("step"); ax1.set_ylabel("min(axis - shell mean) x_HII, x and y")
ax2.set_xlabel("step"); ax2.set_ylabel("rms of x_HII over shells")
ax1.legend(fontsize=8)
fig.savefig(a.out, dpi=110, bbox_inches="tight")
print("wrote", a.out)
