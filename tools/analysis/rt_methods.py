#!/usr/bin/env python3
"""Rays vs the moment RT (M1, OTVET) on four RT-only problems
(docs/rt_moments.md).

usage: rt_methods.py RT_DRIVER WORKDIR [--outdir docs/figures] [--threads n]

Problems (N = 32, z = 10, uniform neutral gas from tests/problems/uniform_1src
unless noted; legacy_rt = legacy_physics = 0, c_reduced = 0.1).

The rays rotate their HEALPix frame for every source and RT substep,
ray_rotate = 3, as in production (the B16 fix, docs/port_notes.md): with a
fixed frame (ray_rotate = 0, the original) the 192 directions of every
emission are the same, and the gaps between them leave grid-aligned holes in
the fronts at large radius. The fixed frame is run as well ("rays, fixed
frame") to show what the rotation changes; the moment methods have no rays
and ignore the key.

  ifront   one source (Luv = 8e28) in the centre, 15 steps of da = 2e-4:
           the I-front radius from the ionised volume, against the finite-c
           solution R = R_inf(t - R/c_hat) built from the c = 1 ray run;
  shadow   a 3x brighter source and a clump of radius 2 cells, 50x the mean
           density, centred 4 cells away along +x: the shadow behind it;
  pair     two sources 12 cells apart on the x axis: the I-fronts where the
           radiation fields of the two sources meet (where M1 fails: its
           closure has one flux direction per cell);
  box      tests/problems/lognormal_16src (16 sources in a lognormal density
           field), 10 steps: global ionised fraction and temperature, and the
           cell-by-cell agreement with the rays.
Writes rt_methods_ifront.png, rt_methods_maps.png and rt_methods.txt (the
numbers quoted in docs/rt_moments.md) to --outdir.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
N = 32
METHODS = ["rays", "rays0", "m1", "otvet"]
COLORS = {"rays": "#2a78d6", "rays0": "#8fb9e8", "m1": "#eb6834", "otvet": "#1baf7a"}
LABELS = {"rays": "rays", "rays0": "rays, fixed frame", "m1": "M1", "otvet": "OTVET"}
SHORT = {"rays": "rays", "rays0": "rays0", "m1": "M1", "otvet": "OTVET"}
# rt_method and the HEALPix frame rotation of each entry of METHODS
RT_METHOD = {"rays": "rays", "rays0": "rays", "m1": "m1", "otvet": "otvet"}
RAY_ROTATE = {"rays": 3, "rays0": 0, "m1": 0, "otvet": 0}
COMMON = "z_init = 10.0\nda_sim = 2e-4\nlegacy_rt = 0\nlegacy_physics = 0\n"


def write_problem(dst, gas, sources, nsteps):
    os.makedirs(dst, exist_ok=True)
    gas.astype(np.float32).tofile(os.path.join(dst, "gas.bin"))
    with open(os.path.join(dst, "params.txt"), "w") as f:
        f.write(f"N = {N}\nnsteps = {nsteps}\n" + COMMON)
    with open(os.path.join(dst, "sources.txt"), "w") as f:
        for (i, j, k, L) in sources:
            f.write(f"{i} {j} {k} {L:e} 1.0\n")


def uniform_gas():
    g = np.fromfile(os.path.join(ROOT, "tests/problems/uniform_1src/gas.bin"), dtype=np.float32)
    return g.reshape(N, N, N, 9).astype(np.float64)  # (k, j, i, component)


def run(driver, prob, out, method, extra=(), threads=0):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    cmd = [driver, prob, out, "--set", f"rt_method={RT_METHOD[method]}", "--set", "c_reduced=0.1",
           "--set", f"ray_rotate={RAY_ROTATE[method]}"]
    if threads:
        cmd += ["--set", f"threads={threads}"]
    for kv in extra:
        cmd += ["--set", kv]
    t0 = time.time()
    with open(out + ".stdout", "w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    return time.time() - t0


def field(out, step, name, ncomp=1):
    a = np.fromfile(os.path.join(out, f"step_{step:04d}", name + ".bin"), dtype=np.float32)
    return a.reshape((N, N, N, ncomp) if ncomp > 1 else (N, N, N)).astype(np.float64)


def stats(out, step, fn="stats.txt"):
    d = {}
    for line in open(os.path.join(out, f"step_{step:04d}", fn)):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = [float(x) for x in v.split()]
    return d


def radius(out, step):
    V = np.sum(1 - field(out, step, "fHI"))
    return (3 * V / (4 * np.pi)) ** (1 / 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("driver")
    ap.add_argument("work")
    ap.add_argument("--outdir", default=os.path.join(ROOT, "docs/figures"))
    ap.add_argument("--threads", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(a.outdir, exist_ok=True)
    rep = []
    sys.path.insert(0, os.path.join(ROOT, "tests"))
    import ifront_rsla_check as irc

    # ---- ifront ---------------------------------------------------------
    gas = uniform_gas()
    p_if = os.path.join(a.work, "p_ifront")
    write_problem(p_if, gas, [(16, 16, 16, 8e28)], 15)
    times = {}
    R, T, CH = {}, {}, {}
    for m in METHODS:
        out = os.path.join(a.work, f"ifront_{m}")
        times[("ifront", m)] = run(a.driver, p_if, out, m, threads=a.threads)
        t, sec = [], 0.0
        for n in range(1, 16):
            sec += stats(out, n, "rt_stats.txt")["dtsec_rt"][0] * stats(out, n)["nstep_rt"][0]
            t.append(sec)
        T[m] = np.array(t)
        R[m] = np.array([radius(out, n) for n in range(1, 16)])
        CH[m] = 1 / stats(out, 1, "rt_stats.txt")["dtsec_rt"][0]
    out1 = os.path.join(a.work, "ifront_rays_c1")
    run(a.driver, p_if, out1, "rays", ["c_reduced=1"], threads=a.threads)
    t1, sec = [], 0.0
    for n in range(1, 16):
        sec += stats(out1, n, "rt_stats.txt")["dtsec_rt"][0] * stats(out1, n)["nstep_rt"][0]
        t1.append(sec)
    t1 = np.array(t1)
    R1 = np.array([radius(out1, n) for n in range(1, 16)])
    c1 = 1 / stats(out1, 1, "rt_stats.txt")["dtsec_rt"][0]
    chat = np.full(15, CH["rays"])
    pred = irc.retarded(T["rays"], None, chat, np.concatenate([[0], t1 - R1 / c1]),
                        np.concatenate([[0], R1]))
    rep.append("I-front (uniform, c_reduced = 0.1): R [cells] and R/R_pred - 1")
    rep.append(" step  c_hat t  R_pred " + " ".join(f"{SHORT[m]:>14s}" for m in METHODS))
    for n in range(15):
        rep.append(f" {n + 1:4d} {chat[n] * T['rays'][n]:8.2f} {pred[n]:7.3f} " + " ".join(
            f"{R[m][n]:6.3f} ({R[m][n] / pred[n] - 1:+.3f})" for m in METHODS))
    # the roughness of the front: the spread over directions of the radius of
    # the partly ionised cells at the last step (the fixed HEALPix frame
    # puts every emission's 192 rays along the same directions)
    kk, jj, ii = np.meshgrid(np.arange(1, N + 1), np.arange(1, N + 1), np.arange(1, N + 1),
                             indexing="ij")
    r_c = np.sqrt((ii - 16) ** 2 + (jj - 16) ** 2 + (kk - 16) ** 2)
    rep.append("")
    rep.append("I-front at step 15: radius of the partly ionised cells (0.3 < x_HII < 0.7), "
               "mean +- rms over the cells")
    for m in METHODS:
        x = 1 - field(os.path.join(a.work, f"ifront_{m}"), 15, "fHI")
        shell = (x > 0.3) & (x < 0.7)
        rr = r_c[shell]
        rep.append(f"  {SHORT[m]:6s} r = {rr.mean():.2f} +- {rr.std():.2f} ({int(shell.sum())} cells)")
    # why the moment fronts run ahead early: the ionised volume beyond the
    # light cone r > c_hat t (+ half a cell diagonal), and the photons in
    # flight (stats.txt ngam_rays_total) relative to the rays'. Early on the
    # recombinations are negligible, so the ionised volume is the photons
    # emitted minus those in flight: fewer in flight means a larger volume.
    rep.append("")
    rep.append("I-front ahead of the finite-c solution: ionised volume beyond c_hat t + 0.87 "
               "[cells] (fraction of the ionised volume), photons in flight / the rays'")
    for n in (1, 2, 3, 5, 10, 15):
        ct = chat[n - 1] * T["rays"][n - 1]
        cols = []
        for m in METHODS:
            out = os.path.join(a.work, f"ifront_{m}")
            x = 1 - field(out, n, "fHI")
            fl = stats(out, n)["ngam_rays_total"][0] / stats(os.path.join(a.work, "ifront_rays"),
                                                            n)["ngam_rays_total"][0]
            b = x[r_c > ct + 0.87].sum()
            cols.append(f"{SHORT[m]} {b:5.1f} ({b / x.sum():4.1%}) {fl:.2f}")
        rep.append(f"  step {n:2d} c_hat t {ct:5.2f}: " + "   ".join(cols))

    # ---- shadow ---------------------------------------------------------
    gs = uniform_gas()
    kk, jj, ii = np.meshgrid(np.arange(1, N + 1), np.arange(1, N + 1), np.arange(1, N + 1),
                             indexing="ij")
    clump = (ii - 20) ** 2 + (jj - 16) ** 2 + (kk - 16) ** 2 <= 2.0 ** 2
    gs[clump] *= 50
    p_sh = os.path.join(a.work, "p_shadow")
    write_problem(p_sh, gs, [(16, 16, 16, 2.4e29)], 15)
    # ---- pair -----------------------------------------------------------
    p_pr = os.path.join(a.work, "p_pair")
    write_problem(p_pr, uniform_gas(), [(10, 16, 16, 4e28), (22, 16, 16, 4e28)], 15)
    maps = {}
    for prob, pdir in (("shadow", p_sh), ("pair", p_pr)):
        for m in METHODS:
            out = os.path.join(a.work, f"{prob}_{m}")
            times[(prob, m)] = run(a.driver, pdir, out, m, threads=a.threads)
            maps[(prob, m)] = field(out, 15, "fHI")[15]  # the k = 16 plane, (j, i)
    # shadow: neutral fraction on the axis behind the clump (i = 23..26)
    rep.append("")
    rep.append("Shadow: mean x_HII on the axis behind the clump (i = 23..26, j = k = 16), "
               "beside it at the same distance (i = 23..26, j = 16 +- 6), and in front (i = 17)")
    for m in METHODS:
        x = 1 - maps[("shadow", m)]
        rep.append(f"  {SHORT[m]:6s} behind {x[15, 22:26].mean():.3f}   beside "
                   f"{0.5 * (x[9, 22:26].mean() + x[21, 22:26].mean()):.3f}   in front {x[15, 16]:.3f}"
                   f"   ionised volume {np.sum(1 - field(os.path.join(a.work, f'shadow_{m}'), 15, 'fHI')):.1f} cells")
    # holes: neutral cells (x_HII < 0.5) well inside the front (r < R - 2,
    # R from the ionised volume) and outside the clump's shadow cone (more
    # than 45 degrees from +x); the fixed frame leaves them in the gaps
    # between its 192 directions
    cosx = (ii - 16) / np.maximum(r_c, 1e-9)
    rep.append("Shadow holes: neutral cells (x_HII < 0.5) at r < R - 2 outside the 45-degree "
               "cone behind the clump")
    for m in METHODS:
        x = 1 - field(os.path.join(a.work, f"shadow_{m}"), 15, "fHI")
        Rv = (3 * np.sum(x) / (4 * np.pi)) ** (1 / 3)
        inner = (r_c < Rv - 2) & (cosx < np.cos(np.pi / 4))
        rep.append(f"  {SHORT[m]:6s} {int(np.sum(inner & (x < 0.5))):4d} of {int(inner.sum())} "
                   f"cells (R = {Rv:.2f})")
    rep.append("")
    rep.append("Pair: extent of the ionised region (x_HII > 0.5) along y through the midpoint "
               "(i = 16) and through a source (i = 10), in cells")
    for m in METHODS:
        x = 1 - maps[("pair", m)]
        rep.append(f"  {SHORT[m]:6s} midpoint {np.sum(x[:, 15] > 0.5):3d}   source {np.sum(x[:, 9] > 0.5):3d}")

    # ---- box ------------------------------------------------------------
    p_box = os.path.join(ROOT, "tests/problems/lognormal_16src")
    rep.append("")
    rep.append("Box (lognormal_16src, legacy_rt = legacy_physics = 0): volume- and mass-weighted "
               "x_HII, fraction of partly ionised cells (0.05 < x < 0.95), median T [K], cells "
               "with |x_HII - x_HII(rays)| > 0.1 (rays: ray_rotate = 3)")
    xr = {}
    hist = {}
    for m in METHODS:
        out = os.path.join(a.work, f"box_{m}")
        times[("box", m)] = run(a.driver, p_box, out, m,
                                ["nsteps=60", "legacy_rt=0", "legacy_physics=0"],
                                threads=a.threads)
        hist[m] = []
        for n in (10, 20, 40, 60):
            fHI = field(out, n, "fHI")
            rho = field(out, n, "rhog")
            temp = field(out, n, "temp")
            x = 1 - fHI
            if m == "rays":
                xr[n] = x
            nbad = int(np.sum(np.abs(x - xr[n]) > 0.1))
            part = float(np.mean((x > 0.05) & (x < 0.95)))
            rep.append(f"  step {n:2d} {SHORT[m]:6s} {x.mean():.4f} {np.sum(x * rho) / rho.sum():.4f} "
                       f"{part:.3f} {np.median(temp):8.1f}   {nbad} of {N ** 3}")
            hist[m].append((stats(out, n)["redshift"][0], x.mean(), np.sum(x * rho) / rho.sum()))
        maps[("box", m)] = x[15]
    rep.append("")
    rep.append("Wall-clock seconds (rt_driver, whole run)")
    for prob in ("ifront", "shadow", "pair", "box"):
        rep.append(f"  {prob:7s} " + "  ".join(f"{SHORT[m]} {times[(prob, m)]:6.2f}" for m in METHODS))
    txt = "\n".join(rep)
    print(txt)
    with open(os.path.join(a.outdir, "rt_methods.txt"), "w") as f:
        f.write(txt + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    my = 3.156e13
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for m in METHODS:
        ls = "o--" if m == "rays0" else "o-"
        ax[0].plot(T[m] / my, R[m], ls, ms=4, lw=2, color=COLORS[m], label=LABELS[m])
        ax[1].plot(chat * T[m], R[m] / pred - 1, ls, ms=4, lw=2, color=COLORS[m], label=LABELS[m])
    ax[0].plot(T["rays"] / my, pred, "k--", lw=1.5, label="finite-c solution")
    ax[0].plot(T["rays"] / my, chat * T["rays"], ":", color="0.5", lw=1.5, label=r"$\hat c\,t$")
    ax[0].set_ylim(0, 1.1 * max(R["rays"].max(), pred.max()))
    ax[0].set_xlabel("t [Myr]")
    ax[0].set_ylabel("I-front radius [cells]")
    ax[0].legend(frameon=False, fontsize=9)
    ax[1].axhline(0, color="0.6", lw=1)
    ax[1].set_xlabel(r"light travel $\hat c\,t$ [cells]")
    ax[1].set_ylabel(r"$R/R_{\rm pred} - 1$")
    ax[1].set_ylim(-0.2, 0.7)
    for x in ax:
        x.spines[["top", "right"]].set_visible(False)
    fig.suptitle("I-front in a uniform medium, c_reduced = 0.1", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "rt_methods_ifront.png"), dpi=120)

    fig, ax = plt.subplots(3, len(METHODS), figsize=(12, 9.2))
    for r, prob in enumerate(("shadow", "pair", "box")):
        for c, m in enumerate(METHODS):
            img = maps[(prob, m)]
            if prob != "box":
                img = 1 - img
            im = ax[r, c].imshow(img, origin="lower", cmap="Blues", vmin=0, vmax=1,
                                 extent=(0.5, N + 0.5, 0.5, N + 0.5))
            ax[r, c].set_xticks([])
            ax[r, c].set_yticks([])
            if r == 0:
                ax[r, c].set_title(LABELS[m])
                circ = plt.Circle((20, 16), 2.0, fill=False, color="#eb6834", lw=1)
                ax[r, c].add_patch(circ)
            if c == 0:
                ax[r, c].set_ylabel({"shadow": "shadow (clump circled)", "pair": "two sources",
                                     "box": "16 sources, lognormal, step 60"}[prob])
    cb = fig.colorbar(im, ax=ax, shrink=0.6, label=r"$x_{\rm HII}$ (plane through the sources)")
    cb.outline.set_visible(False)
    fig.savefig(os.path.join(a.outdir, "rt_methods_maps.png"), dpi=110, bbox_inches="tight")


if __name__ == "__main__":
    main()
