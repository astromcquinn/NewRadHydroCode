#!/usr/bin/env python3
"""Figures of the AMR refinement tests (docs/stage_f.md sections 2 and 3).

usage: amr_figures.py [--runs DIR] [--work DIR] [--out DIR] [--only sedov,slices,errors,fgam]

  sedov   the hydro suite with refinement (tools/hyak/hydro_suite_amr.slurm):
          RUNS/hydro_suite_amr2/quick_unsplit{,_amr1,_amr2}, the profiles and
          errors from docs/data/stage_f/sedov_*.json -> amr_sedov.png
  slices  the 8-source I-fronts (tools/hyak/rt_amr_ifront8.slurm, OUT =
          RUNS/rt_amr_ifront8c): x_HII slices and differences -> amr_ifront8_slices.png
  errors  their error histories, docs/data/stage_f/ifront8_L{1,2}.json
          -> amr_ifront8_errors.png
  fine    the 8-source runs with fields_levels = 1 (WORK/fields_levels/{L1,L2},
          fields_levels_index = 16): x_HII at the finest resolution from
          levels_slice.bin, the levels outlined, beside the uniform 64^3 and
          128^3 runs at the same planes -> amr_ifront8_fine.png
  fgam    the photons in flight before and after the covered-cell fgam fix:
          WORK/f3_8src/static/{p_amr (before), p_amr2 (after), p_k2 (uniform
          64^3, rt_cells_per_substep = 2)} with fields = fHI,ngam -> amr_fgam_fix.png

The refined runs' field files are the level-0 composite (every refined
cell's children averaged down), so those panels show the base resolution;
the uniform runs are shown at their own resolution, and `fine` uses the
plane fields_levels writes at the finest resolution. Defaults: the klone paths,
out = docs/figures.
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
RUNS = "/gscratch/astro/mcquinn/radhydro/runs"
WORK = "/gscratch/astro/mcquinn/radhydro/work"
OUT = os.path.join(REPO, "docs", "figures")
sys.path.insert(0, os.path.join(REPO, "tests/hydro_suite"))
import hs_common as hc  # noqa: E402

# reference palette (dataviz skill, categorical slots 1-3), ink, surface
C_BASE, C_L1, C_L2 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list("seq", ["#f4f8fd", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("div", ["#c2481a", "#eb6834", "#f0efec", "#2a78d6", "#184f95"])

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 10.5, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "lines.linewidth": 2, "legend.frameon": False,
})


def img_axes(ax):
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color(GRID)
    ax.set_xticks([])
    ax.set_yticks([])


def box(ax, lo_i, hi_i, lo_j, hi_j, color=INK, ls="-", label=None):
    ax.add_patch(Rectangle((lo_i - 0.5, lo_j - 0.5), hi_i - lo_i + 1, hi_j - lo_j + 1, fill=False,
                           ec=color, lw=1.4, ls=ls, label=label))


# ---------------------------------------------------------------- 1. Sedov
def fig_sedov():
    H = os.path.join(RUNS, "hydro_suite_amr2")  # after the F2 hierarchy fixes
    runs = [("uniform 32³", "quick_unsplit", 32, C_BASE, "-"),
            ("32³ + 1 level", "quick_unsplit_amr1", 32, C_L1, "--"),
            ("32³ + 2 levels", "quick_unsplit_amr2", 32, C_L2, "-."),
            ("uniform 64³", "quick_unsplit", 64, INK2, ":")]
    fig = plt.figure(figsize=(12, 6.4))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.9], hspace=0.32, wspace=0.08)
    vmax = 0
    slices = []
    for name, tag, N, col, ls in runs:
        run = os.path.join(H, tag, f"sedov_N{N}")
        meta = json.load(open(os.path.join(H, tag, "problems", f"sedov_N{N}", "meta.json")))
        d = hc.outputs(run)[-1]
        rho, vel, p, info = hc.read_state(d, N)
        ctr = np.mod(np.array(meta["centre"]) + np.array(meta["bulk"]) * info["t"], N)
        k = int(ctr[2])
        sl = rho[:, :, k].T  # [j, i]
        # recentre on the blast
        sl = np.roll(np.roll(sl, N // 2 - int(ctr[0]), axis=1), N // 2 - int(ctr[1]), axis=0)
        slices.append((name, sl / meta["rho0"], N))
        vmax = max(vmax, sl.max() / meta["rho0"])
    for n, (name, sl, N) in enumerate(slices):
        ax = fig.add_subplot(gs[0, n])
        im = ax.imshow(sl, origin="lower", cmap=SEQ, vmin=0, vmax=vmax, extent=(0, 32, 0, 32),
                       interpolation="nearest")
        img_axes(ax)
        ax.set_title(name)
    cax = fig.add_axes([0.915, 0.555, 0.008, 0.33])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("ρ / ρ₀", color=INK)
    cb.outline.set_visible(False)
    # profiles from the suite's analysis (binned in r, exact solution alongside)
    axp = fig.add_subplot(gs[1, :2])
    axr = fig.add_subplot(gs[1, 2:])
    exact_done = False
    for name, tag, N, col, ls in runs:
        js = os.path.join(REPO, "docs/data/stage_f",
                          "sedov_quick_unsplit.json" if tag == "quick_unsplit" else f"sedov_{tag}.json")
        o = json.load(open(js))["cases"][f"sedov_N{N}"]["outputs"][-1]
        pr = o["profile"]
        r = np.array(pr["r"]) / (N / 32)  # in 32³ cells
        axp.plot(r, np.array(pr["rho"]) / pr["rho2"], color=col, ls=ls, label=name)
        if not exact_done:
            axp.plot(r, np.array(pr["rho_ex"]) / pr["rho2"], color=INK, lw=1, label="exact (Sedov–Taylor)")
            exact_done = True
        Rd = np.array(o["R_dirs"]) if "R_dirs" in o else None
    axp.set_xlim(0, 15)
    axp.set_xlabel("r  [32³ cells]")
    axp.set_ylabel("ρ / ρ_shock (exact)")
    axp.set_title("Radial density profile at t = 10")
    axp.legend(loc="upper left", fontsize=9)
    # L1 errors and shock-radius errors (table of docs/stage_f.md section 2)
    labels = ["ρ L1 error", "p L1 error", "|R/R_exact − 1|"]
    vals = {}
    for name, tag, N, col, ls in runs:
        js = os.path.join(REPO, "docs/data/stage_f",
                          "sedov_quick_unsplit.json" if tag == "quick_unsplit" else f"sedov_{tag}.json")
        o = json.load(open(js))["cases"][f"sedov_N{N}"]["outputs"][-1]
        vals[name] = [o["l1"]["rho"], o["l1"]["p"], abs(o["R_err_rel"])]
    x = np.arange(len(labels))
    w = 0.19
    for n, (name, tag, N, col, ls) in enumerate(runs):
        b = axr.bar(x + (n - 1.5) * w, vals[name], w * 0.9, color=col, label=name)
        for rect, v in zip(b, vals[name]):
            axr.text(rect.get_x() + rect.get_width() / 2, v + 0.001, f"{v:.3f}", ha="center", va="bottom",
                     fontsize=7.5, color=INK2)
    axr.set_xticks(x)
    axr.set_xticklabels(labels)
    axr.grid(axis="x", visible=False)
    axr.set_title("Errors against the exact solution (lower is better)")
    axr.legend(fontsize=8.5, loc="upper right")
    fig.suptitle("Sedov blast with refinement (hydro, F2): the whole blast refined, t = 10",
                 x=0.06, ha="left", fontsize=12)
    fig.savefig(os.path.join(OUT, "amr_sedov.png"), dpi=140, bbox_inches="tight")


# ------------------------------------------------------- 2. 8-source slices
def load_fhi(run, step):
    x = np.fromfile(os.path.join(run, f"step_{step}", "fHI.bin"), dtype=np.float32).astype(float)
    n = round(x.size ** (1 / 3))
    return x.reshape(n, n, n)  # [k, j, i]


def down(a, N=32):
    f = a.shape[0] // N
    return a.reshape(N, f, N, f, N, f).mean(axis=(1, 3, 5)) if f > 1 else a


def fig_ifront_slices():
    R = os.path.join(RUNS, "rt_amr_ifront8c")
    st = "0010"
    src = np.loadtxt(os.path.join(R, "prob32", "sources_pos.txt"))[:, :3] * 32
    k = 16  # base-cell slice through source 6 (z = 16.5) and the refined boxes
    u32 = 1 - load_fhi(os.path.join(R, "uniform32"), st)
    a1 = 1 - load_fhi(os.path.join(R, "amr32_L1"), st)
    a2 = 1 - load_fhi(os.path.join(R, "amr32_L2"), st)
    u64 = 1 - load_fhi(os.path.join(R, "uniform64"), st)
    u128 = 1 - load_fhi(os.path.join(R, "uniform128"), st)
    fig, axs = plt.subplots(2, 4, figsize=(12.5, 6.6), gridspec_kw={"hspace": 0.18, "wspace": 0.06})
    top = [("uniform 32³", u32[k]), ("32³ + 1 level (composite)", a1[k]),
           ("32³ + 2 levels (composite)", a2[k]), ("uniform 128³ (native)", u128[4 * k + 2])]
    for ax, (name, sl) in zip(axs[0], top):
        ext = (0, 32, 0, 32)
        im = ax.imshow(sl, origin="lower", cmap=SEQ, vmin=0, vmax=1, extent=ext, interpolation="nearest")
        img_axes(ax)
        ax.set_title(name)
    for ax in axs[0, 1:3]:
        box(ax, 0, 15, 4, 19)
    box(axs[0, 2], 2, 13, 6, 17, ls="--")
    for ax in axs[0]:
        for s in src:
            if abs(s[2] - (k + 0.5)) < 6:
                ax.plot(s[0], s[1], marker="*", ms=9, color=C_L1, mec=SURF, mew=0.8)
    cb = fig.colorbar(im, ax=axs[0], fraction=0.015, pad=0.01)
    cb.set_label("x_HII", color=INK)
    cb.outline.set_visible(False)
    # differences against the finest uniform run, averaged to 32³
    r64, r128 = down(u64)[k], down(u128)[k]
    bot = [("uniform 32³ − 64³", u32[k] - r64), ("1 level − 64³", a1[k] - r64),
           ("uniform 32³ − 128³", u32[k] - r128), ("2 levels − 128³", a2[k] - r128)]
    # inside the level-1 box (base cells i 0..15, j 4..19)
    lim = 0.5
    for ax, (name, d) in zip(axs[1], bot):
        im2 = ax.imshow(d, origin="lower", cmap=DIV, vmin=-lim, vmax=lim, extent=(0, 32, 0, 32),
                        interpolation="nearest")
        img_axes(ax)
        inside = np.abs(d[4:20, 0:16]).mean()
        ax.set_title(f"{name}\nmean |Δ| {np.abs(d).mean():.3f};  in box {inside:.3f}", fontsize=9.5)
    box(axs[1, 1], 0, 15, 4, 19)
    box(axs[1, 3], 0, 15, 4, 19)
    box(axs[1, 3], 2, 13, 6, 17, ls="--")
    cb2 = fig.colorbar(im2, ax=axs[1], fraction=0.015, pad=0.01)
    cb2.set_label("Δ x_HII", color=INK)
    cb2.outline.set_visible(False)
    fig.suptitle("8 sources in uniform gas, step 10 (z ≈ 11.6): slice z = 16–17 base cells.  "
                 "Solid box: level 1; dashed: level 2; ★ sources near the slice",
                 x=0.06, ha="left", fontsize=11.5)
    fig.savefig(os.path.join(OUT, "amr_ifront8_slices.png"), dpi=140, bbox_inches="tight")


# ----------------------------------------------------- 3. error histories
def fig_ifront_errors():
    rows = [("1 level vs uniform 64³", "ifront8_L1.json", "amr32_L1", "uniform64", C_L1),
            ("2 levels vs uniform 128³", "ifront8_L2.json", "amr32_L2", "uniform128", C_L2)]
    fig, axs = plt.subplots(2, 3, figsize=(12.5, 6.4), sharex=True,
                            gridspec_kw={"hspace": 0.3, "wspace": 0.28})
    for r, (title, js, new, ref, col) in enumerate(rows):
        d = json.load(open(os.path.join(REPO, "docs/data/stage_f", js)))["runs"]
        runs = {os.path.basename(k): v for k, v in d.items()}
        steps = sorted(runs[ref])
        s = np.array([int(x) for x in steps])

        def err(run, f):
            return np.array([f(runs[run][t], runs[ref][t]) for t in steps])
        panels = [
            ("x_HII,V relative error", lambda a, b: abs(a["xV"] - b["xV"]) / b["xV"]),
            ("photons in flight, relative error", lambda a, b: abs(a["ngam"] - b["ngam"]) / b["ngam"]),
            ("max per-source |ΔR|  [base cells]",
             lambda a, b: max(abs(x - y) for x, y in zip(a["R"], b["R"]))),
        ]
        for c, (lab, f) in enumerate(panels):
            ax = axs[r, c]
            eb, en = err("uniform32", f), err(new, f)
            ax.plot(s, eb, color=C_BASE, marker="o", ms=5, label="uniform 32³")
            ax.plot(s, en, color=col, marker="s", ms=5, ls="--", label=title.split(" vs")[0])
            ax.set_ylim(0, max(eb.max(), en.max()) * 1.18)
            ax.set_title(lab if r == 0 else "")
            if c == 0:
                ax.set_ylabel(title, fontsize=10)
            if r == 1:
                ax.set_xlabel("step  (3 → 10: z ≈ 22 → 11.6)")
            ax.legend(fontsize=8.5, loc="upper right" if c != 2 else "lower right")
    fig.suptitle("8-source I-fronts: errors against the finest uniform run (lower is better; "
                 "no jumps when fronts cross level boundaries)", x=0.06, ha="left", fontsize=11.5)
    fig.savefig(os.path.join(OUT, "amr_ifront8_errors.png"), dpi=140, bbox_inches="tight")


# ------------------------------------- 3b. fields_levels: the fine slices
def read_levels_slice(d):
    t = {}
    for line in open(os.path.join(d, "levels_slice.txt")):
        if " = " in line and not line.startswith("#"):
            k, v = line.split(" = ", 1)
            t[k.strip()] = v.strip()
    n = int(t["n"])
    a = np.fromfile(os.path.join(d, "levels_slice.bin"), dtype=np.float32).reshape(n, n, int(t["ncomp"]))
    return a, t  # a[b, a, comp]: the first in-plane axis fastest


def level_outline(ax, lev, level, **kw):
    """the edges between pixels of level >= `level` and the rest, in base cells"""
    n = lev.shape[0]
    m = lev >= level
    f = 32.0 / n
    segs = []
    for j in range(n):
        for i in range(n):
            if not m[j, i]:
                continue
            if i == 0 or not m[j, i - 1]:
                segs.append(((i * f, j * f), (i * f, (j + 1) * f)))
            if i == n - 1 or not m[j, i + 1]:
                segs.append((((i + 1) * f, j * f), ((i + 1) * f, (j + 1) * f)))
            if j == 0 or not m[j - 1, i]:
                segs.append(((i * f, j * f), ((i + 1) * f, j * f)))
            if j == n - 1 or not m[j + 1, i]:
                segs.append(((i * f, (j + 1) * f), ((i + 1) * f, (j + 1) * f)))
    from matplotlib.collections import LineCollection
    ax.add_collection(LineCollection(segs, **kw))


def fig_fine():
    W = os.path.join(WORK, "fields_levels")
    R = os.path.join(RUNS, "rt_amr_ifront8c")
    st = "0010"
    a1, t1 = read_levels_slice(os.path.join(W, "L1", f"step_{st}"))
    a2, t2 = read_levels_slice(os.path.join(W, "L2", f"step_{st}"))
    u64 = 1 - load_fhi(os.path.join(R, "uniform64"), st)[int(t1["plane"])]
    u128 = 1 - load_fhi(os.path.join(R, "uniform128"), st)[int(t2["plane"])]
    fig, axs = plt.subplots(1, 4, figsize=(15, 4.3), gridspec_kw={"wspace": 0.05})
    panels = [("32³ + 1 level (64³ in the box)", 1 - a1[:, :, 1], a1[:, :, 5], 1),
              ("uniform 64³", u64, None, 0),
              ("32³ + 2 levels (128³ in the inner box)", 1 - a2[:, :, 1], a2[:, :, 5], 2),
              ("uniform 128³", u128, None, 0)]
    for ax, (name, img, lev, L) in zip(axs, panels):
        im = ax.imshow(img, origin="lower", cmap=SEQ, vmin=0, vmax=1, extent=(0, 32, 0, 32), interpolation="nearest")
        img_axes(ax)
        ax.set_title(name)
        if lev is not None:
            for l in range(1, L + 1):
                level_outline(ax, lev, l, colors=C_L1 if l == 1 else C_L2, linewidths=1.6,
                              linestyles="solid" if l == 1 else "dashed")
    cb = fig.colorbar(im, ax=axs, fraction=0.012, pad=0.01)
    cb.set_label("x_HII", color=INK)
    cb.outline.set_visible(False)
    fig.suptitle("fields_levels = 1: the slice at the finest resolution (step 10, z ≈ 11.6, the plane through "
                 "base cell 16).  Outlines: level 1 (orange), level 2 (green)", x=0.06, ha="left", fontsize=11.5)
    fig.savefig(os.path.join(OUT, "amr_ifront8_fine.png"), dpi=140, bbox_inches="tight")


# ----------------------------------------------- 4. the covered-fgam fix
def load_ngam(run, st):
    x = np.fromfile(os.path.join(run, f"step_{st}", "ngam.bin"), dtype=np.float32).astype(float)
    n = round((x.size // 5) ** (1 / 3))
    return down(x.reshape(n, n, n, 5).sum(axis=3))


def fig_fgam_fix():
    W = os.path.join(WORK, "f3_8src", "static")
    st = "0010"
    before, after, ref = (load_ngam(os.path.join(W, r), st) for r in ("p_amr", "p_amr2", "p_k2"))
    k = 12
    fig = plt.figure(figsize=(14, 4.8))
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 0.05, 0.45, 1.25], wspace=0.08)
    for n, (name, a) in enumerate((("before the fix", before), ("after the fix", after))):
        ax = fig.add_subplot(gs[0, n])
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(ref[k] > 1e-3 * ref.max(), a[k] / ref[k], np.nan)
        im = ax.imshow(ratio - 1, origin="lower", cmap=DIV, vmin=-0.3, vmax=0.3, extent=(0, 32, 0, 32),
                       interpolation="nearest")
        img_axes(ax)
        box(ax, 0, 15, 4, 19)
        ax.set_title(name)
    cax = fig.add_subplot(gs[0, 2])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("refined / matched 64³ − 1", color=INK)
    cb.outline.set_visible(False)
    ax = fig.add_subplot(gs[0, 4])
    steps = [f"{i:04d}" for i in range(3, 11)]

    def tot(run):
        out = []
        for t in steps:
            for line in open(os.path.join(W, run, f"step_{t}", "stats.txt")):
                if line.startswith("ngam_rays_total"):
                    out.append(float(line.split("=")[1]))
        return np.array(out)
    rk = tot("p_k2") / 8
    s = np.arange(3, 11)
    ax.plot(s, 100 * (tot("p_amr") / rk - 1), color=C_BASE, marker="o", ms=5, label="before the fix")
    ax.plot(s, 100 * (tot("p_amr2") / rk - 1), color=C_L1, marker="s", ms=5, ls="--", label="after the fix")
    ax.axhline(0, color=INK2, lw=1)
    ax.set_xlabel("step")
    ax.set_ylabel("difference  [%]")
    ax.set_title("Total photons in flight vs matched 64³")
    ax.legend(fontsize=9, loc="lower left")
    fig.suptitle("Photons in flight, refined run against uniform 64³ at k = 2 (step 10, slice z = 12–13; "
                 "box = level 1).  Rays entering the level used a stale covered-cell fgam",
                 x=0.06, ha="left", fontsize=11)
    fig.savefig(os.path.join(OUT, "amr_fgam_fix.png"), dpi=140, bbox_inches="tight")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=RUNS)
    ap.add_argument("--work", default=WORK)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--only", default="sedov,slices,errors,fine,fgam")
    a = ap.parse_args()
    RUNS, WORK, OUT = a.runs, a.work, a.out
    os.makedirs(OUT, exist_ok=True)
    figs = {"sedov": fig_sedov, "slices": fig_ifront_slices, "errors": fig_ifront_errors, "fine": fig_fine,
            "fgam": fig_fgam_fix}
    for name in a.only.split(","):
        figs[name]()
        print("wrote", os.path.join(OUT, "amr_" + {"sedov": "sedov", "slices": "ifront8_slices",
                                                     "errors": "ifront8_errors", "fine": "ifront8_fine", "fgam": "fgam_fix"}[name] + ".png"))
