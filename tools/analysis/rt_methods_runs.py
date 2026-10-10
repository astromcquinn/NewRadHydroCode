#!/usr/bin/env python3
"""Compare apps/radhydro production runs that differ only in rt_method
(rays, M1, OTVET; tools/hyak/rt_methods_runs.sh, docs/rt_moments.md section
6).

usage: rt_methods_runs.py RUN:label [RUN:label ...] --outdir DIR

RUN is a run directory (<problem>/run). The first run is the reference (the
rays). Writes to DIR:
  history.png     x_HII,V and x_HII,M against z (movie/frames.txt), and the
                  difference from the reference;
  maps_zXX.png    slices (the middle plane) of x_HII and T for each run at
                  the zdata1 outputs present in all runs (step_*/fHI.bin);
  compare.txt     z at 10/25/50/75/90/99% ionised (V and M), tau_es, and per
                  common output: <x_HII>_V, <x_HII>_M, the fraction of
                  partly ionised cells, median T in ionised gas, cells
                  differing from the reference by > 0.1 in x_HII, the
                  correlation of log T with the reference, and the
                  power spectrum of x_HII relative to the reference at
                  k = 2, 4, 8, 16 (in units of the fundamental mode).
                  Also the RT time and the peak memory from the logs.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reion_history as rh  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def outputs(run):
    """{redshift rounded to 1e-4: step directory} for steps with fHI.bin"""
    out = {}
    for d in sorted(glob.glob(os.path.join(run, "step_*"))):
        if not os.path.exists(os.path.join(d, "fHI.bin")):
            continue
        z = None
        for line in open(os.path.join(d, "stats.txt")):
            if line.startswith("redshift"):
                z = float(line.split("=")[1])
        # before the RT starts the RT grid fields are not filled (T = 0)
        t = np.fromfile(os.path.join(d, "temp.bin"), dtype=np.float32, count=4096)
        if z is not None and t.size and t.min() > 0:
            out[round(z, 4)] = d
    return out


def frames(run):
    """[(z, {field: path}, xHII_M)] of the movie frames with fHI and temp
    (written at the same movie_dz redshifts in every run, from inside the
    RT substeps, so they line up across runs even when the hydro steps,
    and hence the zdata1 output redshifts, differ)"""
    out = []
    fn = os.path.join(run, "movie", "frames.txt")
    if not os.path.exists(fn):
        return out
    for line in open(fn):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        xm = float(f[5])
        files = {k: os.path.join(run, "movie", "frame_%05d_%s.bin" % (int(f[0]), k))
                 for k in ("fHI", "temp")}
        if np.isfinite(xm) and all(os.path.exists(v) for v in files.values()):
            out.append((float(f[3]), files, xm))
    return out


def load(d, name, N):
    path = d[name] if isinstance(d, dict) else os.path.join(d, name + ".bin")
    return np.fromfile(path, dtype=np.float32).reshape(N, N, N)


def pk(x, kbins):
    N = x.shape[0]
    f = np.fft.rfftn(x - x.mean())
    p = np.abs(f) ** 2
    k = np.sqrt(np.add.outer(np.add.outer(np.fft.fftfreq(N) ** 2, np.fft.fftfreq(N) ** 2),
                             np.fft.rfftfreq(N) ** 2)) * N
    return [p[(k >= kb / 1.41) & (k < kb * 1.41)].mean() for kb in kbins]


def log_info(run):
    """RT seconds (sum of the RT timer phases) and peak RSS from <run>.log"""
    fn = run.rstrip("/") + ".log"
    if not os.path.exists(fn):
        return "", ""
    txt = open(fn, errors="replace").read()
    rt = re.findall(r"RT timers \[s\]: substeps (\d+)(.*)", txt)
    rss = re.findall(r"Maximum resident set size \(kbytes\): (\d+)", txt)
    s = ""
    if rt:
        n, rest = rt[-1]
        tot = sum(float(v) for v in re.findall(r"[a-z_0-9]+ ([0-9.]+)", rest))
        s = f"{int(n)} substeps, RT {tot / 3600:.2f} h (last invocation)"
    return s, (f"{int(rss[-1]) / 1e6:.1f} GB" if rss else "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    runs, labels = [], []
    for r in a.runs:
        p, _, lab = r.partition(":")
        runs.append(p)
        labels.append(lab or os.path.basename(os.path.dirname(p.rstrip("/"))))
    rep = []

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ---- history ----------------------------------------------------------
    hist = {}
    for r, lab in zip(runs, labels):
        try:
            hist[lab] = rh.read(r)
        except OSError:
            print(f"warning: no movie/frames.txt in {r}")
    rep.append("Reionization history (movie/frames.txt)")
    lv = [0.10, 0.25, 0.50, 0.75, 0.90, 0.99]
    rep.append("  run        " + " ".join(f"zV{int(100 * l):02d}" .rjust(7) for l in lv) +
               "  " + " ".join(f"zM{int(100 * l):02d}".rjust(7) for l in lv) + "   tau_es")
    for lab, (z, xv, xm) in hist.items():
        rep.append(f"  {lab:10s} " + " ".join(f"{rh.z_at(z, xv, l):7.3f}" for l in lv) + "  " +
                   " ".join(f"{rh.z_at(z, xm, l):7.3f}" for l in lv) +
                   f"   {rh.tau_es(z, xm):.4f}")
    if hist:
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ref = labels[0]
        for i, (lab, (z, xv, xm)) in enumerate(hist.items()):
            c = COLORS[i % len(COLORS)]
            ax[0].plot(z, xv, "-", color=c, lw=2, label=f"{lab} (V)")
            ax[0].plot(z, xm, "--", color=c, lw=1.5, label=f"{lab} (M)")
            if lab != ref and ref in hist:
                zr, xvr, _ = hist[ref]
                zz = z[(z <= zr.max()) & (z >= zr.min())]
                ax[1].plot(zz, np.interp(zz, z[::-1], xv[::-1]) - np.interp(zz, zr[::-1], xvr[::-1]),
                           "-", color=c, lw=2, label=f"{lab} − {ref}")
        ax[0].set_xlabel("z")
        ax[0].set_ylabel(r"$\langle x_{\rm HII}\rangle$")
        ax[0].invert_xaxis()
        ax[0].legend(frameon=False, fontsize=8)
        ax[1].axhline(0, color="0.6", lw=1)
        ax[1].set_xlabel("z")
        ax[1].set_ylabel(r"$\Delta\langle x_{\rm HII}\rangle_V$")
        ax[1].invert_xaxis()
        ax[1].legend(frameon=False, fontsize=8)
        for x in ax:
            x.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "history.png"), dpi=120)

    # ---- fields at common outputs -------------------------------------------
    # Compare at the movie frames nearest to the reference's zdata1 output
    # redshifts when every run has frames (the zdata1 outputs themselves
    # drift apart once the hydro steps differ); else at step outputs with
    # the same redshift.
    outs = [outputs(r) for r in runs]
    frs = [frames(r) for r in runs]
    sets = []  # (label z, [per-run (z, location, xM or None)])
    if all(frs):
        rep.append("\nFields: movie frames nearest to the reference's zdata1 redshifts; "
                   "<x>_M from frames.txt")
        for zt in sorted(outs[0], reverse=True):
            # the reference frame nearest to zt, then each run's frame
            # nearest to that frame's redshift
            z0 = min(frs[0], key=lambda e: abs(e[0] - zt))[0]
            pick = [min(f, key=lambda e: abs(e[0] - z0)) for f in frs]
            if abs(z0 - zt) <= 0.05 and max(abs(e[0] - z0) for e in pick) <= 0.05:
                sets.append((zt, pick))
    else:
        common = sorted(set.intersection(*[set(o) for o in outs]), reverse=True)
        sets = [(z, [(z, o[z], None) for o in outs]) for z in common]
    if not sets:
        rep.append("\nno common field outputs")
    kb = [2, 4, 8, 16]
    for z, pick in sets:
        d0 = pick[0][1]
        f0 = d0["fHI"] if isinstance(d0, dict) else os.path.join(d0, "fHI.bin")
        N = int(round((os.path.getsize(f0) / 4) ** (1 / 3)))
        rep.append(f"\nz = {z:g} (N = {N}): z_run <x>_V <x>_M partly(0.05..0.95) T_med(x>0.5) "
                   f"|dx|>0.1 corr(logT) P_x(k={','.join(map(str, kb))})/P_ref")
        x0 = t0 = p0 = None
        fig, ax = plt.subplots(2, len(runs), figsize=(3.2 * len(runs), 6.6), squeeze=False)
        for i, ((zr, d, xm), lab) in enumerate(zip(pick, labels)):
            x = 1 - load(d, "fHI", N)
            t = load(d, "temp", N)
            if xm is None:
                rho = load(d, "rhog", N)
                xm = np.sum(x * rho) / rho.sum()
            p = pk(x, kb)
            if i == 0:
                x0, t0, p0 = x, t, p
            ion = x > 0.5
            corr = np.corrcoef(np.log10(t).ravel(), np.log10(t0).ravel())[0, 1]
            rep.append(f"  {lab:10s} {zr:7.4f} {x.mean():.4f} {xm:.4f} "
                       f"{np.mean((x > 0.05) & (x < 0.95)):.3f} "
                       f"{(np.median(t[ion]) if ion.any() else float('nan')):8.0f} "
                       f"{np.mean(np.abs(x - x0) > 0.1):.3f} {corr:.3f} " +
                       " ".join(f"{pi / p0i:.3f}" for pi, p0i in zip(p, p0)))
            ax[0, i].imshow(x[N // 2], origin="lower", cmap="Blues", vmin=0, vmax=1)
            ax[1, i].imshow(np.log10(t[N // 2]), origin="lower", cmap="Oranges", vmin=0, vmax=4.7)
            ax[0, i].set_title(f"{lab} (z = {zr:.3f})")
            for r_ in (0, 1):
                ax[r_, i].set_xticks([])
                ax[r_, i].set_yticks([])
        ax[0, 0].set_ylabel(r"$x_{\rm HII}$")
        ax[1, 0].set_ylabel(r"$\log_{10} T$ [K] (0 – 4.7)")
        fig.suptitle(f"z = {z:g}, middle plane", fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, f"maps_z{z:g}.png"), dpi=100)
        plt.close(fig)

    rep.append("\nCost (from <run>.log)")
    for r, lab in zip(runs, labels):
        rt, rss = log_info(r)
        rep.append(f"  {lab:10s} {rt}  peak RSS {rss}")
    txt = "\n".join(rep)
    print(txt)
    with open(os.path.join(a.outdir, "compare.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
