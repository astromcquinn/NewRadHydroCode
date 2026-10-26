#!/usr/bin/env python3
"""Front profiles and front speeds of a plane-parallel RT-domain run
(docs/plane_parallel.md): the transverse-averaged x_HI and T along the
propagation axis of each lit face, folded over the domains, the front
position per output, its speed, and whether it is limited by the reduced
speed of light.

usage: plane_front.py RUN[:label] [RUN:label ...] --outdir DIR [--axis x|y|z]
                      [--z ...] [--no-fig] [--json FILE]

RUN is a run directory (or the problem directory holding run/ and
params.txt). The fields fHI and temp of every step_NNNN (fields_every /
fields_at_data1) give the profiles; stats.txt gives the time, the RT
substeps and c_reduced of the step. The axis defaults to the first face of
plane_faces. Profiles are averaged over the two transverse axes and over
the domains along the axis (each domain's L cells from its source face,
reversed for a - face), so one profile per output of L cells runs from
the source cells (cell 0) to the far face.

The front position x_f (cells from the source face) is the ionised column
sum(1 - fHI) of the folded profile; the front speed is its change between
outputs over the time between them, in units of the ray speed
c_hat = L_cell / dt_substep of the later step (one cell per substep, times
rt_cells_per_substep): a value near 1 means the front is limited by the
reduced speed of light (DM20 Appendix A: c_sim = 0.1 c until x_HI,V < 1
per cent, then 0.01 c).

Writes DIR/plane_front.txt (per output: z, t [Myr], x_f, v_f/c_hat,
c_reduced, x_HII,V), DIR/plane_profiles.txt (per output: cell, <x_HI>,
<T>), DIR/plane_front.png.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

MYR = 3.156e13


def params_of(run):
    for p in (os.path.join(run, "params.txt"), os.path.join(run, "..", "params.txt")):
        if os.path.exists(p):
            d = {}
            for line in open(p):
                line = line.split("#")[0]
                if "=" in line:
                    k, v = line.split("=", 1)
                    d[k.strip()] = v.strip()
            return d
    sys.exit(f"{run}: no params.txt in it or above it")


def stats(fn):
    d = {}
    for line in open(fn):
        if " = " in line:
            k, v = line.split(" = ", 1)
            try:
                d[k.strip()] = float(v.split()[0])
            except ValueError:
                pass
    return d


def run_dir(spec):
    return os.path.join(spec, "run") if os.path.isdir(os.path.join(spec, "run")) else spec


def fold(f, axis, L, sign):
    """Transverse mean of a [k, j, i] field along the axis, folded over the
    domains from their source faces: a profile of L cells."""
    ax = {"x": 2, "y": 1, "z": 0}[axis]
    prof = f.mean(axis=tuple(a for a in range(3) if a != ax))
    N = prof.size
    ndom = N // L
    out = np.zeros(L)
    for d in range(ndom):
        seg = prof[d * L:(d + 1) * L]
        out += seg if sign > 0 else seg[::-1]
    return out / ndom


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--axis", default=None)
    ap.add_argument("--z", nargs="*", type=float, default=None)
    ap.add_argument("--no-fig", action="store_true")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rows, profs = [], []
    for spec in a.runs:
        path, _, label = spec.partition(":")
        run = run_dir(path)
        P = params_of(run)
        N = int(P["N"])
        ndom = int(P.get("plane_ndom", 1))
        L = N // ndom
        faces = P.get("plane_faces", "+x +y").replace(",", " ").split()
        face = ("+" + a.axis) if a.axis and not a.axis.startswith(("+", "-")) else (a.axis or faces[0])
        axis, sign = face[-1], (-1 if face.startswith("-") else 1)
        k_cells = int(P.get("rt_cells_per_substep", 1))
        label = label or os.path.basename(os.path.dirname(run.rstrip("/"))) or run
        steps = []
        for sd in sorted(glob.glob(os.path.join(run, "step_*"))):
            if os.path.exists(os.path.join(sd, "fHI.bin")) and os.path.exists(os.path.join(sd, "stats.txt")):
                steps.append((stats(os.path.join(sd, "stats.txt")), sd))
        if not steps:
            sys.exit(f"{run}: no step with fHI.bin")
        # the time at the end of every step (year of stats.txt) and the substep length
        prev = None
        for st, sd in steps:
            z = st["redshift"]
            if a.z is not None and not any(abs(z - zz) < 0.05 for zz in a.z) and (st, sd) != steps[-1]:
                pass
            fHI = np.fromfile(os.path.join(sd, "fHI.bin"), dtype=np.float32).reshape(N, N, N).astype(np.float64)
            T = np.fromfile(os.path.join(sd, "temp.bin"), dtype=np.float32).reshape(N, N, N).astype(np.float64)
            pf = fold(fHI, axis, L, sign)
            pT = fold(T, axis, L, sign)
            xf = float(np.sum(1 - pf))
            t = st["year"] / 1e6  # Myr
            nrt = st["nstep_rt"]
            dt_sub = st["dt_yr"] / max(nrt, 1) * 1e6 / 1e6  # yr
            chat_cells = k_cells / max(dt_sub, 1e-300)  # cells per yr at the ray speed
            vf = None
            if prev is not None and t > prev[0] and st.get("plane_emitted", 1) > 0:
                vf = (xf - prev[1]) / ((t - prev[0]) * 1e6) / chat_cells
            rows.append(dict(label=label, z=z, t=t, xf=xf, vf_over_chat=vf, c_reduced=st.get("c_reduced", float(P.get("c_reduced", 1))),
                             xHII_V=st.get("xHII_V", float("nan")), L=L))
            profs.append(dict(label=label, z=z, xHI=pf.tolist(), T=pT.tolist()))
            prev = (t, xf)
            print(f"{label}: z = {z:.4f} t = {t:9.3f} Myr x_f = {xf:7.3f} cells of {L}"
                  + (f" v_f/c_hat = {vf:.3f}" if vf is not None else "") + f" x_HII,V = {st.get('xHII_V', float('nan')):.4f}")
    with open(os.path.join(a.outdir, "plane_front.txt"), "w") as f:
        f.write("# label z t[Myr] x_f[cells] v_f/c_hat c_reduced x_HII,V L\n")
        for r in rows:
            f.write(f"{r['label']} {r['z']:.5f} {r['t']:.5f} {r['xf']:.4f} {r['vf_over_chat'] if r['vf_over_chat'] is not None else float('nan'):.4f} {r['c_reduced']:g} {r['xHII_V']:.5f} {r['L']}\n")
    with open(os.path.join(a.outdir, "plane_profiles.txt"), "w") as f:
        f.write("# label z cell <x_HI> <T>[K]\n")
        for p in profs:
            for c, (x, T) in enumerate(zip(p["xHI"], p["T"])):
                f.write(f"{p['label']} {p['z']:.5f} {c} {x:.6e} {T:.6e}\n")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(dict(front=rows, profiles=profs), f, indent=1)
    if a.no_fig:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: no figure")
        return
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    labels = sorted(set(p["label"] for p in profs))
    zs = sorted(set(round(p["z"], 3) for p in profs), reverse=True)
    cmap = plt.get_cmap("viridis")
    styles = ["-", "--", ":", "-."]
    sel = profs if a.z is None else [p for p in profs if any(abs(p["z"] - zz) < 0.05 for zz in a.z)]
    for p in sel:
        col = cmap(0.1 + 0.8 * zs.index(round(p["z"], 3)) / max(1, len(zs) - 1))
        ls = styles[labels.index(p["label"]) % len(styles)]
        c = np.arange(len(p["xHI"]))
        axs[0].semilogy(c, np.maximum(p["xHI"], 1e-8), ls=ls, color=col, label=f"{p['label']} z={p['z']:.2f}")
        axs[1].semilogy(c, p["T"], ls=ls, color=col)
    axs[0].set_xlabel("cells from the source face")
    axs[0].set_ylabel(r"$\langle x_{\rm HI}\rangle$")
    axs[1].set_xlabel("cells from the source face")
    axs[1].set_ylabel(r"$\langle T\rangle$ [K]")
    axs[0].legend(fontsize=6)
    for lab in labels:
        rr = [r for r in rows if r["label"] == lab]
        axs[2].plot([r["t"] for r in rr], [r["xf"] for r in rr], "o-", ms=3, label=lab)
    axs[2].set_xlabel("t [Myr]")
    axs[2].set_ylabel("front position [cells]")
    axs[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "plane_front.png"), dpi=150)
    print("wrote", os.path.join(a.outdir, "plane_front.png"))


if __name__ == "__main__":
    main()
