#!/usr/bin/env python3
"""Acceptance metrics for the lossy-precision runs (docs/ray_compression.md,
step 6). Compares every variant of run_precision.sh against its baseline.

usage: precision_metrics.py OUTDIR [--json metrics.json] [--fig PREFIX]

Metrics, per problem and variant (worst over all hydro steps unless noted):
  xHII_vol, xHII_mass   max relative error of the volume-/mass-weighted
                        ionised fraction (steps with baseline x > 1e-3)
  T_mass                max relative error of the mass-weighted temperature
  T_ion_cell_rms        rms relative error of T over cells ionised (x > 0.5)
                        in the baseline, last step
  tau_es                relative error of the accumulated electron-scattering
                        optical depth (rt_stats tau), last step
  front_shift           mean I-front displacement in cells: cells whose
                        x > 0.5 state differs / number of ionised-neutral cell
                        faces in the baseline (max over steps)
  R_front               (uniform_1src) |R - R_base| in cells, R from the
                        ionised volume, max over steps
  photons               relative error of the cumulative photoionisations
                        per atom (rt_stats nion, summed) at the last step,
                        and of the photons in flight (ngam_rays_total)
Targets (docs/ROADMAP.md): x_HII within 1-2%, tau_es within 1%, I-front
within half a cell, T within a few percent.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

TARGETS = {"xHII_vol": 0.01, "xHII_mass": 0.01, "tau_es": 0.01, "front_shift": 0.5,
           "T_mass": 0.03}


def read_kv(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = [float(x) for x in v.split()] if len(v.split()) > 1 else float(v)
    return d


def load_run(d, N):
    steps = sorted(Path(d).glob("step_*"))
    out = []
    for s in steps:
        st = read_kv(s / "stats.txt")
        rs = read_kv(s / "rt_stats.txt")
        fHI = np.fromfile(s / "fHI.bin", dtype=np.float32).reshape(N, N, N).astype(np.float64)
        T = np.fromfile(s / "temp.bin", dtype=np.float32).reshape(N, N, N).astype(np.float64)
        rho = np.fromfile(s / "rhog.bin", dtype=np.float32).reshape(N, N, N).astype(np.float64)
        out.append(dict(z=st["redshift"], fHI=fHI, T=T, rho=rho, tau=rs["tau"],
                        nion=rs["nion"][0], ngam_rays=st["ngam_rays_total"]))
    return out


def front_faces(ion):
    n = 0
    for ax in range(3):
        n += np.count_nonzero(ion != np.roll(ion, 1, axis=ax))
    return n


def metrics(base, run):
    m = {"xHII_vol": 0.0, "xHII_mass": 0.0, "T_mass": 0.0, "front_shift": 0.0, "R_front": 0.0}
    hist = {"z": [], "xv_b": [], "xv": [], "xm_b": [], "xm": [], "Tm_b": [], "Tm": [],
            "shift": []}
    nsteps = min(len(base), len(run))
    cum_b = cum = 0.0
    for i in range(nsteps):
        b, r = base[i], run[i]
        xb, xr = 1 - b["fHI"], 1 - r["fHI"]
        xvb, xvr = xb.mean(), xr.mean()
        xmb = (xb * b["rho"]).sum() / b["rho"].sum()
        xmr = (xr * r["rho"]).sum() / r["rho"].sum()
        Tmb = (b["T"] * b["rho"]).sum() / b["rho"].sum()
        Tmr = (r["T"] * r["rho"]).sum() / r["rho"].sum()
        ib, ir_ = xb > 0.5, xr > 0.5
        faces = front_faces(ib)
        shift = np.count_nonzero(ib != ir_) / faces if faces else 0.0
        if xvb > 1e-3:
            m["xHII_vol"] = max(m["xHII_vol"], abs(xvr / xvb - 1))
            m["xHII_mass"] = max(m["xHII_mass"], abs(xmr / xmb - 1))
        m["T_mass"] = max(m["T_mass"], abs(Tmr / Tmb - 1))
        m["front_shift"] = max(m["front_shift"], shift)
        Rb = (3 * ib.sum() / (4 * np.pi)) ** (1 / 3)
        Rr = (3 * ir_.sum() / (4 * np.pi)) ** (1 / 3)
        m["R_front"] = max(m["R_front"], abs(Rr - Rb))
        cum_b += b["nion"]
        cum += r["nion"]
        for k, v in (("z", b["z"]), ("xv_b", xvb), ("xv", xvr), ("xm_b", xmb), ("xm", xmr),
                     ("Tm_b", Tmb), ("Tm", Tmr), ("shift", shift)):
            hist[k].append(float(v))
    b, r = base[nsteps - 1], run[nsteps - 1]
    ion = (1 - b["fHI"]) > 0.5
    if ion.any():
        rel = r["T"][ion] / b["T"][ion] - 1
        m["T_ion_cell_rms"] = float(np.sqrt(np.mean(rel ** 2)))
        xb = 1 - b["fHI"][ion]
        xr = 1 - r["fHI"][ion]
        m["xHII_ion_cell_rms"] = float(np.sqrt(np.mean((xr / xb - 1) ** 2)))
    m["tau_es"] = abs(r["tau"] / b["tau"] - 1) if b["tau"] else 0.0
    m["photons_ionisations"] = abs(cum / cum_b - 1) if cum_b else 0.0
    m["photons_in_rays"] = abs(r["ngam_rays"] / b["ngam_rays"] - 1) if b["ngam_rays"] else 0.0
    m["xHII_vol_final_base"] = hist["xv_b"][-1]
    m["nsteps"] = nsteps
    return {k: float(v) for k, v in m.items()}, hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--json")
    ap.add_argument("--N", type=int, default=32)
    args = ap.parse_args()
    res = {}
    hists = {}
    for pdir in sorted(Path(args.outdir).iterdir()):
        if not (pdir / "base" / "done").exists():
            continue
        base = load_run(pdir / "base", args.N)
        res[pdir.name] = {}
        hists[pdir.name] = {}
        for vdir in sorted(pdir.iterdir()):
            if vdir.name == "base" or not (vdir / "done").exists():
                continue
            run = load_run(vdir, args.N)
            m, h = metrics(base, run)
            res[pdir.name][vdir.name] = m
            hists[pdir.name][vdir.name] = h
        print("\n== %s (baseline x_HII,vol at end = %.3f, %d steps)" % (
            pdir.name, next(iter(res[pdir.name].values()))["xHII_vol_final_base"]
            if res[pdir.name] else float("nan"), len(base)))
        cols = ["xHII_vol", "xHII_mass", "T_mass", "T_ion_cell_rms", "xHII_ion_cell_rms", "tau_es",
                "front_shift", "R_front", "photons_ionisations", "photons_in_rays"]
        print("  %-10s" % "variant" + "".join("%12s" % c[:11] for c in cols))
        for v, m in res[pdir.name].items():
            print("  %-10s" % v + "".join("%12.2e" % m.get(c, float("nan")) for c in cols))
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"metrics": res, "hist": hists}, f, indent=1)


if __name__ == "__main__":
    main()
