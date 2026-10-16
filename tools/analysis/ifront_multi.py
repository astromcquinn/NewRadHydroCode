#!/usr/bin/env python3
"""The I-fronts of several sources in a radhydro run (docs/stage_f.md
section 3; tests/problems_full/make_ifront.py --sources 8src): per step,
x_HII,V, x_HII,M and the photons in flight (ngam_rays_total, per RT-cell
volume, converted to base cells for a finer uniform run) from the step's
stats.txt, and the ionised volume around each source from the fHI field: every
cell goes to its nearest source (periodic), V_s = sum over its cells of
(1 - x_HI), R_s = (3 V_s/4 pi)^(1/3), in base cells (the refined runs'
level-0 field is the composite; a finer uniform run is averaged down to the
base grid first, so that every run has the same partition: with the sources
at base-cell centres, many base cells are equidistant from two sources).

usage: ifront_multi.py --N n --sources sources_pos.txt RUN [RUN ...]
                       [--check REF BASE NEW] [--tol 0.25] [--xtol 0.01]
                       [--gtol 0.03] [--jump 0.15] [--json out.json]

--check REF BASE NEW: the NEW run (refined, base grid N) against REF (the
uniform run at the finest resolution), at every common step:
  * every source: |R_new - R_ref| <= |R_base - R_ref| + tol (base cells);
  * x_HII,V and x_HII,M: |new - ref|/ref <= |base - ref|/ref + xtol;
  * photons in flight: |new - ref|/ref <= |base - ref|/ref + gtol;
  * no jump: for every source, the change of R_new - R_ref from one step to
    the next is at most jump base cells (a front crossing a level boundary
    with the wrong photon normalisation overshot by 0.57).
Exit 1 if any fails.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np


def read_stats(path):
    out = {}
    with open(path) as f:
        for line in f:
            k, _, v = line.partition("=")
            try:
                out[k.strip()] = float(v)
            except ValueError:
                pass
    return out


def nearest(src, n):
    """index of the nearest source of every cell of an n^3 grid (box units)"""
    c = (np.arange(n) + 0.5) / n
    best = np.full((n, n, n), np.inf)
    idx = np.zeros((n, n, n), dtype=np.int32)
    for s, (x, y, z) in enumerate(src):
        d2 = 0
        for m, p in enumerate((x, y, z)):
            d = c - p
            d -= np.round(d)
            shape = [1, 1, 1]
            shape[2 - m] = n  # fields are x fastest: array index [k, j, i]
            d2 = d2 + (d * d).reshape(shape)
        better = d2 < best
        best = np.where(better, d2, best)
        idx = np.where(better, s, idx)
    return idx


def analyse(run, src, N):
    """{step: {xV, xM, ngam, R: [per source], Rtot}}"""
    out = {}
    cache = {}
    for d in sorted(glob.glob(os.path.join(run, "step_*"))):
        f = os.path.join(d, "fHI.bin")
        st = os.path.join(d, "stats.txt")
        if not (os.path.exists(f) and os.path.exists(st)):
            continue
        x = np.fromfile(f, dtype=np.float32)
        n = round(x.size ** (1 / 3))
        if n ** 3 != x.size:
            raise SystemExit(f"{f}: {x.size} values, not a cube")
        ion = 1.0 - x.astype(np.float64).reshape(n, n, n)
        if ion.sum() >= 0.999 * n ** 3:
            continue  # before the RT starts the field is unset (x_HI = 0 everywhere)
        cell = (N / n) ** 3  # a cell's volume in base cells
        if n > N:  # averaged down to the base grid, as the hierarchy's composite: one partition for every run
            f = n // N
            ion = ion.reshape(N, f, N, f, N, f).mean(axis=(1, 3, 5))
        if N not in cache:
            cache[N] = nearest(src, N)
        V = np.bincount(cache[N].ravel(), weights=ion.ravel(), minlength=len(src))
        s = read_stats(st)
        out[os.path.basename(d)[5:]] = {
            "xV": s.get("xHII_V"), "xM": s.get("xHII_M"), "ngam": s.get("ngam_rays_total") * cell,
            "R": [float((3 * v / (4 * np.pi)) ** (1 / 3)) for v in V],
            "Rtot": float((3 * V.sum() / (4 * np.pi)) ** (1 / 3)),
        }
    return out


def rel(a, b):
    return abs(a - b) / abs(b) if b else abs(a - b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--N", type=int, default=32, help="the base grid (R in its cells)")
    ap.add_argument("--sources", required=True, help="sources_pos.txt of make_ifront.py")
    ap.add_argument("--check", nargs=3, metavar=("REF", "BASE", "NEW"))
    ap.add_argument("--tol", type=float, default=0.25)
    ap.add_argument("--xtol", type=float, default=0.01)
    ap.add_argument("--gtol", type=float, default=0.03)
    ap.add_argument("--jump", type=float, default=0.15)
    ap.add_argument("--json")
    a = ap.parse_args()
    src = np.loadtxt(a.sources, ndmin=2)[:, :3]
    res = {}
    for r in a.runs + (list(a.check) if a.check else []):
        if r not in res:
            res[r] = analyse(r, src, a.N)
    for r in a.runs:
        print(f"{r}:")
        print("  step   x_HII,V   x_HII,M   ngam      R per source (base cells)")
        for k, v in res[r].items():
            print(f"  {k}  {v['xV']:.6f}  {v['xM']:.6f}  {v['ngam']:9.3f} " + " ".join(f"{q:6.3f}" for q in v["R"]))
    bad = 0
    summary = {}
    if a.check:
        ref, base, new = (res[r] for r in a.check)
        steps = sorted(set(ref) & set(base) & set(new))
        print(f"check: {a.check[2]} against {a.check[0]}, base {a.check[1]}")
        print("  step  |dxV| base/new     |dxM| base/new     |dngam| base/new   max|dR| base/new  max jump")
        prev = None
        worst = {"dR_base": 0.0, "dR_new": 0.0, "dxV_base": 0.0, "dxV_new": 0.0, "dxM_base": 0.0,
                 "dxM_new": 0.0, "dngam_base": 0.0, "dngam_new": 0.0, "jump_new": 0.0, "jump_base": 0.0}
        for k in steps:
            r, b, n = ref[k], base[k], new[k]
            msg = []
            ex = {}
            for q, tol in (("xV", a.xtol), ("xM", a.xtol), ("ngam", a.gtol)):
                eb, en = rel(b[q], r[q]), rel(n[q], r[q])
                ex[q] = (eb, en)
                worst["d" + q + "_base"] = max(worst["d" + q + "_base"], eb)
                worst["d" + q + "_new"] = max(worst["d" + q + "_new"], en)
                if en > eb + tol:
                    msg.append(f"{q} {en:.4f} > {eb:.4f} + {tol}")
            dRb = [bb - rr for bb, rr in zip(b["R"], r["R"])]
            dRn = [nn - rr for nn, rr in zip(n["R"], r["R"])]
            for s, (eb, en) in enumerate(zip(dRb, dRn)):
                if abs(en) > abs(eb) + a.tol:
                    msg.append(f"source {s + 1}: |dR| {abs(en):.3f} > {abs(eb):.3f} + {a.tol}")
            jump = jumpb = 0.0
            if prev is not None:
                jump = max(abs(x - y) for x, y in zip(dRn, prev[1]))
                jumpb = max(abs(x - y) for x, y in zip(dRb, prev[0]))
                if jump > a.jump:
                    msg.append(f"jump {jump:.3f} > {a.jump}")
            prev = (dRb, dRn)
            mb, mn = max(map(abs, dRb)), max(map(abs, dRn))
            worst["dR_base"] = max(worst["dR_base"], mb)
            worst["dR_new"] = max(worst["dR_new"], mn)
            worst["jump_new"] = max(worst["jump_new"], jump)
            worst["jump_base"] = max(worst["jump_base"], jumpb)
            bad += bool(msg)
            print(f"  {k}  {ex['xV'][0]:.5f}/{ex['xV'][1]:.5f}  {ex['xM'][0]:.5f}/{ex['xM'][1]:.5f}  "
                  f"{ex['ngam'][0]:.4f}/{ex['ngam'][1]:.4f}  {mb:.3f}/{mn:.3f}  {jump:.3f}"
                  + ("  <-- " + "; ".join(msg) if msg else ""))
        summary = {"steps": steps, "worst": worst, "failed_steps": bad}
        print("  worst over the steps: " + ", ".join(f"{k} {v:.4f}" for k, v in worst.items()))
        print("PASS" if not bad else f"FAIL: {bad} steps")
    if a.json:
        with open(a.json, "w") as f:
            json.dump({"N": a.N, "sources": src.tolist(), "runs": res, "check": a.check, "summary": summary}, f,
                      indent=1)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
