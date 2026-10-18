#!/usr/bin/env python3
"""The PM backend's run against the serial code after a restart
(tools/hyak/amrex_pm_scaling.slurm, docs/stage_e_pm.md section 4.1): per
step the gas fields (rms and largest difference over the field's rms and
maximum), the DM and gas statistics of stats.txt, and for the zdata1 output
the particles (periodic position difference in DM-grid cells, velocity
difference over the rms velocity) and P(k). Round-off is expected: the
fixed-point CIC and the per-line FFT against the serial float CIC and 3-D
FFTW plan (section 1.3).

usage: amrex_pm_vs_serial.py SERIAL BACKEND [--ndm NDM]
"""
import argparse
import glob
import os

import numpy as np

STATS = ("amax_dm", "arms_dm", "vmax_dm", "vrms_dm", "cmax_gas", "vmax_gas", "vrms_gas", "Nhalo")


def read_kv(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                d[k.strip()] = float(v.split()[0])
            except ValueError:
                pass
    return d


def field_diff(a, b, ncomp):
    """rms(b - a)/rms(a) and max|b - a|/max|a| per component (memory-mapped, chunked)."""
    x = np.memmap(a, dtype="<f4", mode="r").reshape(-1, ncomp)
    y = np.memmap(b, dtype="<f4", mode="r").reshape(-1, ncomp)
    s2 = np.zeros(ncomp)
    d2 = np.zeros(ncomp)
    dm = np.zeros(ncomp)
    am = np.zeros(ncomp)
    step = 1 << 24
    for i in range(0, x.shape[0], step):
        xa = x[i:i + step].astype(np.float64)
        ya = y[i:i + step].astype(np.float64)
        s2 += (xa ** 2).sum(0)
        d2 += ((ya - xa) ** 2).sum(0)
        dm = np.maximum(dm, np.abs(ya - xa).max(0))
        am = np.maximum(am, np.abs(xa).max(0))
    return np.sqrt(d2 / s2), dm / am


def dm_diff(a, b, ndm):
    x = np.memmap(a, dtype="<f4", mode="r").reshape(-1, 6)
    y = np.memmap(b, dtype="<f4", mode="r").reshape(-1, 6)
    n = x.shape[0]
    dx2 = dv2 = v2 = 0.0
    dxm = dvm = 0.0
    step = 1 << 24
    for i in range(0, n, step):
        xa = x[i:i + step].astype(np.float64)
        ya = y[i:i + step].astype(np.float64)
        d = ya[:, :3] - xa[:, :3]
        d -= ndm * np.round(d / ndm)  # periodic
        dv = ya[:, 3:] - xa[:, 3:]
        dx2 += (d ** 2).sum()
        dv2 += (dv ** 2).sum()
        v2 += (xa[:, 3:] ** 2).sum()
        dxm = max(dxm, np.abs(d).max())
        dvm = max(dvm, np.abs(dv).max())
    vrms = np.sqrt(v2 / n)
    return np.sqrt(dx2 / n), dxm, np.sqrt(dv2 / n) / vrms, dvm / vrms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("serial")
    ap.add_argument("backend")
    ap.add_argument("--ndm", type=float, default=1024)
    a = ap.parse_args()
    steps = sorted(os.path.basename(s) for s in glob.glob(os.path.join(a.serial, "step_*")))
    print("| step | field | rms diff / rms | max diff / max |")
    print("|---|---|---|---|")
    for st in steps:
        for f, nc in (("u_gas", 5), ("e_gas", 1)):
            pa, pb = os.path.join(a.serial, st, f + ".bin"), os.path.join(a.backend, st, f + ".bin")
            if not (os.path.exists(pa) and os.path.exists(pb)):
                continue
            r, m = field_diff(pa, pb, nc)
            print(f"| {st} | {f} | " + ", ".join(f"{v:.1e}" for v in r) + " | " +
                  ", ".join(f"{v:.1e}" for v in m) + " |")
    print()
    print("| step | " + " | ".join(STATS) + " |")
    print("|---|" + "---|" * len(STATS))
    for st in steps:
        sa = read_kv(os.path.join(a.serial, st, "stats.txt"))
        sb = read_kv(os.path.join(a.backend, st, "stats.txt"))
        cells = []
        for k in STATS:
            if k == "Nhalo":
                cells.append(f"{int(sa[k])} / {int(sb[k])}")
            else:
                cells.append(f"{abs(sb[k] / sa[k] - 1):.1e}")
        print(f"| {st} | " + " | ".join(cells) + " |")
    print()
    for fa in sorted(glob.glob(os.path.join(a.serial, "DM", "xv_dm.*"))):
        fb = os.path.join(a.backend, "DM", os.path.basename(fa))
        if os.path.exists(fb):
            rx, mx, rv, mv = dm_diff(fa, fb, a.ndm)
            print(f"{os.path.basename(fa)}: positions rms {rx:.2e} max {mx:.2e} DM cells; "
                  f"velocities rms {rv:.2e} max {mv:.2e} of the rms velocity")
    for fa in sorted(glob.glob(os.path.join(a.serial, "Pk", "pk.*"))):
        fb = os.path.join(a.backend, "Pk", os.path.basename(fa))
        if os.path.exists(fb):
            x, y = np.loadtxt(fa), np.loadtxt(fb)
            if x.ndim == 2 and x.shape == y.shape:
                ok = np.abs(x) > 0
                rel = np.abs(y[ok] - x[ok]) / np.abs(x[ok])
                print(f"{os.path.basename(fa)}: largest relative difference {rel.max():.2e} "
                      f"(median {np.median(rel):.2e}) over {ok.sum()} entries")
            else:
                print(f"{os.path.basename(fa)}: shapes differ or not a table")


if __name__ == "__main__":
    main()
