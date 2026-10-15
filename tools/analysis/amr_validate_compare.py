#!/usr/bin/env python3
"""Compare the refined production-configuration run with the uniform runs
(tools/hyak/amr_validate_runs.sh, docs/stage_f.md section 3).

usage: amr_validate_compare.py ROOT [--ref u256] [--ref-N 256] [--N 128]
                               [--runs u128 amr128] [--json OUT] [--zmax 15]

ROOT holds one problem directory per run (ROOT/<name>/run is the radhydro
output). Every run is compared with the reference (default u256, the uniform
run at the first refined level's resolution):
  * x_HII,V(z), x_HII,M(z) (movie/frames.txt), the redshifts of x_HII,V =
    0.1, 0.5, 0.9, 0.99, tau_es (tools/analysis/reion.py);
  * T_M and T_V of every hydro step (stats.txt), max relative error at
    z <= zmax;
  * at the zdata1 field outputs, the reference block-averaged onto the run's
    grid (rho and fHI by volume, T by mass): the acceptance metrics of
    tools/ray_compression/precision_metrics.py (x_HII,V/M and T_M max
    relative error, the front shift, rms T and x_HII per ionised cell);
  * the cost (timing.txt) and, for the backends, the last hydro timers line
    of run.log (levels, regrid, reflux, pushdown).
The functions are those of tools/analysis/rt_multicell_compare.py.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "ray_compression"))
import precision_metrics as pm  # noqa: E402
from rt_multicell_compare import history_metrics, timing, tm_metrics  # noqa: E402


def coarsen(steps, f):
    """block-average a run's field outputs by f per dimension"""
    if f == 1:
        return steps
    out = []
    for s in steps:
        n = s["rho"].shape[0] // f
        blk = lambda a: a.reshape(n, f, n, f, n, f).mean(axis=(1, 3, 5))  # noqa: E731
        rho = blk(s["rho"])
        c = dict(s)
        c.update(rho=rho, fHI=blk(s["fHI"]), T=blk(s["T"] * s["rho"]) / rho)
        out.append(c)
    return out


def hydro_timers(prob):
    log = Path(prob) / "run.log"
    if not log.exists():
        return ""
    last = ""
    for line in open(log, errors="replace"):
        if "AMReX hydro timers" in line:
            last = line
    m = re.search(r"levels .*?pushdown [0-9.]+", last)
    return m.group(0) if m else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--ref", default="u256")
    ap.add_argument("--ref-N", type=int, default=256)
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--runs", nargs="*", default=None,
                    help="default: every other directory of ROOT with run/movie/frames.txt")
    ap.add_argument("--json")
    ap.add_argument("--zmax", type=float, default=15.0, help="T compared at z <= zmax (RT on)")
    a = ap.parse_args()
    root = Path(a.root)
    ref = root / a.ref / "run"
    names = a.runs if a.runs is not None else sorted(
        p.name for p in root.iterdir()
        if p.is_dir() and p.name != a.ref and (p / "run" / "movie" / "frames.txt").exists())
    if a.ref_N % a.N:
        sys.exit("--ref-N must be a multiple of --N")
    rf = coarsen(pm.load_run(ref, a.ref_N), a.ref_N // a.N)
    res = {a.ref: {"timing": timing(ref), "history": history_metrics(ref, ref)}}
    for n in names:
        run = root / n / "run"
        r = {"history": history_metrics(ref, run), "T": tm_metrics(ref, run, a.zmax),
             "timing": timing(run), "hydro": hydro_timers(root / n)}
        r["fields"], _ = pm.metrics(rf, pm.load_run(run, a.N))
        res[n] = r
    b = res[a.ref]
    print(f"reference {a.ref}: tau_es {b['history']['tau_es_base']:.5f}, z(x_V = 0.5) "
          f"{b['history']['z_x0.5_base']:.3f}, z(0.99) {b['history']['z_x0.99_base']:.3f}; "
          f"wall {b['timing']['wall']:.0f} s, RT {b['timing']['t_rt']:.0f} s, "
          f"{b['timing']['hydro_steps']} hydro steps")
    print("errors against the reference (fields: the reference averaged onto the run's grid)")
    hdr = ("run", "xV(z) rel", "xM(z) rel", "tau_es rel", "dz(0.5)", "dz(0.99)", "T_M rel",
           "fld xV", "fld xM", "fld T_M", "front", "T_ion rms", "x_ion rms", "wall [s]", "RT [s]")
    print(" | ".join(hdr))
    for n in names:
        r = res[n]
        h, f, t, T = r["history"], r["fields"], r["timing"], r["T"]
        row = (n, f"{h['xHII_V_rel_max']:.3e}", f"{h['xHII_M_rel_max']:.3e}", f"{h['tau_es_rel']:+.3e}",
               f"{h['z_x0.5'] - h['z_x0.5_base']:+.3f}", f"{h['z_x0.99'] - h['z_x0.99_base']:+.3f}",
               f"{T['T_M_rel_max']:.3e}", f"{f['xHII_vol']:.3e}", f"{f['xHII_mass']:.3e}",
               f"{f['T_mass']:.3e}", f"{f['front_shift']:.3f}",
               f"{f.get('T_ion_cell_rms', float('nan')):.3e}",
               f"{f.get('xHII_ion_cell_rms', float('nan')):.3e}",
               f"{t['wall']:.0f}", f"{t['t_rt']:.0f}")
        print(" | ".join(row))
    for n in names:
        if res[n]["hydro"]:
            print(f"{n} hydro: {res[n]['hydro']}")
    if a.json:
        with open(a.json, "w") as fjs:
            json.dump(res, fjs, indent=1)


if __name__ == "__main__":
    main()
