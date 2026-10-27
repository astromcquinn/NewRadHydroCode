#!/usr/bin/env python3
"""dual_e_self = 1 in apps/radhydro (docs/hydro_unsplit.md section 7, U7).

usage: dual_e_self_check.py RADHYDRO PROBLEM WORKDIR [--set key=value ...]

Six runs of PROBLEM without RT (zdata1 = 20), hydro_scheme = unsplit,
thermal_floor = 0.01 K (so over-cooled cells are not hidden at the 1 K
floor), u_gas and e_gas every step:
  A  dual_e_self = 1, one thread;
  B  as A with two threads;
  C  no dual_e_self key (today's scheme);
  D  dual_e_self = 0;
  E  thermal_init = recfast, no dual_e_self key;
  F  thermal_init = recfast, dual_e_self = 1.
Checks: B byte-identical to A, D to C, E to F (the key defaults to 1 with
recfast, user decision 2026-10-27); at every step no cell of A colder
than 0.1 x its adiabat temp_min(z) Delta^(2/3) (the defect) and the median
temperature of the mean-density cells (0.9 < Delta < 1.1) within 1 % of
temp_min(z); C shows the defect at the last step (more than 5 % of the cells
below 0.1 x the adiabat), so the problem is sensitive to it; the hottest
cell of A within a factor 2 of C's (shock heating kept). Exit 0 on success.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools", "analysis"))
import thermal_history as th  # noqa: E402
import thermal_check as tc  # noqa: E402  (run, same_tree, params, check)


def history(run, N, box):
    """Per step: z, the fraction of cells below 0.1 x the adiabat, the
    median T/temp_min(z) of the mean-density cells, the largest T."""
    rows = []
    for d in sorted(os.listdir(run)):
        p = os.path.join(run, d)
        if not (d.startswith("step_") and os.path.exists(os.path.join(p, "e_gas.bin"))):
            continue
        st = dict(l.split(" = ", 1) for l in open(os.path.join(p, "stats.txt")) if " = " in l)
        z = float(st["redshift"])
        u = np.fromfile(os.path.join(p, "u_gas.bin"), dtype=np.float32).reshape(-1, 5).astype(np.float64)
        e = np.fromfile(os.path.join(p, "e_gas.bin"), dtype=np.float32).astype(np.float64)
        rho = u[:, 0]
        T = e / rho * th.temp1_unit(1 / (1 + z), box, N)
        D = rho * th.OM / th.OB
        cold = np.mean(T < 0.1 * th.floor(z) * D ** (2 / 3))
        mean = (D > 0.9) & (D < 1.1)
        rows.append((z, cold, np.median(T[mean]) / th.floor(z), T.max()))
    return np.array(rows)


def main():
    exe, prob, work = sys.argv[1:4]
    extra = sys.argv[4:]
    N, box = tc.params(prob)
    base = ["--set", "zdata1=20", "--set", "hydro_scheme=unsplit", "--set", "thermal_floor=0.01",
            "--set", "fields_every=1", "--set", "fields=u_gas,e_gas"] + extra
    A, B, C, D, E, F = (os.path.join(work, n) for n in "ABCDEF")
    tc.run(exe, prob, A, base + ["--set", "dual_e_self=1"], 1)
    tc.run(exe, prob, B, base + ["--set", "dual_e_self=1"], 2)
    tc.run(exe, prob, C, base, 1)
    tc.run(exe, prob, D, base + ["--set", "dual_e_self=0"], 1)
    tc.run(exe, prob, E, base + ["--set", "thermal_init=recfast"], 1)
    tc.run(exe, prob, F, base + ["--set", "thermal_init=recfast", "--set", "dual_e_self=1"], 1)
    tc.same_tree(A, B, "dual_e_self = 1: threads 1 vs 2")
    tc.same_tree(C, D, "dual_e_self = 0 vs no key")
    tc.same_tree(E, F, "thermal_init = recfast: no key vs dual_e_self = 1")
    hA, hC = history(A, N, box), history(C, N, box)
    tc.check(len(hA) >= 10 and len(hA) == len(hC), f"{len(hA)} and {len(hC)} steps")
    print("#      z  cold(A)  cold(C)  Tmean/temp_min(A)  (C)")
    for a, c in zip(hA, hC):
        print(f"{a[0]:8.2f} {a[1]:8.4f} {c[1]:8.4f} {a[2]:10.4f} {c[2]:10.4f}")
    tc.check(hA[:, 1].max() == 0, f"dual_e_self = 1: up to {hA[:, 1].max():.4f} of the cells over-cooled")
    worst = np.abs(hA[:, 2] - 1).max()
    tc.check(worst < 0.01, f"dual_e_self = 1: mean-density T off temp_min(z) by {worst:.4f}")
    tc.check(hC[-1, 1] > 0.05, f"no key: only {hC[-1, 1]:.4f} of the cells over-cooled (problem insensitive)")
    tc.check(hA[-1, 3] > 0.5 * hC[-1, 3],
             f"dual_e_self = 1: hottest cell {hA[-1, 3]:.1f} K against {hC[-1, 3]:.1f} K without")
    print(f"  dual_e_self = 1: max |Tmean/temp_min - 1| = {worst:.4f}; hottest cell {hA[-1, 3]:.1f} K "
          f"(without {hC[-1, 3]:.1f} K); over-cooled without: {hC[-1, 1]:.4f} at z = {hC[-1, 0]:.1f}")
    if tc.bad:
        sys.exit(f"{len(tc.bad)} failure(s)")
    print("OK")


if __name__ == "__main__":
    main()
