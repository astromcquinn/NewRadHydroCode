#!/usr/bin/env python3
"""Timing and memory estimate for the production run (docs/production_run.md §5).

    timing_table.py [--box 10] [--zstart 15] [--zend 5.5] [--z-init 100]

Cost model, fitted to the measured runs (4 threads, deterministic_parallel):
  * hydro step (hydro + PM + summaries + FoF): t_h * (N/128)^3 seconds;
    hydro steps: ra_max steps (da/a = 1/9) from z_init to zend;
  * RT substeps: the light travel at c_hat in RT cells from the first source
    (z_src(N)) to zend: 2c/(H0 sqrt(Om)) (sqrt(a2)-sqrt(a1)) / (box/N) * c_red;
  * per substep: a_cell * N^3 + b_ray * <Nray>, where <Nray> = f_ray *
    Nray_max is the substep-weighted mean ray count (the population grows
    from 0 and saturates near 0.6-0.7 Nray_max with the pool guard);
  * memory: m_cell * N^3 (grids, particles with Ndm = N) + 88 B * peak rays
    (the runs peaked at 0.5-0.65 Nray_max; the pool guard caps it at 0.9).
The fitted constants are set below from the measurements quoted in the doc.
"""
import argparse
import math

C = 2.99792458e10
H0 = 100e5 / 3.0856775807e24  # h/s per (Mpc/h)
OM = 0.27


def light_cells(z1, z2, box, N):
    a1, a2 = 1 / (1 + z1), 1 / (1 + z2)
    d = 2 * C / (H0 * math.sqrt(OM)) * (math.sqrt(a2) - math.sqrt(a1))  # cm per (Mpc/h) units
    d_mpch = d / 3.0856775807e24
    return d_mpch / (box / N)


def hydro_steps(z1, z2):
    a, n = 1 / (1 + z1), 0
    while 1 / a - 1 > z2:
        a += a * 0.1 / 0.9
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--box", type=float, default=10.0)
    ap.add_argument("--z-init", type=float, default=100.0)
    ap.add_argument("--zstart", type=float, default=15.0)
    ap.add_argument("--zend", type=float, default=5.5)
    ap.add_argument("--nrf", type=int, default=16)
    # fitted constants (4 threads), see docs/production_run.md section 5
    ap.add_argument("--t-hydro128", type=float, required=True,
                    help="s per hydro step at N = 128 (hydro+PM+summaries+FoF)")
    ap.add_argument("--a-cell", type=float, required=True, help="s per cell per RT substep")
    ap.add_argument("--b-ray", type=float, required=True, help="s per ray per RT substep")
    ap.add_argument("--f-ray", type=float, default=0.4, help="<Nray>/Nray_max over the RT phase")
    ap.add_argument("--m-cell", type=float, required=True, help="bytes per cell without rays")
    ap.add_argument("--f-peak", type=float, default=0.75, help="peak rays / Nray_max")
    ap.add_argument("--zsrc", default="64:9,128:11,192:12",
                    help="redshift of the first FoF source per N (RT starts there)")
    a = ap.parse_args()
    zsrc = {int(k): float(v) for k, v in (x.split(":") for x in a.zsrc.split(","))}
    print("| N | c_reduced | hydro steps | RT substeps | t hydro | t RT | total | peak memory |")
    print("|---|---|---|---|---|---|---|---|")
    for N in (64, 128, 192):
        nh = hydro_steps(a.z_init, a.zend)
        th = nh * a.t_hydro128 * (N / 128) ** 3
        nray_max = a.nrf * N ** 3
        for c in (0.01, 0.03, 0.1):
            z0 = min(a.zstart, zsrc.get(N, a.zstart))
            ns = light_cells(z0, a.zend, a.box, N) * c
            per = a.a_cell * N ** 3 + a.b_ray * a.f_ray * nray_max
            trt = ns * per
            mem = a.m_cell * N ** 3 + 88 * a.f_peak * nray_max
            fmt = lambda s: f"{s / 3600:.1f} h" if s > 5400 else f"{s / 60:.0f} min"  # noqa: E731
            print(f"| {N} | {c:g} | {nh} | {ns:.0f} | {fmt(th)} | {fmt(trt)} | {fmt(th + trt)} | "
                  f"{mem / 1e9:.1f} GB |")


if __name__ == "__main__":
    main()
