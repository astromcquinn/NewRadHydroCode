#!/usr/bin/env python3
"""Calibrate f_esc of the abundance-matched sources to the Planck tau_es.

    calibrate_fesc.py OUTDIR --N 64 --am-params runs/am/am_params.txt
                      --fesc 0.05 0.1 0.2 [--threads 4] [--target 0.054]
                      [--movie-fields fHI] [--extra key=value ...]

For each f_esc: writes OUTDIR/fesc_<f>/params.txt with
tools/production/make_production.py (the production configuration at N,
abundance-matched sources, xi_ion = 25.38; no checkpoints), links the
shared ICs (OUTDIR/ics, made once), runs build/radhydro, and computes the
reionization history and tau_es (tools/analysis/reion.py) from the movie
log. Prints a table and the f_esc that gives the target tau_es by
interpolation in log f_esc. Runs that already finished are not rerun.

Biases (docs/production_run.md, docs/hmf_validation.md): at N = 64 only
halos above 20 particles = 8e9 Msun emit, and c_reduced = 0.01 and the
missing sub-grid sources delay reionization, so the calibrated f_esc is
biased high relative to the true f_esc of these galaxies; check it with
one run at the production N.
"""
import argparse
import os
import subprocess
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, "..", ".."))
sys.path.insert(0, os.path.join(root, "tools", "analysis"))
import reion  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--N", type=int, default=64)
    ap.add_argument("--am-params", required=True)
    ap.add_argument("--fesc", type=float, nargs="+", required=True)
    ap.add_argument("--xi-ion", type=float, default=25.38)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--target", type=float, default=reion.PLANCK_TAU)
    ap.add_argument("--movie-dz", type=float, default=0.05)
    ap.add_argument("--movie-fields", default="fHI")
    ap.add_argument("--extra", nargs="*", default=["data_files=dm", "fields_at_data1=0"],
                    help="key=value lines added to params.txt (default: lean outputs, about "
                         "0.7 GB per 128^3 run instead of 6 GB)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    ics = os.path.abspath(os.path.join(a.outdir, "ics"))
    if not os.path.exists(os.path.join(ics, "xv_dm.init")):
        subprocess.run([sys.executable, os.path.join(here, "make_production.py"), a.outdir,
                        "--N", str(a.N)], check=True)
    rows = []
    for f in a.fesc:
        d = os.path.join(a.outdir, f"fesc_{f:g}")
        run = os.path.join(d, "run")
        done = os.path.join(run, "movie", "frames.txt")
        if not (os.path.exists(os.path.join(d, "finished"))):
            subprocess.run([sys.executable, os.path.join(here, "make_production.py"), d, "--N",
                            str(a.N), "--no-ics", "--am-params", a.am_params, "--xi-ion",
                            str(a.xi_ion), "--f-esc", str(f), "--movie-dz", str(a.movie_dz),
                            "--movie-fields", a.movie_fields, "--restart-every", "0"], check=True)
            with open(os.path.join(d, "params.txt"), "a") as p:
                for kv in a.extra:
                    p.write(kv.replace("=", " = ", 1) + "\n")
            if not os.path.exists(os.path.join(d, "ics")):
                os.symlink(ics, os.path.join(d, "ics"))
            env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))
            with open(os.path.join(d, "log.txt"), "w") as log:
                subprocess.run([os.path.join(root, "build", "radhydro"), d, run, "--quiet"],
                               check=True, stdout=log, stderr=subprocess.STDOUT, env=env)
            open(os.path.join(d, "finished"), "w").close()
        fr, z, xV, xM = reion.read_frames(run)
        zz, tau = reion.tau_history(z, xM)
        o = np.argsort(-z)
        zs, xs = z[o], xV[o]

        def zcross(lev):
            k = np.where(xs >= lev)[0]
            if not len(k):
                return float("nan")
            k = k[0]
            return zs[k] if k == 0 else float(np.interp(lev, [xs[k - 1], xs[k]], [zs[k - 1], zs[k]]))

        rows.append((f, tau[-1], zcross(0.1), zcross(0.5), zcross(0.99), xM[np.argmin(z)]))
    print("# f_esc  tau_es  z(x_V=0.1)  z(x_V=0.5)  z(x_V=0.99)  x_M(end)")
    for r in rows:
        print(f"{r[0]:.4g}  {r[1]:.5f}  {r[2]:.2f}  {r[3]:.2f}  {r[4]:.2f}  {r[5]:.4f}")
    t = np.array([r[1] for r in rows])
    lf = np.log(np.array([r[0] for r in rows]))
    if len(rows) >= 2 and t.min() <= a.target <= t.max():
        o = np.argsort(t)
        fbest = float(np.exp(np.interp(a.target, t[o], lf[o])))
        print(f"# f_esc for tau_es = {a.target}: {fbest:.4g}")
    else:
        print(f"# target tau_es = {a.target} outside the range of these runs")


if __name__ == "__main__":
    main()
