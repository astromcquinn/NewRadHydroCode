#!/usr/bin/env python3
"""Nightly: the 64^3 production configuration (tools/production/
make_production.py --N 64) with rt_backend = amrex (2 ranks x 2 threads)
against the serial code, both legacy_rt = legacy_physics = 0
(docs/stage_d.md section 4). Compares the history (x_HII,V/M and T of
step_*/stats.txt), tau_es from the movie frames (tools/analysis/reion.py)
and the per-cell fields of the last step. The differences are round-off
(deposit and merge summation order) amplified by the chaotic growth of a
reionisation run; the tolerances are set from the measured values
(docs/stage_d.md section 4) with a margin.

With MODE = state (docs/stage_e_pm.md section 3.1): the production
configuration as make_production.py writes it today (unsplit hydro,
redshift_per_distance) with the three backends and the state on the levels
(rt_backend = hydro_backend = pm_backend = amrex, amrex_state = 1) against
the serial code. There the PM backend's round-off (fixed-point CIC, pencil
FFT) also enters, through the gravity, and the tolerances are TOL_STATE.

usage: amrex_full64_check.py RADHYDRO WORKDIR TOOLS_DIR [rt|state]
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy as np

TOL = {"xHII_V": 2e-2, "xHII_M": 2e-2, "T_V": 2e-2, "T_M": 2e-2, "tau": 1e-2,
       "fHI_rms": 5e-2, "T_rms": 5e-2}
TOL_STATE = dict(TOL, rho_rms=5e-2)


def read_kv(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = float(v.split()[0])
    return d


def run(cmd, log):
    with open(log, "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(open(log).read()[-3000:])
        sys.exit(f"failed: {' '.join(cmd)}")


def main():
    exe, work, tools = sys.argv[1], sys.argv[2], sys.argv[3]
    mode = sys.argv[4] if len(sys.argv) > 4 else "rt"
    if mode not in ("rt", "state"):
        sys.exit(f"unknown mode {mode}")
    tol = TOL if mode == "rt" else TOL_STATE
    launch = os.environ.get("RADRAY_MPIEXEC", "mpirun")
    flags = os.environ.get("RADRAY_MPIEXEC_FLAGS", "--oversubscribe --allow-run-as-root --bind-to none").split()
    prob = os.path.join(work, "prob")
    if not os.path.exists(os.path.join(prob, "params.txt")):
        # hydro_scheme and redshift_per_distance pinned to the configuration the tolerances were
        # measured with (docs/stage_d.md section 4), before make_production.py's production
        # defaults became unsplit and per distance (2026-10-11)
        if mode == "rt":
            run([sys.executable, os.path.join(tools, "production", "make_production.py"), prob, "--N", "64",
                 "--movie-dz", "0", "--hydro-scheme", "legacy", "--redshift-per-distance", "0",
                 "--recomb-case", "B", "--recomb-rad", "0"],  # recomb_rad: serial RT only
                os.path.join(work, "make.log"))
        else:
            # today's production defaults (unsplit, redshift_per_distance, case B
            # without recombination radiation: amrex_state needs recomb_rad = 0)
            run([sys.executable, os.path.join(tools, "production", "make_production.py"), prob, "--N", "64",
                 "--movie-dz", "0", "--recomb-rad", "0"], os.path.join(work, "make.log"))
    common = ["--quiet", "--set", "restart_every=0", "--set", "movie_dz=0.25", "--set", "movie_z0=15"]
    runs = {"serial": [exe, prob, os.path.join(work, "serial")] + common,
            "amrex": [launch] + flags + ["-np", "2", "-x", "OMP_NUM_THREADS=2", exe, prob,
                                           os.path.join(work, "amrex")] + common +
                     ["--set", "rt_backend=amrex", "--set", "deterministic_parallel=0"]}
    if mode == "state":
        runs["amrex"] += ["--set", "hydro_backend=amrex", "--set", "pm_backend=amrex",
                          "--set", "amrex_hydro_max_grid_size=0", "--set", "amrex_state=1"]
    for name, cmd in runs.items():
        d = os.path.join(work, name)
        if os.path.isdir(d):
            shutil.rmtree(d)
        run(cmd, d + ".log")
    bad = []
    # history
    steps = sorted(os.path.basename(s) for s in glob.glob(os.path.join(work, "serial", "step_*")))
    worst = {k: 0.0 for k in ("xHII_V", "xHII_M", "T_V", "T_M")}
    for s in steps:
        a = read_kv(os.path.join(work, "serial", s, "stats.txt"))
        bp = os.path.join(work, "amrex", s, "stats.txt")
        if not os.path.exists(bp):
            bad.append(f"{s}: missing in the amrex run")
            continue
        b = read_kv(bp)
        for k in worst:
            if a[k] > 1e-3:
                worst[k] = max(worst[k], abs(b[k] / a[k] - 1))
    for k, v in worst.items():
        print(f"  {k:7s} max rel diff {v:.3e} (tol {tol[k]})")
        if v > tol[k]:
            bad.append(f"{k} differs by {v:.3e}")
    # tau_es from the frames
    sys.path.insert(0, os.path.join(tools, "analysis"))
    import reion  # noqa: E402
    taus = {}
    for name in ("serial", "amrex"):
        fr = np.loadtxt(os.path.join(work, name, "movie", "frames.txt"))
        z, tau = reion.tau_running(fr[:, 3], fr[:, 5])
        taus[name] = tau[-1]
    dtau = abs(taus["amrex"] / taus["serial"] - 1)
    print(f"  tau_es  serial {taus['serial']:.5f} amrex {taus['amrex']:.5f} rel diff {dtau:.3e} (tol {tol['tau']})")
    if dtau > tol["tau"]:
        bad.append(f"tau_es differs by {dtau:.3e}")
    # per-cell fields of the last step
    last = steps[-1]
    cells = [("fHI", "fHI_rms"), ("temp", "T_rms")] + ([("rhog", "rho_rms")] if mode == "state" else [])
    for f, key in cells:
        x = np.fromfile(os.path.join(work, "serial", last, f + ".bin"), dtype="<f4").astype(float)
        y = np.fromfile(os.path.join(work, "amrex", last, f + ".bin"), dtype="<f4").astype(float)
        rms = np.sqrt(np.mean((y - x) ** 2)) / np.sqrt(np.mean(x ** 2))
        print(f"  {f:5s} last step rms rel diff {rms:.3e} (tol {tol[key]})")
        if rms > tol[key]:
            bad.append(f"{f} rms differs by {rms:.3e}")
    for x in bad:
        print("  FAIL", x)
    print("PASS" if not bad else "FAIL")
    sys.exit(0 if not bad else 1)


if __name__ == "__main__":
    main()
