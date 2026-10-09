#!/usr/bin/env python3
"""Self-test of make_movie.py on a synthetic run directory (docs/movies.md).

    test_make_movie.py [WORKDIR]

Builds a 16^3 run directory with a planar ionisation front moving along +x
at known positions: movie frames carry fHI and temp only (rhog comes from
the step_NNNN field outputs), and the first two frames are pre-RT frames
(fHI = rhog = 0, x_HII,M = nan), which must be skipped. Checks the front
positions found by the `planar` style, the tau_es integral against a
direct evaluation, and that both styles render (MP4 if ffmpeg is
available, PNG stills always). Exit code 77: matplotlib missing (skip).
"""
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

try:
    import matplotlib  # noqa: F401
except ImportError:
    print("matplotlib not installed: skipped")
    sys.exit(77)

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import make_movie as mm  # noqa: E402


def write_cube(path, a):
    np.asarray(a, dtype="<f4").ravel(order="F").tofile(path)


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp()
    run = os.path.join(work, "run")
    if os.path.exists(run):
        shutil.rmtree(run)
    os.makedirs(os.path.join(run, "movie"))
    with open(os.path.join(work, "params.txt"), "w") as f:
        f.write("N = 16\nbox = 8\n")
    N, box = 16, 8.0
    cell = box / N
    x = (np.arange(N) + 0.5) * cell
    zs = [12.0, 11.5, 11.0, 10.5, 10.0, 9.5, 9.0, 8.5]
    fronts = [None, None, 0.6, 1.7, 2.9, 4.0, 5.3, 6.6]  # Mpc/h
    lines = ["# frame istep istep_rt redshift xHII_V xHII_M"]
    rho = 1 + 0.5 * np.sin(2 * np.pi * x / box)[:, None, None] * np.ones((N, N, N))
    for k, (z, R) in enumerate(zip(zs, fronts)):
        if R is None:
            fHI = np.zeros((N, N, N))
            T = np.zeros((N, N, N))
            xv, xm = 1.0, float("nan")
        else:
            # smooth front: x_HII = 0.5 exactly at x = R
            xHII = 1 / (1 + np.exp((x - R) / (0.3 * cell)))
            fHI = np.broadcast_to(1 - xHII[:, None, None], (N, N, N))
            T = 1e2 + 2e4 * (1 - fHI)
            xv = float(np.mean(1 - fHI))
            xm = float(np.sum(rho * (1 - fHI)) / np.sum(rho))
        write_cube(os.path.join(run, "movie", "frame_%05d_fHI.bin" % k), fHI)
        write_cube(os.path.join(run, "movie", "frame_%05d_temp.bin" % k), T)
        lines.append("%6d %6d %6d %.8f %.8e %.8e" % (k, k, 0, z, xv, xm))
    with open(os.path.join(run, "movie", "frames.txt"), "w") as f:
        f.write("\n".join(lines).replace("nan", "-nan") + "\n")
    # field outputs with rhog at two steps
    for istep, z in ((3, 10.6), (6, 9.1)):
        d = os.path.join(run, "step_%04d" % istep)
        os.makedirs(d)
        write_cube(os.path.join(d, "rhog.bin"), rho * (1 + 0.01 * istep))
        with open(os.path.join(d, "stats.txt"), "w") as f:
            f.write("istep = %d\nredshift = %.17g\nxHII_V = 0.5\nxHII_M = 0.5\n" % (istep, z))

    R = mm.Run(run)
    ok = True
    if (R.N, R.box, len(R.frames), R.skipped) != (16, 8.0, 6, 2):
        print("FAIL run index:", R.N, R.box, len(R.frames), R.skipped)
        ok = False
    p, note = R.file_for(R.frames[0], "rhog")  # z = 11.0: nearest field output is step 3
    if not p.endswith("step_0003/rhog.bin") or "step 3" not in note:
        print("FAIL field fallback:", p, note)
        ok = False

    class A:
        ray_axis, source_face = "x", "low"
    got = mm.front_positions(R, A)
    want = np.array(fronts[2:])
    err = np.max(np.abs(got - want))
    print("front positions:", np.round(got, 3), "expected", want, "max error %.3f Mpc/h" % err)
    if err > 0.25 * cell:
        print("FAIL front positions")
        ok = False

    z = np.array([f["redshift"] for f in R.frames])
    xm = np.array([f["xHII_M"] for f in R.frames])
    tau = mm.tau_es_cumulative(z, xm)
    # direct evaluation with a fine grid on the piecewise-linear x_HII(z)
    zz = np.linspace(z.min(), z.max(), 20001)
    xx = np.interp(zz, z[::-1], xm[::-1])
    integ = (mm.SIGMA_T * mm.C_CGS * mm.nH0() * xx * (1 + mm.YHE / (4 * mm.XHY))
             * (1 + zz) ** 2 / mm.hubble(zz))
    direct = mm.trapz(integ, zz)
    rel = abs(tau[-1] - direct) / direct
    print("tau_es over the frames %.6e, direct %.6e (rel. diff %.1e)" % (tau[-1], direct, rel))
    if rel > 2e-2:  # trapezoid on 6 points vs fine grid of the same linear interpolant
        print("FAIL tau_es")
        ok = False
    t = mm.cosmic_time_myr([6.0])[0]
    print("t(z = 6) = %.1f Myr" % t)
    if not 900 < t < 1000:
        print("FAIL cosmic time")
        ok = False

    have_ff = shutil.which("ffmpeg") is not None
    for style in ("cosmo", "planar"):
        cmd = [sys.executable, os.path.join(here, "make_movie.py"), run, "--style", style,
               "-o", os.path.join(work, style + ".mp4"), "--size", "960", "540", "--dpi", "60",
               "--jobs", "2", "--stills", "0,5"]
        if style == "cosmo":
            cmd += ["--axis", "y", "--pos", "2.0", "--thick", "3", "--zoom", "1,7,1,7"]
        if not have_ff:
            cmd.append("--no-mp4")
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(r.stdout.strip())
        if r.returncode != 0:
            print("FAIL render", style, r.stderr)
            ok = False
            continue
        stills = os.listdir(os.path.join(work, style + "_stills"))
        if len(stills) != 2:
            print("FAIL stills", stills)
            ok = False
        if have_ff and os.path.getsize(os.path.join(work, style + ".mp4")) < 1000:
            print("FAIL mp4", style)
            ok = False
    print("PASS" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
