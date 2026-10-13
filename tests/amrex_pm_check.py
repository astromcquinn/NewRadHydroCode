#!/usr/bin/env python3
"""Checks of the AMReX PM backend (pm_backend = amrex, docs/stage_e_pm.md).

usage: amrex_pm_check.py run  OUT LAUNCH "LAUNCH_FLAGS" NPFLAG NP NT DRIVER PROBLEM [--set k=v ...]
       amrex_pm_check.py bitwise  REF NEW
       amrex_pm_check.py roundoff REF NEW [--rtol R] [--atol-frac A] [--stats-rtol S]
       amrex_pm_check.py runcmp bitwise|roundoff REF OUT LAUNCH "LAUNCH_FLAGS" NPFLAG NP NT DRIVER PROBLEM
                         [--rtol R --atol-frac A --stats-rtol S] -- [driver args]

run       runs pm_driver (DRIVER) on PROBLEM under the MPI launcher (NP ranks,
          OMP_NUM_THREADS = NT) into OUT (removed first); NP = 1 runs without
          the launcher when LAUNCH is empty.
bitwise   every step_NNNN/*.bin and stats.txt of REF and NEW byte-identical
          (the rank, thread and box-size independence of the backend).
roundoff  the fields of NEW against REF within rtol (relative to the larger
          magnitude) or atol_frac times the field's maximum magnitude (the
          backend against the serial PM: the deposit order and the FFT);
          the stats: istep/nstep_rt exact; the DM velocity statistics within
          stats-rtol times the velocity scale (the mean velocity of a
          zero-mean field is round-off), the DM position statistics within
          stats-rtol times the box, the other floats within stats-rtol
          relative. The pm_refine > 1 outputs (phi.bin on the
          mesh, Ndm > Ngas) are read by size.
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy as np


def steps(d):
    return sorted(glob.glob(os.path.join(d, "step_*")))


def read_stats(fn):
    out = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = float(v)
    return out


def run(argv):
    out, launch, flags, npflag, np_, nt, drv, prob = argv[:8]
    args = argv[8:]
    flags = [f for f in flags.replace(";", " ").split() if f]
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    np_ = int(np_)
    cmd = ([launch] + flags + [npflag, str(np_)] if launch else []) + [drv, prob, out, "--quiet"] + args
    env = dict(os.environ, OMP_NUM_THREADS=str(nt))
    with open(out + ".log", "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        print(open(out + ".log").read()[-3000:])
        sys.exit(f"pm_driver failed ({np_} ranks, {nt} threads)")
    return 0


def bitwise(ref, new):
    sr, sn = steps(ref), steps(new)
    if not sr or len(sr) != len(sn):
        sys.exit(f"step count differs or no steps: {len(sr)} vs {len(sn)}")
    bad, n = [], 0
    for a, b in zip(sr, sn):
        for fa in sorted(glob.glob(os.path.join(a, "*"))):
            fb = os.path.join(b, os.path.basename(fa))
            n += 1
            if not os.path.exists(fb):
                bad.append(f"{fb} missing")
            elif open(fa, "rb").read() != open(fb, "rb").read():
                bad.append(f"{os.path.relpath(fb, new)} differs")
    for x in bad[:20]:
        print("  DIFF", x)
    print(f"{n} files compared: {'IDENTICAL' if not bad else 'DIFFERENT'}")
    return 0 if not bad else 1


def roundoff(ref, new, rtol, atol_frac, stats_rtol):
    sr, sn = steps(ref), steps(new)
    if not sr or len(sr) != len(sn):
        sys.exit(f"step count differs or no steps: {len(sr)} vs {len(sn)}")
    worst = {}
    fails = []
    for a, b in zip(sr, sn):
        step = os.path.basename(a)
        for name in ("xv_dm", "u_gas", "v_gas", "phi", "w_gas"):
            fa, fb = os.path.join(a, name + ".bin"), os.path.join(b, name + ".bin")
            if not os.path.exists(fa):
                continue
            x = np.fromfile(fa, dtype=np.float32).astype(np.float64)
            y = np.fromfile(fb, dtype=np.float32).astype(np.float64)
            if x.size != y.size:
                fails.append(f"{step} {name}: size {x.size} vs {y.size}")
                continue
            d = np.abs(x - y)
            scale = np.maximum(np.abs(x), np.abs(y))
            amax = np.abs(x).max()
            if name == "xv_dm":
                # positions in DM cells (scale: the box), velocities (scale:
                # the largest velocity)
                x6 = x.reshape(-1, 6)
                sc = np.empty_like(x6)
                sc[:, :3] = np.abs(x6[:, :3]).max()
                sc[:, 3:] = np.abs(x6[:, 3:]).max()
                ok = d <= rtol * scale + atol_frac * sc.reshape(-1)
            else:
                ok = d <= rtol * scale + atol_frac * amax
            rel = np.where(scale > 0, d / np.maximum(scale, 1e-300), 0.0)
            w = worst.get(name, (0.0, "", 0.0))
            if rel.max() > w[0]:
                i = int(rel.argmax())
                worst[name] = (float(rel.max()), f"{step} [{i}] ref={x[i]:.6g} new={y[i]:.6g}", float(d.max() / max(amax, 1e-300)))
            if not ok.all():
                i = int(np.argmax(~ok))
                fails.append(f"{step} {name}: {int((~ok).sum())} values beyond tolerance, e.g. [{i}] ref={x[i]:.6g} new={y[i]:.6g}")
        sa, sb = read_stats(os.path.join(a, "stats.txt")), read_stats(os.path.join(b, "stats.txt"))
        vscale = max(abs(sa.get("dm_vx_max", 0)), abs(sa.get("dm_vy_max", 0)), abs(sa.get("dm_vz_max", 0)),
                     abs(sa.get("dm_vx_min", 0)), abs(sa.get("dm_vy_min", 0)), abs(sa.get("dm_vz_min", 0)), 1e-300)
        xscale = max(abs(sa.get("dm_x_max", 0)), abs(sa.get("dm_y_max", 0)), abs(sa.get("dm_z_max", 0)), 1e-300)
        for k, va in sa.items():
            vb = sb.get(k)
            if vb is None:
                fails.append(f"{step} stat {k} missing")
                continue
            if k in ("istep", "nstep_rt"):
                if va != vb:
                    fails.append(f"{step} stat {k}: {va} vs {vb}")
            elif k.startswith("dm_v"):
                if abs(va - vb) > stats_rtol * vscale:
                    fails.append(f"{step} stat {k}: {va:.9g} vs {vb:.9g} (velocity scale {vscale:.3g})")
            elif k.startswith("dm_"):
                # positions: the minimum is a near-zero value, scale = the box
                if abs(va - vb) > stats_rtol * xscale:
                    fails.append(f"{step} stat {k}: {va:.9g} vs {vb:.9g} (position scale {xscale:.3g})")
            else:
                if abs(va - vb) > stats_rtol * max(abs(va), abs(vb)):
                    fails.append(f"{step} stat {k}: {va:.9g} vs {vb:.9g}")
    for name, (r, where, dfrac) in sorted(worst.items()):
        print(f"  {name}: worst relative difference {r:.3e} ({where}); max |diff| / max |field| {dfrac:.3e}")
    for x in fails[:20]:
        print("  FAIL", x)
    print(f"{len(sr)} steps compared: {'PASS' if not fails else 'FAIL'}")
    return 0 if not fails else 1


def main():
    mode = sys.argv[1]
    if mode == "run":
        return run(sys.argv[2:])
    if mode == "runcmp":
        # runcmp bitwise|roundoff REF OUT LAUNCH FLAGS NPFLAG NP NT DRIVER PROBLEM [tolerance options] -- [driver args]
        cmp_mode, ref = sys.argv[2], sys.argv[3]
        rest = sys.argv[4:]
        dash = rest.index("--") if "--" in rest else len(rest)
        head, drv_args = rest[:dash], rest[dash + 1:]
        run(head[:8] + drv_args)
        if cmp_mode == "bitwise":
            return bitwise(ref, head[0])
        import argparse
        ap = argparse.ArgumentParser()
        ap.add_argument("--rtol", type=float, default=1e-4)
        ap.add_argument("--atol-frac", type=float, default=1e-6)
        ap.add_argument("--stats-rtol", type=float, default=1e-5)
        a = ap.parse_args(head[8:])
        return roundoff(ref, head[0], a.rtol, a.atol_frac, a.stats_rtol)
    if mode == "bitwise":
        return bitwise(sys.argv[2], sys.argv[3])
    if mode == "roundoff":
        import argparse
        ap = argparse.ArgumentParser()
        ap.add_argument("ref")
        ap.add_argument("new")
        ap.add_argument("--rtol", type=float, default=1e-4)
        ap.add_argument("--atol-frac", type=float, default=1e-6)
        ap.add_argument("--stats-rtol", type=float, default=1e-5)
        a = ap.parse_args(sys.argv[2:])
        return roundoff(a.ref, a.new, a.rtol, a.atol_frac, a.stats_rtol)
    sys.exit(__doc__)


if __name__ == "__main__":
    sys.exit(main())
