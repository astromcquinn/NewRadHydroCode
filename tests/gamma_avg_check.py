#!/usr/bin/env python3
"""The substep-averaged Gamma_HI of apps/radhydro (fields / movie_fields
gammaHI_avg; docs/igm_diagnostics.md).

usage: gamma_avg_check.py RADHYDRO PROBLEM WORKDIR [--set key=value ...]

Three runs of PROBLEM with fields every step and a movie frame after every RT
substep (movie_every_substep = 1, movie_slab = 1):
  A  one thread, gammaHI and gammaHI_avg in fields and movie_fields;
  B  as A with two threads;
  C  as A without gammaHI_avg.
Checks, in every hydro step with RT substeps: the gammaHI_avg frame after
substep k is the mean of the gammaHI frames of substeps 1..k (the substep
length is constant within a hydro step); the step's gammaHI_avg.bin equals
the last frame's slabs exactly, and gammaHI.bin the last gammaHI frame; at
least two such steps with two or more substeps were seen. Then every file of
B is byte-identical to A's (thread independence), and every file of C to A's
file of the same name, the gammaHI_avg files being the only ones C lacks (the
field is output only). Exit 0 on success.
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy as np

SKIP = ("timing.txt", "ray.log", "rt.log", "stdout.log")
bad = []


def check(cond, msg):
    if not cond:
        bad.append(msg)
        print("  FAIL", msg)


def run(exe, prob, out, sets, threads):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    cmd = [exe, prob, out, "--quiet"] + sets
    with open(os.path.join(out, "stdout.log"), "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        print(open(os.path.join(out, "stdout.log")).read()[-3000:])
        sys.exit(f"failed ({r.returncode}): {' '.join(cmd)}")


def slab_geometry(run):
    g = {}
    for line in open(os.path.join(run, "movie", "slab.txt")):
        k, _, v = line.partition("=")
        if _ and not line.startswith("#"):
            g[k.strip()] = v.strip()
    return int(g["N"]), int(g["thick"]), g["axes"], int(g["start"])


def cube_slabs(cube, thick, axes, start):
    """The slab file's content for a (N, N, N) [k, j, i] cube: for each axis
    the cells start .. start+thick-1 along it, Fortran order (nx, ny, nz)."""
    parts = []
    for a in axes:
        if a == "x":
            parts.append(cube[:, :, start:start + thick].ravel())
        elif a == "y":
            parts.append(cube[:, start:start + thick, :].ravel())
        else:
            parts.append(cube[start:start + thick, :, :].ravel())
    return np.concatenate(parts)


def frame_rows(run):
    rows = []
    for line in open(os.path.join(run, "movie", "frames.txt")):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        rows.append((int(f[0]), int(f[1]), int(f[2]), float(f[3])))
    return rows


def frame(run, iframe, field):
    return np.fromfile(os.path.join(run, "movie", "frame_%05d_%s.slab" % (iframe, field)),
                       dtype=np.float32)


def check_average(run):
    N, thick, axes, start = slab_geometry(run)
    rows = frame_rows(run)
    steps = {}
    for fr, istep, isub, _z in rows:
        steps.setdefault(istep, []).append((fr, isub))
    nchecked = 0
    for istep in sorted(steps):
        fr_sub = sorted(steps[istep])
        subs = [s for _, s in fr_sub]
        if subs != list(range(1, len(subs) + 1)):
            print(f"  step {istep}: frames not one per substep ({subs}), skipped")
            continue
        acc = None
        last_g = last_avg = None
        for k, (fr, _isub) in enumerate(fr_sub, start=1):
            g = frame(run, fr, "gammaHI").astype(np.float64)
            avg = frame(run, fr, "gammaHI_avg")
            acc = g if acc is None else acc + g
            want = (acc / k).astype(np.float32)
            check(np.allclose(avg, want, rtol=2e-6, atol=0.0) and avg.shape == want.shape,
                  f"step {istep} substep {k}: gammaHI_avg frame is not the mean of the gammaHI frames "
                  f"(max rel diff {np.max(np.abs(avg - want) / np.maximum(np.abs(want), 1e-30)):.3g})")
            last_g, last_avg = frame(run, fr, "gammaHI"), avg
        d = os.path.join(run, "step_%04d" % istep)
        if os.path.isdir(d):
            for name, ref in (("gammaHI_avg", last_avg), ("gammaHI", last_g)):
                p = os.path.join(d, name + ".bin")
                if not os.path.exists(p):
                    check(False, f"{p} missing")
                    continue
                cube = np.fromfile(p, dtype=np.float32).reshape(N, N, N)
                check(np.array_equal(cube_slabs(cube, thick, axes, start), ref),
                      f"step {istep}: {name}.bin differs from the last frame's slabs")
        else:
            check(False, f"{d} missing (fields_every = 1)")
        if len(subs) >= 2:
            nchecked += 1
    check(nchecked >= 2, f"only {nchecked} hydro step(s) with >= 2 RT substeps were checked")
    print(f"  {nchecked} hydro steps with >= 2 substeps, {len(rows)} frames")


def files_of(run):
    out = {}
    for p in glob.glob(os.path.join(run, "**", "*"), recursive=True):
        if os.path.isfile(p) and os.path.basename(p) not in SKIP:
            out[os.path.relpath(p, run)] = p
    return out


def same_bytes(a, b):
    with open(a, "rb") as fa, open(b, "rb") as fb:
        return fa.read() == fb.read()


def main():
    exe, prob, work = sys.argv[1], sys.argv[2], sys.argv[3]
    sets = sys.argv[4:]
    common = sets + ["--set", "fields_every=1", "--set", "movie_dz=0.5", "--set",
                     "movie_every_substep=1", "--set", "movie_slab=1"]
    with_avg = common + ["--set", "fields=fHI,gammaHI,gammaHI_avg", "--set",
                         "movie_fields=gammaHI,gammaHI_avg"]
    without = common + ["--set", "fields=fHI,gammaHI", "--set", "movie_fields=gammaHI"]
    A, B, C = (os.path.join(work, x) for x in ("A", "B", "C"))
    print("run A (1 thread, gammaHI_avg)")
    run(exe, prob, A, with_avg, 1)
    print("run B (2 threads)")
    run(exe, prob, B, with_avg, 2)
    print("run C (without gammaHI_avg)")
    run(exe, prob, C, without, 1)

    print("the average against the per-substep frames")
    check_average(A)

    print("B against A (thread independence)")
    fa, fb = files_of(A), files_of(B)
    check(set(fa) == set(fb), f"file sets differ: {sorted(set(fa) ^ set(fb))[:10]}")
    for rel in sorted(set(fa) & set(fb)):
        check(same_bytes(fa[rel], fb[rel]), f"{rel} differs between 1 and 2 threads")

    print("C against A (output only)")
    fc = files_of(C)
    extra = sorted(set(fa) - set(fc))
    check(all(os.path.basename(r).startswith("gammaHI_avg") or "_gammaHI_avg." in r for r in extra),
          f"A has files besides gammaHI_avg that C lacks: {[r for r in extra if 'gammaHI_avg' not in r][:10]}")
    check(len(extra) > 0, "A wrote no gammaHI_avg files")
    check(not (set(fc) - set(fa)), f"C has files A lacks: {sorted(set(fc) - set(fa))[:10]}")
    for rel in sorted(set(fa) & set(fc)):
        check(same_bytes(fa[rel], fc[rel]), f"{rel} differs with and without gammaHI_avg")

    if bad:
        print(f"{len(bad)} check(s) failed")
        sys.exit(1)
    print("gamma_avg_check: ok")


if __name__ == "__main__":
    main()
