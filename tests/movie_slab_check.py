#!/usr/bin/env python3
"""Slab-only movie frames of apps/radhydro (movie_slab, docs/movies.md).

usage: movie_slab_check.py RADHYDRO PROBLEM WORKDIR [--set key=value ...]

Runs PROBLEM three times with movie frames (fHI, temp, rhog, gamma, gammaHI):
  full   movie_slab = 0 (default): full cubes, frame_NNNNN_<field>.bin;
  slab   movie_slab = 4 (all three axes);
  slab3  movie_slab = 3, movie_slab_axes = zx (odd thickness, axis order);
  slab1  movie_slab = 1 (the make_production.py default).
Checks that the slab runs write no cubes, that every slab file holds exactly
the cells of the full cube that slab.txt describes (gamma: 3 components per
cell), and that frames.txt (whole-box x_HII means) is the same in all runs
apart from its header, and that its Gamma_HI columns (volume mean, median
and mean over x_HII > 0.5, volume fraction above 1e-14 s^-1) and the
gammaHI frames match the full gamma cubes. Then runs the tools on the slab frames and compares
with the full cubes: reion.py (identical output), make_movie.py stills of
the default slice and of a projection over the slab (identical PNGs), the
error for a slice outside the slab, stills.py, and the planar front.
Exit 0 on success.
"""
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MOVIES = os.path.join(ROOT, "tools", "movies")
sys.path.insert(0, MOVIES)
import movie_frames as mf  # noqa: E402

bad = []


def check(cond, msg):
    if not cond:
        bad.append(msg)
        print("  FAIL", msg)


def run(cmd, log, ok_codes=(0,)):
    with open(log, "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode not in ok_codes:
        print(open(log).read()[-3000:])
        sys.exit(f"failed ({r.returncode}): {' '.join(cmd)}")
    return r.returncode


def frame_lines(d):
    return [x for x in open(os.path.join(d, "movie", "frames.txt")).read().splitlines()
            if not x.startswith("#")]


def main():
    exe, prob, work = sys.argv[1], sys.argv[2], sys.argv[3]
    sets = sys.argv[4:]
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    common = sets + ["--set", "movie_fields=fHI,temp,rhog,gamma,gammaHI"]
    cases = {"full": [], "slab": ["--set", "movie_slab=4"],
             "slab3": ["--set", "movie_slab=3", "--set", "movie_slab_axes=zx"],
             "slab1": ["--set", "movie_slab=1"]}
    dirs = {}
    for name, extra in cases.items():
        d = os.path.join(work, name, "run")
        os.makedirs(os.path.dirname(d))
        # params.txt next to the run directory, for the tools (box)
        shutil.copy(os.path.join(prob, "params.txt"), os.path.join(work, name, "params.txt"))
        run([exe, prob, d, "--quiet"] + common + extra, d + ".log")
        dirs[name] = d
    full = dirs["full"]
    N = int([x.split("=")[1] for x in open(os.path.join(prob, "params.txt"))
             if x.split("=")[0].strip() == "N"][0])
    lines = frame_lines(full)
    print(f"N = {N}, {len(lines)} movie frames")
    check(len(lines) > 2, "too few movie frames")
    check(not os.path.exists(os.path.join(full, "movie", "slab.txt")), "full run wrote slab.txt")

    # Gamma_HI columns of frames.txt and the gammaHI frames, from the full cubes
    ng = 0
    for ln in lines:
        p = ln.split()
        check(len(p) == 10, f"frames.txt line has {len(p)} columns, not 10")
        if len(p) != 10:
            break
        fr = int(p[0])
        g3 = np.fromfile(os.path.join(full, "movie", "frame_%05d_gamma.bin" % fr), "<f4").reshape(3, -1, order="F")
        g1 = np.fromfile(os.path.join(full, "movie", "frame_%05d_gammaHI.bin" % fr), "<f4")
        check(np.array_equal(g1, g3[0]), f"frame {fr}: gammaHI is not gamma(1)")
        x = 1 - np.fromfile(os.path.join(full, "movie", "frame_%05d_fHI.bin" % fr), "<f4").astype(np.float64)
        gi = g1[x > 0.5]
        want = [g1.astype(np.float64).sum() / g1.size,
                float(np.median(gi)) if gi.size else 0.0,
                gi.astype(np.float64).mean() if gi.size else 0.0,
                np.count_nonzero(g1.astype(np.float64) > 1e-14) / g1.size]
        for i, w in enumerate(want):
            v = float(p[6 + i])
            if not abs(v - w) <= 2e-6 * abs(w) + 1e-300:
                check(False, f"frame {fr}: column {7 + i} = {v:.8e}, from the cubes {w:.8e}")
        ng += gi.size > 0
    print(f"  Gamma_HI columns checked on {len(lines)} frames ({ng} with ionised cells)")

    for name, t, axes in (("slab", 4, "xyz"), ("slab3", 3, "zx"), ("slab1", 1, "xyz")):
        d = dirs[name]
        print(f"{name}: movie_slab = {t}, axes {axes}")
        info = mf.run_slab_info(d)
        check(info == dict(N=N, thick=t, axes=axes, start=N // 2 - t // 2),
              f"{name}: slab.txt {info}")
        check(frame_lines(d) == lines, f"{name}: frames.txt differs from the full-cube run")
        mv = os.listdir(os.path.join(d, "movie"))
        check(not any(f.endswith(".bin") for f in mv), f"{name}: wrote full cubes")
        n = 0
        for ln in lines:
            fr = int(ln.split()[0])
            for fld, nc in (("fHI", 1), ("temp", 1), ("rhog", 1), ("gamma", 3), ("gammaHI", 1)):
                cube = np.fromfile(os.path.join(full, "movie", "frame_%05d_%s.bin" % (fr, fld)),
                                   "<f4").reshape((nc, N, N, N), order="F")
                want = []
                for a in axes:
                    k = "xyz".index(a) + 1
                    sl = [slice(None)] * 4
                    sl[k] = slice(info["start"], info["start"] + t)
                    want.append(cube[tuple(sl)].ravel(order="F"))
                got = open(os.path.join(d, "movie", "frame_%05d_%s.slab" % (fr, fld)), "rb").read()
                n += 1
                if got != np.concatenate(want).astype("<f4").tobytes():
                    check(False, f"{name}: frame {fr} {fld} is not the slab of the full cube")
        print(f"  {n} slab files compared with the full cubes")

    # reion.py: frames.txt only, the same with slab frames
    outs = {}
    for name in ("full", "slab"):
        o = os.path.join(work, name + "_reion.txt")
        run([sys.executable, os.path.join(ROOT, "tools", "analysis", "reion.py"), dirs[name],
             "--out", o], o + ".log")
        outs[name] = open(o).read()
    check(outs["full"] == outs["slab"], "reion.py output differs between full and slab frames")

    try:
        import matplotlib  # noqa: F401
    except ImportError:
        print("matplotlib not installed: movie tools skipped")
        return finish()

    mm_ = os.path.join(MOVIES, "make_movie.py")

    def movie(name, out, *extra, ok=(0,)):
        return run([sys.executable, mm_, dirs[name], "--no-mp4", "--stills", "0,1,2",
                    "--size", "640", "360", "--dpi", "40", "--still-dir", os.path.join(work, out)]
                   + list(extra), os.path.join(work, out + ".log"), ok)

    def same_pngs(a, b, what):
        fa = sorted(os.listdir(os.path.join(work, a)))
        fb = sorted(os.listdir(os.path.join(work, b)))
        check(fa == fb and len(fa) == 3, f"{what}: stills {fa} vs {fb}")
        for f in fa:
            if f in fb:
                check(open(os.path.join(work, a, f), "rb").read()
                      == open(os.path.join(work, b, f), "rb").read(), f"{what}: {f} differs")

    # default slice (z, cell N//2): the same image from the cube and the slab
    movie("full", "mv_full")
    movie("slab", "mv_slab")
    same_pngs("mv_full", "mv_slab", "default slice")
    # projection over the 4 stored cells: default with slabs, --pos with cubes
    cell = float([x.split("=")[1] for x in open(os.path.join(prob, "params.txt"))
                  if x.split("=")[0].strip() == "box"][0]) / N
    s0 = N // 2 - 2
    movie("full", "mv_full_x4", "--axis", "x", "--thick", "4", "--pos", repr((s0 + 0.25) * cell))
    movie("slab", "mv_slab_x4", "--axis", "x", "--thick", "4")
    same_pngs("mv_full_x4", "mv_slab_x4", "x projection over the slab")
    # outside the stored slab / an axis not stored: a clear error
    rc = movie("slab", "mv_slab_out", "--pos", repr(0.5 * cell), ok=(1,))
    check("choose --pos/--thick inside the slab" in open(os.path.join(work, "mv_slab_out.log")).read(),
          f"slice outside the slab: exit {rc}, no message")
    movie("slab3", "mv_slab3_y", "--axis", "y", ok=(1,))
    check("movie_slab_axes" in open(os.path.join(work, "mv_slab3_y.log")).read(),
          "axis not stored: no message")
    movie("slab3", "mv_slab3_z", "--axis", "z")
    # planar style: front from the slab (runs, prints the note)
    movie("slab", "mv_slab_planar", "--style", "planar")
    check("slab normal to" in open(os.path.join(work, "mv_slab_planar.log")).read(),
          "planar style on slab frames: no note")
    # stills.py on 1-cell slabs (the make_production.py default) with its
    # default --slab 2: falls back to the stored cell
    out = os.path.join(work, "stills_slab1")
    run([sys.executable, os.path.join(MOVIES, "stills.py"), dirs["slab1"], out, "--z", "24,22"],
        out + ".log")
    check(os.path.exists(os.path.join(out, "stills_grid.png")), "stills.py slab1: no grid")
    check("using 1 cell" in open(out + ".log").read(), "stills.py slab1: no fallback note")
    # stills.py: slab frames, default position (most ionised within the slab)
    out = os.path.join(work, "stills_slab")
    run([sys.executable, os.path.join(MOVIES, "stills.py"), dirs["slab"], out, "--z", "24,22"],
        out + ".log")
    check(os.path.exists(os.path.join(out, "stills_grid.png")), "stills.py slab: no grid")
    # the same slab position in both: identical images
    for name in ("full", "slab"):
        out = os.path.join(work, "stills_pos_" + name)
        run([sys.executable, os.path.join(MOVIES, "stills.py"), dirs[name], out, "--z", "24,22",
             "--pos", str(s0 + 1)], out + ".log")
    a, b = (os.path.join(work, "stills_pos_" + n, "stills_grid.png") for n in ("full", "slab"))
    check(open(a, "rb").read() == open(b, "rb").read(), "stills.py --pos: full and slab differ")
    return finish()


def finish():
    print("PASS" if not bad else f"FAILED ({len(bad)})")
    sys.exit(0 if not bad else 1)


if __name__ == "__main__":
    main()
