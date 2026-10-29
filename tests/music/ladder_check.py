#!/usr/bin/env python3
"""The ladder property (docs/zoom_turbulence.md sections 3 and 4.1): a deeper
rung of the zoom, generated with tools/ics/music_ladder.py from the same seed,
base level and region, keeps every coarser level's modes. Needs MUSIC, so it
runs on klone.

usage: ladder_check.py WORKDIR MUSIC [SRC] [--base 5] [--box 4]

Three MUSIC runs, converted by tools/ics/music_convert.py:
  U   the base level alone (uniform), seed[base] = seed;
  A   base .. base + 1, the region given at level base + 2 (--ladder-max):
      A's finest grid is B's level base + 1;
  B   base .. base + 2, the same ladder.
Checks:
 1. White noise (wnoise_LLLL.bin): level base identical in U, A and B; A's
    level base + 1 identical to the overlapping part of B's (exact: MUSIC
    draws each level's noise per cube from (seed, cube) and constrains it
    to its parent).
 2. The converted fields: A's and B's level 1 (base + 1) gas density on A's
    grid averaged onto level-0 cells (the coarser modes) agree, and the
    level-0 gas density of A and B outside the refined region matches U's;
    the level-1 particles common to A and B have the same displacement.
    Tolerances from MUSIC's own solver accuracy (multigrid, 1e-5) and the
    isolated convolutions of a padded patch; the numbers are printed.
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np


def run(cmd, log, cwd=None, env=None):
    with open(log, "w") as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=cwd, env=env)
    if r.returncode != 0:
        print(open(log).read()[-3000:])
        sys.exit(f"failed: {' '.join(cmd)}")


def noise(path):
    with open(path, "rb") as f:
        dims = np.fromfile(f, "<u4", 3)
        return np.fromfile(f, "<f8").reshape(tuple(dims))


def levels(path):
    out, cur = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith("[level"):
            cur = int(line.split()[1].rstrip("]"))
            out[cur] = {}
        elif cur is not None and "=" in line:
            k, v = (x.strip() for x in line.split("=", 1))
            out[cur][k] = v
    return out


def ivec(s):
    return np.array([int(x) for x in s.split()])


def gas_rho(d, name, size):
    return np.fromfile(os.path.join(d, name), "<f4").reshape((5,) + tuple(size), order="F")[0].astype(float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("music")
    ap.add_argument("src", nargs="?", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    ap.add_argument("--base", type=int, default=5)
    ap.add_argument("--box", type=float, default=4.0)
    ap.add_argument("--center", type=float, nargs=3, default=(0.3, 0.6, 0.5))
    ap.add_argument("--extent", type=float, default=0.15)
    args = ap.parse_args()
    src = os.path.abspath(args.src)
    sys.path.insert(0, os.path.join(src, "tools", "ics"))
    import music_convert as mc

    work = os.path.abspath(args.work)
    if os.path.isdir(work):
        shutil.rmtree(work)
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    lad = os.path.join(src, "tools", "ics", "music_ladder.py")
    b = args.base
    common = ["--box", str(args.box), "--base-level", str(b), "--transfer", "eisenstein"]
    rungs = {"U": ["--uniform"], "A": ["--levelmax", str(b + 1), "--ladder-max", str(b + 2)],
             "B": ["--levelmax", str(b + 2), "--ladder-max", str(b + 2)]}
    for name, extra in rungs.items():
        d = os.path.join(work, name)
        os.makedirs(d)
        if name != "U":
            extra = extra + ["--center"] + [str(c) for c in args.center] + ["--extent", str(args.extent)]
        run([sys.executable, lad, os.path.join(d, "r.conf")] + common + extra, os.path.join(d, "ladder.log"))
        run([args.music, "r.conf"], os.path.join(d, "music.out"), cwd=d, env=env)
        mc.main([os.path.join(d, "r.hdf5"), os.path.join(d, "ics"), "--conf", os.path.join(d, "r.conf")])
    fails = []

    # 1. white noise
    def shift(x):
        for line in open(os.path.join(work, x, "ics", "ics_info.txt")):
            if line.startswith("shift ="):
                return tuple(int(v) for v in line.split("=")[1].split())

    # the base noise is stored in MUSIC's shifted frame; the refined levels'
    # files are in the box's own coordinates
    nU, nA, nB = (np.roll(noise(os.path.join(work, x, f"wnoise_{b:04d}.bin")),
                          tuple(-v for v in shift(x)), axis=(0, 1, 2)) for x in "UAB")
    e0 = max(np.abs(nU - nA).max(), np.abs(nU - nB).max())
    print(f"level {b} white noise, U against A and B (shifts {shift('A')}, {shift('B')} undone): "
          f"max difference {e0:.1e}")
    if not e0 < 1e-12:
        fails.append(f"level {b} noise differs: {e0:.1e}")
    a1, b1 = (noise(os.path.join(work, x, f"wnoise_{b + 1:04d}.bin")) for x in "AB")
    e1 = float(np.abs(a1 - b1).max()) if a1.shape == b1.shape else float("inf")
    print(f"level {b + 1} white noise: A's {a1.shape} against B's {b1.shape}: max difference {e1:.1e}")
    if not e1 < 1e-12:
        fails.append(f"A's level {b + 1} noise differs from B's: {e1:.1e}")

    # 2. converted fields
    LA, LB = (levels(os.path.join(work, x, "ics", "ics_levels.txt")) for x in "AB")
    oA, sA = ivec(LA[1]["offset"]), ivec(LA[1]["size"])
    oB, sB = ivec(LB[1]["offset"]), ivec(LB[1]["size"])
    rA = gas_rho(os.path.join(work, "A", "ics"), "u_gas.init.L1", sA)
    rB = gas_rho(os.path.join(work, "B", "ics"), "u_gas.init.L1", sB)
    d = oA - oB
    if np.any(d < 0) or np.any(d + sA > sB):
        sys.exit(f"A's level 1 ({oA}, {sA}) is not inside B's ({oB}, {sB})")
    rBa = rB[d[0]:d[0] + sA[0], d[1]:d[1] + sA[1], d[2]:d[2] + sA[2]]
    fb = rA.mean() / (1 + 0)  # ~ f_b
    dA, dB = rA / fb - 1, rBa / fb - 1
    sig = dA.std()
    blk = lambda x: x.reshape(x.shape[0] // 2, 2, x.shape[1] // 2, 2, x.shape[2] // 2, 2).mean(axis=(1, 3, 5))  # noqa
    e_coarse = float(np.sqrt(np.mean((blk(dA) - blk(dB)) ** 2)) / blk(dA).std())
    e_fine = float(np.sqrt(np.mean((dA - dB) ** 2)) / sig)
    print(f"level 1 gas delta on A's grid, A against B: rms difference {e_fine:.3e} of the rms "
          f"(all modes), {e_coarse:.3e} averaged onto level-0 cells (the coarser modes)")
    cov = np.zeros(tuple(sA), dtype=bool)
    no, ns = ivec(LB[1]["next_offset"]) - d, ivec(LB[1]["next_size"])
    cov[no[0]:no[0] + ns[0], no[1]:no[1] + ns[1], no[2]:no[2] + ns[2]] = True
    for lab, m in (("under B's level 2", cov), ("outside it", ~cov)):
        print(f"   {lab}: rms difference {np.sqrt(np.mean((dA[m] - dB[m]) ** 2)) / sig:.3e} of the rms")
    if not e_coarse < 2e-2:
        fails.append(f"level-1 coarse modes differ by {e_coarse:.3e} >= 2e-2")

    N = 1 << b
    rU = gas_rho(os.path.join(work, "U", "ics"), "u_gas.init", (N, N, N))
    outside = np.ones((N, N, N), dtype=bool)
    o0, s0 = ivec(LA[0]["next_offset"]), ivec(LA[0]["next_size"])
    oB0, sB0 = ivec(LB[0]["next_offset"]), ivec(LB[0]["next_size"])
    for o, s in ((o0, s0), (oB0, sB0)):
        idx = [np.mod(np.arange(o[c] - 2, o[c] + s[c] + 2), N) for c in range(3)]
        outside[np.ix_(*idx)] = False
    for name in "AB":
        r0 = gas_rho(os.path.join(work, name, "ics"), "u_gas.init", (N, N, N))
        e = float(np.sqrt(np.mean((r0[outside] / rU[outside] - 1) ** 2)) / (rU[outside] / fb - 1).std())
        print(f"level 0 gas delta of {name} against U outside the refined region (+2 cells): rms difference "
              f"{e:.3e} of the rms")
        if not e < 2e-2:
            fails.append(f"{name}'s base differs from U's outside the region: {e:.3e}")

    # particles of level 1 common to A and B
    xa = np.fromfile(os.path.join(work, "A", "ics", "xv_dm.L1"), "<f4").reshape(-1, 6).astype(float)
    xb = np.fromfile(os.path.join(work, "B", "ics", "xv_dm.L1"), "<f4").reshape(-1, 6).astype(float)

    def lattice(L, xv):
        o, s = ivec(L[1]["offset"]), ivec(L[1]["size"])
        cov = np.zeros(tuple(s), dtype=bool)
        if "next_offset" in L[1]:
            no, ns = ivec(L[1]["next_offset"]), ivec(L[1]["next_size"])
            cov[no[0]:no[0] + ns[0], no[1]:no[1] + ns[1], no[2]:no[2] + ns[2]] = True
        idx = np.argwhere(~cov.transpose(2, 1, 0))[:, ::-1]      # i fastest
        key = [(o[c] + idx[:, c]) for c in range(3)]
        return {k: n for n, k in enumerate(zip(*key))}, xv

    ka, _ = lattice(LA, xa)
    kb, _ = lattice(LB, xb)
    common_keys = [k for k in ka if k in kb]
    ia = np.array([ka[k] for k in common_keys])
    ib = np.array([kb[k] for k in common_keys])
    dx = xa[ia, :3] - xb[ib, :3]
    dx -= N * np.round(dx / N)
    q = (np.array(common_keys) + 0.5) / 2.0
    psi = xa[ia, :3] - q
    psi -= N * np.round(psi / N)
    e_p = float(np.sqrt(np.mean(dx ** 2)) / np.sqrt(np.mean(psi ** 2)))
    print(f"level-1 particles common to A and B ({len(common_keys)}): displacement rms difference "
          f"{e_p:.3e} of the rms displacement")
    if not e_p < 5e-2:
        fails.append(f"level-1 displacements differ by {e_p:.3e}")
    if fails:
        print("FAIL:\n  " + "\n  ".join(fails))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
