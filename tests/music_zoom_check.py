#!/usr/bin/env python3
"""The MUSIC zoom converter on a small committed two-level zoom
(tools/ics/music_convert.py, tools/ics/music_ladder.py; docs/zoom_turbulence.md
section 4.1; tools/ics/music/README.md). Numpy only: the MUSIC output is
tests/data/music/zoom2.npz (music_convert.py --pack --mantissa-bits 14 of a
MUSIC 2.0a run: 4 Mpc/h, base 32^3 (level 5), levels 6 and 7 over a 0.5
Mpc/h region, padding 4, Eisenstein-Hu, the original cosmology, z = 100).

usage: music_zoom_check.py WORKDIR SRC [RADHYDRO]

1. The generator: music_ladder.py with the arguments recorded in the packed
   .conf writes that .conf again, and its grid of every level is the one
   MUSIC built (the converter's unshifted offsets and sizes).
2. The converter: refuses the Planck target cosmology for this run; with
   the original one writes the files; the particle classes tile the box
   once (the summed class masses are the box's DM mass); every particle sits
   at its lattice cell plus a displacement of under 1.5 base cells; v/psi is
   the code's (dD/dtau)/D in every class; on the intermediate level the gas
   velocity under the finer level is the mean of the finer level's (MUSIC's
   restriction), and the gas density's mean there is the finer level's.
3. With RADHYDRO: the code runs from the base files (xv_dm.init, u_gas.init)
   to z = 50.
"""
import os
import re
import shlex
import shutil
import subprocess
import sys

import numpy as np


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


def main():
    work, src = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    exe = sys.argv[3] if len(sys.argv) > 3 else None
    sys.path.insert(0, os.path.join(src, "tools", "ics"))
    import music_convert as mc
    import music_ladder as ml
    import pmcosmo as pc

    data = os.path.join(src, "tests", "data", "music", "zoom2.npz")
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    fails = []
    md = mc.MusicData(data)
    conf = md.conf_text

    # 1. the generator reproduces the conf and MUSIC's grids
    cmd = re.match(r"# written by tools/ics/music_ladder.py: (.*)", conf.splitlines()[0]).group(1)
    argv = shlex.split(cmd)
    argv[0] = os.path.join(work, "zoom2.conf")
    ml.main(argv)
    again = open(argv[0]).read()
    if again.splitlines()[1:] != conf.splitlines()[1:]:
        fails.append("music_ladder.py no longer writes the packed .conf")
    chain = {int(m.group(1)): (np.array(eval(m.group(2))), np.array(eval(m.group(3))))
             for m in re.finditer(r"#   level (\d+): lo (\[[^]]*\]) hi (\[[^]]*\])", conf)}
    # the region from particle ids, across the periodic boundary: the cells
    # i = 30, 31, 0, 1 (x), j = 4..6, k = 10 of a 32^3 lattice
    n = 32
    cells = [(i, j, 10) for i in (30, 31, 0, 1) for j in (4, 5, 6)]
    ids = [1 + i + n * j + n * n * k for i, j, k in cells]
    np.savetxt(os.path.join(work, "ids.txt"), ids, fmt="%d")

    class A:
        pass
    ar = A()
    ar.ids, ar.ndm, ar.pad, ar.box = os.path.join(work, "ids.txt"), n, 0.0, 4.0
    lo, w = ml.region(ar)
    want_lo, want_w = [-2 / n, 4 / n, 10 / n], [4 / n, 3 / n, 1 / n]
    print(f"region from ids: lower corner {np.round(np.array(lo) * n, 6)} cells, extent "
          f"{np.round(np.array(w) * n, 6)} cells")
    if not (np.allclose(np.mod(lo, 1), np.mod(want_lo, 1)) and np.allclose(w, want_w)):
        fails.append(f"region from ids: {lo}, {w} (want {want_lo}, {want_w})")

    # 2. the converter
    out = os.path.join(work, "ics")
    r = subprocess.run([sys.executable, os.path.join(src, "tools", "ics", "music_convert.py"), data,
                        os.path.join(work, "wrong"), "--cosmology", "planck_c26"], capture_output=True, text=True)
    if r.returncode == 0 or "differs from the target" not in r.stdout + r.stderr:
        fails.append("the converter accepted a cosmology mismatch")
    else:
        print("cosmology mismatch refused:", (r.stdout + r.stderr).strip().splitlines()[0])
    info, _ = mc.convert(md, out, pc.ORIGINAL)
    L = levels(os.path.join(out, "ics_levels.txt"))
    N = int(info["N"])
    nlev = len(L)
    for lv in range(1, nlev):
        lev = md.levelmin + lv
        lo, hi = chain[lev]
        n = 1 << lev
        if not (np.array_equal(ivec(L[lv]["offset"]), np.mod(lo, n)) and np.array_equal(ivec(L[lv]["size"]), hi - lo)):
            fails.append(f"level {lv}: MUSIC's grid {L[lv]['offset']} / {L[lv]['size']} is not the "
                         f"generator's {lo} / {hi - lo}")
    print(f"generator: the .conf and the grids of levels 1..{nlev - 1} reproduced")

    # tiling: base particles outside level 1 plus every finer class
    nb = ivec(L[0]["next_size"]).prod()
    mass = (N ** 3 - nb) * float(L[0]["mass"]) + sum(int(L[lv]["npart"]) * float(L[lv]["mass"])
                                                     for lv in range(1, nlev))
    total = (1 - pc.FB) * N ** 3
    print(f"class masses sum to {mass / total:.15f} of the box's DM mass")
    if abs(mass / total - 1) > 1e-12:
        fails.append(f"the classes do not tile the box: {mass / total}")

    g = pc.Growth(1.0)
    a = 1.0 / (1.0 + md.zstart)
    vfac = float(g.dtau(a) / g(a))
    for lv in range(nlev):
        o, s = ivec(L[lv]["offset"]), ivec(L[lv]["size"])
        sc = 1 << lv
        cov = np.zeros(tuple(s), dtype=bool)
        if "next_offset" in L[lv] and lv > 0:
            no, ns = ivec(L[lv]["next_offset"]), ivec(L[lv]["next_size"])
            cov[no[0]:no[0] + ns[0], no[1]:no[1] + ns[1], no[2]:no[2] + ns[2]] = True
        idx = np.argwhere(~cov.transpose(2, 1, 0))[:, ::-1]        # the files' order: i fastest
        q = (o + idx + 0.5) / sc
        xv = np.fromfile(os.path.join(out, L[lv]["dm_file"]), "<f4").reshape(-1, 6).astype(float)
        if len(xv) != int(L[lv]["npart"]) or len(xv) != len(q):
            fails.append(f"level {lv}: {len(xv)} particles, header {L[lv]['npart']}, lattice {len(q)}")
            continue
        psi = xv[:, :3] - q
        psi -= N * np.round(psi / N)
        ratio = np.sum(xv[:, 3:] * psi) / np.sum(psi * psi) / vfac
        print(f"level {lv}: {len(xv)} particles, max |psi| {np.abs(psi).max():.3f} base cells, "
              f"v/psi = {ratio:.6f} x (dD/dtau)/D")
        if np.abs(psi).max() > 1.5 or abs(ratio - 1) > 1e-3:
            fails.append(f"level {lv}: displacement {np.abs(psi).max():.3f} or v/psi {ratio:.6f}")

    # gas: the intermediate level under the finest one
    for lv in range(1, nlev - 1):
        s, s1 = ivec(L[lv]["size"]), ivec(L[lv + 1]["size"])
        u = np.fromfile(os.path.join(out, L[lv]["gas_file"]), "<f4").reshape((5,) + tuple(s), order="F").astype(float)
        u1 = np.fromfile(os.path.join(out, L[lv + 1]["gas_file"]), "<f4").reshape((5,) + tuple(s1),
                                                                                  order="F").astype(float)
        no, ns = ivec(L[lv]["next_offset"]), ivec(L[lv]["next_size"])
        sub = (slice(no[0], no[0] + ns[0]), slice(no[1], no[1] + ns[1]), slice(no[2], no[2] + ns[2]))
        blk = lambda x: x.reshape(ns[0], 2, ns[1], 2, ns[2], 2).mean(axis=(1, 3, 5))  # noqa: E731
        ev = max(float(np.abs(blk(u1[1 + c] / u1[0]) - u[1 + c][sub] / u[0][sub]).max()) for c in range(3))
        vs = max(float(np.abs(u1[1 + c] / u1[0]).max()) for c in range(3))
        er = abs(blk(u1[0]).mean() / u[0][sub].mean() - 1)
        print(f"level {lv} under level {lv + 1}: gas velocity against the children's mean {ev / vs:.2e} "
              f"(of the max), density mean {er:.2e}")
        if ev / vs > 1e-4 or er > 1e-5:
            fails.append(f"level {lv}: restriction {ev / vs:.2e}, {er:.2e}")

    # 3. the code runs from the base files
    if exe:
        p = os.path.join(work, "problem")
        os.makedirs(p)
        shutil.copytree(out, os.path.join(p, "ics"))
        txt = open(os.path.join(src, "tests", "problems_full", "smoke", "params.txt")).read()
        for k, v in (("N", N), ("Ndm", N), ("box", md.boxlength), ("z_init", md.zstart), ("zdata1", 50)):
            txt = re.sub(rf"(?m)^{k} = .*$", f"{k} = {v}", txt)
        open(os.path.join(p, "params.txt"), "w").write(txt)
        o = os.path.join(work, "run")
        with open(o + ".log", "w") as f:
            rr = subprocess.run([exe, p, o, "--quiet"], stdout=f, stderr=subprocess.STDOUT)
        if rr.returncode != 0:
            print(open(o + ".log").read()[-2000:])
            fails.append("radhydro failed from the converted files")
        else:
            gas = sorted(x for x in os.listdir(os.path.join(o, "Gas")) if x.startswith("u_gas"))[-1]
            rho = np.fromfile(os.path.join(o, "Gas", gas), "<f4")[: 5 * N ** 3].reshape(5, N, N, N, order="F")[0]
            print(f"radhydro ran to {gas}: mean gas density {rho.mean():.6f} (f_b = {pc.FB:.6f})")
    if fails:
        print("FAIL:\n  " + "\n  ".join(fails))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
