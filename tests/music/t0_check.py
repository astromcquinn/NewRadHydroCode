#!/usr/bin/env python3
"""T0 (docs/zoom_turbulence.md section 5): a uniform MUSIC box converted by
tools/ics/music_convert.py against tools/ics/make_ics.py at the same P(k) and
the same white noise. Needs MUSIC (tools/hyak/setup_music.sh), so it runs on
klone, not in the fast suite.

usage: t0_check.py WORKDIR MUSIC [RADHYDRO SRC] [--n 32] [--box 4] [--z 300]

1. P(k): make_ics.py's own Eisenstein-Hu no-wiggle transfer function
   (pmcosmo, the original cosmology) is written as a CAMB-format table
   (transfer = camb_file; the velocity columns are f(a_init) times the
   density's, CAMB's convention, so the converter's velocity-TF path is
   used), MUSIC runs a uniform N^3 box from it, the converter writes the
   code's files.
2. The reference: make_ics.py --kind zeldovich --Ndm N --Ngas N with MUSIC's
   white noise (wnoise_LLLL.bin) in place of its own draws.
3. Mode by mode (FFTs of the gas overdensity and of the particle
   displacements; the Nyquist planes, which make_ics.py zeroes and MUSIC
   keeps, excluded): the overall amplitude ratio (the sigma_8 integrals of
   MUSIC and pmcosmo, and the growth factors), and the scatter of the
   per-mode ratio about it, which must be float round-off; the
   velocity-to-displacement ratio of the particles and of the gas in both
   against (dD/dtau)/D of the code.
4. With RADHYDRO and SRC: the code runs from both sets of files (the
   tests/problems_full/smoke problem at N^3, gravity and hydro) and the
   final gas density and particles are compared.
"""
import argparse
import os
import re
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


def nyquist_mask(n):
    m = np.fft.fftfreq(n, 1.0 / n)
    mx, my, mz = np.meshgrid(m, m, m, indexing="ij")
    nyq = (np.abs(mx) == n // 2) | (np.abs(my) == n // 2) | (np.abs(mz) == n // 2)
    nyq[0, 0, 0] = True
    return ~nyq


def half_cell_phase(n):
    """make_ics.py puts the white noise at the integer points and evaluates
    the particles and the gas at n + 1/2 (eval_on_lattice, offset 0.5); MUSIC
    has the noise at the cell centres. The same noise therefore gives
    make_ics's modes = MUSIC's times exp(i pi (mx + my + mz)/n) (Nyquist
    planes excluded)."""
    m = np.fft.fftfreq(n, 1.0 / n)
    ph = np.exp(1j * np.pi * m / n)
    return ph[:, None, None] * ph[None, :, None] * ph[None, None, :]


def compare_modes(a, b, keep, label):
    """Amplitude ratio <b a*>/<|a|^2> and the rms scatter of b/a about it,
    weighted by |a|^2, over the kept modes; a from make_ics.py, b from the
    converter (MUSIC's noise at the cell centres, shifted half a cell)."""
    fa = np.fft.fftn(a)[keep]
    fb = (np.fft.fftn(b) * half_cell_phase(a.shape[0]))[keep]
    w = np.abs(fa) ** 2
    amp = float(np.real(np.sum(fb * np.conj(fa)) / np.sum(w)))
    scat = float(np.sqrt(np.sum(np.abs(fb - amp * fa) ** 2) / np.sum(w)))
    print(f"{label}: amplitude ratio converter/make_ics = {amp:.9f}, per-mode scatter {scat:.3e}")
    return amp, scat


def read_xv(path, n):
    return np.fromfile(path, "<f4").reshape(-1, 6).astype(np.float64).reshape(n, n, n, 6, order="F")


def read_u(path, n):
    return np.fromfile(path, "<f4").reshape(5, n, n, n, order="F").astype(np.float64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("music")
    ap.add_argument("radhydro", nargs="?")
    ap.add_argument("src", nargs="?", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    ap.add_argument("--n", type=int, default=32)
    ap.add_argument("--box", type=float, default=4.0)
    ap.add_argument("--z", type=float, default=300.0)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--run-z", type=float, default=30.0, help="radhydro: run to this redshift")
    args = ap.parse_args()
    src = os.path.abspath(args.src)
    sys.path.insert(0, os.path.join(src, "tools", "ics"))
    import make_ics
    import music_convert as mc
    import pmcosmo as pc

    work = os.path.abspath(args.work)
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    n, lev = args.n, int(round(np.log2(args.n)))
    a_init = 1.0 / (1.0 + args.z)
    g = pc.Growth(1.0)
    fails = []

    # 1. the table and MUSIC
    k = np.logspace(-5, 4.5, 4000)
    T = pc.transfer_eh_nowiggle(k)
    fT = float(g.f(a_init)) * T
    tab = np.column_stack([k, T, T, T, T, np.zeros_like(k), T, T, T, np.zeros_like(k), fT, fT, np.zeros_like(k)])
    np.savetxt(os.path.join(work, "eh_nowiggle.dat"), tab, fmt="%.12e",
               header="pmcosmo.transfer_eh_nowiggle (original cosmology) in CAMB's 13 columns; "
                      "v columns = f(a_init) T")
    conf = os.path.join(work, "t0.conf")
    run([sys.executable, os.path.join(src, "tools", "ics", "music_ladder.py"), conf, "--box", str(args.box),
         "--base-level", str(lev), "--uniform", "--seed", str(args.seed), "--zstart", str(args.z),
         "--cosmology", "original", "--transfer", "camb_file", "--transfer-file", "eh_nowiggle.dat"],
        os.path.join(work, "ladder.log"))
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    run([args.music, "t0.conf"], os.path.join(work, "music.out"), cwd=work, env=env)
    mc.main([os.path.join(work, "t0.hdf5"), os.path.join(work, "music_ics"), "--conf", conf,
             "--cosmology", "original"])

    # 2. make_ics.py with MUSIC's white noise
    with open(os.path.join(work, f"wnoise_{lev:04d}.bin"), "rb") as f:
        dims = np.fromfile(f, "<u4", 3)
        noise = np.fromfile(f, "<f8").reshape(tuple(dims))

    class NoiseRng:
        def standard_normal(self, shape):
            assert tuple(shape) == noise.shape
            return noise.copy()

    real_rng = np.random.default_rng
    np.random.default_rng = lambda seed=None: NoiseRng()
    try:
        make_ics.main([os.path.join(work, "make_ics"), "--kind", "zeldovich", "--Ngas", str(n), "--Ndm", str(n),
                       "--z-init", str(args.z), "--box", str(args.box)])
    finally:
        np.random.default_rng = real_rng

    # 3. mode by mode
    keep = nyquist_mask(n)
    fb = pc.FB
    vfac = float(g.dtau(a_init) / g(a_init))
    q = (np.arange(n) + 0.5)
    Q = np.meshgrid(q, q, q, indexing="ij")
    res = {}
    for name in ("music_ics", "make_ics"):
        u = read_u(os.path.join(work, name, "u_gas.init"), n)
        xv = read_xv(os.path.join(work, name, "xv_dm.init"), n)
        psi = [xv[..., c] - Q[c] for c in range(3)]
        psi = [p - n * np.round(p / n) for p in psi]
        vel = [xv[..., 3 + c] for c in range(3)]
        vg = [u[1 + c] / u[0] for c in range(3)]
        res[name] = dict(delta=u[0] / fb - 1, psi=psi, vel=vel, vg=vg)
    m, r = res["music_ics"], res["make_ics"]
    amp, scat = compare_modes(r["delta"], m["delta"], keep, "gas overdensity")
    amps, scat_psi = [amp], 0.0
    for c in range(3):
        a_, s_ = compare_modes(r["vel"][c], m["vel"][c], keep, "particle velocity " + "xyz"[c])
        amps.append(a_)
        scat = max(scat, s_)
        a_, s_ = compare_modes(r["vg"][c], m["vg"][c], keep, "gas velocity " + "xyz"[c])
        amps.append(a_)
        scat = max(scat, s_)
        # psi = x - q from float32 positions x ~ n: ~6e-8 n absolute round-off
        a_, s_ = compare_modes(r["psi"][c], m["psi"][c], keep, "particle displacement " + "xyz"[c])
        amps.append(a_)
        scat_psi = max(scat_psi, s_)
    for name, d in res.items():
        num = sum(np.sum(d["vel"][c] * d["psi"][c]) for c in range(3))
        den = sum(np.sum(d["psi"][c] ** 2) for c in range(3))
        ratio = num / den / vfac
        rms = np.sqrt(sum(np.sum((d["vel"][c] - vfac * d["psi"][c]) ** 2) for c in range(3)) /
                      sum(np.sum((vfac * d["psi"][c]) ** 2) for c in range(3)))
        print(f"{name}: particles v/psi = {ratio:.9f} x (dD/dtau)/D, rms deviation {rms:.2e}")
        gg = sum(np.sum(d["vg"][c] * d["psi"][c]) for c in range(3)) / den / vfac
        print(f"{name}: gas v / particle psi (same lattice) = {gg:.9f} x (dD/dtau)/D")
        # positions are float32 (x ~ n): psi carries ~1e-7 n absolute round-off
        if not abs(ratio - 1) < 1e-5 or not rms < 1e-3:
            fails.append(f"{name}: v/psi {ratio:.9f}, rms {rms:.2e}")
    if not scat < 1e-5:
        fails.append(f"per-mode scatter {scat:.3e} >= 1e-5")
    psi_round = 6e-8 * n / np.sqrt(np.mean(r["psi"][0] ** 2))
    if not scat_psi < 4 * psi_round:
        fails.append(f"displacement per-mode scatter {scat_psi:.3e} >= 4 x float32 position round-off {psi_round:.1e}")
    if max(abs(x - amps[0]) for x in amps) > 1e-6 * 4:
        fails.append(f"amplitudes differ between fields: {amps}")
    if not abs(amps[0] - 1) < 1e-4:
        fails.append(f"amplitude ratio {amps[0]:.6f} differs from 1 by > 1e-4")

    # 4. the code from both, and from the converter's files moved half a cell
    # onto make_ics.py's sampling (a spectral shift; the Nyquist planes, which
    # make_ics.py zeroes and a half-cell shift cannot carry, dropped): that
    # run must match make_ics's
    if args.radhydro:
        nyq = os.path.join(work, "music_ics_shifted")
        shutil.copytree(os.path.join(work, "music_ics"), nyq)
        keep_dc = keep.copy()
        keep_dc[0, 0, 0] = True
        ph = half_cell_phase(n) * keep_dc
        drop = lambda f: np.real(np.fft.ifftn(np.fft.fftn(f) * ph))  # noqa: E731
        u = read_u(os.path.join(nyq, "u_gas.init"), n)
        rho = drop(u[0])
        vg = [drop(u[1 + c] / u[0]) for c in range(3)]
        u = np.array([rho] + [rho * v for v in vg] + [0.5 * rho * sum(v * v for v in vg)])
        mc.f32_fortran(os.path.join(nyq, "u_gas.init"), u)
        xv = np.empty((6, n, n, n))
        for c in range(3):
            xv[c] = np.mod(Q[c] + drop(m["psi"][c]), n)
            xv[3 + c] = drop(m["vel"][c])
        mc.f32_fortran(os.path.join(nyq, "xv_dm.init"), xv.reshape(6, -1, order="F"))
        smoke = open(os.path.join(src, "tests", "problems_full", "smoke", "params.txt")).read()
        out = {}
        for name in ("music_ics", "music_ics_shifted", "make_ics"):
            p = os.path.join(work, "p_" + name)
            os.makedirs(p)
            shutil.copytree(os.path.join(work, name), os.path.join(p, "ics"))
            txt = re.sub(r"(?m)^N = .*$", f"N = {n}", smoke)
            txt = re.sub(r"(?m)^Ndm = .*$", f"Ndm = {n}", txt)
            txt = re.sub(r"(?m)^box = .*$", f"box = {args.box}", txt)
            txt = re.sub(r"(?m)^z_init = .*$", f"z_init = {args.z}", txt)
            txt = re.sub(r"(?m)^zdata1 = .*$", f"zdata1 = {args.run_z}", txt)
            open(os.path.join(p, "params.txt"), "w").write(txt)
            o = os.path.join(work, "run_" + name)
            run([args.radhydro, p, o, "--quiet"], o + ".log")
            gas = sorted(f for f in os.listdir(os.path.join(o, "Gas")) if f.startswith("u_gas"))[-1]
            dm = sorted(f for f in os.listdir(os.path.join(o, "DM")) if f.startswith("xv_dm"))[-1]
            out[name] = (np.fromfile(os.path.join(o, "Gas", gas), "<f4")[: 5 * n ** 3].reshape(5, n, n, n, order="F")[0],
                         np.fromfile(os.path.join(o, "DM", dm), "<f4")[: 6 * n ** 3].reshape(-1, 6))
        rb, xb = out["make_ics"]
        for name in ("music_ics", "music_ics_shifted"):
            ra, xa = out[name]
            e_rho = float(np.sqrt(np.mean((ra.astype(float) / rb - 1) ** 2)))
            dx = xa[:, :3].astype(float) - xb[:, :3]
            dx -= n * np.round(dx / n)
            e_x = float(np.sqrt(np.mean(dx ** 2)))
            print(f"radhydro to z = {args.run_z}, {name} against make_ics: gas density rms relative "
                  f"difference {e_rho:.3e}, particle positions rms difference {e_x:.3e} cells")
            if name.endswith("shifted") and not (e_rho < 1e-4 and e_x < 1e-4):
                fails.append(f"the run from the converter's files, shifted, differs from "
                             f"make_ics's: {e_rho:.2e}, {e_x:.2e}")
    if fails:
        print("FAIL:\n  " + "\n  ".join(fails))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
