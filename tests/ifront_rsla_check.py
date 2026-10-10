#!/usr/bin/env python3
"""I-front in a uniform medium with a reduced speed of light (c_reduced).

usage: ifront_rsla_check.py RT_DRIVER UNIFORM_PROBLEM_DIR WORKDIR [--plot PNG]
                            [--set key=value ...]

--set is passed on to both rt_driver runs, e.g. --set rt_method=m1 for the
moment RT (docs/rt_moments.md). --rmin R (default 2) sets the smallest
front radius in cells at which the match to (1) is checked: the first-order
moment solver spreads the photons over a few cells (numerical diffusion), so
its fronts run ahead of (1) while they are small (+16% at R ~ 4, < 10% from
R ~ 5); the moment tests use --rmin 5.

Builds a problem from UNIFORM_PROBLEM_DIR (uniform neutral gas at the mean
density, N = 32, z = 10): one bright source in the centre (Luv = 8e28),
da_sim = 2e-4, 15 steps. Runs rt_driver with c_reduced = 1 and 0.1
(legacy_rt = legacy_physics = 0) and measures the I-front radius
R = (3 V_ion / 4 pi)^(1/3) in RT cells, V_ion = sum over cells of x_HII,
after every hydro step.

Finite-c I-front (Shapiro et al. 2006, AJ 131, 1314; their eq. for the
"relativistic" I-front). Photons emitted at t' reach the front at radius R at
t' + R/c_hat, so the front at time t is where the infinite-c front was at the
retarded time:
        R(t) = R_inf(t - R(t)/c_hat).                                (1)
Without recombinations R_inf^3 = 3 Q t and (1) is exactly Shapiro et al.'s
solution, R^3/3 + (Q/c_hat) R = Q t, with the limits R -> c_hat t (early)
and R -> R_inf (late). Recombinations enter (1) only through R_inf.

R_inf(t) is taken from the c_reduced = 1 run, itself corrected with (1)
(its samples (t - R/c, R) are samples of R_inf); c_hat is what the code
actually does, one RT cell per substep: c_hat = 1/dtsec_rt cells/s (nstep_rt
is rounded, so c_hat is 0.083 c here, not exactly 0.1 c). There is no free
parameter.

Checks (exit 1 on failure):
  * causality: R(t) <= c_hat t + 1.5 cells (rays start at the cell centre)
    for both runs;
  * the c_reduced = 0.1 front matches (1) to 12% (rms over the steps <= 6%)
    for R >= 2 cells;
  * the c_reduced = 0.1 front lags the c = 1 front by > 10% early on (the
    test resolves the finite-c effect), and the lag shrinks with time;
  * total photon budget: the emitted photons are the same in both runs (to
    1e-6; the source rate per unit time does not depend on c_reduced).
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

N = 32
NSTEPS = 15
LUV = 8e28
DA = 2e-4


def make_problem(src, dst):
    os.makedirs(dst, exist_ok=True)
    shutil.copy(os.path.join(src, "gas.bin"), os.path.join(dst, "gas.bin"))
    with open(os.path.join(dst, "params.txt"), "w") as f:
        f.write(f"N = {N}\nz_init = 10.0\nda_sim = {DA}\nnsteps = {NSTEPS}\n"
                "legacy_rt = 0\nlegacy_physics = 0\n")
    with open(os.path.join(dst, "sources.txt"), "w") as f:
        f.write(f"16 16 16 {LUV:e} 1.0\n")


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = [float(x) for x in v.split()]
    return d


EXTRA = []


def run(driver, prob, out, c):
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)
    with open(out + ".stdout", "w") as log:
        subprocess.run([driver, prob, out, "--set", f"c_reduced={c}"] + EXTRA, check=True,
                       stdout=log, stderr=subprocess.STDOUT)
    t, R, chat, ngam = [], [], [], 0.0
    sec = 0.0
    for n in range(1, NSTEPS + 1):
        sd = os.path.join(out, f"step_{n:04d}")
        st = parse(os.path.join(sd, "rt_stats.txt"))
        dtsec_rt = st["dtsec_rt"][0]
        sec += dtsec_rt * parse(os.path.join(sd, "stats.txt"))["nstep_rt"][0]
        fHI = np.fromfile(os.path.join(sd, "fHI.bin"), dtype=np.float32).astype(np.float64)
        V = np.sum(1 - fHI)
        t.append(sec)
        R.append((3 * V / (4 * np.pi)) ** (1 / 3))
        chat.append(1 / dtsec_rt)  # cells per second
        ngam += st["ngam"][0]  # photons emitted (> 13.6 eV), per H atom of the box
    return np.array(t), np.array(R), np.array(chat), ngam


def retarded(t, R, chat, Rinf_t, Rinf_R):
    """Solve R = R_inf(t - R/chat) by bisection (R_inf^3 interpolated in t)."""
    out = []
    for ti, ci in zip(t, chat):
        lo, hi = 0.0, ci * ti
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            tr = ti - mid / ci
            # interpolate the ionised volume R^3, which grows ~linearly in t
            Rr = np.interp(tr, Rinf_t, Rinf_R ** 3, left=0.0) ** (1 / 3)
            if mid > Rr:
                hi = mid
            else:
                lo = mid
        out.append(0.5 * (lo + hi))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("driver")
    ap.add_argument("uniform")
    ap.add_argument("work")
    ap.add_argument("--plot")
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--rmin", type=float, default=2.0)
    a = ap.parse_args()
    for kv in a.set:
        EXTRA.extend(["--set", kv])
    prob = os.path.join(a.work, "ifront_uniform")
    make_problem(a.uniform, prob)
    t1, R1, c1, g1 = run(a.driver, prob, os.path.join(a.work, "c1"), 1.0)
    t2, R2, c2, g2 = run(a.driver, prob, os.path.join(a.work, "c0.1"), 0.1)

    ok = True
    # R_inf samples from the c = 1 run: (t - R/c, R), plus the origin
    Rinf_t = np.concatenate([[0.0], t1 - R1 / c1])
    Rinf_R = np.concatenate([[0.0], R1])
    pred = retarded(t2, R2, c2, Rinf_t, Rinf_R)
    print(" step   t[Myr]  c_hat t  R(c=1)  R(c=0.1)  R_pred(0.1)  rel.err")
    errs = []
    for n in range(NSTEPS):
        e = R2[n] / pred[n] - 1
        print(f" {n + 1:4d} {t2[n] / 3.156e13:8.3f} {c2[n] * t2[n]:8.3f} {R1[n]:7.3f} "
              f"{R2[n]:9.3f} {pred[n]:12.3f} {e:+8.3f}")
        if R2[n] >= a.rmin:
            errs.append(e)
            if abs(e) > 0.12:
                print(f"  FAIL: step {n + 1}: |R/R_pred - 1| > 0.12")
                ok = False
    rms = float(np.sqrt(np.mean(np.square(errs))))
    print(f" rms error (R >= {a.rmin:g}): {rms:.4f}; c_hat(c=0.1)/c_hat(c=1) = {c2[0] / c1[0]:.4f}")
    if rms > 0.06:
        print("  FAIL: rms error > 0.06")
        ok = False
    for (t, R, c, lab) in ((t1, R1, c1, "c=1"), (t2, R2, c2, "c=0.1")):
        bad = R > c * t + 1.5
        if np.any(bad):
            print(f"  FAIL: {lab}: front faster than c_hat at steps {np.nonzero(bad)[0] + 1}")
            ok = False
    lag = 1 - R2 / R1
    print(f" lag 1 - R(0.1)/R(1): step 1 {lag[0]:.3f}, step 3 {lag[2]:.3f}, last {lag[-1]:.3f}")
    if not (lag[2] > 0.10 and lag[-1] < lag[2]):
        print("  FAIL: the finite-c lag is not resolved or does not shrink")
        ok = False
    print(f" photons emitted (sum over steps, per H atom): c=1 {g1:.9g}, c=0.1 {g2:.9g}")
    if abs(g2 / g1 - 1) > 1e-6:
        print("  FAIL: emitted photons depend on c_reduced")
        ok = False
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.5, 4))
        my = 3.156e13
        ax.plot(t1 / my, R1, "o-", label="c_reduced = 1")
        ax.plot(t2 / my, R2, "s-", label="c_reduced = 0.1")
        ax.plot(t2 / my, pred, "k--", label="R = R_inf(t - R/c_hat)")
        ax.plot(t2 / my, c2 * t2, ":", color="gray", label="c_hat t (0.1)")
        ax.set_ylim(0, 1.15 * max(R1))
        ax.set_xlabel("t [Myr]")
        ax.set_ylabel("I-front radius [RT cells]")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=120)
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
