#!/usr/bin/env python3
"""Stromgren sphere with case A, case B and case A + recombination radiation
(docs/recomb_radiation.md section 3).

usage: recomb_stromgren.py RT_DRIVER UNIFORM_PROBLEM_DIR WORKDIR [--plot PNG]
                           [--check] [--threads n] [--reuse]
                           [--method rays|m1|otvet]

Builds a uniform problem from UNIFORM_PROBLEM_DIR (tests/problems/
uniform_1src, mean density, N = 32) at 30x the mean density and z = 30
(n_H ~ 0.17 cm^-3, recombination time ~ 0.7 Myr at 1e4 K, below
rho_source = 50 so that every photoionisation uses a photon), with one
source of Luv = 6e31 in the centre, c_reduced = 1, legacy_rt =
legacy_physics = 0, 40 steps of da = 4e-5 (about 12 recombination times).
Runs rt_driver with

    B      recomb_case = B (the original)
    A      recomb_case = A
    A+RR   recomb_case = A, recomb_rad = 1 (block 1, and block 2)

and measures the ionised radius R = (3 V_HII / 4 pi)^(1/3) in cells after
every step. At equilibrium a Stromgren sphere has R^3 ~ 1/alpha, so
R_A / R_B = (alpha_B / alpha_A)^(1/3) ~ 0.86 at 1e4 K. With the
recombination radiation transported, the diffuse photons are absorbed in the
optically thick sphere and R should come back to about R_B (the on-the-spot
approximation is exact for the total in this geometry).

--check (the ctest): exit 1 unless, at the last step (T = the median
temperature of each run's ionised gas),
  * R_A is within 4% of R_B (alpha_B(T_B) / alpha_A(T_A))^(1/3),
  * R(A+RR) is within [0.97, 1.06] of R_B (alpha_B(T_B)/alpha_B(T_RR))^(1/3)
    for both blocks. It comes out a few per cent above 1: the diffuse
    photons widen the soft edge, and they fill the holes that the stellar
    rays' HEALPix sampling leaves inside the case B sphere (B16).
--method m1 | otvet runs the three cases with the moment RT instead (A+RR
once: every cell emits into the moments, no blocks). The radii are compared
with that method's own case B.
The fractions of the ground-state recombinations launched as rays and kept
on the spot (case B: optically thick cells, or no room in the ray pool) are
printed from ray.log.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

import numpy as np

N = 32
NSTEPS = 40
DA = 4e-5
Z0 = 30.0
DENS = 30.0
LUV = 6e31


def make_problem(src, dst):
    os.makedirs(dst, exist_ok=True)
    g = np.fromfile(os.path.join(src, "gas.bin"), dtype=np.float32).reshape(-1, 9)
    g = g.copy()
    # density, momenta (zero), energy and the species / metal densities all
    # scale with the density at fixed temperature and abundances
    g[:, 0] *= DENS
    g[:, 1:4] *= DENS
    g[:, 4] *= DENS
    g[:, 5:9] *= DENS
    g.astype(np.float32).tofile(os.path.join(dst, "gas.bin"))
    with open(os.path.join(dst, "params.txt"), "w") as f:
        f.write(f"N = {N}\nz_init = {Z0}\nda_sim = {DA}\nnsteps = {NSTEPS}\n"
                "legacy_rt = 0\nlegacy_physics = 0\nc_reduced = 1\n")
    with open(os.path.join(dst, "sources.txt"), "w") as f:
        f.write(f"16 16 16 {LUV:e} 1.0\n")


def parse(fn):
    d = {}
    for line in open(fn):
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                d[k.strip()] = [float(x) for x in v.split()]
            except ValueError:
                pass
    return d


def run(driver, prob, out, sets, threads, reuse=False):
    if not (reuse and os.path.exists(os.path.join(out, f"step_{NSTEPS:04d}", "fHI.bin"))):
        if os.path.isdir(out):
            shutil.rmtree(out)
        os.makedirs(out)
        cmd = [driver, prob, out]
        for kv in sets:
            cmd += ["--set", kv]
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        with open(out + ".stdout", "w") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT, env=env)
    R, T = [], []
    for n in range(1, NSTEPS + 1):
        sd = os.path.join(out, f"step_{n:04d}")
        fHI = np.fromfile(os.path.join(sd, "fHI.bin"), dtype=np.float32).astype(np.float64)
        R.append((3 * np.sum(1 - fHI) / (4 * np.pi)) ** (1 / 3))
        tp = os.path.join(sd, "temp.bin")
        if os.path.exists(tp):
            t = np.fromfile(tp, dtype=np.float32)
            ion = fHI < 0.5
            T.append(float(np.median(t[ion])) if ion.any() else np.nan)
    budget = None
    rl = os.path.join(out, "ray.log")
    if os.path.exists(rl):
        last = [l for l in open(rl) if l.startswith(" Recomb")]
        if last:
            v = dict(re.findall(r"(launched|dropped|onspot|emitted) ([0-9.eE+-]+)", last[-1]))
            budget = {k: float(x) for k, x in v.items()}
            if "emitted" in budget:  # moment RT: everything is emitted
                budget = {"launched": budget["emitted"], "dropped": 0.0, "onspot": 0.0}
    return np.array(R), np.array(T), budget


def alpha_A(T):  # Hui & Gnedin 1997, HII case A
    lam = 2 * 157807.0 / T
    return 1.269e-13 * lam ** 1.503 / (1 + (lam / 0.522) ** 0.470) ** 1.923


def alpha_B(T):  # Hui & Gnedin 1997, HII case B
    lam = 2 * 157807.0 / T
    return 2.753e-14 * lam ** 1.500 / (1 + (lam / 2.740) ** 0.407) ** 2.242


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("driver")
    ap.add_argument("uniform")
    ap.add_argument("work")
    ap.add_argument("--plot")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--reuse", action="store_true", help="keep finished runs in WORKDIR")
    ap.add_argument("--method", default="rays", choices=["rays", "m1", "otvet"])
    a = ap.parse_args()
    prob = os.path.join(a.work, "stromgren")
    make_problem(a.uniform, prob)
    cases = [("B", []), ("A", ["recomb_case=A"]),
             ("A+RR (block 1)", ["recomb_case=A", "recomb_rad=1", "recomb_rad_block=1"]),
             ("A+RR (block 2)", ["recomb_case=A", "recomb_rad=1", "recomb_rad_block=2"])]
    if a.method != "rays":
        cases = [(lab, sets + [f"rt_method={a.method}"]) for lab, sets in cases[:3]]
        cases[2] = ("A+RR", cases[2][1])
    res = {}
    for lab, sets in cases:
        tag = lab.replace(" ", "").replace("(", "_").replace(")", "").replace("+", "p")
        if a.method != "rays":
            tag = f"{a.method}_{tag}"
        res[lab] = run(a.driver, prob, os.path.join(a.work, tag), sets, a.threads, a.reuse)
    RB = res["B"][0]
    TB = res["B"][1][-1]
    print(f"step  " + "  ".join(f"{lab:>14s}" for lab, _ in cases))
    for n in range(NSTEPS):
        print(f"{n + 1:4d}  " + "  ".join(f"{res[lab][0][n]:14.3f}" for lab, _ in cases))
    # Equilibrium R^3 ~ 1/alpha(T), with each run's own temperature (the
    # case A sphere is hotter); the case A + RR run is compared with case B
    # corrected to its temperature.
    print(f"\nlast step (R in cells; T = median of the ionised gas):")
    ok = True
    for lab, _ in cases:
        R, T, bud = res[lab]
        Tl = T[-1]
        if lab == "B":
            pred = RB[-1]
        elif lab == "A":
            pred = RB[-1] * (alpha_B(TB) / alpha_A(Tl)) ** (1 / 3)
        else:
            pred = RB[-1] * (alpha_B(TB) / alpha_B(Tl)) ** (1 / 3)
        rat = R[-1] / pred
        line = f"  {lab:16s} R = {R[-1]:6.3f}  T = {Tl:6.0f} K  expected {pred:6.3f}  ratio {rat:.3f}"
        if bud:
            tot = bud["launched"] + bud["dropped"] + bud["onspot"]
            line += (f"\n      ground-state recombinations per H atom {tot:.4g}: launched"
                     f" {bud['launched'] / tot:.1%}, on the spot {bud['onspot'] / tot:.1%},"
                     f" dropped {bud['dropped'] / tot:.2%}")
        print(line)
        if lab == "A" and abs(rat - 1) > 0.04:
            ok = False
            print("    FAIL: case A radius off (alpha_B/alpha_A)^(1/3) R_B by more than 4%")
        if lab.startswith("A+RR") and not (0.97 <= rat <= 1.06):
            ok = False
            print("    FAIL: R(A+RR) outside [0.97, 1.06] of the temperature-corrected R_B")
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4))
        cols = ["#2a78d6", "#eb6834", "#1baf7a", "#8c5bd6"]
        for (lab, _), c in zip(cases, cols):
            ax.plot(np.arange(1, NSTEPS + 1), res[lab][0], color=c, lw=2, label=lab)
        ax.set_xlabel("step (da = 4e-5 each, z = 30)")
        ax.set_ylabel("ionised radius [cells]")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=120)
    if a.check:
        print("PASS" if ok else "FAIL")
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
