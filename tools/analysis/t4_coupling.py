#!/usr/bin/env python3
"""T4: how short must the RT-hydro coupling interval be while a plane-parallel
front crosses the box (docs/zoom_turbulence.md section 4.4)?

usage: t4_coupling.py make    ROOT [--N 64 --box 0.25 --z-re 7 --gamma12 0.3 --seed 7]
       t4_coupling.py run     ROOT --exe build/radhydro [--threads 4] [--variants base,dt100,...]
       t4_coupling.py analyse ROOT [--outdir ROOT/t4] [--ref sub]

The problem: the C1 set-up of docs/zoom_turbulence.md section 3b at N^3
(tools/production/make_production.py --plane-parallel: 0.25 Mpc/h, 8^3 RT
domains lit from +x and +y at Gamma_-12 from z_re, DM20's speed-of-light
schedule 0.1 c -> 0.01 c at x_HI,V = 1 per cent, cooling on everywhere,
Zel'dovich initial conditions of the existing generator), with the gas
state (u_gas) written at the outputs z_re, and 25, 50, 100, 200 and 400 Myr
after it.

The variants, the same problem with one params.txt override each:
  base    the hydro step as it is: ra_max (da/a <= 1/9) and the CFL, so the
          RT runs its 10^3-10^4 substeps inside hydro steps of 50-100 Myr and
          the gas feels the heating once per step;
  dt100, dt10, dt1, dt0.1   dt_max_myr = 100, 10, 1, 0.1: the heating is
          coupled to the hydro every that many Myr;
  sub     rt_substeps_max = 1: the hydro step is one RT substep, the light
          crossing of a cell at c_hat (0.02 Myr at 0.1 c, 0.2 Myr at 0.01 c
          here): coupling every substep, as C26's uniform boxes.

analyse: the folded front profiles (tools/analysis/plane_front.py) and the
kinetic energy spectra E(k) (tools/analysis/velocity_power.py) at every
output of every variant, and the comparison with `sub`: the largest
relative difference of E(k) over 0.5 < k < k_Nyquist/2 (where the
turbulence of C26 lives; from 4 k_f at low N) and of the folded x_HI and T profiles, per
output. Writes OUTDIR/t4_coupling.json, t4_coupling.txt (the table) and
t4_coupling.png (E(k) and the profiles of every variant at the last two
outputs).
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
VARIANTS = {
    "base": [],
    "dt100": ["dt_max_myr=100"],
    "dt10": ["dt_max_myr=10"],
    "dt1": ["dt_max_myr=1"],
    "dt0.1": ["dt_max_myr=0.1"],
    "sub": ["rt_substeps_max=1"],
}
# outputs: z_re and 25, 50, 100, 200, 400 Myr after it (z(t) for the code's
# EdS-like early-time growth: t ~ a^1.5, computed in `make`)
DT_MYR = [25, 50, 100, 200, 400]

OM, H0 = float(np.float32(0.27)), float(np.float32(0.70))
H0_CGS = 3.24086e-18
MYR = 3.156e13


def t_of_a(a):
    # Omega_m + Omega_Lambda flat: t = 2/(3 H0 sqrt(OL)) asinh(sqrt(OL/OM) a^1.5)
    OL = 1 - OM
    return 2 / (3 * H0_CGS * H0 * np.sqrt(OL)) * np.arcsinh(np.sqrt(OL / OM) * a ** 1.5)


def a_of_t(t):
    OL = 1 - OM
    return (np.sinh(1.5 * H0_CGS * H0 * np.sqrt(OL) * t) * np.sqrt(OM / OL)) ** (2 / 3)


def zdata_after(z_re):
    t0 = t_of_a(1 / (1 + z_re))
    zs = [z_re] + [float(1 / a_of_t(t0 + dt * MYR) - 1) for dt in DT_MYR]
    return " ".join(f"{z:.4f}" for z in zs)


def make(a):
    os.makedirs(a.root, exist_ok=True)
    prob = os.path.join(a.root, "prob")
    cmd = [sys.executable, os.path.join(ROOT, "tools/production/make_production.py"), prob,
           "--N", str(a.N), "--ndm", str(a.N), "--box", str(a.box), "--z-init", "100", "--seed", str(a.seed),
           "--plane-parallel", "--plane-ndom", str(a.ndom), "--plane-z-re", str(a.z_re),
           "--plane-gamma12", str(a.gamma12), "--zdata1", zdata_after(a.z_re),
           "--fields", "fHI,temp,rhog,u_gas", "--data-files", "gas", "--movie-dz", "0",
           "--restart-every", "0", "--nray-max-factor", "8"]
    subprocess.check_call(cmd)
    print("problem in", prob, "outputs at z =", zdata_after(a.z_re))


def run(a):
    prob = os.path.join(a.root, "prob")
    for v in a.variants.split(","):
        out = os.path.join(a.root, v)
        if os.path.exists(os.path.join(out, "DONE")) and not a.force:
            print(v, "done")
            continue
        cmd = [a.exe, prob, out, "--quiet", "--set", f"threads={a.threads}"]
        for kv in VARIANTS[v]:
            cmd += ["--set", kv]
        print("running", v, "...", flush=True)
        with open(out + ".log", "w") as log:
            r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            sys.exit(f"{v} failed: see {out}.log")
        open(os.path.join(out, "DONE"), "w").close()


def analyse(a):
    outdir = a.outdir or os.path.join(a.root, "t4")
    os.makedirs(outdir, exist_ok=True)
    variants = [v for v in a.variants.split(",") if os.path.isdir(os.path.join(a.root, v))]
    ref = a.ref
    if ref not in variants:
        sys.exit(f"analyse needs the reference variant `{ref}` (--ref)")
    specs = [f"{os.path.join(a.root, v)}:{v}" for v in variants]
    subprocess.check_call([sys.executable, os.path.join(HERE, "velocity_power.py")] + specs +
                          ["--outdir", outdir, "--json", os.path.join(outdir, "ek.json"), "--no-fig"])
    subprocess.check_call([sys.executable, os.path.join(HERE, "plane_front.py")] + specs +
                          ["--outdir", outdir, "--json", os.path.join(outdir, "front.json"), "--no-fig"])
    ek = json.load(open(os.path.join(outdir, "ek.json")))
    fr = json.load(open(os.path.join(outdir, "front.json")))
    ref_ek = {round(r["z"], 2): r for r in ek if r["label"] == ref}
    ref_pr = {round(p["z"], 2): p for p in fr["profiles"] if p["label"] == ref}
    table = []
    for v in variants:
        for r in [x for x in ek if x["label"] == v]:
            zk = min(ref_ek, key=lambda z: abs(z - r["z"]))
            if abs(zk - r["z"]) > 0.05:
                continue
            R = ref_ek[zk]
            k = np.array(r["k"])
            E, E0 = np.array(r["E"]), np.array(R["E"])
            kf, kN = 2 * np.pi / (r["box"] * 1e3), np.pi * r["N"] / (r["box"] * 1e3)
            sel = (k > min(0.5, 4 * kf)) & (k < 0.5 * kN) & (E0 > 0)  # the inertial range of C26 (0.5 < k < kN/2), or 4 k_f .. kN/2 at low N
            dE = float(np.max(np.abs(E[sel] / E0[sel] - 1))) if sel.any() else float("nan")
            dxHI = dT = float("nan")
            pz = min(ref_pr, key=lambda z: abs(z - r["z"]))
            pv = [p for p in fr["profiles"] if p["label"] == v and abs(p["z"] - r["z"]) < 0.05]
            if pv and abs(pz - r["z"]) < 0.05:
                x0, x1 = np.array(ref_pr[pz]["xHI"]), np.array(pv[0]["xHI"])
                T0, T1 = np.array(ref_pr[pz]["T"]), np.array(pv[0]["T"])
                dxHI = float(np.max(np.abs(x1 - x0)))
                dT = float(np.max(np.abs(T1[T0 > 0] / T0[T0 > 0] - 1))) if np.any(T0 > 0) else float("nan")  # T = 0 at the turn-on output
            table.append(dict(variant=v, z=r["z"], v2=r["v2"], v2_ref=R["v2"], dE_max=dE, dxHI_max=dxHI, dT_max=dT))
    with open(os.path.join(outdir, "t4_coupling.json"), "w") as f:
        json.dump(table, f, indent=1)
    lines = [f"# variant z <v^2>[(km/s)^2] <v^2>_{ref} max|E/E_{ref}-1| (0.5<k<kN/2) max|xHI-xHI_{ref}| max|T/T_{ref}-1|"]
    for t in table:
        lines.append(f"{t['variant']:6s} {t['z']:.3f} {t['v2']:10.4g} {t['v2_ref']:10.4g} {t['dE_max']:8.3f} {t['dxHI_max']:8.3f} {t['dT_max']:8.3f}")
    open(os.path.join(outdir, "t4_coupling.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    zs = sorted(set(round(r["z"], 2) for r in ek))[:2]  # the last two outputs
    fig, axs = plt.subplots(2, len(zs), figsize=(5 * len(zs), 8), squeeze=False)
    styles = {"base": "-", "dt100": "--", "dt10": "-.", "dt1": ":", "dt0.1": (0, (5, 1)), "sub": "-"}
    for j, z in enumerate(zs):
        for v in variants:
            rr = [r for r in ek if r["label"] == v and abs(r["z"] - z) < 0.05]
            pp = [p for p in fr["profiles"] if p["label"] == v and abs(p["z"] - z) < 0.05]
            if rr:
                k, E = np.array(rr[0]["k"]), np.array(rr[0]["E"])
                axs[0, j].loglog(k[E > 0], E[E > 0], ls=styles.get(v, "-"), lw=2 if v == ref else 1, label=v)
            if pp:
                axs[1, j].semilogy(np.maximum(pp[0]["T"], 1), ls=styles.get(v, "-"), lw=2 if v == ref else 1, label=v)
        axs[0, j].set_title(f"z = {z}")
        axs[0, j].set_xlabel("k [h/kpc]")
        axs[0, j].set_ylabel("E(k) [(km/s)^2 kpc/h]")
        axs[1, j].set_xlabel("cells from the source face")
        axs[1, j].set_ylabel("<T> [K]")
        axs[0, j].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "t4_coupling.png"), dpi=140)
    print("wrote", os.path.join(outdir, "t4_coupling.png"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("phase", choices=("make", "run", "analyse"))
    ap.add_argument("root")
    ap.add_argument("--N", type=int, default=64)
    ap.add_argument("--box", type=float, default=0.25)
    ap.add_argument("--ndom", type=int, default=8)
    ap.add_argument("--z-re", type=float, default=7.0)
    ap.add_argument("--gamma12", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--exe", default=os.path.join(ROOT, "build", "radhydro"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--ref", default="sub", help="the reference variant of analyse (default sub)")
    a = ap.parse_args()
    {"make": make, "run": run, "analyse": analyse}[a.phase](a)


if __name__ == "__main__":
    main()
