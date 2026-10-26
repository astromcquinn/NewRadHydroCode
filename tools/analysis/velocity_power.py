#!/usr/bin/env python3
"""Kinetic energy spectrum E(k) of the gas velocity field (C26 Fig. 3;
docs/plane_parallel.md, docs/zoom_turbulence.md section 3b).

usage: velocity_power.py RUN[:label] [RUN:label ...] --outdir DIR [--z 6.5 5.5 ...]
                         [--steps ...] [--nbins 32] [--weight none|mass|sqrtrho]
                         [--ref FILE] [--no-fig] [--json FILE]

RUN is a run directory (the `run/` of a radhydro problem, or the problem
directory holding run/ and params.txt). The velocity comes from the gas
state u_gas (rho, rho v, E) of each step_NNNN/u_gas.bin (fields = ...,u_gas)
or of run/Gas/u_gas.<ztag> (the zdata1/zdata2 outputs; data_files = gas or
all), converted to proper km/s with vel_unit = dl/time_unit
(docs/igm_diagnostics.md section 1). For every requested redshift the
nearest output is used.

E(k): the velocity field v(x) (km/s, peculiar, proper) is Fourier
transformed on the N^3 grid; E(k) dk is the kinetic energy per unit mass in
the shell [k, k + dk], normalised so that int E(k) dk = <v^2>/2 (the
specific kinetic energy of the box, (km/s)^2). --weight mass uses
sqrt(rho/<rho>) v (the kinetic energy per unit volume over the mean
density, int E dk = <rho v^2>/(2 <rho>)); none (default) is the plain
velocity field, as a volume-weighted spectrum. k is comoving, in h/kpc
(k = 2 pi / lambda), the box's fundamental k_f = 2 pi / L to the Nyquist
k_N = pi N / L, in --nbins logarithmic bins of |k|; a k^-5/3 line is drawn
through the bin nearest k = 2 h/kpc for the eye. --ref FILE overlays a
two-column table (k [h/kpc], E(k) [(km/s)^2 kpc/h]) such as digitised
points of C26 Fig. 3.

Writes DIR/velocity_power.txt (label, z, k, E(k), the shell mode count),
DIR/velocity_power.png (unless --no-fig) and, with --json, the same table.
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np

# cosmo_parameters.f90 / physical_constants.f90 (include/radray/constants.hpp)
OM, H0 = float(np.float32(0.27)), float(np.float32(0.70))
H0_CGS, MPC = 3.24086e-18, 3.08560e24


def params_of(run):
    """params.txt of the problem: run/../params.txt or run/params.txt."""
    for p in (os.path.join(run, "params.txt"), os.path.join(run, "..", "params.txt")):
        if os.path.exists(p):
            d = {}
            for line in open(p):
                line = line.split("#")[0]
                if "=" in line:
                    k, v = line.split("=", 1)
                    d[k.strip()] = v.strip()
            return d
    sys.exit(f"{run}: no params.txt in it or above it")


def ztag(z):
    return f"z=0{z:6.4f}" if z < 10 else f"z={z:7.4f}"


def run_dir(spec):
    if os.path.isdir(os.path.join(spec, "run")):
        return os.path.join(spec, "run")
    return spec


def outputs(run):
    """[(z, path of a u_gas file)] of the run, by redshift."""
    out = []
    for sd in sorted(glob.glob(os.path.join(run, "step_*"))):
        ug = os.path.join(sd, "u_gas.bin")
        st = os.path.join(sd, "stats.txt")
        if os.path.exists(ug) and os.path.exists(st):
            for line in open(st):
                if line.startswith("redshift"):
                    out.append((float(line.split("=")[1]), ug))
    for p in glob.glob(os.path.join(run, "Gas", "u_gas.z=*")):
        m = re.search(r"z=([0-9.]+)", os.path.basename(p))
        if m:
            out.append((float(m.group(1)), p))
    if not out:
        sys.exit(f"{run}: no u_gas output (fields = ...,u_gas or data_files = gas|all)")
    return sorted(out, key=lambda t: -t[0])


def velocity(path, N, z, box, weight):
    """v (km/s, proper peculiar), [3, N, N, N] in [k, j, i]; optionally
    weighted by sqrt(rho/<rho>)."""
    u = np.fromfile(path, dtype=np.float32)
    if u.size != 5 * N ** 3:
        sys.exit(f"{path}: {u.size} floats, not 5 N^3 with N = {N}")
    u = u.reshape(N, N, N, 5).astype(np.float64)
    a = 1 / (1 + z)
    dl = a * box / H0 / N * MPC
    time_unit = 2 / (3 * H0_CGS * H0 * np.sqrt(OM)) * a ** 2
    vel_unit = dl / time_unit / 1e5  # km/s per grid velocity
    rho = u[..., 0]
    v = u[..., 1:4] / rho[..., None] * vel_unit
    if weight == "mass":
        v = v * np.sqrt(rho / rho.mean())[..., None]
    return np.moveaxis(v, -1, 0)


class Shells:
    def __init__(self, N, box_kpc, nbins):
        kf = 2 * np.pi / box_kpc
        n = np.fft.fftfreq(N, 1.0 / N)
        nz = np.arange(N // 2 + 1)
        k2 = n[:, None, None] ** 2 + n[None, :, None] ** 2 + nz[None, None, :] ** 2
        kmag = np.sqrt(k2.astype(np.float64)) * kf
        self.edges = np.logspace(np.log10(kf * 0.999), np.log10(np.pi * N / box_kpc * 1.001), nbins + 1)
        self.idx = np.digitize(kmag.ravel(), self.edges) - 1
        w = np.full((N, N, N // 2 + 1), 2.0)
        w[:, :, 0] = 1.0
        if N % 2 == 0:
            w[:, :, N // 2] = 1.0
        self.w = w.ravel()
        self.ok = (self.idx >= 0) & (self.idx < nbins)
        self.nb = nbins
        self.cnt = np.bincount(self.idx[self.ok], weights=self.w[self.ok], minlength=nbins)
        self.k = np.bincount(self.idx[self.ok], weights=(self.w * kmag.ravel())[self.ok], minlength=nbins) / np.maximum(self.cnt, 1)
        self.dk = np.diff(self.edges)
        self.N = N

    def energy(self, v):
        """E(k) per bin: sum over the shell of |v_hat|^2 / (2 N^6) / dk."""
        p = np.zeros(self.N * self.N * (self.N // 2 + 1))
        for c in range(3):
            vh = np.fft.rfftn(v[c])
            p += (vh.real ** 2 + vh.imag ** 2).ravel()
        s = np.bincount(self.idx[self.ok], weights=(self.w * p)[self.ok], minlength=self.nb)
        return s / (2.0 * float(self.N) ** 6) / self.dk


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--z", nargs="*", type=float, default=None, help="redshifts (default: every output)")
    ap.add_argument("--nbins", type=int, default=32)
    ap.add_argument("--weight", choices=("none", "mass"), default="none")
    ap.add_argument("--ref", default=None, help="k [h/kpc], E(k) table to overlay (e.g. C26 Fig. 3)")
    ap.add_argument("--no-fig", action="store_true")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    for spec in a.runs:
        path, _, label = spec.partition(":")
        run = run_dir(path)
        P = params_of(run)
        N, box = int(P["N"]), float(P.get("box", 5))
        label = label or os.path.basename(os.path.dirname(run.rstrip("/"))) or run
        outs = outputs(run)
        want = outs if a.z is None else [min(outs, key=lambda t: abs(t[0] - z)) for z in a.z]
        shells = Shells(N, box * 1e3, a.nbins)
        for z, ug in want:
            E = shells.energy(velocity(ug, N, z, box, a.weight))
            v2 = 2 * np.sum(E * shells.dk)
            rows.append(dict(label=label, z=z, N=N, box=box, k=shells.k.tolist(), E=E.tolist(),
                             count=shells.cnt.tolist(), v2=float(v2)))
            print(f"{label}: z = {z:.3f}, N = {N}, <v^2> = {v2:.3g} (km/s)^2, sqrt = {np.sqrt(v2):.3g} km/s")
    with open(os.path.join(a.outdir, "velocity_power.txt"), "w") as f:
        f.write("# label z k[h/kpc] E(k)[(km/s)^2 kpc/h] modes\n")
        for r in rows:
            for k, E, c in zip(r["k"], r["E"], r["count"]):
                f.write(f"{r['label']} {r['z']:.4f} {k:.6e} {E:.6e} {c:.0f}\n")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(rows, f, indent=1)
    if a.no_fig:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: no figure")
        return
    fig, ax = plt.subplots(figsize=(6, 4.5))
    labels = sorted(set(r["label"] for r in rows))
    zs = sorted(set(round(r["z"], 3) for r in rows), reverse=True)
    cmap = plt.get_cmap("viridis")
    styles = ["-", "--", ":", "-.", (0, (5, 1)), (0, (3, 1, 1, 1))]
    for r in rows:
        col = cmap(0.1 + 0.8 * zs.index(round(r["z"], 3)) / max(1, len(zs) - 1))
        ls = styles[labels.index(r["label"]) % len(styles)]
        k, E = np.array(r["k"]), np.array(r["E"])
        ok = E > 0
        ax.loglog(k[ok], E[ok], ls=ls, color=col, label=f"{r['label']} z = {r['z']:.2f}")
    if rows:
        k = np.array(rows[-1]["k"])
        E = np.array(rows[-1]["E"])
        i0 = int(np.argmin(np.abs(k - 2.0)))
        if E[i0] > 0:
            kk = np.array([k[i0] / 4, k[i0] * 4])
            ax.loglog(kk, E[i0] * (kk / k[i0]) ** (-5 / 3), color="0.5", lw=0.8, label=r"$k^{-5/3}$")
    if a.ref:
        ref = np.loadtxt(a.ref)
        ax.loglog(ref[:, 0], ref[:, 1], "ko", ms=3, label=os.path.basename(a.ref))
    ax.set_xlabel(r"$k$ [$h$ kpc$^{-1}$]")
    ax.set_ylabel(r"$E(k)$ [(km/s)$^2$ kpc/$h$]")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "velocity_power.png"), dpi=150)
    print("wrote", os.path.join(a.outdir, "velocity_power.png"))


if __name__ == "__main__":
    main()
