#!/usr/bin/env python3
"""T_GammaL(k) = [P_Gamma(k)/P_L(k)]^1/2 of runs against their own source
field, with arctan(k lambda)/(k lambda) (Gnedin 2025, eq. 2) for each
snapshot's mean free path (docs/rt_moments.md 8.6, A4).

usage: gamma_source_transfer.py RUN:label:IGM_JSON [...] --outdir DIR
                                [--steps 69,84,93,99] [--box 40] [--nbins 24]
                                [--field gammaHI_avg] [--peak-delta 10] [--workers N]

RUN is a run directory (<problem>/run); IGM_JSON its igm_stats.py output
(the mean free path per snapshot: mfp.ionized.lambda_fit in proper Mpc, from
sight lines started in ionised cells, converted to comoving Mpc/h with the
run's hubble0; mfp.random is printed beside it). The first run is the
reference for the steps; every run's snapshot nearest in redshift is used.

The source field. The snapshots hold no source list (no Luv.bin), and the
checkpoints with the particles are gone, so L is a stand-in from the gas:
the local maxima (3^3, periodic) of the run's own rhog above --peak-delta
times the mean (the halo centres' stand-in; the run.log 'Nsource' counts
source halos, not cells, and exceeds the number of density maxima in the
box, so it cannot set the count), each with the weight
  'mass'  the gas mass of its 3^3 neighbourhood (sum of rhog), a proxy for
          the halo mass that the abundance-matched luminosity rises with;
  'num'   1 (the number density of sources).
delta_L = L/<L> - 1 and delta_Gamma = Gamma/<Gamma>_V - 1 on the same grid.
The stand-in's own errors (missed and spurious sources, the weights) enter
P_L and the cross spectrum; r(k) of Gamma with L measures how well they
agree.

Writes DIR/gamma_source_transfer.txt: per step, per k bin and run,
T = (P_Gamma/P_L)^1/2 for both weights, T normalised to its value in the
lowest k bin, r(k), and arctan(k lambda)/(k lambda).
"""
import argparse
import json
import os
import sys

import numpy as np
import scipy.fft as sfft
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gamma_power as gp  # noqa: E402


def hubble(run):
    p = os.path.join(os.path.dirname(run.rstrip("/")), "params.txt")
    for line in open(p):
        k, _, v = line.partition("=")
        if k.strip() == "hubble0":
            return float(v.split("#")[0])
    return 0.70  # constants.hpp hubble0 (real(4) 0.70)


def mfp(js, z):
    d = json.load(open(js))
    s = min(d["snapshots"], key=lambda e: abs(e["z"] - z))
    m = s["mfp"]
    return s["z"], m["ionized"]["lambda_fit"], m["random"]["lambda_fit"], s["xHII_V"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--steps", default="69,84,93,99")
    ap.add_argument("--box", type=float, default=40.0)
    ap.add_argument("--nbins", type=int, default=24)
    ap.add_argument("--field", default="gammaHI_avg")
    ap.add_argument("--peak-delta", type=float, default=10.0)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    a = ap.parse_args()
    gp.FIELD = a.field
    runs = [tuple(s.split(":")) for s in a.runs]
    os.makedirs(a.outdir, exist_ok=True)
    out = open(os.path.join(a.outdir, "gamma_source_transfer.txt"), "w")

    def pr(s):
        print(s, flush=True)
        out.write(s + "\n")
        out.flush()

    pr("# gamma_source_transfer.py " + " ".join(sys.argv[1:]))
    snaps = {lab: gp.snapshots(r, 0.0, 99.0) for r, lab, _ in runs}
    B = None
    for st in [int(x) for x in a.steps.split(",")]:
        d0 = os.path.join(runs[0][0], f"step_{st:04d}")
        zr = gp.read_z(d0)
        res = {}
        for r, lab, js in runs:
            d = snaps[lab][min(snaps[lab], key=lambda z: abs(z - zr))]
            z = gp.read_z(d)
            g = gp.load_gamma(d).astype(np.float64)
            N = g.shape[0]
            if B is None:
                B = gp.Binner(N, a.box, a.nbins)
            rho = np.fromfile(os.path.join(d, "rhog.bin"), dtype=np.float32).reshape(N, N, N)
            mx = ndimage.maximum_filter(rho, size=3, mode="wrap")
            cand = np.flatnonzero((rho == mx).ravel())
            src = cand[rho.ravel()[cand] > a.peak_delta]
            ns = src.size
            mass = ndimage.uniform_filter(rho.astype(np.float64), size=3, mode="wrap").ravel()[src] * 27
            del mx, cand
            fg = sfft.rfftn(g / g.mean() - 1.0, workers=a.workers)
            Pg = B.power(fg)
            zj, lam_ion, lam_rnd, xv = mfp(js, z)
            h = hubble(r)
            lam_c = lam_ion * (1 + zj) * h  # comoving Mpc/h
            res[lab] = dict(z=z, xv=xv, lam=lam_c, lam_rnd=lam_rnd * (1 + zj) * h, Pg=Pg)
            for w in ("mass", "num"):
                L = np.zeros(N ** 3, np.float64)
                np.add.at(L, src, mass if w == "mass" else 1.0)
                L = L.reshape(N, N, N)
                fl = sfft.rfftn(L / L.mean() - 1.0, workers=a.workers)
                Pl = B.power(fl)
                X = B.power(fg, fl)
                res[lab][w] = (np.sqrt(Pg / Pl), X / np.sqrt(Pg * Pl))
                del L, fl
            del g, rho, fg
            pr(f"# step {st} {lab}: {d} z {z:.4f} x_V {xv:.4f}; N_src {ns}; lambda_mfp (ionised starts) "
               f"{lam_c:.2f} cMpc/h ({lam_c / a.box:.3f} box), random starts {res[lab]['lam_rnd']:.2f} cMpc/h")
        labs = [lab for _, lab, _ in runs]
        pr(f"\n## step {st} z {zr:.3f}: T = (P_Gamma/P_L)^1/2 (L: mass-weighted | number), "
           f"T/T(k_min), r(Gamma, L), arctan(k lam)/(k lam)")
        pr("#  k[h/Mpc] " + " ".join(f"{'T_m:' + l:>9s} {'T_n:' + l:>9s} {'Tm/T0:' + l:>10s} {'r_m:' + l:>8s} "
                                     f"{'r_n:' + l:>8s} {'exact:' + l:>9s}" for l in labs))
        i0 = int(np.flatnonzero(B.cnt > 0)[0])
        for i, k in enumerate(B.k):
            if B.cnt[i] == 0:
                continue
            row = []
            for l in labs:
                e = res[l]
                tm, rm = e["mass"][0][i], e["mass"][1][i]
                tn, rn = e["num"][0][i], e["num"][1][i]
                kl = k * e["lam"]
                row.append(f"{tm:9.4f} {tn:9.4f} {tm / e['mass'][0][i0]:10.4f} {rm:8.3f} {rn:8.3f} "
                           f"{np.arctan(kl) / kl:9.4f}")
            pr(f"{k:10.3f} " + " ".join(row))
    out.close()


if __name__ == "__main__":
    main()
