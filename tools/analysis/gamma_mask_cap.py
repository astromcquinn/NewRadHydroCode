#!/usr/bin/env python3
"""Is the Gamma_HI power excess of the moment closures at k ~ 5 h/Mpc the
source cells, the proximity zones, or neither? (docs/rt_moments.md 8.6, A1)

usage: gamma_mask_cap.py REF_RUN:label RUN:label [...] --outdir DIR
                         [--steps 69,84,93,99] [--box 40] [--nbins 24]
                         [--field gammaHI_avg] [--cap 3] [--spike 2] [--peak-delta 10]
                         [--workers N]

REF_RUN (the rays) and RUN are run directories (<problem>/run, the layout of
gamma_power.py). For each reference step of --steps, every other run's
snapshot nearest in redshift is used (the runs write the same zdata1
outputs). P(k) of delta_Gamma = Gamma/<Gamma>_V - 1 with the binning of
gamma_power.py / rt_suite_figures.py (the paper's Fig. 2), three ways, each
applied alike to every run including the reference:

  asis   the field as written (--field, by default the substep-averaged rate)
  spike  every single-cell spike replaced by the mean of its 6 face
         neighbours: the local maxima (3^3, periodic) of the run's own Gamma
         above --spike times that mean. The photons a moment method injects
         into a source cell in the current substep (the injection floor,
         docs/rt_moments.md section 7) make such spikes; the rays' source
         cells, which their rays leave half a cell from the centre, less so
  peaks  the halo stand-ins replaced by the mean of their non-masked face
         neighbours (as $RUNS/cone_blend_256/gamma_masked.py): the local
         maxima (3^3, periodic) of the run's own rhog above --peak-delta
         times the mean
  cap    Gamma capped at --cap times its volume mean <Gamma>_V before
         delta is formed (rt_suite_figures.py --clip): the proximity zones'
         peaks drop out

Why stand-ins: the snapshots hold no source list (no Luv.bin, no
sources.txt) and the checkpoints with the halos are gone. The run.log
'Nsource' counts source halos, not cells (halo_sources.cpp), and at 512^3
in 40 Mpc/h it exceeds the number of density maxima in the box, so the
source cells cannot be taken as "the N_src densest maxima". Printed per run
and mask: the masked cells, their share of the volume and of sum Gamma, and
the median Gamma(cell)/(mean of the 6 neighbours) over them.

Writes DIR/gamma_mask_cap.txt: per step and variant, per k bin, the ratio
P_Gamma(run)/P_Gamma(ref) and r(k) = P_ab/sqrt(P_aa P_bb) for every run.
"""
import argparse
import os
import sys

import numpy as np
import scipy.fft as sfft
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gamma_power as gp  # noqa: E402


def load(d, name):
    a = np.fromfile(os.path.join(d, name + ".bin"), dtype=np.float32)
    N = round(a.size ** (1 / 3))
    return a.reshape(N, N, N)


def neighbour_mean(g, m):
    s = np.zeros(g.shape, np.float64)
    c = np.zeros(g.shape, np.float64)
    for ax in range(3):
        for sh in (-1, 1):
            mn = np.roll(m, sh, axis=ax)
            s += np.where(mn, 0.0, np.roll(g, sh, axis=ax))
            c += ~mn
    return s / np.maximum(c, 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--steps", default="69,84,93,99")
    ap.add_argument("--box", type=float, default=40.0)
    ap.add_argument("--nbins", type=int, default=24)
    ap.add_argument("--field", default="gammaHI_avg")
    ap.add_argument("--cap", type=float, default=3.0)
    ap.add_argument("--spike", type=float, default=2.0)
    ap.add_argument("--peak-delta", type=float, default=10.0)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    a = ap.parse_args()
    gp.FIELD = a.field
    runs = [(s.rsplit(":", 1)[0], s.rsplit(":", 1)[1]) for s in a.runs]
    os.makedirs(a.outdir, exist_ok=True)
    out = open(os.path.join(a.outdir, "gamma_mask_cap.txt"), "w")

    def pr(s):
        print(s, flush=True)
        out.write(s + "\n")
        out.flush()

    pr("# gamma_mask_cap.py " + " ".join(sys.argv[1:]))
    pr(f"# P_Gamma(run)/P_Gamma({runs[0][1]}) and r(k); variants asis | spike (single-cell Gamma maxima "
       f"> {a.spike}x the 6-neighbour mean -> that mean) | peaks (rhog maxima > {a.peak_delta} -> mean of "
       f"non-masked neighbours) | cap (Gamma <= {a.cap} <Gamma>_V)")
    snaps = {lab: gp.snapshots(r, 0.0, 99.0) for r, lab in runs}
    B = None
    for st in [int(x) for x in a.steps.split(",")]:
        dref = os.path.join(runs[0][0], f"step_{st:04d}")
        zr = gp.read_z(dref)
        VAR = ("asis", "spike", "peaks", "cap")
        F = {v: {} for v in VAR}
        for r, lab in runs:
            d = dref if lab == runs[0][1] else snaps[lab][min(snaps[lab], key=lambda z: abs(z - zr))]
            z = gp.read_z(d)
            g = gp.load_gamma(d).astype(np.float64)
            N = g.shape[0]
            if B is None:
                B = gp.Binner(N, a.box, a.nbins)
            mean = g.mean()
            none = np.zeros(g.shape, bool)
            nb6 = neighbour_mean(g, none)
            con = g / np.maximum(nb6, 1e-300)
            ms = (g == ndimage.maximum_filter(g, size=3, mode="wrap")) & (con > a.spike)
            gs = np.where(ms, nb6, g)
            del nb6
            rho = load(d, "rhog")
            mp = (rho == ndimage.maximum_filter(rho, size=3, mode="wrap")) & (rho > a.peak_delta)
            del rho
            gpk = np.where(mp, neighbour_mean(g, mp), g)
            gc = np.minimum(g, a.cap * mean)
            info = []
            for nm, m in (("spike", ms), ("peaks", mp)):
                info.append(f"{nm}: {int(m.sum())} cells, share of sum Gamma {g[m].sum() / g.sum():.4f}, "
                            f"median Gamma/6-neighbour mean {np.median(con[m]) if m.any() else float('nan'):.2f}")
            pr(f"# step {st} (z_ref {zr:.4f}) {lab}: {d} z {z:.4f}; <Gamma>_V {mean:.4e}; " + "; ".join(info)
               + f"; overlap spike&peaks {int((ms & mp).sum())}; cap: {(g > a.cap * mean).mean():.4f} of the volume, "
               f"{np.maximum(g - a.cap * mean, 0).sum() / g.sum():.4f} of sum Gamma")
            for v, f in (("asis", g), ("spike", gs), ("peaks", gpk), ("cap", gc)):
                F[v][lab] = sfft.rfftn(f / f.mean() - 1.0, workers=a.workers).astype(np.complex64)
            del g, gs, gpk, gc, con, ms, mp
        ref = runs[0][1]
        labs = [lab for _, lab in runs[1:]]
        for v in VAR:
            P = {lab: B.power(F[v][lab]) for lab in F[v]}
            X = {lab: B.power(F[v][lab], F[v][ref]) for lab in labs}
            pr(f"\n## step {st} z {zr:.3f} variant {v}")
            pr("#  k[h/Mpc] " + " ".join(f"{'P/P_' + ref + ':' + l:>14s}" for l in labs) + " | "
               + " ".join(f"{'r:' + l:>8s}" for l in labs) + f" | Delta2_{ref}")
            for i, k in enumerate(B.k):
                if B.cnt[i] == 0:
                    continue
                rat = [P[l][i] / P[ref][i] for l in labs]
                rr = [X[l][i] / np.sqrt(P[l][i] * P[ref][i]) for l in labs]
                d2 = k ** 3 * P[ref][i] / (2 * np.pi ** 2)
                pr(f"{k:10.3f} " + " ".join(f"{x:14.3f}" for x in rat) + " | "
                   + " ".join(f"{x:8.3f}" for x in rr) + f" | {d2:.3e}")
        del F
    out.close()


if __name__ == "__main__":
    main()
