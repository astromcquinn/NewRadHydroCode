#!/usr/bin/env python3
"""Spectra of the convergence series (tools/hyak/rt_convergence.sh,
docs/rt_moments.md section 8) at matched ionised fraction.

usage: rt_convergence_cmp.py --runs DIR --out FILE.json [--res 128,256,512]
                             [--rayf 8,15,32,64,120] [--box 20] [--workers 16]

gamma_power.py pairs snapshots by redshift; here every run's snapshot is
paired with the reference's snapshot nearest in x_HII,V (the runs differ in
their reionization history, between resolutions by a lot); after overlap
(x_HII,V > 0.995 in both) with the one nearest in redshift.

Groups: for each N the runs rtconv_N<N>_rays (reference), _rays_m1,
_rays_m1_inj0, _rays_otvet and _rays_m1_o2 (moment_order = 2); the ray budgets rtconv_N256_rays_f<f> (reference: the largest
f present; f = 15 is rtconv_N256_rays).

Per snapshot (step_NNNN with gammaHI_avg.bin, fHI.bin, rhog.bin; z < 7.5):
z, x_HII,V, x_HII,M, <Gamma>_V, the median Gamma of ionised cells
(x_HII > 0.5), the statistics of the cells above 100x the mean gas density
(count, their volume- and mass-weighted neutral fraction, mean and median
Gamma), and P(k) of delta_Gamma = Gamma/<Gamma>_V - 1 and of x_HII (the field
itself) in k bins that are the same for every N (log bins from the box's
fundamental, 4 per factor of 2), so spectra of different N
compare at fixed k.

Per pair (run snapshot, reference snapshot of the same grid): the cross
spectra, hence r(k), of both fields.
"""
import argparse
import glob
import json
import os

import numpy as np
import scipy.fft as sfft


XLATE = 0.995  # x_HII,V above which snapshots are paired by redshift


def stat(d, key):
    for line in open(os.path.join(d, "stats.txt")):
        k, _, v = line.partition("=")
        if k.strip() == key:
            return float(v)
    return None


def snapshots(run, field, zmax):
    out = []
    for d in sorted(glob.glob(os.path.join(run, "step_[0-9][0-9][0-9][0-9]"))):
        if all(os.path.exists(os.path.join(d, f + ".bin")) for f in (field, "fHI", "rhog")):
            z = stat(d, "redshift")
            if z is not None and z < zmax:
                out.append({"dir": d, "z": z, "xV_stats": stat(d, "xHII_V")})
    return out


class Binner:
    """Fixed |k| bins [h/Mpc] of the rfftn modes of an N^3 grid."""

    def __init__(self, N, box, edges):
        kf = 2 * np.pi / box
        n = np.fft.fftfreq(N, 1.0 / N)
        nz = np.arange(N // 2 + 1)
        k2 = (n[:, None, None] ** 2 + n[None, :, None] ** 2) + nz[None, None, :] ** 2
        kmag = np.sqrt(k2.astype(np.float64)).ravel() * kf
        self.nb = len(edges) - 1
        idx = np.digitize(kmag, edges) - 1
        w = np.full((N, N, N // 2 + 1), 2.0)
        w[:, :, 0] = 1.0
        if N % 2 == 0:
            w[:, :, N // 2] = 1.0
        w = w.ravel()
        ok = (idx >= 0) & (idx < self.nb) & (kmag > 0)
        self.idx, self.w, self.ok = idx[ok], w[ok], ok
        self.cnt = np.bincount(self.idx, weights=self.w, minlength=self.nb)
        self.k = np.bincount(self.idx, weights=self.w * kmag[ok], minlength=self.nb) / np.maximum(self.cnt, 1)
        self.norm = box ** 3 / float(N) ** 6

    def power(self, fa, fb=None):
        p = (fa * np.conj(fb)).real if fb is not None else (fa.real ** 2 + fa.imag ** 2)
        p = p.ravel()[self.ok].astype(np.float64)
        return np.bincount(self.idx, weights=self.w * p, minlength=self.nb) / np.maximum(self.cnt, 1) * self.norm


def load(d, name):
    g = np.fromfile(os.path.join(d, name + ".bin"), dtype=np.float32)
    N = round(g.size ** (1 / 3))
    return g.reshape(N, N, N)


def measure(sn, field, binner, workers):
    """Scalars and auto spectra of one snapshot; returns the two transforms."""
    d = sn["dir"]
    g = load(d, field).astype(np.float64)
    x = 1.0 - load(d, "fHI").astype(np.float64)
    rho = load(d, "rhog").astype(np.float64)
    rho /= rho.mean()
    gm = g.mean()
    ion = x > 0.5
    dense = rho > 100.0
    sn.update(xV=x.mean(), xM=(x * rho).mean(), gamma_V=gm,
              gamma_med_ion=float(np.median(g[ion])) if ion.any() else None,
              gamma_mean_ion=float(g[ion].mean()) if ion.any() else None,
              N=g.shape[0])
    nd = int(dense.sum())
    sn["dense"] = {"n": nd}
    if nd:
        sn["dense"].update(
            xHI_V=float((1 - x[dense]).mean()),
            xHI_M=float(((1 - x[dense]) * rho[dense]).sum() / rho[dense].sum()),
            xHI_med=float(np.median(1 - x[dense])),
            gamma_mean=float(g[dense].mean()), gamma_med=float(np.median(g[dense])),
            mass_frac=float(rho[dense].sum() / rho.size))
    fg = sfft.rfftn(g / gm - 1.0, workers=workers)
    fx = sfft.rfftn(x, workers=workers)
    sn["k"] = binner.k.tolist()
    sn["nmodes"] = binner.cnt.tolist()
    sn["P_gamma"] = binner.power(fg).tolist()
    sn["P_x"] = binner.power(fx).tolist()
    return fg, fx


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--res", default="128,256,512")
    ap.add_argument("--rayf", default="8,15,32,64,120")
    ap.add_argument("--box", type=float, default=20.0)
    ap.add_argument("--field", default="gammaHI_avg")
    ap.add_argument("--zmax", type=float, default=7.5)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    a = ap.parse_args()
    kf = 2 * np.pi / a.box
    # 4 bins per factor of 2 from the fundamental to beyond the 512^3 corner
    edges = kf * 0.999 * 2.0 ** (np.arange(0, 41) / 4.0)

    groups = []
    for n in [int(t) for t in a.res.split(",") if t]:
        b = os.path.join(a.runs, "rtconv_N%d_rays" % n)
        runs = [("rays", b), ("M1", b + "_m1"), ("M1_inject0", b + "_m1_inj0"), ("OTVET", b + "_otvet"),
                ("M1_order2", b + "_m1_o2")]
        # a finished run with the largest ray pool is the reference where the
        # series' f = 15 pool was full (128^3; at 256^3 see the group "rays")
        if n != 256 and os.path.exists(os.path.join(b + "_f120", "run", "DONE")):
            runs.insert(0, ("rays_f120", b + "_f120"))
        groups.append(("N%d" % n, runs))
    fs = sorted([int(t) for t in a.rayf.split(",") if t], reverse=True)
    rr = []
    for f in fs:
        p = os.path.join(a.runs, "rtconv_N256_rays" + ("" if f == 15 else "_f%d" % f))
        if os.path.exists(os.path.join(p, "run", "DONE")):
            rr.append(("f%d" % f, p))
    if len(rr) > 1:
        groups.append(("rays", rr))

    res = {"box": a.box, "field": a.field, "edges": edges.tolist(), "groups": {}}
    for gname, runs in groups:
        runs = [(lab, p) for lab, p in runs if glob.glob(os.path.join(p, "run", "step_*", a.field + ".bin"))]
        if not runs:
            continue
        out = {"reference": runs[0][0], "runs": {}, "pairs": []}
        snaps = {lab: snapshots(os.path.join(p, "run"), a.field, a.zmax) for lab, p in runs}
        binner = None
        ref = snaps[runs[0][0]]
        # the reference first; its transforms stay in memory one at a time:
        # loop over reference snapshots, pairing each other run's snapshot
        # nearest in x_HII,V (from stats.txt) with it
        done = {lab: set() for lab in snaps}
        for i, sr in enumerate(ref):
            if binner is None:
                N = round((os.path.getsize(os.path.join(sr["dir"], "fHI.bin")) // 4) ** (1 / 3))
                binner = Binner(N, a.box, edges)
            fgr, fxr = measure(sr, a.field, binner, a.workers)
            for lab, _ in runs[1:]:
                for j, so in enumerate(snaps[lab]):
                    # so is paired with the reference snapshot nearest in x_V
                    # (by redshift among the reference's snapshots after
                    # overlap once this one is past it too)
                    late = [q for q in range(len(ref)) if ref[q]["xV_stats"] > XLATE]
                    if so["xV_stats"] > XLATE and late:
                        near = min(late, key=lambda q: abs(ref[q]["z"] - so["z"]))
                    else:
                        near = min(range(len(ref)), key=lambda q: abs(ref[q]["xV_stats"] - so["xV_stats"]))
                    if near != i:
                        continue
                    fg, fx = measure(so, a.field, binner, a.workers)
                    done[lab].add(j)
                    out["pairs"].append({"run": lab, "i": j, "ref_i": i,
                                         "C_gamma": binner.power(fg, fgr).tolist(),
                                         "C_x": binner.power(fx, fxr).tolist()})
                    del fg, fx
            del fgr, fxr
            print(gname, "reference snapshot z = %.3f done" % sr["z"], flush=True)
        for lab in snaps:
            for sn in snaps[lab]:
                sn["dir"] = os.path.relpath(sn["dir"], a.runs)
            out["runs"][lab] = snaps[lab]
        res["groups"][gname] = out
        json.dump(res, open(a.out, "w"))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
