#!/usr/bin/env python3
"""Real-space tau_eff scatter from the movie slabs, finer in z than the
snapshots (Part A4 of docs/taueff_21cmfast.md).

    taueff_slabs.py RUN[:label] [...] --out DIR [--zmin 4.95] [--zmax 6.25]
                    [--seg 50] [--jobs 4] [--calib]

The movie slabs (movie/frame_NNNNN_<field>.slab, movie/slab.txt;
docs/movies.md "Slab frames") hold, every dz = 0.025, `thick` cells through
the box centre normal to each axis, with x_HI, T and rho_gas but no
velocities. Lya skewers run along both in-plane axes of each slab (3 slabs x
2 directions x N x thick skewers, pixel = cell), in real space: thermal
broadening only (igm_stats.lya_tau with v = 0). Segments as igm_stats.py
--taueff (box length, and seg Mpc/h joined from a second skewer), tau
rescaled by A to the Bosman+22 mean flux at the frame's z.

The slabs see only 3 x thick planes: the skewers of a slab share their
large-scale structure (they are not independent sightlines), and real space
lacks the peculiar-velocity redistribution. Hence --calib: at each snapshot
the same slab skewers are cut from the 3D cubes and their tau_eff computed
in real space and in redshift space (Gas/u_gas velocities), and against the
redshift-space values of all 3 x nskew random skewers of igm_stats.py
--taueff (that JSON, if given with --rs-json). The ratio of the widths gives
the real-space offset, reported before any use of the slab curve.

Writes DIR/<label>_slabs.txt (per frame: z, x_V, A, the 5/16/50/84/95 %
points of the rescaled seg tau_eff, its 84-16 width, the same raw) and with
--calib DIR/<label>_calib.txt.
"""
import argparse
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import igm_stats as ig  # noqa: E402

QS = (5, 16, 50, 84, 95)


def slab_geom(mdir):
    g = {}
    for line in open(os.path.join(mdir, "slab.txt")):
        if "=" in line and not line.startswith("#"):
            k, v = [x.strip() for x in line.split("=", 1)]
            g[k] = v
    return int(g["N"]), int(g["thick"]), g["axes"], int(g["start"])


def frames(mdir):
    out = []
    for line in open(os.path.join(mdir, "frames.txt")):
        if line.startswith("#"):
            continue
        t = line.split()
        out.append((int(t[0]), int(t[1]), float(t[3]), float(t[4])))
    return out


def slab_blocks(path, N, t, axes):
    """{axis: array [k, j, i]} with the slab axis of length t."""
    a = np.fromfile(path, dtype=np.float32)
    out, off = {}, 0
    n = t * N * N
    for ax in axes:
        blk = a[off:off + n]
        off += n
        shape = {"x": (N, N, t), "y": (N, t, N), "z": (t, N, N)}[ax]  # C order of Fortran (nx,ny,nz)
        out[ax] = blk.reshape(shape)
    return out


def skewers_from_blocks(B):
    """All in-plane skewers, as rows (M, N): for the x slab along y and z, for
    the y slab along x and z, for the z slab along x and y."""
    rows = []
    for ax, b in B.items():
        if ax == "x":    # [k, j, i], i in slab
            rows += [b.transpose(0, 2, 1).reshape(-1, b.shape[1]),   # along j
                     b.transpose(1, 2, 0).reshape(-1, b.shape[0])]   # along k
        elif ax == "y":  # j in slab
            rows += [b.reshape(-1, b.shape[2]),                       # along i
                     b.transpose(1, 2, 0).reshape(-1, b.shape[0])]   # along k
        else:            # k in slab
            rows += [b.reshape(-1, b.shape[2]),                       # along i
                     b.transpose(0, 2, 1).reshape(-1, b.shape[1])]   # along j
    return np.concatenate(rows)


def cube_slab_blocks(cube, N, t, axes, s0):
    """The movie geometry cut from a 3D cube [k, j, i]."""
    out = {}
    for ax in axes:
        if ax == "x":
            out[ax] = cube[:, :, s0:s0 + t]
        elif ax == "y":
            out[ax] = cube[:, s0:s0 + t, :]
        else:
            out[ax] = cube[s0:s0 + t]
    return out


def velocity_rows(gas, N, t, axes, s0, vel_unit):
    """Peculiar velocity along each skewer of skewers_from_blocks()."""
    rows = []
    for ax in axes:
        if ax == "x":
            sl = gas[:, :, s0:s0 + t]
            vy = sl[..., 2] / sl[..., 0]
            vz = sl[..., 3] / sl[..., 0]
            rows += [vy.transpose(0, 2, 1).reshape(-1, N), vz.transpose(1, 2, 0).reshape(-1, N)]
        elif ax == "y":
            sl = gas[:, s0:s0 + t, :]
            vx = sl[..., 1] / sl[..., 0]
            vz = sl[..., 3] / sl[..., 0]
            rows += [vx.reshape(-1, N), vz.transpose(1, 2, 0).reshape(-1, N)]
        else:
            sl = gas[s0:s0 + t]
            vx = sl[..., 1] / sl[..., 0]
            vy = sl[..., 2] / sl[..., 0]
            rows += [vx.reshape(-1, N), vy.transpose(0, 2, 1).reshape(-1, N)]
    return np.concatenate(rows).astype(np.float64) * vel_unit


class Geo:
    """igm_stats.Snapshot's unit conversions at a redshift z."""

    def __init__(self, run, z):
        a = 1.0 / (1.0 + z)
        h = ig.HUBBLE0
        self.dl = a * run.box / h / run.N * ig.MPC_CM
        self.nH_avg = ig.RHOC_CGS * ig.OMEGAM * h * h / a ** 3 * ig.OMEGAB / ig.OMEGAM * ig.XHY / ig.MH_CGS
        time_unit = 2 / (3 * ig.H0_CGS * h * math.sqrt(ig.OMEGAM)) * a * a
        self.vel_unit = self.dl / time_unit
        self.H = ig.H0_CGS * h * math.sqrt(ig.OMEGAM / a ** 3 + (1 - ig.OMEGAM) + ig.OMEGAR / a ** 4)


def tau_rows(x, rho, T, v, geo, chunk=2048):
    out = []
    for c in range(0, x.shape[0], chunk):
        s = slice(c, c + chunk)
        nHI = x[s].astype(np.float64) * rho[s] * geo.nH_avg
        vv = np.zeros(nHI.shape) if v is None else v[s]
        out.append(ig.lya_tau(nHI, T[s], vv, geo.dl, geo.H)[0])
    return np.concatenate(out)


def seg_stats(tau, box, seg, Fo, seed):
    A = ig.rescale_A(tau, Fo)
    S = ig.taueff_segments(tau, box, seg, A, np.random.default_rng(seed))
    res = {"A": A, "p_raw": np.full(len(QS), np.nan), "p_A": np.full(len(QS), np.nan)}
    for tag in ("raw", "A"):
        k = "F_seg_" + tag
        if k in S:
            F = np.asarray(S[k])
            te = np.where(F > 0, -np.log(np.maximum(F, 1e-300)), np.inf)
            res["p_" + tag] = np.percentile(te, QS)
    return res


def do_frame(args):
    run_path, fr, z, xv, seg, data = args
    run = ig.Run(run_path)
    mdir = os.path.join(run.out, "movie")
    N, t, axes, s0 = slab_geom(mdir)
    f = {n: skewers_from_blocks(slab_blocks(os.path.join(mdir, f"frame_{fr:05d}_{n}.slab"), N, t, axes))
         for n in ("fHI", "rhog", "temp")}
    geo = Geo(run, z)
    tau = tau_rows(f["fHI"], f["rhog"], f["temp"], None, geo)
    Fo = ig.obs_interp(ig.read_obs(data), "mean_flux", z)
    r = seg_stats(tau, run.box, seg, Fo, 5 + fr)
    return fr, z, xv, r, tau.shape[0]


def do_calib(args):
    run_path, step, seg, data = args
    run = ig.Run(run_path)
    snap = ig.Snapshot(run, step, run.stats(step))
    mdir = os.path.join(run.out, "movie")
    N, t, axes, s0 = slab_geom(mdir)
    rows = {}
    for n in ("fHI", "rhog", "temp"):
        rows[n] = skewers_from_blocks(cube_slab_blocks(np.asarray(snap.field(n)), N, t, axes, s0))
    v = velocity_rows(snap.gas(), N, t, axes, s0, snap.vel_unit)
    geo = Geo(run, snap.z)
    Fo = ig.obs_interp(ig.read_obs(data), "mean_flux", snap.z)
    out = {}
    for name, vv in (("real", None), ("redshift", v)):
        tau = tau_rows(rows["fHI"], rows["rhog"], rows["temp"], vv, geo)
        out[name] = seg_stats(tau, run.box, seg, Fo, 5 + step)
    return step, snap.z, out


def fmt(p):
    return " ".join("nan" if not np.isfinite(x) else f"{x:.4g}" for x in p)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--zmin", type=float, default=4.95)
    ap.add_argument("--zmax", type=float, default=6.25)
    ap.add_argument("--seg", type=float, default=50.0)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--calib", action="store_true")
    ap.add_argument("--data", default=ig.DATA_DEFAULT)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for spec in a.runs:
        path, _, label = spec.partition(":")
        run = ig.Run(path)
        label = label or run.name
        mdir = os.path.join(run.out, "movie")
        fl = [f for f in frames(mdir) if a.zmin <= f[2] <= a.zmax]
        # one frame per redshift (a restart can repeat one)
        seen, todo = set(), []
        for fr, ist, z, xv in fl:
            if round(z, 3) not in seen:
                seen.add(round(z, 3))
                todo.append((run.top, fr, z, xv, a.seg, a.data))
        with ProcessPoolExecutor(max_workers=a.jobs) as ex:
            res = sorted(ex.map(do_frame, todo), key=lambda r: -r[1])
        with open(os.path.join(a.out, f"{label}_slabs.txt"), "w") as f:
            f.write(f"# {run.top}: real-space tau_eff over {a.seg} Mpc/h segments from the movie slabs "
                    "(taueff_slabs.py, docs/taueff_21cmfast.md)\n")
            f.write("# frame z x_V A nskew then tau_eff rescaled at " + ",".join(map(str, QS)) +
                    " %, width(84-16), then raw at the same points\n")
            for fr, z, xv, r, m in res:
                pA = r["p_A"]
                f.write(f"{fr} {z:.4f} {xv:.5f} {r['A']:.4g} {m} {fmt(pA)} {pA[3] - pA[1]:.4g} "
                        f"{fmt(r['p_raw'])}\n")
        print(f"{label}: {len(res)} frames", flush=True)
        if a.calib:
            steps = [s.step for s in run.snapshots(a.zmin - 0.05, a.zmax + 0.05)]
            with ProcessPoolExecutor(max_workers=a.jobs) as ex:
                cal = sorted(ex.map(do_calib, [(run.top, s, a.seg, a.data) for s in steps]),
                             key=lambda r: -r[1])
            with open(os.path.join(a.out, f"{label}_calib.txt"), "w") as f:
                f.write("# the movie-slab skewers cut from the snapshot cubes: tau_eff (rescaled, "
                        f"{a.seg} Mpc/h) in real space and in redshift space\n")
                f.write("# step z A_real A_rs then real 5..95 %, width; redshift-space 5..95 %, width; "
                        "width ratio rs/real; median difference rs - real\n")
                for st, z, o in cal:
                    pr, ps = o["real"]["p_A"], o["redshift"]["p_A"]
                    wr, ws = pr[3] - pr[1], ps[3] - ps[1]
                    f.write(f"{st} {z:.4f} {o['real']['A']:.4g} {o['redshift']['A']:.4g} {fmt(pr)} {wr:.4g} "
                            f"{fmt(ps)} {ws:.4g} {ws / wr if wr > 0 else float('nan'):.4g} "
                            f"{ps[2] - pr[2]:.4g}\n")
            print(f"{label}: calibration at {len(cal)} snapshots", flush=True)


if __name__ == "__main__":
    main()
