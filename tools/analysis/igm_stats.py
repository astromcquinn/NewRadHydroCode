#!/usr/bin/env python3
"""IGM statistics of apps/radhydro runs, to compare with the z = 5-6 data of
data/igm/igm_z5_6.txt (docs/igm_diagnostics.md).

    igm_stats.py RUN [--out FILE.json] [--zmin 4.5] [--zmax 7.1] [--steps 43,48]
                 [--nlos 4096] [--rmax 40] [--nskew 4096] [--jobs 4] [--seed 1]

RUN is a run directory (params.txt and run/) or its run/ subdirectory. For
every step with fields (run/step_NNNN/{fHI,rhog,temp,gamma}.bin, written at
the zdata1 outputs; params `fields = fHI,fHeI,fHeII,temp,rhog,gamma`) between
zmin and zmax it measures:

1. Gamma_HI, the HI photoionization rate [s^-1] of the RT's last substep of
   the hydro step: the first of the three components of rt_grid%gamma
   (HI, HeI, HeII), summed over the five frequency groups (rt_tools.cpp,
   photo_cell: gamma = n_gam sigma_HI(nu) len_rt/dtsec_rt, reset every
   substep in rt_step_1). Volume average, and averages (volume, mass,
   median) over the cells with x_HI below 0.5, 0.1, 1e-2, 1e-3.
2. The mean free path of 912 A photons: tau_912 = sigma_912 int n_HI dl
   along random sightlines (random start points and isotropic directions,
   periodic, longer than the box), <F(r)> = <exp(-tau(r))>; lambda_e is the
   r of <F> = 1/e, lambda_fit the best fit of exp(-r/lambda) to <F(r)>.
   Three sets of start points: random (volume weighted), ionized cells
   (x_HI < 0.01), and source cells (rt_grid%Luv > 0, from RT/rt_grid.<z>,
   luminosity weighted; the source's own cell is skipped).
3. The Lya forest mean flux along skewers parallel to the three axes, in
   redshift space: Hubble flow plus the line-of-sight peculiar velocity
   (Gas/u_gas.<z>: momentum / density), the gas of each cell spread
   between its faces in redshift space (face velocities interpolated),
   thermal Doppler profile (no damping wings), integrated over each pixel
   (pixel = cell).
   <F> = <exp(-tau)> over all pixels, the factor A for which
   <exp(-A tau)> equals the Bosman et al. 2022 mean flux at the same z, and
   the transmission-weighted Gamma_HI of the gas (each cell weighted by the
   flux of the pixel its gas lands in).
4. ndot_ion, the escaped ionizing photon rate of the sources at every hydro
   step, from step_NNNN/stats.txt: Luv is sum(source L)/V_box [erg/s/Hz per
   comoving Mpc^3] (rt_tools.cpp calc_halo), and the RT emits
   fResc*sum(pop2_init) photons per erg/Hz (ray_source; halo_sources.hpp:
   the source L already carries xi_ion and f_esc/fResc), so
   ndot = 0.2 * 10^25.2 * 1.0001 * Luv  [photons s^-1 cMpc^-3].

Units follow cosmology.cpp calc_sim_var at the snapshot's scale factor:
proper cell size dl = a box/h/N Mpc; n_H = rhog <n_H>, <n_H> =
rho_crit Omega_b h^2 X / m_H / a^3 (rhog is the density in units of the mean
baryon density); velocity unit = dl / time_unit, time_unit =
2 a^2 / (3 H0 h sqrt(Omega_m)) (u_gas(2:4)/u_gas(1) are peculiar velocities
in this unit, component 1+d along grid axis d). Fields are float32 in
Fortran order, read in slabs of k (memory ~ 1.5 x 4 N^3 bytes per job).

Output: a JSON file (default igm_<run name>.json) read by igm_plot.py.
"""
import argparse
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from box_lf import read_params  # noqa: E402


def f32(x):
    return float(np.float32(x))


# include/radray/constants.hpp (real(4) literals widened, as in the code)
OMEGAM, OMEGAB, HUBBLE0, XHY = f32(0.27), f32(0.045), f32(0.70), f32(0.76)
OMEGAR = float.fromhex("0x1.6559c34a8bc79p-14")
RHOC_CGS, MH_CGS, MPC_CM = 1.87890e-29, 1.67223e-24, 3.08560e24
H0_CGS, K_CGS = 3.24086e-18, 1.38066e-16
POP2_INIT = [float.fromhex(x) for x in ("0x1.5f23100000000p+82", "0x1.edb0780000000p+81",
                                         "0x1.5120020000000p+81", "0x1.102f000000000p+80",
                                         "0x1.21f89e0000000p+76")]
FRESC = f32(0.2)
PHOT_PER_LUV = FRESC * sum(POP2_INIT)  # photons per (erg/s/Hz) emitted by the RT

# Standard atomic constants (not taken from the code)
SIGMA_912 = 6.30e-18          # cm^2, HI photoionization cross-section at 13.6 eV
C_CGS = 2.99792458e10         # cm/s
E_ESU = 4.80320471e-10        # electron charge [esu]
ME_G = 9.1093837e-28          # electron mass [g]
F_LYA = 0.4164                # Lya oscillator strength
LAMBDA_LYA = 1215.67e-8       # Lya wavelength [cm]
# tau = K_LYA int n_HI phi_v dl with phi_v the line profile in velocity
# [s/cm]: K_LYA = pi e^2 f lambda / (m_e c) = sigma_Lya c = 1.343e-7 cm^3/s
# (sigma_Lya = 4.48e-18 cm^2)
K_LYA = math.pi * E_ESU ** 2 * F_LYA * LAMBDA_LYA / (ME_G * C_CGS)

XHI_THRESHOLDS = (0.5, 0.1, 1e-2, 1e-3)
DATA_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data",
                            "igm", "igm_z5_6.txt")


def ztag(z):
    """pm.cpp ztag: the suffix of the Gas/, RT/ and DM/ files."""
    return f"z=0{z:6.4f}" if z < 10 else f"z={z:7.4f}"


# ---------------------------------------------------------------------------
# Observations
# ---------------------------------------------------------------------------

def read_obs(path=DATA_DEFAULT):
    """Blocks of data/igm/igm_z5_6.txt: name -> dict of arrays z, v, ep, em, ref."""
    blocks, cur = {}, None
    for line in open(path):
        s = line.strip()
        if s.startswith("# block"):
            cur = s.split()[2]
            blocks[cur] = []
        elif s and not s.startswith("#") and cur is not None:
            t = s.split()
            blocks[cur].append((float(t[0]), float(t[1]), float(t[2]), float(t[3]), t[4]))
    out = {}
    for k, rows in blocks.items():
        out[k] = {"z": np.array([r[0] for r in rows]), "v": np.array([r[1] for r in rows]),
                  "ep": np.array([r[2] for r in rows]), "em": np.array([r[3] for r in rows]),
                  "ref": [r[4] for r in rows]}
    return out


def obs_interp(obs, name, z, ref=None):
    """Linear interpolation of a block (only the rows of reference ref, if
    given) at z; nan outside its range."""
    b = obs.get(name)
    if b is None:
        return float("nan")
    m = np.ones(len(b["z"]), bool) if ref is None else np.array([r == ref for r in b["ref"]])
    zz, vv = b["z"][m], b["v"][m]
    if len(zz) == 0 or z < zz.min() or z > zz.max():
        return float("nan")
    o = np.argsort(zz)
    return float(np.interp(z, zz[o], vv[o]))


# ---------------------------------------------------------------------------
# Run and snapshot readers
# ---------------------------------------------------------------------------

def read_stats(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = [x.strip() for x in line.split("=", 1)]
            try:
                d[k] = float(v)
            except ValueError:
                d[k] = v
    return d


class Run:
    def __init__(self, path):
        path = os.path.abspath(path.rstrip("/"))
        if os.path.exists(os.path.join(path, "params.txt")) and os.path.isdir(
                os.path.join(path, "run")):
            self.top, self.out = path, os.path.join(path, "run")
        elif os.path.exists(os.path.join(os.path.dirname(path), "params.txt")):
            self.top, self.out = os.path.dirname(path), path
        else:
            raise SystemExit(f"{path}: no params.txt (expected RUN/params.txt and RUN/run/)")
        self.name = os.path.basename(self.top)
        self.p = read_params(os.path.join(self.top, "params.txt"))
        self.N = int(self.p.get("N", 32))
        self.box = float(self.p.get("box", 5))
        self.c_reduced = float(self.p.get("c_reduced", 1))
        self.f_esc = float(self.p["f_esc"]) if "f_esc" in self.p else None

    def steps(self):
        out = []
        for d in sorted(os.listdir(self.out)):
            if d.startswith("step_") and os.path.exists(os.path.join(self.out, d, "stats.txt")):
                out.append(int(d[5:]))
        return out

    def stats(self, step):
        return read_stats(os.path.join(self.out, f"step_{step:04d}", "stats.txt"))

    def snapshots(self, zmin=-1.0, zmax=1e9, steps=None):
        need = ("fHI", "rhog", "temp", "gamma")
        snaps = []
        for st in self.steps():
            if steps and st not in steps:
                continue
            d = os.path.join(self.out, f"step_{st:04d}")
            if not all(os.path.exists(os.path.join(d, n + ".bin")) for n in need):
                continue
            s = Snapshot(self, st, self.stats(st))
            if zmin <= s.z <= zmax:
                snaps.append(s)
        return snaps

    def ndot(self):
        """Escaped ionizing photon rate [s^-1 cMpc^-3] at every hydro step with sources."""
        z, nd, xv, ns = [], [], [], []
        for st in self.steps():
            s = self.stats(st)
            if s.get("Luv", 0) > 0:
                z.append(s["redshift"])
                nd.append(PHOT_PER_LUV * s["Luv"])
                xv.append(s.get("xHII_V", float("nan")))
                ns.append(s.get("Nsource", 0))
        return {"z": z, "ndot": nd, "xHII_V": xv, "Nsource": ns}


class Snapshot:
    def __init__(self, run, step, stats):
        self.run, self.step, self.st = run, step, stats
        self.N = run.N
        self.dir = os.path.join(run.out, f"step_{step:04d}")
        self.z = stats["redshift"]
        self.a = stats["scalefactor"]
        self.tag = ztag(self.z)
        a, h = self.a, HUBBLE0
        # cosmology.cpp calc_sim_var
        self.dl = a * run.box / h / self.N * MPC_CM            # proper cell [cm]
        self.dl_pmpc = self.dl / MPC_CM
        self.nH_avg = RHOC_CGS * OMEGAM * h * h / a ** 3 * OMEGAB / OMEGAM * XHY / MH_CGS
        self.time_unit = 2 / (3 * H0_CGS * h * math.sqrt(OMEGAM)) * a * a
        self.vel_unit = self.dl / self.time_unit               # cm/s
        self.H = H0_CGS * h * math.sqrt(OMEGAM / a ** 3 + (1 - OMEGAM) + OMEGAR / a ** 4)

    def field(self, name, ncomp=1):
        """Memory map of step_NNNN/<name>.bin as [k, j, i(, comp)]."""
        shape = (self.N,) * 3 + ((ncomp,) if ncomp > 1 else ())
        return np.memmap(os.path.join(self.dir, name + ".bin"), np.float32, "r", shape=shape)

    def gas(self):
        """Gas/u_gas.<z> as [k, j, i, 5]: rho, rho v (3), E (code units)."""
        p = os.path.join(self.run.out, "Gas", "u_gas." + self.tag)
        return np.memmap(p, np.float32, "r", shape=(self.N,) * 3 + (5,))

    def rt_grid(self):
        """RT/rt_grid.<z> as [k, j, i, 36] (apps/radhydro.cpp write_data_files)."""
        p = os.path.join(self.run.out, "RT", "rt_grid." + self.tag)
        return np.memmap(p, np.float32, "r", shape=(self.N,) * 3 + (36,))

    def slabs(self, bytes_per_cell=4):
        """k ranges of about 64 MB."""
        nk = max(1, int(64e6 // (bytes_per_cell * self.N * self.N)))
        for k0 in range(0, self.N, nk):
            yield k0, min(self.N, k0 + nk)


# ---------------------------------------------------------------------------
# 1. Gamma_HI
# ---------------------------------------------------------------------------

LG_BINS = np.linspace(-20, -9, 2201)  # log10 Gamma histogram for medians (0.005 dex)


def _hist_median(h):
    c = np.cumsum(h)
    if c[-1] == 0:
        return float("nan")
    i = int(np.searchsorted(c, 0.5 * c[-1]))
    return float(10 ** (0.5 * (LG_BINS[i] + LG_BINS[i + 1])))


def gamma_stats(snap, thresholds=XHI_THRESHOLDS):
    fHI, rho, gam = snap.field("fHI"), snap.field("rhog"), snap.field("gamma", 3)
    n = 0
    sg = 0.0
    acc = {t: {"n": 0, "g": 0.0, "gm": 0.0, "m": 0.0, "h": np.zeros(len(LG_BINS) - 1)}
           for t in thresholds}
    for k0, k1 in snap.slabs(20):
        x = np.asarray(fHI[k0:k1]).ravel()
        r = np.asarray(rho[k0:k1], dtype=np.float64).ravel()
        g = np.asarray(gam[k0:k1, :, :, 0], dtype=np.float64).ravel()
        n += x.size
        sg += g.sum()
        lg = np.log10(np.maximum(g, 1e-30))
        for t in thresholds:
            m = x < t
            a = acc[t]
            a["n"] += int(m.sum())
            a["g"] += g[m].sum()
            a["gm"] += (g[m] * r[m]).sum()
            a["m"] += r[m].sum()
            a["h"] += np.histogram(lg[m], LG_BINS)[0]
    out = {"V": sg / n, "thresholds": list(thresholds), "fV": [], "ion_V": [], "ion_M": [],
           "ion_med": []}
    for t in thresholds:
        a = acc[t]
        out["fV"].append(a["n"] / n)
        out["ion_V"].append(a["g"] / a["n"] if a["n"] else float("nan"))
        out["ion_M"].append(a["gm"] / a["m"] if a["m"] else float("nan"))
        out["ion_med"].append(_hist_median(a["h"]))
    return out


# ---------------------------------------------------------------------------
# 2. Mean free path at 912 A
# ---------------------------------------------------------------------------

def kappa_field(snap):
    """tau_912 per cell length, sigma_912 n_HI dl, as float32 [k, j, i]."""
    fHI, rho = snap.field("fHI"), snap.field("rhog")
    c = SIGMA_912 * snap.nH_avg * snap.dl
    kap = np.empty((snap.N,) * 3, np.float32)
    for k0, k1 in snap.slabs(12):
        kap[k0:k1] = (c * np.asarray(fHI[k0:k1], np.float64)
                      * np.asarray(rho[k0:k1], np.float64)).astype(np.float32)
    return kap


def source_cells(snap):
    """Linear indices ([k, j, i] order) and Luv of the cells with sources."""
    rt = snap.rt_grid()
    idx, lum = [], []
    for k0, k1 in snap.slabs(144):
        L = np.asarray(rt[k0:k1, :, :, 5]).ravel()
        w = np.nonzero(L > 0)[0]
        idx.append(w + k0 * snap.N * snap.N)
        lum.append(L[w])
    return np.concatenate(idx), np.concatenate(lum).astype(np.float64)


def random_dirs(rng, m):
    v = rng.normal(size=(m, 3))
    return v / np.linalg.norm(v, axis=1)[:, None]


def sightlines(kap, x0, dirs, rmax_cells, ds=0.5, skip=0.0, nblock=16, fmin=1e-7):
    """<exp(-tau(r))> along rays (x0, dirs in cell units, x = i + 0.5 at a
    cell centre, coordinates (x, y, z) = (i, j, k)), periodic, NGP sampling
    at the midpoints of steps of ds cells; tau starts at r = skip. Returns
    r (end of each step, cells) and the block sums of F (nblock, S), with
    the rays split into nblock interleaved groups for jackknife errors."""
    N = kap.shape[0]
    flat = kap.ravel()
    m = len(x0) - len(x0) % nblock
    x0, dirs = x0[:m], dirs[:m]
    S = int(math.ceil(rmax_cells / ds))
    r = (np.arange(S) + 1) * ds
    F = np.zeros((nblock, S))
    tau = np.zeros(m)
    for s in range(S):
        rm = (s + 0.5) * ds
        if rm >= skip:
            p = np.floor(x0 + dirs * rm).astype(np.int64) % N
            tau += flat[(p[:, 2] * N + p[:, 1]) * N + p[:, 0]] * ds
        f = np.exp(-tau)
        F[:, s] = f.reshape(-1, nblock).sum(axis=0)
        if f.max() < fmin:
            break
    return r, F, m


def _lambda_e(r, f):
    w = np.nonzero(f < math.exp(-1))[0]
    if len(w) == 0:
        return float("nan")
    i = w[0]
    if i == 0:
        return float(r[0])
    l1, l2 = math.log(max(f[i - 1], 1e-300)), math.log(max(f[i], 1e-300))
    return float(r[i - 1] + (r[i] - r[i - 1]) * (-1 - l1) / (l2 - l1))


def _lambda_fit(r, f, fmin=0.05):
    """Least-squares exp(-r/lambda) to <F(r)> where <F> > fmin."""
    from scipy.optimize import minimize_scalar
    w = f > fmin
    if w.sum() < 3:
        return float("nan")
    rr, ff = r[w], f[w]
    res = minimize_scalar(lambda q: np.sum((ff - np.exp(-rr / math.exp(q))) ** 2),
                          bounds=(math.log(rr[0] * 0.01), math.log(rr[-1] * 100)),
                          method="bounded")
    return float(math.exp(res.x))


def mfp_from_sums(r, F, m, dl_pmpc):
    """lambda_e and lambda_fit [pMpc] with jackknife errors over the blocks."""
    nb = F.shape[0]
    per = m / nb
    fall = F.sum(axis=0) / m
    le, lf = _lambda_e(r, fall), _lambda_fit(r, fall)
    jk_e, jk_f = [], []
    for b in range(nb):
        fb = (F.sum(axis=0) - F[b]) / (m - per)
        jk_e.append(_lambda_e(r, fb))
        jk_f.append(_lambda_fit(r, fb))

    def jkerr(v):
        v = np.array(v)
        if not np.all(np.isfinite(v)):
            return float("nan")
        return float(math.sqrt((nb - 1) / nb * np.sum((v - v.mean()) ** 2)))

    return {"lambda_e": le * dl_pmpc, "lambda_e_err": jkerr(jk_e) * dl_pmpc,
            "lambda_fit": lf * dl_pmpc, "lambda_fit_err": jkerr(jk_f) * dl_pmpc,
            "rmax": float(r[-1] * dl_pmpc)}


def mfp_stats(snap, nlos=4096, rmax_pmpc=40.0, ds=0.5, seed=1, kap=None, fHI_ion=0.01,
              src_skip=1.0, keep_curve=True):
    rng = np.random.default_rng(seed + 1000 * snap.step)
    N = snap.N
    if kap is None:
        kap = kappa_field(snap)
    rmax = rmax_pmpc / snap.dl_pmpc
    out = {"box_pmpc": N * snap.dl_pmpc, "nlos": nlos, "ds_cells": ds}

    def run(name, x0, skip=0.0):
        r, F, m = sightlines(kap, x0, random_dirs(rng, len(x0)), rmax, ds=ds, skip=skip)
        res = mfp_from_sums(r, F, m, snap.dl_pmpc)
        if keep_curve:  # <F(r)> on a coarse grid of r for plots
            fm = F.sum(axis=0) / m
            sel = np.unique(np.geomspace(1, len(r), 200).astype(int) - 1)
            res["r_pmpc"] = (r[sel] * snap.dl_pmpc).tolist()
            res["F"] = fm[sel].tolist()
        out[name] = res

    # random (volume-weighted) starts
    run("random", rng.uniform(0, N, size=(nlos, 3)))
    # starts in ionized cells (x_HI < fHI_ion)
    fl = snap.field("fHI").reshape(-1)
    nion = sum(int((np.asarray(fl[k0 * N * N:k1 * N * N]) < fHI_ion).sum())
               for k0, k1 in snap.slabs(5))
    out["f_ion"] = nion / fl.size
    if nion:
        if out["f_ion"] > 1e-2:  # rejection sampling: no N^3 index array
            c = np.empty(0, np.int64)
            while len(c) < nlos:
                t = np.sort(rng.integers(0, fl.size, 4 * nlos + int(4 * nlos / out["f_ion"])))
                c = np.concatenate([c, t[np.asarray(fl[t]) < fHI_ion]])
            c = rng.permutation(c)[:nlos]
        else:
            ion = np.nonzero(np.asarray(fl) < fHI_ion)[0]
            c = ion[rng.integers(0, len(ion), nlos)]
        x0 = np.stack([c % N, (c // N) % N, c // (N * N)], axis=1) + rng.uniform(0, 1, (nlos, 3))
        run("ionized", x0)
    # starts at sources (luminosity weighted), skipping their own cell
    if os.path.exists(os.path.join(snap.run.out, "RT", "rt_grid." + snap.tag)):
        idx, lum = source_cells(snap)
        out["nsource_cells"] = int(len(idx))
        if len(idx):
            c = idx[rng.choice(len(idx), nlos, p=lum / lum.sum())]
            x0 = np.stack([c % N, (c // N) % N, c // (N * N)], axis=1) + 0.5
            run("sources", x0.astype(np.float64), skip=src_skip)
    return out


# ---------------------------------------------------------------------------
# 3. Lya forest
# ---------------------------------------------------------------------------

def choose_pencils(N, nskew, rng):
    """nskew distinct (a, b) per axis: axis 0 runs along i at (k, j), axis 1
    along j at (k, i), axis 2 along k at (j, i)."""
    sel = {}
    for ax in range(3):
        p = rng.choice(N * N, size=min(nskew, N * N), replace=False)
        sel[ax] = (p // N, p % N)
    return sel


def read_pencils(snap, sel):
    """fHI, rhog, temp, Gamma_HI and v_los [cm/s] on the pencils (one pass
    over each file, in k slabs)."""
    N = snap.N
    names = ("fHI", "rhog", "temp", "gamma", "v")
    out = {ax: {n: np.empty((len(sel[ax][0]), N), np.float32) for n in names} for ax in sel}
    srcs = {"fHI": snap.field("fHI"), "rhog": snap.field("rhog"), "temp": snap.field("temp"),
            "gamma": snap.field("gamma", 3), "v": snap.gas()}
    for name, mm in srcs.items():
        bpc = mm.itemsize * (mm.shape[3] if mm.ndim == 4 else 1)
        for k0, k1 in snap.slabs(bpc):
            blk = np.asarray(mm[k0:k1])
            for ax in sel:
                a, b = sel[ax]
                if name == "v":  # momentum along the axis / density
                    arr = blk[..., 1 + ax] / blk[..., 0] * np.float32(snap.vel_unit)
                elif name == "gamma":
                    arr = blk[..., 0]
                else:
                    arr = blk
                if ax == 2:
                    out[ax][name][:, k0:k1] = arr[:, a, b].T
                else:
                    m = (a >= k0) & (a < k1)
                    if m.any():
                        out[ax][name][m] = arr[a[m] - k0, b[m], :] if ax == 0 else \
                            arr[a[m] - k0, :, b[m]]
    return out


def lya_tau(nHI, T, v, dl, H, cut=1e-6, profile="tophat"):
    """Lya optical depth on pixels = cells along periodic skewers (M, N).
    nHI [cm^-3], T [K], v [cm/s] (line-of-sight peculiar), dl proper cell
    [cm], H [s^-1]; pixel units: dv = H dl.

    profile = "tophat" (default): the gas of cell i is spread uniformly in
    redshift space between its faces, y_lo = i + v_{i-1/2}/dv and y_hi =
    i + 1 + v_{i+1/2}/dv (face velocities: the mean of the two cells), and
    every element has a Doppler profile, b = sqrt(2kT/m_H); the result is
    integrated over each pixel analytically. "point": all the gas of the
    cell at its centre, x_i = i + 1/2 + v_i/dv. The point form leaves
    empty pixels between cold (T ~ 1 K, b ~ 0.1 km/s) cells whose redshift-
    space spacing the velocity field stretches beyond a pixel, i.e. spurious
    transmission through neutral voids; the top-hat closes these gaps.
    Returns tau (M, N) and the pixel of each cell's centre."""
    from scipy.special import erf
    M, N = nHI.shape
    dv = H * dl
    A = K_LYA * nHI.astype(np.float64) * dl / dv
    beta = np.sqrt(2 * K_CGS * np.maximum(T.astype(np.float64), 1.0) / MH_CGS) / dv
    w = v.astype(np.float64) / dv
    i = np.arange(N)[None, :]
    xc = i + 0.5 + w
    if profile == "point":
        ylo = yhi = xc
    else:
        wf = 0.5 * (w + np.roll(w, -1, axis=1))          # face i+1/2
        ylo, yhi = i + np.roll(wf, 1, axis=1), i + 1 + wf
        ylo, yhi = np.minimum(ylo, yhi), np.maximum(ylo, yhi)
    delta = yhi - ylo
    point = delta < 1e-4 * beta
    base = np.floor(ylo)
    sig = A > cut
    if sig.any():
        w5 = math.ceil(5 * beta[sig].max()) + 1
        wlo, whi = min(N, w5), min(N, w5 + int(math.ceil(delta[sig].max())))
    else:
        wlo = whi = 1
    sq = 1 / math.sqrt(math.pi)

    def G(u):
        return u * erf(u) + np.exp(-u * u) * sq

    def cdf(q):  # fraction of the cell's gas at y < q
        with np.errstate(divide="ignore", invalid="ignore"):
            c = 0.5 * beta * (G((q - ylo) / beta) - G((q - yhi) / beta)) / delta
        if point.any():
            c = np.where(point, 0.5 * erf((q - 0.5 * (ylo + yhi)) / beta), c)
        return c

    rows = (np.arange(M) * N)[:, None]
    bi = base.astype(np.int64)
    tau = np.zeros(M * N)
    elo = cdf(base - wlo)
    for m in range(-wlo, whi + 1):
        ehi = cdf(base + m + 1)
        tau += np.bincount((rows + (bi + m) % N).ravel(), weights=(A * (ehi - elo)).ravel(),
                           minlength=M * N)
        elo = ehi
    return tau.reshape(M, N), np.floor(xc).astype(np.int64) % N


def rescale_A(tau, F_target):
    """A such that <exp(-A tau)> = F_target."""
    from scipy.optimize import brentq
    if not np.isfinite(F_target):
        return float("nan")
    t = tau.ravel()
    g = lambda q: np.mean(np.exp(-math.exp(q) * t)) - F_target
    lo, hi = math.log(1e-4), math.log(1e4)
    if g(lo) < 0 or g(hi) > 0:
        return float("nan")
    return float(math.exp(brentq(g, lo, hi, xtol=1e-6)))


def lya_stats(snap, nskew=4096, seed=1, obs=None, chunk=1024, profile="tophat"):
    rng = np.random.default_rng(seed + 7 + 1000 * snap.step)
    sel = choose_pencils(snap.N, nskew, rng)
    P = read_pencils(snap, sel)
    taus, gw, fw = [], 0.0, 0.0
    for ax in sel:
        d = P[ax]
        for c0 in range(0, d["fHI"].shape[0], chunk):
            s = slice(c0, c0 + chunk)
            nHI = d["fHI"][s].astype(np.float64) * d["rhog"][s] * snap.nH_avg
            tau, pix = lya_tau(nHI, d["temp"][s], d["v"][s], snap.dl, snap.H, profile=profile)
            taus.append(tau)
            # transmission-weighted Gamma: each cell weighted by the flux of
            # the pixel its gas falls in
            Fc = np.exp(-np.take_along_axis(tau, pix, axis=1))
            gw += float((Fc * d["gamma"][s]).sum())
            fw += float(Fc.sum())
    tau = np.concatenate(taus)
    F = np.exp(-tau)
    out = {"nskew": int(tau.shape[0]), "profile": profile, "F_mean": float(F.mean()),
           "F_mean_axis": [float(np.exp(-t).mean()) for t in np.split(tau, 3)],
           "tau_eff": float(-math.log(max(F.mean(), 1e-300))),
           "tau_median": float(np.median(tau)),
           "gamma_Fweighted": gw / fw if fw > 0 else float("nan"),
           "dv_kms": snap.H * snap.dl / 1e5}
    if obs is not None:
        Fo = obs_interp(obs, "mean_flux", snap.z)
        out["F_obs"] = Fo
        out["A"] = rescale_A(tau, Fo)
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def measure(args):
    run_path, step, opt = args
    run = Run(run_path)
    snap = Snapshot(run, step, run.stats(step))
    obs = read_obs(opt["data"])
    res = {"step": step, "z": snap.z, "a": snap.a, "xHII_V": snap.st.get("xHII_V"),
           "xHII_M": snap.st.get("xHII_M"), "T_V": snap.st.get("T_V"),
           "dl_pmpc": snap.dl_pmpc, "box_pmpc": snap.N * snap.dl_pmpc}
    res["gamma"] = gamma_stats(snap)
    res["mfp"] = mfp_stats(snap, nlos=opt["nlos"], rmax_pmpc=opt["rmax"], ds=opt["ds"],
                           seed=opt["seed"])
    res["lya"] = lya_stats(snap, nskew=opt["nskew"], seed=opt["seed"], obs=obs,
                           profile=opt.get("profile", "tophat"))
    g = res["gamma"]
    i2 = g["thresholds"].index(1e-2)
    res["gamma"]["req_from_F"] = g["ion_V"][i2] / res["lya"]["A"] if res["lya"].get("A") else \
        float("nan")
    print(f"# {run.name} step {step} z = {snap.z:.3f}: Gamma_V = {g['V']:.3g}, "
          f"Gamma(x<1e-2) = {g['ion_V'][i2]:.3g}, lambda_e(random) = "
          f"{res['mfp']['random']['lambda_e']:.3g} pMpc, <F> = {res['lya']['F_mean']:.4f}, "
          f"A = {res['lya'].get('A', float('nan')):.3g}", flush=True)
    return res


def selftest():
    """Analytic checks: a uniform field without velocities has the
    Gunn-Peterson tau = K_LYA n_HI / H in every pixel (any T); a single
    cell's line integrates to K_LYA n dl / dv over the skewer; sightlines
    through a uniform kappa give F = exp(-kappa r), lambda = 1/kappa."""
    N, M, dl, H = 64, 3, 3e22, 2e-17
    n = np.full((M, N), 1e-10)
    T = np.array([[1e2], [1e4], [1e5]]) * np.ones((1, N))
    tau, _ = lya_tau(n, T, np.zeros((M, N)), dl, H)
    gp = K_LYA * 1e-10 / H
    e1 = np.abs(tau / gp - 1).max()
    n1 = np.zeros((M, N))
    n1[:, 10] = 1e-10
    v = np.array([[0], [3.3 * H * dl], [-40.7 * H * dl]]) * np.ones((1, N))  # uniform shifts,
    e2 = 0.0                                                             # also across the edge
    for prof in ("tophat", "point"):
        tau1, pix = lya_tau(n1, T, v, dl, H, profile=prof)
        e2 = max(e2, np.abs(tau1.sum(axis=1) / (K_LYA * 1e-10 * dl / (H * dl)) - 1).max())
    # a stretched cold void: the top-hat fills every pixel, the point form leaves gaps
    vs = (np.arange(N) - N / 2)[None, :] * 0.8 * H * dl * np.ones((M, 1))
    vs[:, :8] = vs[:, -8:] = 0
    tt, _ = lya_tau(n, np.ones((M, N)), vs, dl, H)
    tp, _ = lya_tau(n, np.ones((M, N)), vs, dl, H, profile="point")
    e5 = (tt[:, 12:-12] / gp).min()
    gaps = int((tp[:, 12:-12] < 1e-3 * gp).sum())
    e1 = max(e1, np.abs(lya_tau(n, T, np.zeros((M, N)), dl, H, profile="point")[0] / gp - 1).max())
    kap = np.full((16, 16, 16), np.float32(0.01))
    rng = np.random.default_rng(0)
    r, F, m = sightlines(kap, rng.uniform(0, 16, (64, 3)), random_dirs(rng, 64), 400)
    res = mfp_from_sums(r, F, m, 1.0)
    e3 = abs(res["lambda_e"] / 100 - 1)
    e4 = abs(res["lambda_fit"] / 100 - 1)
    print(f"GP tau: max rel. error {e1:.2e}; single-line integral {e2:.2e}; "
          f"lambda_e {e3:.2e}; lambda_fit {e4:.2e}; peak pixel of the shifted lines "
          f"{tau1.argmax(axis=1).tolist()} (expect 10, 13, 33); stretched cold void: "
          f"min tau/tau_GP(top-hat) {e5:.3f} (expect 1/1.8 = 0.556), empty pixels (point) {gaps}")
    ok = e1 < 1e-6 and e2 < 1e-6 and e3 < 1e-3 and e4 < 1e-3 and \
        tau1.argmax(axis=1).tolist() == [10, 13, 33] and abs(e5 - 1 / 1.8) < 1e-3 and gaps > 0
    print("selftest", "passed" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run", nargs="?")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--zmin", type=float, default=4.5)
    ap.add_argument("--zmax", type=float, default=7.1)
    ap.add_argument("--steps", default="", help="comma-separated steps (default: all with fields)")
    ap.add_argument("--nlos", type=int, default=4096, help="sightlines per start set (912 A)")
    ap.add_argument("--rmax", type=float, default=40.0, help="sightline length [pMpc]")
    ap.add_argument("--ds", type=float, default=0.5, help="sightline step [cells]")
    ap.add_argument("--nskew", type=int, default=4096, help="Lya skewers per axis")
    ap.add_argument("--lya-profile", default="tophat", choices=("tophat", "point"),
                    help="gas of a cell spread between its faces in redshift space, or at its centre")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--jobs", type=int, default=1, help="snapshots measured in parallel")
    ap.add_argument("--data", default=DATA_DEFAULT)
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if not a.run:
        ap.error("RUN is required")
    run = Run(a.run)
    steps = [int(s) for s in a.steps.split(",") if s]
    snaps = run.snapshots(a.zmin, a.zmax, steps or None)
    opt = {"nlos": a.nlos, "rmax": a.rmax, "ds": a.ds, "nskew": a.nskew, "seed": a.seed,
           "data": a.data, "profile": a.lya_profile}
    jobs = [(run.top, s.step, opt) for s in snaps]
    if a.jobs > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=a.jobs) as ex:
            res = list(ex.map(measure, jobs))
    else:
        res = [measure(j) for j in jobs]
    out = {"run": run.top, "name": run.name, "N": run.N, "box": run.box,
           "c_reduced": run.c_reduced, "f_esc": run.f_esc, "options": opt,
           "snapshots": sorted(res, key=lambda d: -d["z"]), "ndot": run.ndot()}
    path = a.out or f"igm_{run.name}.json"
    with open(path, "w") as fp:
        json.dump(out, fp, indent=1)
    print("#", path)


if __name__ == "__main__":
    main()
