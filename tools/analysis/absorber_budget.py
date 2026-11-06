#!/usr/bin/env python3
"""Where the ionizing photons go, by gas density, in apps/radhydro runs that
differ in the RT method (docs/absorbers_closures.md): the photon sink of the
dense, self-shielding gas against the diffuse IGM and the redshift losses,
and how resolved the absorbers are.

usage: absorber_budget.py RUN:label [RUN:label ...] --outdir DIR
                          [--zlist 5.98,5.49,5.25,4.99] [--tol 0.03]
                          [--field gammaHI_avg] [--jobs 4]
                          [--absorbers 5.49,4.99] [--no-absorbers]
       absorber_budget.py --selftest

RUN is a run directory (params.txt and run/). The first run is the
reference (the rays): for each of its snapshots nearest the redshifts of
--zlist, every other run's snapshot nearest in redshift (within --tol) is
analysed. A snapshot is run/step_NNNN with fHI, fHeI, fHeII, rhog, temp and
the Gamma_HI field (--field, default gammaHI_avg: the substep-averaged
rate; gammaHI is the last substep's).

Per cell, with the code's units (cosmology.cpp calc_sim_var at the step's
scale factor; constants.hpp) and rate functions (rt_rates.cpp,
kernels/rt_cell.hpp):

  n_H = rhog <n_H>, n_HI = fHI n_H, n_HII = n_H - n_HI,
  n_He = rhog <n_He>, n_e = n_HII + n_He (fHeII + 2 (1 - fHeI - fHeII)),
  photoionisations  P = Gamma_HI n_HI V          [1/s]
  removed photons   X = x(delta) P,  x = min(exp(-(delta - 50)/50), 1)
                    (photo_cell's rho_source rule, the original and the
                    suite's; P11 `photon_consume_dense = 1` makes x = 1)
  recombinations    R = alpha_B(T) n_e n_HII V   [1/s]
                    (the HII case-B fit of rt_rates.cpp, tabulated at the
                    code's 0.01 dex nodes and interpolated log-linearly as
                    with legacy_physics = 0, P9)

binned in log10 delta (delta = rhog, the density in units of the mean
baryon density; 0.1 dex) and log10 x_HI (0.25 dex), with 2D histograms of
log Gamma and log x_HI against delta for the medians. The box totals are
compared with the sources' emission rate ndot = 0.2 10^25.2 Luv V_box
(igm_stats.py) and, for the moment methods, with the cumulative photon
budget of Log/ray.log (" Moments : in flight I injected J absorbed A
lost L", photons per hydrogen atom of a mean-density RT cell: inject_cell
divides by vol_rt <n_H>); for the rays, with ngam_rays_total of stats.txt
(the photons in flight, same units).

The absorbers (--absorbers, on the reference run only): connected regions
(6-connectivity, periodic) of x_HI > 0.5, their sizes in cells and the
largest HI column through each along the three grid axes
(N_HI = sum n_HI dl over the region's cells on one grid line), the number
above the Lyman limit N_HI > 1.6e17 cm^-2, and the share of the
volume-averaged 912 A opacity <sigma n_HI> in cells with tau_912 =
sigma_912 n_HI dl > 1. Needs scipy (scipy.ndimage.label).

Output: DIR/budget.json (every accumulator), DIR/budget.txt (the tables),
DIR/absorbers.txt, DIR/absorber_budget.png.
"""
import argparse
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igm_stats import (HUBBLE0, MPC_CM, OMEGAB, OMEGAM, OMEGAR, PHOT_PER_LUV,  # noqa: E402
                       RHOC_CGS, SIGMA_912, XHY, Run, Snapshot, f32)

# constants.hpp
MH_CGS = 1.67223e-24
YHE = f32(0.24)
EV_CGS, K_CGS = 1.60218e-12, 1.38066e-16
HI_K = f32(13.60) / (K_CGS / EV_CGS)   # 13.6 eV in K
RHO_SOURCE = 50.0
YR2SEC = 365 * 86400
N_LL = 1.6e17                            # Lyman-limit column [cm^-2]

# bins
LGD_EDGES = np.arange(-2.0, 4.0001, 0.1)          # log10 delta, 0.1 dex
LGX_EDGES = np.arange(-8.0, 0.0001, 0.25)         # log10 x_HI, 0.25 dex
LGG_EDGES = np.arange(-16.0, -9.0 + 1e-9, 0.02)   # log10 Gamma for medians
LGXF_EDGES = np.arange(-8.0, 0.0001, 0.05)        # log10 x_HI for medians
DELTA_CUTS = (10.0, 50.0, 100.0)
XHI_CUTS = (0.1, 0.5)


def alpha_b_formula(T):
    """rt_rates.cpp: HII recombination case B at temperature T [K],
    x = 2 HI_K/T (real(4) literals widened)."""
    x = 2 * HI_K / np.asarray(T, dtype=np.float64)
    return 2.753e-14 * x ** f32(1.5) / (1.0 + (x / f32(2.740)) ** f32(0.407)) ** f32(2.242)


_T_NODES = 10 ** (np.arange(1000) * f32(0.01))     # lgT_min = 0, lgT_del = 0.01, Ntemp = 1000
_LN_ALPHA_NODES = np.log(alpha_b_formula(_T_NODES))


def alpha_b(T):
    """The code's lookup with legacy_physics = 0 (P9): the fit at the table's
    0.01 dex nodes, log R linear in log T between them, the edge node
    outside."""
    lgT = np.log10(np.maximum(np.asarray(T, dtype=np.float64), 1e-30))
    u = np.clip(lgT / f32(0.01), 0.0, 999.0)
    i = np.minimum(u.astype(np.int64), 998)
    w = u - i
    return np.exp(_LN_ALPHA_NODES[i] * (1 - w) + _LN_ALPHA_NODES[i + 1] * w)


def consume_factor(delta):
    """photo_cell: the fraction of a cell's photoionisations taken out of the
    radiation, min(exp(-(delta - rho_source)/rho_source), 1)."""
    return np.minimum(np.exp(-(delta - RHO_SOURCE) / RHO_SOURCE), 1.0)


# ---------------------------------------------------------------------------
# the accumulators of one snapshot
# ---------------------------------------------------------------------------

class Acc:
    """Sums per (log delta, log x_HI) bin and the 2D histograms for the
    medians; the first index is the delta bin (with an under- and an
    overflow bin at each end), as np.digitize gives it."""

    def __init__(self):
        nd, nx = len(LGD_EDGES) + 1, len(LGX_EDGES) + 1
        self.n = np.zeros((nd, nx))
        self.P = np.zeros((nd, nx))     # photoionisations [1/s]
        self.X = np.zeros((nd, nx))     # removed photons, x(delta) P
        self.R = np.zeros((nd, nx))     # case-B recombinations [1/s]
        self.hg = np.zeros((nd, len(LGG_EDGES) + 1))    # log Gamma
        self.hx = np.zeros((nd, len(LGXF_EDGES) + 1))   # log x_HI (fine)
        self.ncell = 0
        self.tau_gt1_nHI = 0.0     # sum n_HI over cells with tau_912 > 1
        self.sum_nHI = 0.0
        self.n_tau_gt1 = 0
        self.gamma_sum = 0.0

    def add(self, delta, xhi, gam, P, X, R, tau):
        ib = np.digitize(np.log10(np.maximum(delta, 1e-30)), LGD_EDGES)
        jb = np.digitize(np.log10(np.maximum(xhi, 1e-30)), LGX_EDGES)
        nd, nx = self.n.shape
        key = ib * nx + jb
        self.n += np.bincount(key, minlength=nd * nx).reshape(nd, nx)
        self.P += np.bincount(key, weights=P, minlength=nd * nx).reshape(nd, nx)
        self.X += np.bincount(key, weights=X, minlength=nd * nx).reshape(nd, nx)
        self.R += np.bincount(key, weights=R, minlength=nd * nx).reshape(nd, nx)
        ng = self.hg.shape[1]
        kg = np.digitize(np.log10(np.maximum(gam, 1e-30)), LGG_EDGES)
        self.hg += np.bincount(ib * ng + kg, minlength=nd * ng).reshape(nd, ng)
        nxf = self.hx.shape[1]
        kx = np.digitize(np.log10(np.maximum(xhi, 1e-30)), LGXF_EDGES)
        self.hx += np.bincount(ib * nxf + kx, minlength=nd * nxf).reshape(nd, nxf)
        self.ncell += delta.size
        nHI_rel = delta * xhi
        self.sum_nHI += nHI_rel.sum()
        m = tau > 1
        self.tau_gt1_nHI += nHI_rel[m].sum()
        self.n_tau_gt1 += int(m.sum())
        self.gamma_sum += gam.sum()


def _median_from_hist(h, edges):
    """Median of the binned values (bin centres; the under/overflow bins at
    the edges)."""
    c = np.cumsum(h)
    if c[-1] == 0:
        return float("nan")
    i = int(np.searchsorted(c, 0.5 * c[-1]))
    i = min(max(i, 1), len(edges) - 1)
    return float(10 ** (0.5 * (edges[i - 1] + edges[i])))


def snapshot_budget(snap, field):
    """The accumulators of one snapshot (slab by slab)."""
    fHI, fHeI, fHeII = snap.field("fHI"), snap.field("fHeI"), snap.field("fHeII")
    rho, temp = snap.field("rhog"), snap.field("temp")
    gam = snap.field(field)
    V = snap.dl ** 3
    nH_avg = snap.nH_avg
    nHe_avg = nH_avg * (YHE / XHY) / 4.0          # <n_He> = rhob YHe / (4 mH)
    acc = Acc()
    for k0, k1 in snap.slabs(40):
        d = np.asarray(rho[k0:k1], np.float64).ravel()
        x = np.asarray(fHI[k0:k1], np.float64).ravel()
        he1 = np.asarray(fHeI[k0:k1], np.float64).ravel()
        he2 = np.asarray(fHeII[k0:k1], np.float64).ravel()
        T = np.asarray(temp[k0:k1], np.float64).ravel()
        g = np.asarray(gam[k0:k1], np.float64).ravel()
        nH = d * nH_avg
        nHI = x * nH
        nHII = nH - nHI
        nHe = d * nHe_avg
        ne = nHII + nHe * (he2 + 2 * (1 - he1 - he2))
        P = g * nHI * V
        X = consume_factor(d) * P
        R = alpha_b(T) * ne * nHII * V
        tau = SIGMA_912 * nHI * snap.dl
        acc.add(d, x, g, P, X, R, tau)
    return acc


def budget_summary(acc, snap):
    """Fractions and totals of one snapshot's accumulators."""
    nd, nx = acc.n.shape
    lgd_lo = np.concatenate([[-np.inf], LGD_EDGES])          # lower edge of delta bin ib
    lgx_lo = np.concatenate([[-np.inf], LGX_EDGES])
    P_tot, X_tot, R_tot = acc.P.sum(), acc.X.sum(), acc.R.sum()
    out = {"ncell": int(acc.ncell), "P_tot": P_tot, "X_tot": X_tot, "R_tot": R_tot,
           "gamma_V": acc.gamma_sum / acc.ncell,
           "frac_opacity_tau_gt1": acc.tau_gt1_nHI / acc.sum_nHI if acc.sum_nHI > 0 else 0.0,
           "frac_cells_tau_gt1": acc.n_tau_gt1 / acc.ncell,
           "mean_nHI_rel": acc.sum_nHI / acc.ncell}
    # by delta cut
    for dc in DELTA_CUTS:
        m = lgd_lo >= np.log10(dc) - 1e-9
        out[f"delta>{dc:g}"] = {"fV": acc.n[m].sum() / acc.ncell,
                                "P": acc.P[m].sum() / P_tot if P_tot else 0.0,
                                "X": acc.X[m].sum() / X_tot if X_tot else 0.0,
                                "R": acc.R[m].sum() / R_tot if R_tot else 0.0,
                                "X_abs": acc.X[m].sum(), "P_abs": acc.P[m].sum(),
                                "R_abs": acc.R[m].sum()}
    for xc in XHI_CUTS:
        m = lgx_lo >= np.log10(xc) - 1e-9
        out[f"xHI>{xc:g}"] = {"fV": acc.n[:, m].sum() / acc.ncell,
                              "P": acc.P[:, m].sum() / P_tot if P_tot else 0.0,
                              "X": acc.X[:, m].sum() / X_tot if X_tot else 0.0,
                              "R": acc.R[:, m].sum() / R_tot if R_tot else 0.0}
    # the emission rate of the sources at this step
    V_box = (snap.run.box / HUBBLE0) ** 3           # cMpc^3
    ndot = PHOT_PER_LUV * float(snap.st.get("Luv", 0.0)) * V_box
    out["ndot_box"] = ndot
    out["dt_s"] = float(snap.st.get("dt_yr", 0.0)) * YR2SEC
    out["P_over_ndot"] = P_tot / ndot if ndot else float("nan")
    out["X_over_ndot"] = X_tot / ndot if ndot else float("nan")
    out["R_over_ndot"] = R_tot / ndot if ndot else float("nan")
    # the curves against delta: medians of Gamma and x_HI, the sums per bin
    lgd_c = np.concatenate([[LGD_EDGES[0] - 0.05], 0.5 * (LGD_EDGES[1:] + LGD_EDGES[:-1]),
                            [LGD_EDGES[-1] + 0.05]])
    out["lgd_centres"] = lgd_c.tolist()
    out["n_delta"] = acc.n.sum(axis=1).tolist()
    out["P_delta"] = acc.P.sum(axis=1).tolist()
    out["X_delta"] = acc.X.sum(axis=1).tolist()
    out["R_delta"] = acc.R.sum(axis=1).tolist()
    out["gamma_med_delta"] = [_median_from_hist(acc.hg[i], LGG_EDGES) for i in range(nd)]
    out["xhi_med_delta"] = [_median_from_hist(acc.hx[i], LGXF_EDGES) for i in range(nd)]
    # the self-shielding density: the first delta bin (rising) whose median
    # x_HI reaches 0.5, interpolated in log delta
    xm = np.array(out["xhi_med_delta"])
    ok = np.isfinite(xm) & (np.array(out["n_delta"]) > 100)
    out["delta_ss"] = float("nan")
    for i in range(1, nd):
        if ok[i] and ok[i - 1] and xm[i] >= 0.5 > xm[i - 1]:
            f = (np.log10(0.5) - np.log10(xm[i - 1])) / (np.log10(xm[i]) - np.log10(xm[i - 1]))
            out["delta_ss"] = float(10 ** (lgd_c[i - 1] + f * (lgd_c[i] - lgd_c[i - 1])))
            break
    return out


# ---------------------------------------------------------------------------
# the photon budget of Log/ray.log (moment methods) and stats.txt (rays)
# ---------------------------------------------------------------------------

def read_moment_budget(run):
    """{step: (in_flight, injected, absorbed, lost)} from the ' Moments :'
    lines of run/Log/ray.log, in photons per hydrogen atom of a mean-density
    RT cell (cumulative over the run; the log is rewritten at a restart, so
    only the steps it holds appear)."""
    path = os.path.join(run.out, "Log", "ray.log")
    out = {}
    if not os.path.exists(path):
        return out
    step = None
    with open(path) as f:
        for line in f:
            if line.startswith(" Hydro step :"):
                step = int(line.split()[3])
            elif line.startswith(" Moments") and step is not None:
                t = line.split()
                # ' Moments : in flight I injected J absorbed A lost L'
                out[step] = (float(t[4]), float(t[6]), float(t[8]), float(t[10]))
    return out


def step_budget(run, snap):
    """The photon budget of the snapshot's hydro step: emitted (ndot dt),
    the change of the photons in flight, absorbed and lost (moments), all
    in absolute photons; the rays give emitted and the in-flight change
    only (ngam_rays_total of stats.txt), and their absorbed + lost is the
    remainder."""
    st = snap.step
    V_box = (run.box / HUBBLE0) ** 3
    dt = float(snap.st.get("dt_yr", 0.0)) * YR2SEC
    emitted = PHOT_PER_LUV * float(snap.st.get("Luv", 0.0)) * V_box * dt
    per_atom = snap.nH_avg * snap.dl ** 3 * snap.N ** 3   # photons per (atom per mean cell) x N^3
    out = {"emitted": emitted, "dt_s": dt}
    mb = read_moment_budget(run)
    if st in mb and (st - 1) in mb:
        a, b = mb[st - 1], mb[st]
        out["method"] = "moments"
        out["in_flight"] = b[0] * per_atom
        out["d_in_flight"] = (b[0] - a[0]) * per_atom
        out["injected"] = (b[1] - a[1]) * per_atom
        out["absorbed"] = (b[2] - a[2]) * per_atom
        out["lost"] = (b[3] - a[3]) * per_atom
        return out
    # rays: ngam_rays_total of this and the previous step's stats.txt
    try:
        prev = run.stats(st - 1)
    except (OSError, KeyError):
        prev = None
    g1 = float(snap.st.get("ngam_rays_total", float("nan")))
    g0 = float(prev.get("ngam_rays_total", float("nan"))) if prev else float("nan")
    out["method"] = "rays"
    out["in_flight"] = g1 * snap.nH_avg * snap.dl ** 3
    out["d_in_flight"] = (g1 - g0) * snap.nH_avg * snap.dl ** 3
    out["absorbed_plus_lost"] = emitted - out["d_in_flight"]
    return out


# ---------------------------------------------------------------------------
# the absorbers of the reference run
# ---------------------------------------------------------------------------

def _union_find_merge(labels, pairs):
    """Relabel so that every pair (a, b) shares a label; returns the new
    label array (0 kept) and the number of components."""
    nl = int(labels.max()) + 1
    parent = np.arange(nl)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in pairs:
        ra, rb = find(int(a)), find(int(b))
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    root = np.array([find(i) for i in range(nl)])
    root[0] = 0
    uniq, new = np.unique(root, return_inverse=True)
    new = new.astype(np.int64)
    if uniq[0] != 0:
        new += 1
    return new[labels], int(len(uniq) - (1 if uniq[0] == 0 else 0))


def periodic_components(mask):
    """6-connected components of a periodic boolean cube: labels (0 =
    outside) and their number. scipy.ndimage.label, then the components
    touching across each pair of opposite faces merged."""
    from scipy import ndimage
    labels, n = ndimage.label(mask)
    labels = labels.astype(np.int64)
    pairs = set()
    for ax in range(3):
        a = np.take(labels, 0, axis=ax).ravel()
        b = np.take(labels, -1, axis=ax).ravel()
        m = (a > 0) & (b > 0) & (a != b)
        for pa, pb in zip(a[m], b[m]):
            pairs.add((int(pa), int(pb)))
    if not pairs:
        return labels, int(n)
    return _union_find_merge(labels, sorted(pairs))


def absorber_stats(nHI_dl, mask):
    """Sizes and the largest grid-line HI column per axis of the connected
    regions of mask; nHI_dl = n_HI dl per cell [cm^-2]."""
    labels, ncomp = periodic_components(mask)
    idx = np.nonzero(labels)
    lab = labels[idx]
    col = nHI_dl[idx]
    sizes = np.bincount(lab, minlength=ncomp + 1)[1:]
    N = mask.shape[0]
    peak = np.zeros((3, ncomp + 1))
    for ax in range(3):
        others = [d for d in range(3) if d != ax]
        line = idx[others[0]].astype(np.int64) * N + idx[others[1]]
        key = lab * (N * N) + line
        uk, inv = np.unique(key, return_inverse=True)
        s = np.bincount(inv, weights=col)
        klab = uk // (N * N)
        np.maximum.at(peak[ax], klab, s)
    peak = peak[:, 1:]
    return {"ncomp": int(ncomp), "sizes": sizes, "peak_col": peak, "labels": labels}


def snapshot_absorbers(snap):
    fHI = np.asarray(snap.field("fHI"), np.float64)
    rho = np.asarray(snap.field("rhog"), np.float64)
    nHI_dl = fHI * rho * snap.nH_avg * snap.dl
    mask = fHI > 0.5
    r = absorber_stats(nHI_dl, mask)
    sizes, peak = r["sizes"], r["peak_col"]
    tau = SIGMA_912 * nHI_dl
    out = {"z": snap.z, "step": snap.step, "n_neutral_cells": int(mask.sum()),
           "f_neutral_V": float(mask.mean()), "ncomp": r["ncomp"],
           "frac_opacity_tau_gt1": float(nHI_dl[tau > 1].sum() / nHI_dl.sum()),
           "frac_opacity_neutral": float(nHI_dl[mask].sum() / nHI_dl.sum()),
           "n_cells_tau_gt1": int((tau > 1).sum()),
           "mean_nHI": float(nHI_dl.mean() / snap.dl),
           "lambda_from_mean_nHI_pMpc": float(1.0 / (SIGMA_912 * nHI_dl.mean() / snap.dl) / MPC_CM),
           "dl_pkpc": snap.dl / MPC_CM * 1e3}
    if r["ncomp"] > 0:
        pmax = peak.max(axis=0)
        pmin = peak.min(axis=0)
        out["size_quantiles"] = {q: float(np.percentile(sizes, q)) for q in (50, 90, 99)}
        out["size_max"] = int(sizes.max())
        out["n_size_1"] = int((sizes == 1).sum())
        out["n_size_ge8"] = int((sizes >= 8).sum())
        out["n_size_ge27"] = int((sizes >= 27).sum())
        out["frac_cells_in_ge8"] = float(sizes[sizes >= 8].sum() / sizes.sum())
        out["n_LL_any_axis"] = int((pmax > N_LL).sum())
        out["n_LL_all_axes"] = int((pmin > N_LL).sum())
        out["frac_cells_in_LL_any"] = float(sizes[pmax > N_LL].sum() / sizes.sum())
        out["peak_col_quantiles"] = {q: float(np.percentile(pmax, q)) for q in (50, 90, 99)}
        out["peak_col_max"] = float(pmax.max())
        # the size-column relation: median peak column per size class
        cls = {}
        for lo, hi in ((1, 1), (2, 7), (8, 26), (27, 124), (125, 10 ** 9)):
            m = (sizes >= lo) & (sizes <= hi)
            if m.any():
                cls[f"{lo}-{hi}"] = {"n": int(m.sum()), "col_med": float(np.median(pmax[m])),
                                     "n_LL": int((pmax[m] > N_LL).sum())}
        out["by_size"] = cls
    return out


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def nearest_snapshot(run, z, tol):
    best = None
    for s in run.snapshots():
        if best is None or abs(s.z - z) < abs(best.z - z):
            best = s
    if best is None or abs(best.z - z) > tol:
        return None
    return best


def _work(args):
    runpath, step, field, do_abs = args
    run = Run(runpath)
    snap = Snapshot(run, step, run.stats(step))
    acc = snapshot_budget(snap, field)
    res = {"z": snap.z, "step": step, "budget": budget_summary(acc, snap),
           "step_budget": step_budget(run, snap)}
    if do_abs:
        res["absorbers"] = snapshot_absorbers(snap)
    return runpath, step, res


def fmt_frac(x):
    return "nan" if not np.isfinite(x) else f"{100 * x:.1f} %"


def write_tables(results, labels, zlist, outdir):
    L = []
    L.append("# absorber_budget.py: where the photons go, by density\n")
    L.append("P = photoionisations Gamma n_HI V, X = photons removed x(delta) P (the suite's rule),")
    L.append("R = case-B recombinations alpha_B(T) n_e n_HII V; ndot = the sources' emission rate.")
    L.append("Fractions are of the box total of each quantity; fV the volume fraction.\n")
    for z in zlist:
        L.append(f"## z ~ {z}\n")
        L.append("| run | z | X/ndot | P/ndot | R/ndot | X(d>10) | X(d>50) | X(d>100) | "
                 "P(d>50) | P(d>100) | R(d>10) | R(d>50) | R(d>100) | fV(d>50) |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for lab in labels:
            r = results[lab].get(z)
            if r is None:
                L.append(f"| {lab} | - |")
                continue
            b = r["budget"]
            L.append(f"| {lab} | {r['z']:.3f} | {b['X_over_ndot']:.3f} | {b['P_over_ndot']:.3f} | "
                     f"{b['R_over_ndot']:.3f} | {fmt_frac(b['delta>10']['X'])} | "
                     f"{fmt_frac(b['delta>50']['X'])} | {fmt_frac(b['delta>100']['X'])} | "
                     f"{fmt_frac(b['delta>50']['P'])} | {fmt_frac(b['delta>100']['P'])} | "
                     f"{fmt_frac(b['delta>10']['R'])} | {fmt_frac(b['delta>50']['R'])} | "
                     f"{fmt_frac(b['delta>100']['R'])} | {100 * b['delta>50']['fV']:.3f} % |")
        L.append("")
        L.append("| run | X(xHI>0.1) | X(xHI>0.5) | P(xHI>0.1) | P(xHI>0.5) | R(xHI>0.1) | "
                 "R(xHI>0.5) | fV(xHI>0.5) | delta_ss | <Gamma>_V | opacity in tau>1 cells |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for lab in labels:
            r = results[lab].get(z)
            if r is None:
                continue
            b = r["budget"]
            L.append(f"| {lab} | {fmt_frac(b['xHI>0.1']['X'])} | {fmt_frac(b['xHI>0.5']['X'])} | "
                     f"{fmt_frac(b['xHI>0.1']['P'])} | {fmt_frac(b['xHI>0.5']['P'])} | "
                     f"{fmt_frac(b['xHI>0.1']['R'])} | {fmt_frac(b['xHI>0.5']['R'])} | "
                     f"{100 * b['xHI>0.5']['fV']:.3f} % | {b['delta_ss']:.0f} | "
                     f"{b['gamma_V']:.3e} | {fmt_frac(b['frac_opacity_tau_gt1'])} |")
        L.append("")
        L.append("Photon budget of the hydro step (photons; the moment methods from Log/ray.log,"
                 " the rays from ngam_rays_total): emitted = ndot dt")
        L.append("| run | emitted | X dt | P dt | R dt | absorbed (log) | lost (log) | "
                 "d(in flight) | in flight | absorbed/emitted | lost/emitted | X dt/absorbed |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for lab in labels:
            r = results[lab].get(z)
            if r is None:
                continue
            b, s = r["budget"], r["step_budget"]
            dt = s["dt_s"]
            e = s["emitted"]
            ab = s.get("absorbed", float("nan"))
            lo = s.get("lost", float("nan"))
            L.append(f"| {lab} | {e:.3e} | {b['X_tot'] * dt:.3e} | {b['P_tot'] * dt:.3e} | "
                     f"{b['R_tot'] * dt:.3e} | {ab:.3e} | {lo:.3e} | {s['d_in_flight']:.3e} | "
                     f"{s['in_flight']:.3e} | {ab / e:.3f} | {lo / e:.4f} | "
                     f"{b['X_tot'] * dt / ab if np.isfinite(ab) and ab else float('nan'):.3f} |")
        L.append("")
        # Wu's mechanism: the recombination rate in delta bins relative to the reference
        ref = results[labels[0]].get(z)
        if ref is not None:
            L.append("Recombination rate in delta bins relative to the reference (Wu+21: 1.3-1.4x "
                     "in resolved r^-2 absorbers), and the reference's share of all recombinations "
                     "in the bin:")
            bins = [(10, 30), (30, 100), (100, 300), (300, 1000), (1000, 1e9)]
            head = "| run | " + " | ".join(f"{lo:g}-{hi:g}" if hi < 1e8 else f">{lo:g}"
                                            for lo, hi in bins) + " |"
            L.append(head)
            L.append("|---|" + "---|" * len(bins))
            lgd = np.array(ref["budget"]["lgd_centres"])
            Rref = np.array(ref["budget"]["R_delta"])
            row = "| share (ref) |"
            for lo, hi in bins:
                m = (lgd >= np.log10(lo) - 1e-9) & (lgd < np.log10(hi) - 1e-9)
                row += f" {100 * Rref[m].sum() / Rref.sum():.2f} % |"
            L.append(row)
            for lab in labels[1:]:
                r = results[lab].get(z)
                if r is None:
                    continue
                R = np.array(r["budget"]["R_delta"])
                row = f"| {lab} |"
                for lo, hi in bins:
                    m = (lgd >= np.log10(lo) - 1e-9) & (lgd < np.log10(hi) - 1e-9)
                    row += f" {R[m].sum() / Rref[m].sum():.3f} |" if Rref[m].sum() > 0 else " - |"
                L.append(row)
            # the same for the photoionisations
            Pref = np.array(ref["budget"]["P_delta"])
            L.append("")
            L.append("Photoionisation rate (P) in the same bins relative to the reference:")
            L.append(head)
            L.append("|---|" + "---|" * len(bins))
            for lab in labels[1:]:
                r = results[lab].get(z)
                if r is None:
                    continue
                P = np.array(r["budget"]["P_delta"])
                row = f"| {lab} |"
                for lo, hi in bins:
                    m = (lgd >= np.log10(lo) - 1e-9) & (lgd < np.log10(hi) - 1e-9)
                    row += f" {P[m].sum() / Pref[m].sum():.3f} |" if Pref[m].sum() > 0 else " - |"
                L.append(row)
        L.append("")
        # self-shielding curves
        L.append("Median Gamma_HI [1e-12 /s] and x_HI against delta (bins of 0.1 dex, selected):")
        sel = [1, 3, 10, 30, 100, 300, 1000, 3000]
        L.append("| run | " + " | ".join(f"d={d}" for d in sel) + " |")
        L.append("|---|" + "---|" * len(sel))
        for lab in labels:
            r = results[lab].get(z)
            if r is None:
                continue
            lgd = np.array(r["budget"]["lgd_centres"])
            gm = np.array(r["budget"]["gamma_med_delta"])
            xm = np.array(r["budget"]["xhi_med_delta"])
            row = f"| {lab} |"
            for d in sel:
                i = int(np.argmin(np.abs(lgd - np.log10(d))))
                row += f" {gm[i] * 1e12:.2f} / {xm[i]:.1e} |"
            L.append(row)
        L.append("")
    with open(os.path.join(outdir, "budget.txt"), "w") as f:
        f.write("\n".join(L) + "\n")


def write_absorbers(results, labels, outdir):
    L = ["# absorbers of the reference run: connected regions of x_HI > 0.5 (periodic, "
         "6-connected)\n"]
    for lab in labels:
        for z, r in sorted(results[lab].items(), reverse=True):
            a = r.get("absorbers")
            if a is None:
                continue
            L.append(f"## {lab}, z = {a['z']:.3f} (step {a['step']}), cell {a['dl_pkpc']:.1f} pkpc")
            L.append(f"neutral cells (x_HI > 0.5): {a['n_neutral_cells']} "
                     f"({100 * a['f_neutral_V']:.3f} % of the volume), components: {a['ncomp']}")
            L.append(f"share of the volume-averaged 912 A opacity: in neutral cells "
                     f"{100 * a['frac_opacity_neutral']:.1f} %, in cells with tau_912 > 1 "
                     f"{100 * a['frac_opacity_tau_gt1']:.1f} % ({a['n_cells_tau_gt1']} cells); "
                     f"1/(sigma <n_HI>) = {a['lambda_from_mean_nHI_pMpc']:.2f} pMpc")
            if a["ncomp"]:
                sq, cq = a["size_quantiles"], a["peak_col_quantiles"]
                L.append(f"sizes [cells]: median {sq[50]:.0f}, 90 % {sq[90]:.0f}, 99 % "
                         f"{sq[99]:.0f}, max {a['size_max']}; single cells {a['n_size_1']}, "
                         f">= 8 cells {a['n_size_ge8']} ({100 * a['frac_cells_in_ge8']:.1f} % of "
                         f"the neutral cells), >= 27 cells {a['n_size_ge27']}")
                L.append(f"largest grid-line column N_HI [cm^-2] per component (max over axes): "
                         f"median {cq[50]:.2e}, 90 % {cq[90]:.2e}, 99 % {cq[99]:.2e}, max "
                         f"{a['peak_col_max']:.2e}")
                L.append(f"Lyman-limit systems (N_HI > 1.6e17): {a['n_LL_any_axis']} on at least "
                         f"one axis, {a['n_LL_all_axes']} on all three "
                         f"({100 * a['frac_cells_in_LL_any']:.1f} % of the neutral cells)")
                L.append("| size [cells] | n | median peak N_HI | n with N_HI > 1.6e17 |")
                L.append("|---|---|---|---|")
                for k, v in a["by_size"].items():
                    L.append(f"| {k} | {v['n']} | {v['col_med']:.2e} | {v['n_LL']} |")
            L.append("")
    with open(os.path.join(outdir, "absorbers.txt"), "w") as f:
        f.write("\n".join(L) + "\n")


def make_figure(results, labels, zlist, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = ["black", "#2d7bd6", "#d9541e", "#2a9d6b", "#8e5bc2", "#c99a1b", "#d64278"]
    nz = len(zlist)
    fig, axes = plt.subplots(3, nz, figsize=(3.6 * nz, 9.5), squeeze=False)
    for j, z in enumerate(zlist):
        ax0, ax1, ax2 = axes[0, j], axes[1, j], axes[2, j]
        for i, lab in enumerate(labels):
            r = results[lab].get(z)
            if r is None:
                continue
            b = r["budget"]
            lgd = np.array(b["lgd_centres"])[1:-1]
            n = np.array(b["n_delta"])[1:-1]
            X = np.array(b["X_delta"])[1:-1]
            R = np.array(b["R_delta"])[1:-1]
            c = colors[i % len(colors)]
            # cumulative share of the removed photons and the recombinations above delta
            cX = np.cumsum(X[::-1])[::-1] / max(X.sum(), 1e-300)
            cR = np.cumsum(R[::-1])[::-1] / max(R.sum(), 1e-300)
            ax0.plot(10 ** lgd, cX, color=c, lw=1.6, label=lab if j == 0 else None)
            ax0.plot(10 ** lgd, cR, color=c, lw=1.0, ls="--")
            gm = np.array(b["gamma_med_delta"])[1:-1]
            xm = np.array(b["xhi_med_delta"])[1:-1]
            ok = n > 100
            ax1.plot(10 ** lgd[ok], gm[ok], color=c, lw=1.6)
            ax2.plot(10 ** lgd[ok], xm[ok], color=c, lw=1.6)
        ax0.set_xscale("log")
        ax0.set_ylim(0, 1)
        ax0.set_title(f"z = {z}")
        ax0.set_ylabel("share above $\\delta$: removed photons (solid), recombinations (dashed)"
                       if j == 0 else "")
        ax0.axvline(50, color="0.6", lw=0.8)
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_ylabel("median $\\Gamma_{\\rm HI}$ [s$^{-1}$]" if j == 0 else "")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_ylabel("median $x_{\\rm HI}$" if j == 0 else "")
        ax2.set_xlabel("$\\delta = \\rho/\\bar\\rho_b$")
        for ax in (ax0, ax1, ax2):
            ax.set_xlim(0.03, 1e4)
    axes[0, 0].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "absorber_budget.png"), dpi=130)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--outdir")
    ap.add_argument("--zlist", default="5.98,5.49,5.25,4.99")
    ap.add_argument("--tol", type=float, default=0.03)
    ap.add_argument("--field", default="gammaHI_avg")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--absorbers", default="5.49,4.99",
                    help="redshifts (of --zlist) at which the reference run's absorbers are "
                         "analysed")
    ap.add_argument("--no-absorbers", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if not a.runs or not a.outdir:
        ap.error("runs and --outdir are required")
    os.makedirs(a.outdir, exist_ok=True)
    zlist = [float(x) for x in a.zlist.split(",")]
    zabs = [] if a.no_absorbers else [float(x) for x in a.absorbers.split(",") if x]
    labels, paths = [], {}
    for spec in a.runs:
        p, _, lab = spec.rpartition(":")
        if not p:
            p, lab = lab, os.path.basename(lab.rstrip("/"))
        labels.append(lab)
        paths[lab] = p
    ref = Run(paths[labels[0]])
    tasks = []
    for z in zlist:
        s = nearest_snapshot(ref, z, a.tol)
        if s is None:
            print(f"reference: no snapshot within {a.tol} of z = {z}", file=sys.stderr)
            continue
        for lab in labels:
            run = Run(paths[lab])
            sn = nearest_snapshot(run, s.z, a.tol)
            if sn is None:
                print(f"{lab}: no snapshot within {a.tol} of z = {s.z:.3f}", file=sys.stderr)
                continue
            do_abs = lab == labels[0] and any(abs(z - za) < 0.05 for za in zabs)
            tasks.append((paths[lab], sn.step, a.field, do_abs, lab, z))
    results = {lab: {} for lab in labels}
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for (path, step, res), t in zip(ex.map(_work, [t[:4] for t in tasks]), tasks):
            results[t[4]][t[5]] = res
            print(f"{t[4]} z = {res['z']:.3f} step {step}: X/ndot {res['budget']['X_over_ndot']:.3f}"
                  f" R/ndot {res['budget']['R_over_ndot']:.3f}", flush=True)
    with open(os.path.join(a.outdir, "budget.json"), "w") as f:
        json.dump({"labels": labels, "zlist": zlist, "field": a.field,
                   "results": {lab: {str(z): r for z, r in d.items()}
                               for lab, d in results.items()}}, f, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    write_tables(results, labels, zlist, a.outdir)
    write_absorbers(results, labels, a.outdir)
    make_figure(results, labels, zlist, a.outdir)
    return 0


# ---------------------------------------------------------------------------
# self-test on a synthetic cube
# ---------------------------------------------------------------------------

class _FakeRun:
    box = 10.0


def selftest():
    """A 16^3 cube: ionised background at the mean density with a clump of
    2x2x2 cells at delta = 100 across the periodic corner, one cell at
    delta = 200 and one neutral single cell at delta = 1. The sums, the
    rho_source factor, the cuts and (with scipy) the periodic components and
    their columns are checked against direct arithmetic."""
    N = 16
    a = 1.0 / 6.5
    dl = a * _FakeRun.box / HUBBLE0 / N * MPC_CM
    nH_avg = RHOC_CGS * OMEGAM * HUBBLE0 ** 2 / a ** 3 * OMEGAB / OMEGAM * XHY / MH_CGS
    nHe_avg = nH_avg * (YHE / XHY) / 4
    V = dl ** 3
    delta = np.ones((N, N, N))
    xhi = np.full((N, N, N), 1e-4)
    T = np.full((N, N, N), 1e4)
    gam = np.full((N, N, N), 1e-12)
    fHeI = np.full((N, N, N), 0.0)
    fHeII = np.full((N, N, N), 1.0)
    # the clump across the corner: cells {0, N-1}^3
    ci = [0, N - 1]
    for k in ci:
        for j in ci:
            for i in ci:
                delta[k, j, i] = 100.0
                xhi[k, j, i] = 0.9
                gam[k, j, i] = 1e-14
                T[k, j, i] = 2e4
    delta[8, 8, 8] = 200.0
    xhi[8, 8, 8] = 0.3
    gam[8, 8, 8] = 1e-13
    delta[4, 4, 4] = 1.0
    xhi[4, 4, 4] = 0.99
    gam[4, 4, 4] = 0.0

    nH = delta * nH_avg
    nHI = xhi * nH
    nHII = nH - nHI
    ne = nHII + delta * nHe_avg * (fHeII + 2 * (1 - fHeI - fHeII))
    P = gam * nHI * V
    X = consume_factor(delta) * P
    R = alpha_b(T) * ne * nHII * V
    acc = Acc()
    acc.add(delta.ravel(), xhi.ravel(), gam.ravel(), P.ravel(), X.ravel(), R.ravel(),
            (SIGMA_912 * nHI * dl).ravel())
    errs = []

    def check(name, got, want, rtol=1e-9):
        if not (abs(got - want) <= rtol * max(abs(want), 1e-300)):
            errs.append(f"{name}: got {got!r}, want {want!r}")

    # the rate function against the formula at a node and between nodes
    check("alpha_B node", float(alpha_b(1e4)), float(alpha_b_formula(1e4)), 1e-6)  # node at 400 lgT_del (float)
    T_mid = 10 ** (4.005)
    want = math.exp(0.5 * (math.log(alpha_b_formula(1e4)) + math.log(alpha_b_formula(10 ** 4.01))))
    check("alpha_B interp", float(alpha_b(T_mid)), want, 1e-9)
    check("alpha_B(1e4) value", float(alpha_b_formula(1e4)), 2.59e-13, 0.03)  # HG97 case B
    # the rho_source factor
    check("x(100)", float(consume_factor(100.0)), math.exp(-1.0))
    check("x(1)", float(consume_factor(1.0)), 1.0)
    # totals
    check("P_tot", acc.P.sum(), P.sum())
    check("X_tot", acc.X.sum(), X.sum())
    check("R_tot", acc.R.sum(), R.sum())
    check("ncell", acc.ncell, N ** 3)
    # the cuts: delta > 50 holds the 8 clump cells and the delta = 200 cell
    m50 = delta > 50
    nd = len(LGD_EDGES) + 1
    lgd_lo = np.concatenate([[-np.inf], LGD_EDGES])
    sel = lgd_lo >= np.log10(50) - 1e-9
    check("n(d>50)", acc.n[sel].sum(), m50.sum())
    check("X(d>50)", acc.X[sel].sum(), X[m50].sum())
    check("P(d>50)", acc.P[sel].sum(), P[m50].sum())
    check("R(d>50)", acc.R[sel].sum(), R[m50].sum())
    # x_HI > 0.5: the 8 clump cells and the neutral single cell
    mx = xhi > 0.5
    lgx_lo = np.concatenate([[-np.inf], LGX_EDGES])
    selx = lgx_lo >= np.log10(0.5) - 1e-9
    check("n(x>0.5)", acc.n[:, selx].sum(), mx.sum())
    check("P(x>0.5)", acc.P[:, selx].sum(), P[mx].sum())
    # removed with x = 1 against the rule: the clump cells remove e^-1 of theirs
    clump = delta == 100
    check("X clump", X[clump].sum(), math.exp(-1.0) * P[clump].sum())
    check("X d=200", float(X[8, 8, 8]), math.exp(-3.0) * float(P[8, 8, 8]))
    # the medians: at delta = 1 the median Gamma is 1e-12 and x_HI 1e-4
    i1 = int(np.argmin(np.abs(np.array(
        [LGD_EDGES[0] - 0.05] + list(0.5 * (LGD_EDGES[1:] + LGD_EDGES[:-1])) + [LGD_EDGES[-1] + 0.05]))))
    gm = _median_from_hist(acc.hg[i1], LGG_EDGES)
    xm = _median_from_hist(acc.hx[i1], LGXF_EDGES)
    check("median Gamma(d=1)", gm, 1e-12, 0.05)
    check("median x_HI(d=1)", xm, 1e-4, 0.13)
    # the opacity share of the cells with tau_912 > 1 (the clump and the
    # neutral cell, here)
    tau = SIGMA_912 * nHI * dl
    check("tau>1 share", acc.tau_gt1_nHI, (delta * xhi)[tau > 1].sum())
    check("n tau>1", acc.n_tau_gt1, int((tau > 1).sum()))
    # the components (scipy)
    try:
        import scipy  # noqa: F401
        have_scipy = True
    except ImportError:
        have_scipy = False
    if have_scipy:
        nHI_dl = nHI * dl
        r = absorber_stats(nHI_dl, xhi > 0.5)
        check("ncomp", r["ncomp"], 2)
        sizes = sorted(r["sizes"].tolist())
        check("sizes", sizes[0], 1)
        check("sizes clump", sizes[1], 8)
        # the clump's peak column along any axis: two cells of n_HI dl
        big = int(np.argmax(r["sizes"]))
        want_col = 2 * 0.9 * 100 * nH_avg * dl
        for ax in range(3):
            check(f"peak column axis {ax}", float(r["peak_col"][ax, big]), want_col, 1e-9)
        small = 1 - big
        check("single-cell column", float(r["peak_col"][0, small]), 0.99 * nH_avg * dl, 1e-9)
        # a component touching only one face is not merged with another
        mask2 = np.zeros((N, N, N), bool)
        mask2[0, 3, 3] = True
        mask2[N - 1, 7, 7] = True
        lab2, n2 = periodic_components(mask2)
        check("no false merge", n2, 2)
        mask2[N - 1, 3, 3] = True
        lab3, n3 = periodic_components(mask2)
        check("merge across face", n3, 2)
    else:
        print("scipy not available: the component checks were skipped")
    if errs:
        print("absorber_budget selftest FAILED:\n  " + "\n  ".join(errs))
        return 1
    print("absorber_budget selftest OK" + ("" if have_scipy else " (components skipped)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
