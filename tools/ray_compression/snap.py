"""Load ray snapshots written by the RT's optional ray dump
(include/radray/ray_dump.hpp) and derive the quantities used by the
compression study (docs/ray_compression.md).

    a, meta = load("dir/rays_s0000120")       # .npy + .json
    d = derive(a, meta)                         # dict of numpy arrays
"""
import json
import math

import numpy as np

try:
    import healpy as hp
except ImportError:  # pragma: no cover
    hp = None

NFREQ = 5
HI, HeI, HeII = 13.6, 24.59, 54.42  # eV, as pc::HI_eV ... (band edges only)
FREQ_EDGES = np.array([HI, math.sqrt(HI * HeI), HeI, math.sqrt(HeI * HeII), HeII, 100.0])
FREQ_INIT = np.array([15.71, 21.06, 29.64, 42.88, 74.12], dtype=np.float32)
RAY_PER_CELL = 1.0
LSIDE_MIN, LSIDE_MAX = 2, 13
FLAG_MERGED, FLAG_MIXED, FLAG_CLAMPED = 1, 2, 4


def load(stem):
    stem = str(stem)
    if stem.endswith(".npy"):
        stem = stem[:-4]
    a = np.load(stem + ".npy")
    with open(stem + ".json") as f:
        meta = json.load(f)
    assert len(a) == meta["nrays"]
    return a, meta


def rmax_of_nside(nside):
    return nside * math.sqrt(3 / (math.pi * RAY_PER_CELL))


def rmin_of_nside(nside):
    nside = np.asarray(nside)
    return np.where(nside == 4, 0.0, (nside // 2) * math.sqrt(3 / (math.pi * RAY_PER_CELL)))


def lside_of_r(r):
    """ray_tracing_tools lside_from_radius (in double; ties may differ from
    the real(4) log by one level, which is what we want to measure)."""
    r = np.asarray(r, dtype=np.float64)
    with np.errstate(divide="ignore"):
        l = np.ceil(np.log(math.pi * RAY_PER_CELL / 3.0 * r * r) / np.log(4.0))
    l = np.where(r < 1, LSIDE_MIN, l)
    return np.clip(l, LSIDE_MIN, LSIDE_MAX).astype(np.int64)


def tangent_basis(c):
    """Unit (e_theta, e_phi) at unit vectors c (n,3)."""
    z = c[:, 2]
    st = np.sqrt(np.maximum(1 - z * z, 0))
    phi = np.arctan2(c[:, 1], c[:, 0])
    cp, sp = np.cos(phi), np.sin(phi)
    e_th = np.stack([z * cp, z * sp, -st], axis=1)
    e_ph = np.stack([-sp, cp, np.zeros_like(sp)], axis=1)
    return e_th, e_ph


def pixel_radius(nside):
    """sqrt(pixel solid angle): the angular size of one pixel (rad)."""
    return np.sqrt(4 * np.pi / (12.0 * np.asarray(nside, dtype=np.float64) ** 2))


def derive(a, meta):
    N = meta["N"]
    d = {}
    cell = a["cell"].astype(np.int64)
    ijk = np.stack([cell % N, (cell // N) % N, cell // (N * N)], axis=1)
    x = a["x"].astype(np.float64)
    off = x - ijk
    # A ray is chained in the cell of its stored float x; off is in [0,1)
    d["ijk"] = ijk
    d["off"] = off
    nside = a["nside"].astype(np.int64)
    lside = np.round(np.log2(nside)).astype(np.int64)
    d["lside"] = lside
    d["uniq"] = 4 * nside * nside + a["ipix"].astype(np.int64)
    # Direction relative to the pixel centre, in the local tangent basis,
    # in units of the pixel size
    cen = np.array(hp.pix2vec(nside, a["ipix"].astype(np.int64), nest=True)).T
    dx = a["dx"].astype(np.float64)
    dev = dx - cen
    e_th, e_ph = tangent_basis(cen)
    pr = pixel_radius(nside)
    d["cen"] = cen
    d["dev_th"] = (dev * e_th).sum(1) / pr
    d["dev_ph"] = (dev * e_ph).sum(1) / pr
    d["dev_abs"] = np.sqrt((dev * dev).sum(1))
    d["dx_norm_err"] = np.sqrt((dx * dx).sum(1)) - 1
    # r and dr
    r = a["r"].astype(np.float64)
    d["r"] = r
    rmax = nside * math.sqrt(3 / (math.pi * RAY_PER_CELL))
    d["dr_pred"] = rmax - r
    d["dr_resid"] = a["dr"].astype(np.float64) - (rmax - r)
    d["lside_r"] = lside_of_r(r)
    # photons
    ng = a["ngam"].astype(np.float64)
    ntot = ng.sum(1)
    d["ntot"] = ntot
    with np.errstate(divide="ignore", invalid="ignore"):
        d["shape"] = ng / ntot[:, None]
    # nu position in band [0,1]
    lo, hi = FREQ_EDGES[:-1], FREQ_EDGES[1:]
    d["nu_pos"] = (a["nu"].astype(np.float64) - lo) / (hi - lo)
    d["nu_is_init"] = a["nu"] == FREQ_INIT
    # provenance
    flags = a["flags"]
    d["merged"] = (flags & FLAG_MERGED) != 0
    d["clamped"] = (flags & FLAG_CLAMPED) != 0
    d["never_merged"] = ~d["merged"]
    # origin-based prediction of position for never-merged rays
    src = a["src"].astype(np.int64)
    ok = src >= 0
    sijk = np.stack([src % N, (src // N) % N, src // (N * N)], axis=1).astype(np.float32)
    xs = sijk + np.float32(0.5)  # (irt,jrt,krt)-0.5 in 1-based = ijk+0.5
    d["src_x"] = xs.astype(np.float64)
    xp = (xs.astype(np.float64) + r[:, None] * cen) % N
    res = x - xp
    res = (res + N / 2) % N - N / 2
    res[~ok] = np.nan
    d["x_resid_cen"] = res  # x - (x_src + r*centre(ipix))
    xp2 = (xs.astype(np.float64) + r[:, None] * dx) % N
    res2 = (x - xp2 + N / 2) % N - N / 2
    res2[~ok] = np.nan
    d["x_resid_dx"] = res2  # x - (x_src + r*dx)
    return d


def tau_ctx(a):
    """Cell optical-depth context: min over groups of tau, log10."""
    t = a["c_tau"].astype(np.float64)
    return np.log10(np.maximum(t.min(1), 1e-12))
