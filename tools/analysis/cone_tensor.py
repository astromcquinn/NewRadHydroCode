#!/usr/bin/env python3
"""Cone-traced, opacity-attenuated Eddington tensor (prototype for the a
priori comparison of docs/otvet_closure.md; numpy, not the solver).

The direct radiation of the sources at a receiver x, per HEALPix direction
n_p (cone of solid angle Omega = 4 pi / N_pix):

    I_p(x) = sum over steps of  eta(x + s n_p) Omega ds / (4 pi) * T(N_HI(<s))
    P(x)   = P_near(x) + sum_p I_p n_p n_p^T

with eta and the HI column read from mip-maps (level with cells ~ the cone
diameter, trilinear, periodic, so the periodic images out to s_max are
included), and T(N_HI) = sum_g W_g exp(-sigma_g N_HI) / sum_g W_g the
Gamma-weighted transmission of the source spectrum. Receivers are the cells
of a grid --coarse times coarser; the cone part is prolonged to the RT grid.
The near field (r < r_c) is the thin minimum-image sum on the RT grid, with
a smooth hand-over between r_c and 2 r_c (weights w and 1 - w).

Units: rho = Q/(4 pi) per source cell (Q photons per substep, Gamma-
weighted, as in tensor_compare.py), so P is in the rays' S0 units.
"""
import numpy as np
from scipy import ndimage


def healpix_dirs(nside):
    """Pixel centres of a HEALPix nested-free (ring) scheme, unit vectors."""
    npix = 12 * nside * nside
    out = np.empty((npix, 3))
    p = 0
    # polar caps and equatorial belt (ring scheme)
    for i in range(1, 4 * nside):
        if i < nside:
            z = 1 - i * i / (3.0 * nside * nside)
            nr = 4 * i
            shift = 0.5
        elif i <= 3 * nside:
            z = 4.0 / 3 - 2.0 * i / (3.0 * nside)
            nr = 4 * nside
            shift = 0.5 * ((i - nside + 1) % 2)
        else:
            ii = 4 * nside - i
            z = -(1 - ii * ii / (3.0 * nside * nside))
            nr = 4 * ii
            shift = 0.5
        st = np.sqrt(max(0.0, 1 - z * z))
        for j in range(nr):
            phi = 2 * np.pi * (j + shift) / nr
            out[p] = (st * np.cos(phi), st * np.sin(phi), z)
            p += 1
    assert p == npix
    return out


def pyramid(a, nlev):
    """a (N,N,N) -> [a, mean 2^3, ...] down to nlev levels."""
    lev = [a]
    for _ in range(nlev - 1):
        b = lev[-1]
        n = b.shape[0] // 2
        if n < 1:
            break
        lev.append(b.reshape(n, 2, n, 2, n, 2).mean(axis=(1, 3, 5)))
    return lev


def smoothstep(r, r0, r1):
    t = np.clip((r - r0) / (r1 - r0), 0, 1)
    return t * t * (3 - 2 * t)


def cone_tensor(rho, col, sigma18, W, nside=2, coarse=4, rc=12.0, smax_boxes=8.0,
                step=0.5, tmin=1e-4, near_kernel=None, workers=8, lam_extra=(None,), sbins=None):
    """rho: (N,N,N) Q/(4 pi) per cell; col: (N,N,N) HI column per cell in
    1e-18 cm^-2 (tau_g = sigma18_g * col); W: group weights. lam_extra: a
    list of extra screening lengths in cells (None: none), each multiplying
    the transmission by exp(-s / lam): the attenuation the HI column does not
    give (light horizon, redshifting out of the band). Returns, per
    lam_extra, the cone part on the RT grid, shape (N,N,N,6), and on the
    coarse grid. sbins: distance bin edges (cells); if given, also returns
    the unscreened cone part per bin on the coarse grid, shape
    (nbins, nc, nc, nc, 6), so that any extra screening can be applied
    afterwards bin by bin (budget_screen)."""
    N = rho.shape[0]
    W = np.asarray(W, float) / np.sum(W)
    s18 = np.asarray(sigma18, float)
    dirs = healpix_dirs(nside)
    npix = len(dirs)
    omega = 4 * np.pi / npix
    theta = np.arccos(1 - omega / (2 * np.pi))  # cone half angle
    nlev = int(np.log2(N)) + 1
    R = pyramid(rho, nlev)
    C = pyramid(col, nlev)
    nc = N // coarse
    g = (np.arange(nc) + 0.5) * coarse  # receiver centres, fine-cell units
    zz, yy, xx = np.meshgrid(g, g, g, indexing="ij")
    X = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], 1)
    lam_extra = list(lam_extra)
    Pc = np.zeros((len(lam_extra), X.shape[0], 6))
    smax = smax_boxes * N
    nb = 0 if sbins is None else len(sbins) - 1
    Pb = np.zeros((nb, X.shape[0], 6))

    def sample(levels, l, pos):
        a = levels[l]
        f = 2 ** l
        # array index order [z, y, x]; cell centres at (i + 0.5) f
        coords = np.stack([pos[:, 2] / f - 0.5, pos[:, 1] / f - 0.5, pos[:, 0] / f - 0.5])
        return ndimage.map_coordinates(a, coords, order=1, mode="grid-wrap")

    def trans(c):
        return (W[None, :] * np.exp(-np.outer(c, s18))).sum(1)

    for p in range(npix):
        n = dirs[p]
        nn = [n[0] * n[0], n[1] * n[1], n[2] * n[2], n[0] * n[1], n[0] * n[2], n[1] * n[2]]
        s = rc
        colacc = np.zeros(X.shape[0])
        I = np.zeros((len(lam_extra), X.shape[0]))
        Ib = np.zeros((nb, X.shape[0]))
        alive = np.ones(X.shape[0], bool)
        while s < smax:
            d = max(2 * s * np.tan(theta), 1.0)
            ds = step * d
            sm = s + ds / 2
            l = int(np.clip(np.floor(np.log2(d)), 0, nlev - 1))
            pos = X + sm * n[None, :]
            r = sample(R, l, pos)
            c = sample(C, l, pos)
            w = smoothstep(np.array([sm]), rc, 2 * rc)[0]  # hand-over from the near field
            T = trans(colacc + c * ds / 2)
            for j, lx in enumerate(lam_extra):
                I[j] += w * r * omega * ds * T * (np.exp(-sm / lx) if lx else 1.0)
            if nb:
                Ib[min(max(np.searchsorted(sbins, sm, "right") - 1, 0), nb - 1)] += w * r * omega * ds * T
            colacc += c * ds
            s += ds
            if T.max() < tmin:
                break
        for m in range(6):
            Pc[:, :, m] += I * nn[m]
            if nb:
                Pb[:, :, m] += Ib * nn[m]
    Pc = Pc.reshape(len(lam_extra), nc, nc, nc, 6)
    Ps = [prolong(Pc[j], N) for j in range(len(lam_extra))]
    if sbins is not None:
        return Ps, Pc, Pb.reshape(nb, nc, nc, nc, 6)
    return Ps, Pc


def prolong(Pc, N):
    """Coarse (nc,nc,nc,6) to the RT grid (N,N,N,6), trilinear, periodic."""
    coarse = N // Pc.shape[0]
    f = np.arange(N)
    cz = (f + 0.5) / coarse - 0.5
    Z, Y, Xg = np.meshgrid(cz, cz, cz, indexing="ij")
    coords = np.stack([Z.ravel(), Y.ravel(), Xg.ravel()])
    P = np.empty((N, N, N, 6))
    for m in range(6):
        P[..., m] = ndimage.map_coordinates(Pc[..., m], coords, order=1, mode="grid-wrap").reshape(N, N, N)
    return P


def budget_screen(Pb, sbins, near_trace, target, coarse):
    """The extra screening length lam_x (cells) for which the direct light
    sum_cells tr P_dir = near_trace + coarse^3 sum_b tr Pb_b exp(-s_b/lam_x)
    equals target (the photons the solution holds; here the rays' sum S0),
    with s_b the geometric centre of bin b. None (no screening) if the
    unscreened light is already <= target. Returns (lam_x, Pc screened)."""
    sb = np.sqrt(sbins[:-1] * sbins[1:])
    tb = coarse ** 3 * (Pb[..., 0] + Pb[..., 1] + Pb[..., 2]).reshape(len(sb), -1).sum(1)

    def light(lx):
        return near_trace + (tb * (np.exp(-sb / lx) if lx else 1.0)).sum()

    if light(None) <= target:
        lx = None
    elif near_trace >= target:
        lx = 1.0  # the unscreened near field alone exceeds the budget
    else:
        lo, hi = 0.0, 12.0  # log10 lam_x in cells; light(lam) increases with lam
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if light(10 ** mid) > target:
                hi = mid
            else:
                lo = mid
        lx = 10 ** (0.5 * (lo + hi))
    fac = np.ones(len(sb)) if lx is None else np.exp(-sb / lx)
    return lx, np.tensordot(fac, Pb, axes=(0, 0))
