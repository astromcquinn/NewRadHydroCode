#!/usr/bin/env python3
"""Eddington tensors from the rays against closures (docs/otvet_closure.md).

usage: tensor_compare.py STEP [STEP ...] --out DIR [--block 8] [--lam L ...]
                         [--moment RUNSTEP:label ...] [--workers 16]

STEP is a step_NNNN directory of a ray run made with tensor_dump = 1
(tensor.bin, tensor_odd.bin, otvet_h.bin, sources.txt, fHI.bin, rhog.bin,
gammaHI.bin). The true tensor is h = S2/S0 of the rays (Gamma_HI-weighted,
averaged over the substeps of the snapshot step); its sampling noise is
estimated from the odd and even substeps separately.

The closures are evaluated a priori, on the rays' own radiation field (S0,
and S1 for M1):
  otvet     the code's optically thin tensor (otvet_h.bin)
  m1        Levermore's M1 closure of the rays' own f = S1/S0
  iso       I/3
  split     P = N I/3 + min(N, N_loc) (h_loc - I/3): the direct light of the
            sources, N_loc = sum_s w_s / (4 pi r^2) (minimum image, in the
            rays' units: sources.txt qfac and sigma18), carries the thin
            tensor h_loc and the rest of the radiation is isotropic
  atten     the thin sum with exp(-r/lambda) over all periodic images (27
            explicit, the rest as an isotropic continuum), lambda in cells
            (--lam; 'auto' = lambda_eff = sum S0 / sum_s w_s, the photons'
            mean path, from the rays' own budget; 'auto*f' = f lambda_eff)
  split_atten  split with the attenuated, image-summed N_loc and h_loc
  coneNcC   (--cone N:C[:cap[:rc]]) the direct light attenuated by the HI along
            HEALPix cones (nside N) from receivers on a C-times coarser grid
            through mip-maps of the emissivity and the HI column, with the
            periodic images out to 8 boxes, plus the thin near field
            (r < rc, default 2C cells, handed over to the cones by 2 rc;
            cone_tensor.py); cap limits each cell's group-1 tau (clumps
            saturate: a crude covering-fraction correction of the mip average);
            suffix x<lam>: also screened by exp(-s/lam) (--cone-lamx), the
            attenuation the HI does not give (light horizon, redshifting)
  split_cone  split with the cone-traced direct light
  blend_A_B_s (--blend A:B:s) h_A + G_s * (h_B - h_A): the small scales of
            closure A and the large ones of B (Gaussian of s cells)

Metrics, S0-weighted over the cells of each regime:
  err   the Frobenius distance |h - h_true| (off-diagonals counted twice)
  A     the anisotropy |h - I/3| / sqrt(2/3) (0 isotropic, 1 one beam)
  noise the same distance between the odd- and even-substep h_true / 2
The block metrics use the tensors summed over --block^3 cells (S0 h) before
normalising: the large-scale tensor that sets P_Gamma at low k, nearly free
of ray noise. Scale by scale: the radiation force D = div(S0 h), the only
way the closure enters the flux equation; per k bin (factors of 2 from the
fundamental) the error power |D_c - D_true|^2 / |D_true|^2 and the
cross-correlation with the truth.

--moment RUNSTEP:label adds a moment run's own tensor.bin at the nearest
redshift (a posteriori; its own N, so cell-by-cell only where the
morphology agrees).

Writes DIR/tensor_compare.txt and DIR/tensor_zXX.png per step.
"""
import argparse
import os

import sys

import numpy as np
import scipy.fft as sfft

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

COMP = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]  # xx yy zz xy xz yz
WOFF = np.array([1, 1, 1, 2, 2, 2], dtype=float)  # Frobenius weights


def read_sources(d):
    hdr = {}
    rows = []
    for line in open(os.path.join(d, "sources.txt")):
        if line.startswith("#"):
            t = line[1:].split()
            if t and t[0] in ("qfac", "sigma18"):
                hdr[t[0]] = np.array([float(x) for x in t[1:]])
            elif t and t[0] == "redshift":
                hdr["z"] = float(t[1])
                hdr["N"] = int(t[3])
            elif "substeps" in t:
                hdr["nsub"] = int(t[t.index("substeps") + 1])
                hdr["nodd"] = int(t[t.index("odd") + 1]) if "odd" in t else 0
            continue
        rows.append([float(x) for x in line.split()])
    src = np.array(rows).reshape(-1, 4)
    return hdr, src


def load(d, N, name, ncomp):
    a = np.fromfile(os.path.join(d, name), dtype=np.float32)
    return a.reshape(N, N, N, ncomp).astype(np.float64) if ncomp > 1 else a.reshape(N, N, N).astype(np.float64)


def tensor_of(S):
    """(S0, F (3), P (6)) from a 10-component dump."""
    return S[..., 0], S[..., 1:4], S[..., 4:10]


def hi_column(d, N, z, fHI, box=None):
    """HI column per RT cell, 1e-18 cm^-2 (tau_g = sigma18_g * col), from
    rhog and fHI (cosmology of include/radray/constants.hpp)."""
    if box is None:
        pf = os.path.join(d, "..", "..", "params.txt")
        for line in open(pf):
            t = line.split("=")
            if len(t) == 2 and t[0].strip() == "box":
                box = float(t[1])
    rhoc, mpc, h0, ob, XHy, mH = 1.87890e-29, 3.08560e24, 0.70, 0.045, 0.76, 1.67223e-24
    a = 1 / (1 + z)
    nH = rhoc * ob * h0 ** 2 / a ** 3 * XHy / mH
    len_rt = box * mpc / h0 / N * a
    rhog = load(d, N, "rhog.bin", 1)
    return nH * rhog * fHI * len_rt * 1e-18


def m1_D(f3):
    fn = np.sqrt((f3 ** 2).sum(-1))
    f = np.minimum(fn, 1.0)
    chi = (3 + 4 * f * f) / (5 + 2 * np.sqrt(4 - 3 * f * f))
    a = (1 - chi) / 2
    b = (3 * chi - 1) / 2
    with np.errstate(invalid="ignore", divide="ignore"):
        u = np.where(fn[..., None] > 0, f3 / fn[..., None], 0.0)
    D = np.empty(f3.shape[:-1] + (6,))
    for m, (p, q) in enumerate(COMP):
        D[..., m] = b * u[..., p] * u[..., q] + (a if p == q else 0.0)
    iso = fn == 0
    D[iso] = [1 / 3, 1 / 3, 1 / 3, 0, 0, 0]
    return D


def thin_kernel(N, lam=None, images=False):
    """K_m(d) = sum exp(-r/lam) d_p d_q / r^4 (self cell 4/3 I), and the
    trace kernel 1/r^2; arrays [z, y, x] with d = minimum image (or all
    images with images=True)."""
    ax = np.arange(N)
    ax = np.where(ax > N // 2, ax - N, ax).astype(np.float64)
    z, y, x = np.meshgrid(ax, ax, ax, indexing="ij")
    shifts = [(0, 0, 0)]
    if images:
        shifts = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    K = np.zeros((6, N, N, N))
    half = N / 2
    for (sa, sb, sc) in shifts:
        dx, dy, dz = x + sa * N, y + sb * N, z + sc * N
        r2 = dx * dx + dy * dy + dz * dz
        with np.errstate(divide="ignore", invalid="ignore"):
            w = 1.0 / (r2 * r2)
            if lam:
                w = w * np.exp(-np.sqrt(r2) / lam)
        w[r2 == 0] = 0
        d = (dx, dy, dz)
        for m, (p, q) in enumerate(COMP):
            t = w * d[p] * d[q]
            if not images and p != q:
                # the code: a component at exactly N/2 has two images of
                # opposite sign
                t[(np.abs(d[p]) == half) | (np.abs(d[q]) == half)] = 0
            K[m] += t
    # self cell: isotropic 4 (1/r0^2, r0 = half a cell)
    K[0, 0, 0, 0] += 4 / 3
    K[1, 0, 0, 0] += 4 / 3
    K[2, 0, 0, 0] += 4 / 3
    if images and lam:
        # images beyond the 3^3 block: isotropic continuum with the source
        # density 1/N^3, from R (the radius of a sphere of volume (3N)^3)
        R = 3 * N * (3 / (4 * np.pi)) ** (1 / 3)
        bg = 4 * np.pi * lam * np.exp(-R / lam) / N ** 3
        for m in range(3):
            K[m] += bg / 3
    return K


def convolve(rho, K, workers):
    fr = sfft.rfftn(rho, workers=workers)
    out = np.empty(K.shape)
    for m in range(K.shape[0]):
        out[m] = sfft.irfftn(fr * sfft.rfftn(K[m], workers=workers), s=rho.shape, workers=workers)
    return out


def normalise(P):
    tr = P[..., 0] + P[..., 1] + P[..., 2]
    with np.errstate(invalid="ignore", divide="ignore"):
        h = P / tr[..., None]
    h[~(tr > 0)] = [1 / 3, 1 / 3, 1 / 3, 0, 0, 0]
    return h


def fro(a, b):
    return np.sqrt(((a - b) ** 2 * WOFF).sum(-1))


def aniso(h):
    iso = np.array([1 / 3, 1 / 3, 1 / 3, 0, 0, 0])
    return fro(h, iso) / np.sqrt(2 / 3)


def block_sum(a, b):
    N = a.shape[0]
    n = N // b
    sh = (n, b, n, b, n, b) + a.shape[3:]
    return a.reshape(sh).sum(axis=(1, 3, 5))


def wmean(v, w, m):
    ww = w[m]
    return float((v[m] * ww).sum() / ww.sum()) if ww.sum() > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("steps", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--block", type=int, default=8)
    ap.add_argument("--lam", nargs="*", default=["auto"])
    ap.add_argument("--moment", nargs="*", default=[])
    ap.add_argument("--blend", nargs="*", default=[],
                    help="A:B:sigma -> h_A + G_sigma * (h_B - h_A): A's small scales, B's large ones "
                         "(xauto in a name: x<lambda_eff>)")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--cone", nargs="*", default=[],
                    help="cone-traced attenuated tensors nside[:coarse[:tau_cap[:r_c]]] (cone_tensor.py)")
    ap.add_argument("--cone-lamx", nargs="*", default=["inf"],
                    help="extra screening lengths for the cones, exp(-s/lam): inf (none), cells, "
                         "auto or auto*f (lambda_eff times f); one tensor each, label suffix x<lam>")
    ap.add_argument("--box", type=float, default=None, help="Mpc/h (default: the run's params.txt)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    txt = open(os.path.join(a.out, "tensor_compare.txt"), "w")

    def out(s=""):
        print(s)
        txt.write(s + "\n")
        txt.flush()

    for d in a.steps:
        hdr, src = read_sources(d)
        N, z = hdr["N"], hdr["z"]
        S = load(d, N, "tensor.bin", 10)
        S0, F, P = tensor_of(S)
        h_true = normalise(P)
        # noise: odd vs even substeps
        n, no = hdr["nsub"], hdr["nodd"]
        So = load(d, N, "tensor_odd.bin", 10)
        Se = (S.astype(np.float64) * n - So * no) / max(n - no, 1)
        # float32 cancellation in Se can leave a cell's half nearly empty:
        # clip to the largest distance of two trace-1 tensors (sqrt 2) / 2
        noise = np.minimum(fro(normalise(So[..., 4:10]), normalise(Se[..., 4:10])) / 2, np.sqrt(0.5))
        del So, Se
        fHI = load(d, N, "fHI.bin", 1)
        x = 1 - fHI
        gam = load(d, N, "gammaHI.bin", 1) if os.path.exists(os.path.join(d, "gammaHI.bin")) else None
        with np.errstate(invalid="ignore", divide="ignore"):
            f3 = np.where(S0[..., None] > 0, F / S0[..., None], 0.0)
        # sources: weights in the rays' units (per substep, Gamma-weighted)
        wfac = float((hdr["qfac"] * hdr["sigma18"]).sum())
        rho = np.zeros((N, N, N))
        ii, jj, kk = (src[:, 0].astype(int) - 1, src[:, 1].astype(int) - 1, src[:, 2].astype(int) - 1)
        np.add.at(rho, (kk, jj, ii), src[:, 3] * wfac / (4 * np.pi))
        lam_eff = S0.sum() / (rho.sum() * 4 * np.pi) if rho.sum() > 0 else float("nan")

        cl = {}
        cl["otvet"] = load(d, N, "otvet_h.bin", 6)
        cl["m1"] = m1_D(f3)
        cl["iso"] = np.broadcast_to(np.array([1 / 3, 1 / 3, 1 / 3, 0, 0, 0]), h_true.shape)
        Ploc = np.moveaxis(convolve(rho, thin_kernel(N), a.workers), 0, -1)
        Nloc = Ploc[..., 0] + Ploc[..., 1] + Ploc[..., 2]
        hloc = normalise(Ploc)
        iso = np.array([1 / 3, 1 / 3, 1 / 3, 0, 0, 0])

        def split(N_, hl):
            m = np.minimum(S0, N_)
            Pc = S0[..., None] * iso + m[..., None] * (hl - iso)
            return normalise(Pc)

        cl["split"] = split(Nloc, hloc)
        lams = []
        for l in a.lam:
            if l.startswith("auto"):
                lams.append(lam_eff * (float(l[5:]) if l.startswith("auto*") else 1.0))
            else:
                lams.append(float(l))
        for l in lams:
            Pa = np.moveaxis(convolve(rho, thin_kernel(N, lam=l, images=True), a.workers), 0, -1)
            Na = Pa[..., 0] + Pa[..., 1] + Pa[..., 2]
            cl["atten%.0f" % l] = normalise(Pa)
            cl["split_atten%.0f" % l] = split(Na, normalise(Pa))
            del Pa
        if a.cone:
            from cone_tensor import cone_tensor, smoothstep
            col = hi_column(d, N, z, fHI, a.box)
            s18 = hdr["sigma18"]
            Wg = hdr["qfac"] * s18
            ax_ = np.arange(N)
            ax_ = np.where(ax_ > N // 2, ax_ - N, ax_).astype(np.float64)
            zr, yr, xr = np.meshgrid(ax_, ax_, ax_, indexing="ij")
            rr = np.sqrt(xr * xr + yr * yr + zr * zr)
            del zr, yr, xr
        for spec in a.cone:
            t = spec.split(":")
            ns = int(t[0])
            co = int(t[1]) if len(t) > 1 else 4
            cap = float(t[2]) if len(t) > 2 else float("inf")
            # near field r < 2 rc (unattenuated), cones beyond rc
            rc = float(t[3]) if len(t) > 3 else 2.0 * co
            colc = col if not np.isfinite(cap) else np.minimum(col, cap / s18[0])
            Kn = thin_kernel(N) * (1 - smoothstep(rr, rc, 2 * rc))[None]
            Pn = np.moveaxis(convolve(rho, Kn, a.workers), 0, -1)
            del Kn
            lamx = [None if l == "inf" else lam_eff * (float(l[5:]) if l.startswith("auto*") else 1.0)
                    if l.startswith("auto") else float(l) for l in a.cone_lamx]
            Pks, _ = cone_tensor(rho, colc, s18, Wg, nside=ns, coarse=co, rc=rc, workers=a.workers,
                                 lam_extra=lamx)
            for lx, Pk in zip(lamx, Pks):
                Pd = Pn + Pk
                Nd = Pd[..., 0] + Pd[..., 1] + Pd[..., 2]
                lab = (f"cone{ns}c{co}" + ("" if not np.isfinite(cap) else f"t{cap:g}")
                       + ("" if len(t) < 4 else f"r{rc:g}") + ("" if lx is None else f"x{lx:.0f}"))
                cl[lab] = normalise(Pd)
                cl["split_" + lab] = split(Nd, cl[lab])
                out(f"# {lab}: sum N_dir / sum S0 = {Nd.sum() / S0.sum():.3f}")
                del Pd
            del Pn, Pks
        for spec in a.blend:
            # small scales of closure A, large scales (Gaussian sigma cells) of B
            ka, kb, sg = (t_.replace("xauto", "x%.0f" % lam_eff) for t_ in spec.split(":"))
            if ka in cl and kb in cl:
                from scipy.ndimage import gaussian_filter
                dh = cl[kb] - cl[ka]
                for m in range(6):
                    dh[..., m] = gaussian_filter(dh[..., m], float(sg), mode="wrap")
                cl[f"blend_{ka}_{kb}_{sg}"] = cl[ka] + dh
                del dh
            else:
                out(f"# blend {spec}: no closure {ka if ka not in cl else kb}")
        for spec in a.moment:
            path, lab = spec.rsplit(":", 1)
            cl["post_" + lab] = normalise(load(path, N, "tensor.bin", 10)[..., 4:10])

        # regimes
        gi = gam[x > 0.9].mean() if gam is not None and (x > 0.9).any() else 1.0
        reg = {
            "all": S0 > 0,
            "front 0.1<x<0.9": (S0 > 0) & (x > 0.1) & (x < 0.9),
            "ionised x>0.9": (S0 > 0) & (x > 0.9),
        }
        if gam is not None:
            g = gam / gi
            reg["ion G<0.5<G>"] = reg["ionised x>0.9"] & (g < 0.5)
            reg["ion 0.5-2"] = reg["ionised x>0.9"] & (g >= 0.5) & (g < 2)
            reg["ion 2-10"] = reg["ionised x>0.9"] & (g >= 2) & (g < 10)
            reg["ion >10 (sources)"] = reg["ionised x>0.9"] & (g >= 10)
        with np.errstate(invalid="ignore", divide="ignore"):
            q = np.where(Nloc > 0, S0 / Nloc, np.inf)
        reg["N<N_loc/2 (shadowed)"] = (S0 > 0) & (q < 0.5)
        reg["N>2N_loc (diffuse)"] = (S0 > 0) & (q > 2)

        out(f"## {d}")
        out(f"z = {z:.3f}  N = {N}  sources = {len(src)}  x_V = {x.mean():.4f}  "
            f"lambda_eff = {lam_eff:.1f} cells  <N/N_loc>_ion = {np.median(q[reg['ionised x>0.9']]) if reg['ionised x>0.9'].any() else float('nan'):.3f} (median)")
        names = list(cl)
        # the same numbers, untruncated, for tensor_summary.py
        M = {"z": z, "xv": float(x.mean()), "lam": float(lam_eff), "cell": {}, "Atrue": {}, "noise": {},
             "A": {}, "force": {}, "r": {}, "block": {}, "blockA": {}}
        out("cell level: S0-weighted err |h - h_true| (A = anisotropy; truth A and noise first)")
        out("%-24s %7s %7s %7s " % ("regime", "frac", "A_true", "noise") + " ".join("%10s" % k[:10] for k in names))
        for r, m in reg.items():
            if not m.any():
                continue
            row = "%-24s %7.4f %7.3f %7.3f " % (r, (S0 * m).sum() / (S0 * (S0 > 0)).sum(), wmean(aniso(h_true), S0, m),
                                             wmean(noise, S0, m))
            M["cell"][r] = {k: float(wmean(fro(cl[k], h_true), S0, m)) for k in names}
            M["Atrue"][r] = float(wmean(aniso(h_true), S0, m))
            M["noise"][r] = float(wmean(noise, S0, m))
            row += " ".join("%10.3f" % M["cell"][r][k] for k in names)
            out(row)
        out("anisotropy A of each closure (S0-weighted)")
        for r, m in reg.items():
            if not m.any():
                continue
            M["A"][r] = {k: float(wmean(aniso(cl[k]), S0, m)) for k in names}
            out("%-24s " % r + " ".join("%10.3f" % M["A"][r][k] for k in names))
        # scale by scale: the radiation force div(N h) that drives the flux
        # in the moment equations; error spectrum |D_c - D_true|^2 / |D_true|^2
        # and the cross-correlation r(k), in k bins [h/Mpc] (box from params)
        Pt_hat = [sfft.rfftn(P[..., m], workers=a.workers) for m in range(6)]
        kf = 2 * np.pi * np.fft.fftfreq(N)  # rad per cell
        kr = 2 * np.pi * np.fft.rfftfreq(N)
        KZ, KY, KX = np.meshgrid(kf, kf, kr, indexing="ij")
        kmag = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
        idx = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (0, 1): 3, (1, 0): 3, (0, 2): 4, (2, 0): 4, (1, 2): 5, (2, 1): 5}
        kv = (KX, KY, KZ)

        def divP(Ph):
            return [sum(1j * kv[j] * Ph[idx[(i, j)]] for j in range(3)) for i in range(3)]
        Dt = divP(Pt_hat)
        del Pt_hat
        boxl = a.box
        if boxl is None:
            for line in open(os.path.join(d, "..", "..", "params.txt")):
                t = line.split("=")
                if len(t) == 2 and t[0].strip() == "box":
                    boxl = float(t[1])
        kunit = N / boxl  # grid k (rad/cell) -> h/Mpc
        kb = np.array([0.5, 1, 2, 4, 8, 16, 32, 64]) * 2 * np.pi / N * 2
        kb = kb[kb < np.pi * np.sqrt(3)]
        mult = np.where(KX == 0, 1.0, 2.0)
        out("radiation force div(S0 h): error power / true power by k (and r(k) with the truth)")
        hdrk = "%-16s" % "k [h/Mpc]" + " ".join("%9.2f" % (0.5 * (kb[i] + kb[i + 1]) * kunit) for i in range(len(kb) - 1))
        out(hdrk)
        M["k"] = [float(0.5 * (kb[i] + kb[i + 1]) * kunit) for i in range(len(kb) - 1)]
        tpow = []
        for i in range(len(kb) - 1):
            mk = (kmag >= kb[i]) & (kmag < kb[i + 1])
            tpow.append(sum((mult * np.abs(Dt[q]) ** 2)[mk].sum() for q in range(3)))
        for k_ in names:
            Pc_hat = [sfft.rfftn(S0 * cl[k_][..., m], workers=a.workers) for m in range(6)]
            Dc = divP(Pc_hat)
            del Pc_hat
            row_e, row_r = [], []
            for i in range(len(kb) - 1):
                mk = (kmag >= kb[i]) & (kmag < kb[i + 1])
                e = sum((mult * np.abs(Dc[q] - Dt[q]) ** 2)[mk].sum() for q in range(3))
                cc = sum((mult * (Dc[q] * np.conj(Dt[q])).real)[mk].sum() for q in range(3))
                pc = sum((mult * np.abs(Dc[q]) ** 2)[mk].sum() for q in range(3))
                row_e.append(e / tpow[i] if tpow[i] > 0 else np.nan)
                row_r.append(cc / np.sqrt(pc * tpow[i]) if pc * tpow[i] > 0 else np.nan)
            M["force"][k_], M["r"][k_] = [float(v) for v in row_e], [float(v) for v in row_r]
            out("%-16s" % k_[:16] + " ".join("%9.3f" % v for v in row_e))
            out("%-16s" % ("  r") + " ".join("%9.3f" % v for v in row_r))
            del Dc
        del Dt, KX, KY, KZ, kmag, mult
        # block level
        b = a.block
        S0b = block_sum(S0, b)
        Pb_true = normalise(block_sum(P, b))
        xb = block_sum(x, b) / b ** 3
        # anisotropy against the distance to the nearest source (periodic),
        # ionised cells: the isotropisation of the true field away from sources
        from scipy.spatial import cKDTree
        tree = cKDTree(np.stack([ii, jj, kk], 1) + 0.5, boxsize=N)
        g1 = np.arange(N) + 0.5
        dist = np.empty((N, N, N))
        for kz in range(N):
            yy, xx = np.meshgrid(g1, g1, indexing="ij")
            pts = np.stack([xx.ravel(), yy.ravel(), np.full(xx.size, g1[kz])], 1)
            dist[kz] = tree.query(pts, workers=a.workers)[0].reshape(N, N)
        edges = [0, 1.5, 3, 6, 12, 24, 48, 96, 1e9]
        out("anisotropy A vs distance to the nearest source [cells], ionised cells (S0-weighted); err in brackets")
        out("%-12s %7s %7s %7s " % ("d", "frac", "A_true", "noise") + " ".join("%14s" % k[:14] for k in names))
        mi = reg["ionised x>0.9"]
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = mi & (dist >= lo) & (dist < hi)
            if not m.any():
                continue
            row = "%-12s %7.4f %7.3f %7.3f " % (f"{lo:g}-{hi:g}", m.sum() / mi.sum(), wmean(aniso(h_true), S0, m),
                                                wmean(noise, S0, m))
            row += " ".join("%6.3f (%5.3f)" % (wmean(aniso(cl[k]), S0, m), wmean(fro(cl[k], h_true), S0, m))
                            for k in names)
            out(row)
        out(f"block {b}^3: S0-weighted err of the block tensors (sum S0 h over the block)")
        rowh = "%-24s " % "regime" + " ".join("%10s" % k[:10] for k in names)
        out(rowh)
        blk = {k: normalise(block_sum(S0[..., None] * cl[k], b)) for k in names}
        for r, m in {"all": S0b > 0, "ionised blocks x>0.9": (S0b > 0) & (xb > 0.9),
                     "partly ionised": (S0b > 0) & (xb <= 0.9)}.items():
            if not m.any():
                continue
            M["block"][r] = {k: float(wmean(fro(blk[k], Pb_true), S0b, m)) for k in names}
            out("%-24s " % r + " ".join("%10.3f" % M["block"][r][k] for k in names))
        M["blockA"] = {k: float(wmean(aniso(blk[k]), S0b, S0b > 0)) for k in names}
        M["blockA"]["truth"] = float(wmean(aniso(Pb_true), S0b, S0b > 0))
        out("block anisotropy: truth %.3f | " % M["blockA"]["truth"] +
            " ".join("%s %.3f" % (k, M["blockA"][k]) for k in names))
        out()
        import json
        with open(os.path.join(a.out, f"metrics_z{z:.2f}.json"), "w") as fj:
            json.dump(M, fj, indent=1)

        if not a.no_plot:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            k0 = N // 2
            show = ["truth"] + [k for k in names if k in ("otvet", "split") or k.startswith("atten")
                                or k.startswith("cone") or k.startswith("split_cone")][:8]
            fig, axs = plt.subplots(2, len(show), figsize=(3.2 * len(show), 6.4), squeeze=False)
            for i, k in enumerate(show):
                h = h_true if k == "truth" else cl[k]
                axs[0, i].imshow(aniso(h)[k0], origin="lower", vmin=0, vmax=1, cmap="magma")
                axs[0, i].set_title(f"A, {k}")
                if k != "truth":
                    axs[1, i].imshow(fro(h, h_true)[k0], origin="lower", vmin=0, vmax=0.8, cmap="viridis")
                    axs[1, i].set_title(f"|h-h_true|, {k}")
                else:
                    axs[1, i].imshow(np.log10(S0[k0] + 1e-30), origin="lower", cmap="inferno",
                                     vmin=np.log10(np.median(S0[S0 > 0])) - 2,
                                     vmax=np.log10(np.median(S0[S0 > 0])) + 2)
                    axs[1, i].set_title("log S0 (rays)")
            for ax in axs.flat:
                ax.set_xticks([])
                ax.set_yticks([])
            fig.suptitle(f"z = {z:.2f}, x_V = {x.mean():.3f}")
            fig.tight_layout()
            fig.savefig(os.path.join(a.out, f"tensor_z{z:.2f}.png"), dpi=90)
            plt.close(fig)


if __name__ == "__main__":
    main()
