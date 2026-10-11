#!/usr/bin/env python3
"""Escaped ionizing emissivity of the simulations vs redshift.

    emissivity.py --halos CATALOG_DIR [--run DIR:LABEL ...] [--moments DIR:LABEL ...]
                  [--box 10] [--nmin 20] [--beta 2.8] [--f6 F] [--out fig.png]

--run (sfr_model = abundance): applies the run's L_UV(M_halo, z) (the
am_* lines of DIR/params.txt, as box_lf.py and the C++ LuvModel), xi_ion
and f_esc to the FoF halos (N >= nmin) of every catalog in CATALOG_DIR,
  ndot_esc = f_esc 10^xi_ion sum_i L_UV,i / V   [photons s^-1 cMpc^-3].
That is the rate the code hands to the RT (halo_sources.cpp, assign_sources:
L_rt = L_UV 10^xi_ion / xi_pop2 * f_esc / fResc). The catalogs need not come
from the run itself: hmf256 (DM + gas, the reion256 ICs) has 24 redshifts,
the production runs find their halos on the fly and do not write them.

--moments (rt_method = m1 | otvet, any source model): the photons the run
injected, from the cumulative `Moments : ... injected` line of
DIR/run/Log/ray.log (photons per H atom; includes recombination photons when
recomb_rad = 1), differenced between hydro steps.

Also drawn: the median of the Lopez et al. (in prep.) joint forest + CMB
emissivity posterior, read by eye from Fig. 1 of D'Aloisio et al.
("Models of the declining ionizing emissivity near the end of
reionization"); and, for the first --run, f_esc(z) = f_6 ((1+z)/7)^beta
(their toy model with a halo-mass independent f_esc: beta = 2.6-3.0), with
f_6 from --f6 or, by default, chosen so that the curve meets the Lopez
median at z = 6.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from box_lf import log_luv, read_params  # noqa: E402
from hmf_compare import catalogs, read_catalog  # noqa: E402

# include/radray/constants.hpp (real(4) literals)
OMEGAM, OMEGAB, H0, YHE = 0.27, 0.045, 0.70, 0.24
MPC_CM, MH_G, SEC_PER_YR = 3.0857e24, 1.6726e-24, 3.156e7
RHO_CRIT100 = 1.8785e-29  # g cm^-3 for h = 1
trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz  # numpy < 2.0

# Lopez et al. in prep. median, read by eye from D'Aloisio et al. Fig. 1
LOPEZ_Z = [4, 5, 6, 7, 8, 9, 10, 11, 12]
LOPEZ_LOG = [50.40, 50.64, 50.71, 50.73, 50.70, 50.62, 50.51, 50.40, 50.24]


def read_obs(path):
    """Blocks of data/igm/igm_z5_6.txt: {name: [(z, v, +err, -err, ref)]}."""
    out, cur = {}, None
    for l in open(path):
        t = l.split()
        if l.startswith("# block"):
            cur = t[2]
            out[cur] = []
        elif t and not l.startswith("#") and cur:
            out[cur].append((float(t[0]), float(t[1]), float(t[2]), float(t[3]), t[4]))
    return out


def nH_comoving():
    """Mean comoving hydrogen number density [cm^-3]."""
    return (1 - YHE) * OMEGAB * RHO_CRIT100 * H0 ** 2 / MH_G


def dtdz_s(z):
    Hz = 100 * H0 * 1e5 / MPC_CM * np.sqrt(OMEGAM * (1 + z) ** 3 + 1 - OMEGAM)
    return 1 / ((1 + z) * Hz)


def abundance_run(path, cats, box, nmin):
    p = read_params(os.path.join(path, "params.txt"))
    if p.get("sfr_model") != "abundance":
        raise SystemExit(f"{path}: sfr_model is not abundance; use --moments")
    xi, f0 = float(p["xi_ion"]), float(p["f_esc"])
    beta = float(p.get("f_esc_beta", 0))
    z0, fmax = float(p.get("f_esc_z0", 6)), float(p.get("f_esc_max", 1))
    # halo_sources.cpp f_esc_at
    fesc = np.array([f0 if beta == 0 else min(fmax, f0 * ((1 + z) / (1 + z0)) ** beta)
                     for z, _ in cats])
    desc = f"f_esc = {f0:g}" if beta == 0 else \
        f"f_esc = {f0:g} ((1+z)/{1 + z0:g})$^{{{beta:g}}}$"
    vol = (box / H0) ** 3
    zs, nint = [], []
    for z, cat in cats:
        h = read_catalog(cat)
        h = h[h["N"] >= nmin]
        L = np.sum(10.0 ** log_luv(p, np.log10(h["M"].astype(float)), z)) if len(h) else 0.0
        zs.append(z)
        nint.append(10 ** xi * L / vol)  # f_esc = 1
    return np.array(zs), np.array(nint), fesc, desc


def moment_run(path):
    zs, inj = [], []
    zcur = None
    for l in open(os.path.join(path, "run", "Log", "ray.log")):
        if l.startswith(" Hydro step"):
            zcur = float(l.split()[-1])
        elif l.startswith(" Moments") and zcur is not None:
            t = l.split()
            inj.append(float(t[t.index("injected") + 1]))
            zs.append(zcur)
    zs, inj = np.array(zs), np.array(inj)
    # the Moments line closes the step that started at zs[k]; it ends at zs[k+1]
    za, zb = zs[:-1], zs[1:]
    dN = np.diff(np.concatenate([[0.0], inj]))[:-1]
    dt = np.array([trapz(dtdz_s(np.linspace(b, a, 50)), np.linspace(b, a, 50))
                   for a, b in zip(za, zb)])
    ndot = dN * nH_comoving() * MPC_CM ** 3 / dt
    return np.sqrt((1 + za) * (1 + zb)) - 1, ndot


def history(path):
    f = os.path.join(path, "run", "reion.txt")
    if not os.path.exists(f):
        return None
    a = np.loadtxt(f, comments="#")
    return a[:, 1], a[:, 2]  # z, x_HII,V


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--halos", required=True)
    ap.add_argument("--run", action="append", default=[])
    ap.add_argument("--moments", action="append", default=[])
    ap.add_argument("--box", type=float, default=10.0, help="Mpc/h")
    ap.add_argument("--nmin", type=int, default=20)
    ap.add_argument("--beta", type=float, default=2.8)
    ap.add_argument("--f6", type=float, default=None)
    ap.add_argument("--obs", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "..", "..", "data", "igm", "igm_z5_6.txt"))
    ap.add_argument("--out", default="emissivity.png")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cats = catalogs(a.halos)
    fig, (ax, axf, axx) = plt.subplots(3, 1, figsize=(6.4, 9.6), sharex=True,
                                       gridspec_kw={"height_ratios": [2.2, 1.2, 1]},
                                  constrained_layout=True)
    ax.plot(LOPEZ_Z, 10 ** np.array(LOPEZ_LOG), "-", color="0.6", lw=1.0,
            label="Lopez+ in prep. model median (by eye, D'Aloisio+ Fig. 1)")
    if a.obs and os.path.exists(a.obs):
        obs = read_obs(a.obs)
        for blk, mk, lab in (("ndot", "s", "Becker+21 (Γ, λ$_{mfp}$)"),
                             ("ndot_gz", "D", "Gaikwad+23 Γ + Zhu+23 λ$_{mfp}$")):
            r = np.array([o[:4] for o in obs.get(blk, [])])
            if len(r):
                v = 10 ** r[:, 1]
                ax.errorbar(r[:, 0], v, yerr=[v - 10 ** (r[:, 1] - r[:, 3]), 10 ** (r[:, 1] + r[:, 2]) - v],
                            fmt=mk, color="k", ms=5, capsize=2, zorder=5, label=lab)
    rows = []
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    first = None
    for i, spec in enumerate(a.run):
        path, label = spec.rsplit(":", 1)
        z, nint, fesc, desc = abundance_run(path, cats, a.box, a.nmin)
        c = colors[i % len(colors)]
        ax.plot(z, fesc * nint, "-", color=c, lw=2, label=f"{label} ({desc})")
        rows.append((label, z, fesc * nint))
        if first is None:
            first = (label, z, nint, fesc, desc, c)
        hs = history(path)
        if hs is not None:
            axx.plot(hs[0], hs[1], "-", color=c, lw=1.5, label=label)
    if first is not None:
        label, z, nint, fesc, desc, c = first
        ax.plot(z, nint, ":", color="0.5", lw=1.2, label="intrinsic (f_esc = 1)")
        # f_esc(z) that turns the intrinsic emissivity into the Lopez median,
        # and a power-law fit to it over the redshifts both cover
        ok = (z >= LOPEZ_Z[0]) & (z <= LOPEZ_Z[-1]) & (nint > 0)
        freq = 10 ** np.interp(z[ok], LOPEZ_Z, LOPEZ_LOG) / nint[ok]
        lx = np.log10((1 + z[ok]) / 7.0)
        bfit, lf6 = np.polyfit(lx, np.log10(freq), 1)
        axf.plot(z[ok], freq, "o", color="0.6", ms=3, label="Lopez model median / intrinsic")
        axf.plot(z, fesc, "-", color=c, lw=2, label=f"{label}: {desc}")
        zz = np.linspace(5, 13.5, 100)
        curves = [(10 ** lf6, bfit, "-", "fit"), (None, a.beta, "--", "D'Aloisio+ toy")]
        for f6, beta, ls, what in curves:
            if f6 is None:
                f6 = a.f6 if a.f6 is not None else \
                    10 ** np.interp(6.0, LOPEZ_Z, LOPEZ_LOG) / 10 ** np.interp(6.0, z, np.log10(nint))
            lab = f"{f6:.3g} ((1+z)/7)$^{{{beta:.2g}}}$ ({what})"
            axf.plot(zz, np.minimum(1, f6 * ((1 + zz) / 7) ** beta), ls, color="0.35", lw=1.3,
                     label=lab)
            fz = np.minimum(1.0, f6 * ((1 + z) / 7.0) ** beta)
            ax.plot(z, fz * nint, ls, color=c, lw=1.3, alpha=0.8 if what == "fit" else 0.5,
                    label=f"f_esc = {lab}")
            print(f"# f_esc(z) = {f6:.4g} ((1+z)/7)^{beta:.3g} ({what}): f_esc(z=6, 8, 10, 12) = "
                  + ", ".join(f"{min(1, f6 * ((1 + q) / 7) ** beta):.3f}" for q in (6, 8, 10, 12)))
            rows.append((f"{label}, f_esc(z) {what}", z, fz * nint))
        if a.obs and os.path.exists(a.obs):
            for blk, mk, lab in (("ndot", "s", "Becker+21 / intrinsic"),
                                 ("ndot_gz", "D", "Gaikwad+23 Γ + Zhu+23 λ / intrinsic")):
                r = np.array([o[:4] for o in read_obs(a.obs).get(blk, [])])
                if len(r):
                    li = np.interp(r[:, 0], z[::-1], np.log10(nint[::-1])) if z[0] > z[-1] else \
                        np.interp(r[:, 0], z, np.log10(nint))
                    fv = 10 ** (r[:, 1] - li)
                    axf.errorbar(r[:, 0], fv, yerr=[fv - 10 ** (r[:, 1] - r[:, 3] - li),
                                                    10 ** (r[:, 1] + r[:, 2] - li) - fv],
                                 fmt=mk, color="k", mfc="w", ms=5, capsize=2, label=lab)
                    for q, f_ in zip(r[:, 0], fv):
                        print(f"# f_esc needed for {blk} at z = {q:.2f}: {f_:.3f}")
        axf.set_yscale("log")
        axf.set_ylim(0.03, 1.1)
        axf.set_ylabel(r"$f_{\rm esc}$")
        axf.legend(frameon=False, fontsize=7.5, loc="upper left")
        axf.grid(alpha=0.25, lw=0.5, which="both")
    for j, spec in enumerate(a.moments):
        path, label = spec.rsplit(":", 1)
        z, nd = moment_run(path)
        c = colors[(len(a.run) + j) % len(colors)]
        ax.plot(z, nd, "-", color=c, lw=1.5, alpha=0.8, label=f"{label} (injected)")
        rows.append((label, z, nd))
        hs = history(path)
        if hs is not None:
            axx.plot(hs[0], hs[1], "-", color=c, lw=1.2, label=label)
    ax.set_yscale("log")
    ax.set_ylim(1e49, 3e52)
    ax.set_ylabel(r"$\dot n_{\rm ion}$ [s$^{-1}$ cMpc$^{-3}$]")
    ax.legend(frameon=False, fontsize=7.5, loc="lower center", ncol=1)
    ax.grid(alpha=0.25, lw=0.5)
    axx.set_xlim(4.5, 14)
    axx.invert_xaxis()
    axx.set_ylim(0, 1.02)
    axx.set_xlabel("redshift")
    axx.set_ylabel(r"$x_{\rm HII,V}$")
    axx.grid(alpha=0.25, lw=0.5)
    fig.savefig(a.out, dpi=140)
    print("#", a.out)
    print("# escaped emissivity [log10 s^-1 cMpc^-3]")
    for label, z, nd in rows:
        print(f"## {label}")
        for zz, v in zip(z, nd):
            if v > 0:
                print(f"  z = {zz:6.3f}  {np.log10(v):6.3f}")


if __name__ == "__main__":
    main()
