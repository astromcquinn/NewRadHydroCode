#!/usr/bin/env python3
"""Write a CAMB transfer-function table in the 13-column format that MUSIC's
`transfer = camb_file` reader expects (src/plugins/transfer_camb.cc: k/h, CDM,
baryon, photon, massless nu, massive nu, total, no-nu, total-de, Weyl,
v_CDM, v_b, v_b - v_CDM), evaluated at the starting redshift of the ICs
(MUSIC asks for the table at zstart when baryons = yes).

Defaults: the cosmology of Cain et al. 2026 (Planck-like: Omega_m = 0.305,
Omega_b = 0.048, h = 0.68, n_s = 0.9667, sigma_8 = 0.82), z = 300.
The transfer functions do not depend on A_s; MUSIC renormalises to the
sigma_8 of its .conf, so pass the same --sigma8 here (it is printed for
the record). One CAMB call.

usage: camb_transfer.py OUT.dat [--zstart 300] [--omega-m 0.305] ...
needs: pip install camb   (or conda install -c conda-forge camb)
"""
import argparse
import numpy as np
import camb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--zstart", type=float, default=300.0)
    ap.add_argument("--omega-m", type=float, default=0.305)
    ap.add_argument("--omega-b", type=float, default=0.048)
    ap.add_argument("--h", type=float, default=0.68)
    ap.add_argument("--ns", type=float, default=0.9667)
    ap.add_argument("--sigma8", type=float, default=0.82)
    ap.add_argument("--kmax", type=float, default=5000.0, help="h/Mpc; the 0.3 kpc/h cells need k ~ 1e4")
    ap.add_argument("--k-per-logint", type=int, default=30, help="CAMB's k points per decade")
    a = ap.parse_args()

    # one CAMB call: the transfer functions do not depend on A_s (MUSIC
    # renormalises to its own sigma_8), so A_s only sets the sigma_8 printed
    p = camb.CAMBparams()
    p.set_cosmology(H0=100 * a.h, ombh2=a.omega_b * a.h ** 2, omch2=(a.omega_m - a.omega_b) * a.h ** 2,
                    mnu=0.0, omk=0.0, tau=0.054)
    p.InitPower.set_params(ns=a.ns, As=2.1e-9)
    p.set_matter_power(redshifts=[a.zstart, 0.0], kmax=a.kmax, k_per_logint=a.k_per_logint)
    p.NonLinear = camb.model.NonLinear_none
    p.Transfer.accurate_massive_neutrinos = False
    r = camb.get_results(p)
    s8 = r.get_sigma8_0()
    print(f"sigma_8(z=0) = {s8:.5f} at A_s = 2.1e-9; A_s for sigma_8 = {a.sigma8}: {2.1e-9 * (a.sigma8 / s8) ** 2:.4e} "
          "(MUSIC renormalises to the sigma_8 of its .conf; this is for the record)")
    td = r.get_matter_transfer_data()
    # transfer_data: (13 columns, nk, nz); redshifts are stored in the order given, highest first
    iz = list(r.transfer_redshifts).index(a.zstart) if hasattr(r, "transfer_redshifts") else 0
    tab = td.transfer_data[:, :, iz].T  # (nk, 13)
    k = tab[:, 0]
    assert np.all(np.diff(k) > 0), "k not increasing"
    hdr = (f"CAMB {camb.__version__} transfer functions at z = {a.zstart} for MUSIC (transfer = camb_file); "
           f"Omega_m = {a.omega_m}, Omega_b = {a.omega_b}, h = {a.h}, n_s = {a.ns}, sigma_8 = {a.sigma8}\n"
           "k/h CDM baryon photon nu mass_nu total no_nu total_de Weyl v_CDM v_b v_b-v_CDM")
    np.savetxt(a.out, tab, fmt="%.8e", header=hdr)
    print(f"wrote {a.out}: {len(k)} rows, k = {k[0]:.3e} .. {k[-1]:.3e} h/Mpc")


if __name__ == "__main__":
    main()
