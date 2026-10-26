// The X-ray background of xray_heating = 1 (xray_background.hpp,
// docs/xray_heating.md): the SFR density against the slow mass function,
// the emitting-halo threshold against the halo finder's virial
// temperature, the transfer in the optically thin limit against a
// brute-force quadrature and in the opaque limit against energy
// conservation, the secondary-electron fractions, and the cell physics
// (photo_gas_cell) with the X-ray terms: bit-identical without them, the
// heating and ionisation they add in a neutral cell.
#include <cmath>
#include <cstdio>
#include <vector>

#include "check.hpp"
#include "radray/config.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/halo_sources.hpp"
#include "radray/halos.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/rt_rates.hpp"
#include "radray/xray_background.hpp"

using namespace radray;

namespace {

// data/uv_lf/am_params_st_fof.txt (the production source model)
SourceModelParams production_sources() {
  HaloSourceParams hp;
  const char* kv[][2] = {{"sfr_model", "abundance"},     {"am_l0", "29.096732"},
                         {"am_l1", "-0.1246293432"},     {"am_l2", "-0.05156652063"},
                         {"am_a0", "1.486136421"},       {"am_a1", "0.0255368079"},
                         {"am_m0", "11.38823653"},       {"am_m1", "-0.1404545307"},
                         {"am_b", "0.2810735695"},       {"am_z0", "8"},
                         {"am_zlo", "5.9"},              {"am_zhi", "9"},
                         {"am_norm_z", "10,11,12.5"},
                         {"am_norm_dl", "0.3205722651,0.5840562069,0.7969103106"}};
  for (auto& p : kv) CHECK(set_halo_param(hp, p[0], p[1]));
  return hp.model;
}

double hubble(double z) {
  using namespace cp;
  const double a = 1 / (1 + z);
  return pc::H0_cgs * hubble0 * std::sqrt(omegam / (a * a * a) + omegal + omegar / (a * a * a * a));
}

}  // namespace

int main() {
  using namespace pc;
  const SourceModelParams src = production_sources();
  XrayParams xp;
  xp.on = true;
  xp.dlog = 0.004;
  const XrayBackground bg(xp, src);

  // --- parameters ---
  {
    XrayParams p;
    CHECK(set_xray_param(p, "xray_st", "0.3,0.7,0.25"));
    CHECK(p.st.A == 0.3 && p.st.a == 0.7 && p.st.p == 0.25);
    CHECK(set_xray_param(p, "xray_opacity", "neutral"));
    CHECK(p.opacity == XrayParams::Opacity::Neutral);
    CHECK(!set_xray_param(p, "am_l0", "1"));
    bool threw = false;
    try {
      set_xray_param(p, "xray_heating", "2");
    } catch (const std::exception&) {
      threw = true;
    }
    CHECK(threw);
    SourceModelParams g;  // sfr_model = global: refused
    threw = false;
    try {
      XrayBackground b(xp, g);
    } catch (const std::exception&) {
      threw = true;
    }
    CHECK(threw);
  }

  // --- the emitting-halo threshold: T_vir(M_min) = tvir_min with the halo
  // finder's virial radius and temperature ---
  for (double z : {6.0, 10.0, 20.0}) {
    const double M = bg.m_min_msun(z);
    HaloUnits hu;
    hu.redshift = z;
    hu.hubble0 = cp::hubble0;
    const double H = hubble(z);
    hu.rho_crit = 3 * H * H / (8 * pi * G_cgs);
    const double a = 1 / (1 + z);
    const double om_z = cp::omegam / (a * a * a) / (H * H / (H0_cgs * H0_cgs * cp::hubble0 * cp::hubble0));
    const double R = virial_radius_mpch(M, delta_vir_bryan_norman(om_z), hu);
    CHECK_CLOSE(virial_temperature(M, R, 0.59, hu), 1e4, 1e-9);
  }

  // --- rho_SFR: the table against the direct integral, and the direct
  // integral against dndlnM of the mass function (independent sigma) ---
  for (double z : {6.0, 8.3, 12.7}) CHECK_CLOSE(bg.sfrd(z), bg.sfrd_direct(z), 2e-3);
  {
    const double z = 8.0, h = cp::hubble0;
    MassFunction mf;
    mf.st = xp.st;
    const double l0 = std::log(bg.m_min_msun(z)), l1 = std::log(1e15);
    const int n = 400;
    double s = 0;
    for (int i = 0; i < n; ++i) {
      const double M = std::exp(l0 + (i + 0.5) * (l1 - l0) / n);  // Msun
      s += mf.dndlnM(M * h, z, MassFunction::Kind::ShethTormen) * h * h * h *
           src.am.luv(M, z) / src.L_per_sfr;
    }
    s *= (l1 - l0) / n;
    CHECK_CLOSE(bg.sfrd_direct(z), s, 5e-3);
    std::printf("rho_SFR(z = 8) = %.4e Msun/yr/Mpc^3 (direct %.4e)\n", bg.sfrd_direct(z), s);
  }

  // --- optically thin: J_E against a brute-force quadrature in z ---
  {
    const double z0 = 8.0, x0[3] = {0, 0, 0};
    std::vector<double> E, J;
    bg.spectrum(z0, x0, E, J);
    CHECK(!E.empty());
    for (double Et : {0.3, 1.0, 4.0}) {
      std::size_t i = 0;
      while (i + 1 < E.size() && E[i] < Et) ++i;
      const double Ei = E[i];
      const int n = 400000;
      const double s0 = std::log(1 + z0), s1 = std::log(1 + xp.z_max), ds = (s1 - s0) / n;
      double sum = 0;
      for (int k = 0; k < n; ++k) {
        const double z = std::exp(s0 + (k + 0.5) * ds) - 1;
        sum += bg.emissivity(z) * bg.shape(Ei * (1 + z) / (1 + z0)) / hubble(z);
      }
      const double Jbf = c_cgs / (4 * pi) * std::pow(1 + z0, 3) * sum * ds /
                         std::pow(Mpc2cm, 3);
      std::printf("thin J(%.3f keV) = %.5e, brute force %.5e\n", Ei, J[i], Jbf);
      CHECK_CLOSE(J[i], Jbf, 1e-2);
    }
  }

  // --- opaque: the absorbed power per proper volume equals the proper
  // emissivity (every photon absorbed where it is emitted) ---
  {
    const double z0 = 9.0, xf[3] = {1e6, 1e6, 0};
    XrayParams p = xp;
    p.opacity = XrayParams::Opacity::Box;
    const XrayBackground b(p, src);
    const XrayRates r = b.rates(z0, xf);
    const double nH = cp::omegab * rhoc_cgs * cp::hubble0 * cp::hubble0 * cp::XHy / mH_cgs *
                      std::pow(1 + z0, 3) * xf[0];
    const double nHe = cp::omegab * rhoc_cgs * cp::hubble0 * cp::hubble0 * cp::YHe / mHe_cgs *
                       std::pow(1 + z0, 3) * xf[1];
    const double absorbed = nH * (r.phi[0] + HI_eV * eV_cgs * r.gamma[0]) +
                            nHe * (r.phi[1] + HeI_eV * eV_cgs * r.gamma[1]);
    // emitted: eps over [e0, emax] in units of the [e0, eband] band
    const int n = 100000;
    double band = 0;
    const double l0 = std::log(p.e0_keV), l1 = std::log(p.emax_keV), dl = (l1 - l0) / n;
    for (int i = 0; i < n; ++i) {
      const double E = std::exp(l0 + (i + 0.5) * dl);
      band += b.shape(E) * E * dl;
    }
    const double emitted = b.emissivity(z0) * band * std::pow(1 + z0, 3) / std::pow(Mpc2cm, 3);
    std::printf("opaque: absorbed %.5e, emitted %.5e erg/s/cm^3\n", absorbed, emitted);
    CHECK_CLOSE(absorbed, emitted, 2e-2);
  }

  // --- neutral IGM: rates positive, heating per H atom of order 1e-29 erg/s
  // at z = 8 for the default L_X/SFR ---
  {
    const double xf[3] = {1, 1, 0};
    const XrayRates r = bg.rates(8.0, xf);
    CHECK(r.gamma[0] > 0 && r.gamma[1] > 0 && r.phi[0] > 0 && r.phi[1] > 0);
    CHECK(r.phi[0] > 1e-31 && r.phi[0] < 1e-27);
  }

  // --- secondary electrons ---
  {
    double fh, fi, fe;
    rrk::xray_secondary(1.0, fh, fi, fe);
    CHECK_CLOSE(fh, 0.9971, 1e-12);
    CHECK(fi == 0 && fe == 0);
    rrk::xray_secondary(0.0, fh, fi, fe);  // clamped to 1e-4
    double fh4, fi4, fe4;
    rrk::xray_secondary(1e-4, fh4, fi4, fe4);
    CHECK(fh == fh4 && fi == fi4 && fe == fe4);
    CHECK(fh > 0.1 && fh < 0.12 && fi > 0.35 && fi < 0.40 && fe > 0.04 && fe < 0.06);
    double prev = 0;
    for (double x = 1e-4; x <= 1; x *= 1.5) {  // heat fraction rises with x
      rrk::xray_secondary(x, fh, fi, fe);
      CHECK(fh >= prev);
      CHECK(fh + fi + fe <= 1.0);
      prev = fh;
    }
  }

  // --- the cell physics ---
  {
    const Config cfg = Config::make(16, 4, 8, 4.0);
    SimClock clk;
    clk.redshift = 10;
    clk.scalefactor = 1 / 11.0;
    Units u;
    calc_sim_var(cfg, clk, u);
    AtomicRates ar;
    calc_atomic_rates(ar, false);
    rrk::PhysCtx X;
    X.u = u;
    X.rates.cool = ar.coef_cool_.data();
    X.rates.alpha = ar.alpha_rec_.data();
    X.rates.gcol = ar.gamma_col_.data();
    X.rates.legacy_physics = false;
    X.dtsec_rt = 3e13;  // 1 Myr
    double nu[5], sig[3][5];
    for (int q = 0; q < 5; ++q) nu[q] = cp::freq_init[q];
    rrk::photoion_sigma(nu, sig);
    const float gng[5] = {0, 0, 0, 0, 0};
    const double T0 = 20;
    auto cell = [&](float uu[5], float ww[4]) {
      const double D = u.Dgas_avg;
      const double nH = D * u.rho_unit * cp::XHy / mH_cgs, nHe = D * u.rho_unit * cp::YHe / mHe_cgs;
      uu[0] = static_cast<float>(D);
      uu[1] = uu[2] = uu[3] = 0;
      uu[4] = static_cast<float>((nH + nHe) * k_cgs * T0 / (cp::gamma_ideal - 1) / u.eden_unit);
      ww[0] = static_cast<float>(D * cp::XHy);
      ww[1] = static_cast<float>(D * cp::YHe);
      ww[2] = 0;
      ww[3] = 0;
    };
    // (a) xray on with zero rates: bit-identical to off
    float u1[5], w1[4], u2[5], w2[4];
    cell(u1, w1);
    cell(u2, w2);
    double nr1[5] = {0, 0, 0, 0, 0}, nr2[5] = {0, 0, 0, 0, 0};
    rrk::PhotoStats s1, s2;
    CHECK(rrk::photo_gas_cell(u1, w1, nullptr, gng, nu, sig, nr1, X, s1) == rrk::kPhotoOK);
    rrk::PhysCtx Xz = X;
    Xz.xray = true;
    CHECK(rrk::photo_gas_cell(u2, w2, nullptr, gng, nu, sig, nr2, Xz, s2) == rrk::kPhotoOK);
    for (int m = 0; m < 5; ++m) CHECK(u1[m] == u2[m]);
    for (int m = 0; m < 4; ++m) CHECK(w1[m] == w2[m]);
    // a neutral cell at T0 > temp_min is left alone by the original
    // physics (no electrons: no cooling)
    float u0[5], w0[4];
    cell(u0, w0);
    CHECK(u1[4] == u0[4] && w1[0] == w0[0]);

    // (b) heating and ionisation: dT = (gamma-1) f_heat edep dt/(k n_tot),
    // dn_HI = -(gamma_HI + f_ion,HI edep/(I_HI n_HI)) n_HI dt, to first order
    rrk::PhysCtx Xr = X;
    Xr.xray = true;
    Xr.xray_gamma[0] = 1e-20;
    Xr.xray_gamma[1] = 1e-20;
    Xr.xray_phi[0] = 2e-29;
    Xr.xray_phi[1] = 2e-29;
    Xr.dtsec_rt = 3e14;  // 10 Myr
    float u3[5], w3[4];
    cell(u3, w3);
    double nr3[5] = {0, 0, 0, 0, 0};
    rrk::PhotoStats s3;
    CHECK(rrk::photo_gas_cell(u3, w3, nullptr, gng, nu, sig, nr3, Xr, s3) == rrk::kPhotoOK);
    const double D = u.Dgas_avg;
    const double nH = D * u.rho_unit * cp::XHy / mH_cgs, nHe = D * u.rho_unit * cp::YHe / mHe_cgs;
    double fh, fi, fe;
    rrk::xray_secondary(0.0, fh, fi, fe);
    const double edep = Xr.xray_phi[0] * nH + Xr.xray_phi[1] * nHe;
    const double dT = (cp::gamma_ideal - 1) * fh * edep * Xr.dtsec_rt / (k_cgs * (nH + nHe));
    const double ntot1 = nH + nHe + (nH - static_cast<double>(w3[0]) * u.H_unit) +
                         (nHe - static_cast<double>(w3[1]) * u.He_unit);
    const double T1 = (cp::gamma_ideal - 1) * static_cast<double>(u3[4]) * u.eden_unit / (k_cgs * ntot1);
    const double dxHI = (Xr.xray_gamma[0] + fi * edep / (HI_eV * eV_cgs * nH)) * Xr.dtsec_rt;
    const double xHI = static_cast<double>(w3[0]) / (D * cp::XHy);
    std::printf("cell: T %.3f -> %.3f K (expected +%.3f), 1 - x_HI = %.4e (expected %.4e)\n", T0,
                T1 * ntot1 / (nH + nHe), dT, 1 - xHI, dxHI);
    // the thermal energy grows by the heat (the per-particle T falls as the
    // particle number grows); Compton heating by the few electrons is ~1e-3
    CHECK_CLOSE(T1 * ntot1 / (nH + nHe) - T0, dT, 2e-2);
    CHECK_CLOSE(1 - xHI, dxHI, 2e-2);
  }

  TEST_MAIN_END();
}
