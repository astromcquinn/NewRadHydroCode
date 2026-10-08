// Framework-independent per-cell RT kernels for the spike: the gas -> RT
// grid average (prep_rt/calc_rt_i), the per-substep reset and optical depth
// (calc_rt_1, legacy_rt = 0), radiative_physics/photo_cell, and the final
// optical depth (post_rt/calc_rt_f). Line-by-line copies of src/rt_tools.cpp
// and src/rt_rates.cpp for N_rt_gas = 1, operating on one cell's values, so
// that every framework calls the same arithmetic. The RT statistics
// (rt_stat: nion, ncol, ...) are diagnostics of ray.log and are not
// accumulated here.
#pragma once

#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "consts_hd.hpp"
#include "rr_hd.hpp"

namespace rrk {

// One RT cell (rt_grid_type fields) plus its gas cell u_gas(5), w_gas(4).
struct CellData {
  float rhog, fHI, fHeI, fHeII, temp, Luv, sfr;
  float dx[3], gamma[3], lambda[3];
  float nu[5], ngam[5], fgam[5], tau[5];
  float u[5], w[4];
};

// rt_rates.f90: photoion_sigma, sigma[species][freq]
RR_HD void photoion_sigma(const double nu[5], double sig[3][5]) {
  using namespace radray;
  using namespace radray::pc;
  auto s = [](float v) { return static_cast<double>(v); };
  for (int i = 0; i < 5; ++i) {
    const double fq = nu[i];
    double x, y;
    if (fq < HI_eV) {
      sig[0][i] = 0;
    } else {
      x = fq / s(0.4298f);
      y = x;
      sig[0][i] = 5.475e-14 * rrf::pow2(x - 1) / std::pow(y, s(4.0185f)) /
                  std::pow(1 + std::sqrt(y / s(32.88f)), s(2.963f));
    }
    if (fq < HeI_eV) {
      sig[1][i] = 0;
    } else {
      x = fq / s(13.61f) - s(0.4434f);
      y = std::sqrt(rrf::pow2(x) + s(2.136f * 2.136f));
      sig[1][i] = 9.492e-16 * (rrf::pow2(x - 1) + s(2.039f * 2.039f)) /
                  std::pow(y, s(3.9060f)) / std::pow(1 + std::sqrt(y / s(1.469f)), s(3.188f));
    }
    if (fq < HeII_eV) {
      sig[2][i] = 0;
    } else {
      x = fq / s(1.720f);
      y = x;
      sig[2][i] = 1.369e-14 * rrf::pow2(x - 1) / std::pow(y, s(4.0185f)) /
                  std::pow(1 + std::sqrt(y / s(32.88f)), s(2.963f));
    }
  }
}

// Atomic-rate tables (AtomicRates::coef_cool_ etc., Fortran layout)
struct RatesView {
  const double* cool = nullptr;   // (16, Ntemp)
  const double* alpha = nullptr;  // (6, Ntemp)
  const double* gcol = nullptr;   // (3, Ntemp)
  bool legacy_physics = true;
};
struct RateSetP {
  double alpha[6], gcol[3], cool[16];
};
RR_HD RateSetP rates_node(const RatesView& R, int iT) {
  RateSetP r;
  for (int n = 0; n < 6; ++n) r.alpha[n] = R.alpha[n + 6 * (iT - 1)];
  for (int n = 0; n < 3; ++n) r.gcol[n] = R.gcol[n + 3 * (iT - 1)];
  for (int n = 0; n < 16; ++n) r.cool[n] = R.cool[n + 16 * (iT - 1)];
  return r;
}
RR_HD double rate_lerp(double a, double b, double w) {
  if (a > 0 && b > 0) return a * std::exp(w * std::log(b / a));
  return a + w * (b - a);
}
RR_HD RateSetP rates_interp(const RatesView& R, double tmp) {
  using namespace radray;
  const double u = (std::log10(tmp) - cp::lgT_min) / cp::lgT_del;
  if (u > 0 && u < cp::Ntemp - 1) {
    const int i = 1 + static_cast<int>(u);
    const double w = u - (i - 1);
    RateSetP r;
    for (int n = 0; n < 6; ++n) r.alpha[n] = rate_lerp(R.alpha[n + 6 * (i - 1)], R.alpha[n + 6 * i], w);
    for (int n = 0; n < 3; ++n) r.gcol[n] = rate_lerp(R.gcol[n + 3 * (i - 1)], R.gcol[n + 3 * i], w);
    for (int n = 0; n < 16; ++n) r.cool[n] = rate_lerp(R.cool[n + 16 * (i - 1)], R.cool[n + 16 * i], w);
    return r;
  }
  return rates_node(R, u > 0 ? cp::Ntemp : 1);
}

struct GasStateP {
  double rhog, nH, nHe, nHI, nHeI, nHeII, nHII, nHeIII, ne, neff, tmp;
};
RR_HD GasStateP gas_state_from_avg(const double us_in[5], const double ws[4],
                                   const radray::Units& u) {
  using namespace radray;
  using namespace radray::pc;
  using namespace radray::cp;
  GasStateP g;
  const double us1 = us_in[0];
  g.rhog = us1 * u.rho_unit;
  g.nH = g.rhog * XHy / mH_cgs;
  g.nHe = g.rhog * YHe / mHe_cgs;
  g.nHI = rrf::min(ws[0] * u.H_unit, g.nH);
  g.nHeI = rrf::min(ws[1] * u.He_unit, g.nHe);
  g.nHeII = rrf::min(ws[2] * u.He_unit, g.nHe - g.nHeI);
  g.nHII = g.nH - g.nHI;
  g.nHeIII = g.nHe - g.nHeI - g.nHeII;
  g.ne = g.nHII + g.nHeII + 2 * g.nHeIII;
  g.neff = g.nH + g.nHe + g.ne;
  const double ksum = ((0.0 + us_in[1] * us_in[1]) + us_in[2] * us_in[2]) + us_in[3] * us_in[3];
  const double eth = (us_in[4] - ksum / us1 / 2) * u.eden_unit;
  g.tmp = rrf::max((gamma_ideal - 1) * eth / k_cgs / g.neff, u.temp_min);
  return g;
}
RR_HD void store_gas_state(CellData& c, const GasStateP& g, const radray::Units& u) {
  c.rhog = static_cast<float>(g.rhog / u.rhob_avg);
  c.fHI = static_cast<float>(rrf::min(rrf::max(g.nHI / g.nH, 0.0), 1.0));
  c.fHeI = static_cast<float>(rrf::min(rrf::max(g.nHeI / g.nHe, 0.0), 1.0));
  c.fHeII = static_cast<float>(rrf::min(rrf::max(g.nHeII / g.nHe, 0.0), 1.0));
  c.temp = static_cast<float>(g.tmp);
}

// average_gas for N_rt_gas = 1 (normalised by N_rt_gas_c = 1, clamped)
RR_HD void average_gas_1(const CellData& c, const radray::Units& u, double us[5], double ws[4]) {
  for (int m = 0; m < 5; ++m) us[m] = (0.0 + static_cast<double>(c.u[m])) / 1;
  for (int m = 0; m < 4; ++m) ws[m] = (0.0 + static_cast<double>(c.w[m])) / 1;
  us[0] = rrf::max(us[0], u.Dgas_min);
}

// prep_rt: calc_rt_i
RR_HD void prep_cell(CellData& c, const radray::Units& u) {
  double us[5], ws[4];
  average_gas_1(c, u, us, ws);
  store_gas_state(c, gas_state_from_avg(us, ws, u), u);
}

// rt_step_1: calc_rt_1 with legacy_rt = 0 (B3, B4 fixed). sig_init =
// photoion_sigma(freq_init).
RR_HD void calc_rt_1_cell(CellData& c, const radray::Units& u, const double sig[3][5]) {
  using namespace radray::cp;
  using namespace radray::pc;
  for (int n = 0; n < 5; ++n) {
    c.nu[n] = 0;
    c.ngam[n] = 0;
    c.fgam[n] = 0;
  }
  for (int d = 0; d < 3; ++d) {
    c.gamma[d] = 0;
    c.lambda[d] = 0;
    c.dx[d] = 0;
  }
  const double rhog = static_cast<double>(c.rhog) * u.rhob_avg;
  const double nH = rhog * XHy / mH_cgs;
  const double nHe = rhog * YHe / mHe_cgs;
  const double nHI = nH * static_cast<double>(c.fHI);
  const double nHeI = nHe * static_cast<double>(c.fHeI);
  const double nHeII = nHe * static_cast<double>(c.fHeII);
  for (int i = 0; i < 5; ++i) {
    const double sm = ((0.0 + nHI * sig[0][i]) + nHeI * sig[1][i]) + nHeII * sig[2][i];
    c.tau[i] = static_cast<float>(sm * u.len_rt);
  }
}

// post_rt: calc_rt_f (gas -> RT grid, then tau from the gas)
RR_HD void post_cell(CellData& c, const radray::Units& u, const double sig[3][5]) {
  double us[5], ws[4];
  average_gas_1(c, u, us, ws);
  const GasStateP g = gas_state_from_avg(us, ws, u);
  store_gas_state(c, g, u);
  for (int i = 0; i < 5; ++i) {
    const double sm = ((0.0 + g.nHI * sig[0][i]) + g.nHeI * sig[1][i]) + g.nHeII * sig[2][i];
    c.tau[i] = static_cast<float>(sm * u.len_rt);
  }
}

struct PhysCtx {
  radray::Units u;
  RatesView rates;
  double dtsec_rt = 0;
};

// rt_tools.f90: radiative_physics for one RT cell (photo_source stats
// skipped, photo_cell for N_rt_gas = 1). Before the call c.ngam/c.nu/c.dx
// hold the deposited photon sums. Returns false if the cell was skipped.
RR_HD bool physics_cell(CellData& c, const PhysCtx& X) {
  using namespace radray;
  using namespace radray::pc;
  using namespace radray::cp;
  const Units& u = X.u;
  const double g = gamma_ideal;
  const double nu0[3] = {HI_eV, HeI_eV, HeII_eV};
  const double dtsec_rt = X.dtsec_rt;
  const double Nc = 1;
  float* gng = c.ngam;
  float* gnu = c.nu;
  float* gfg = c.fgam;
  float* gdx = c.dx;

  for (int n = 0; n < 5; ++n) {
    if (static_cast<double>(gng[n]) < u.ngam_min) {
      gnu[n] = static_cast<float>(rrk::freq_init(n));
      gng[n] = 0;
    } else {
      gnu[n] = gnu[n] / gng[n];
    }
  }
  const double ngam_tot = static_cast<double>(rrf::sum4<5>(gng));
  if (ngam_tot < u.ngam_min) {
    gdx[0] = gdx[1] = gdx[2] = 0;
  } else {
    const float sng = rrf::sum4<5>(gng);
    double dx[3];
    for (int d = 0; d < 3; ++d) dx[d] = static_cast<double>(gdx[d] / sng);
    const double dl = std::sqrt(((0.0 + dx[0] * dx[0]) + dx[1] * dx[1]) + dx[2] * dx[2]);
    if (dl > static_cast<double>(1e-9f)) {
      for (int d = 0; d < 3; ++d) gdx[d] = static_cast<float>(dx[d] / dl);
    } else {
      gdx[0] = gdx[1] = gdx[2] = 0;
    }
  }
  if (!(ngam_tot > u.ngam_min || c.fHI < 0.99f)) return false;

  // ---- photo_cell ----
  double us[5] = {0, 0, 0, 0, 0}, ws[4] = {0, 0, 0, 0};
  double nu[5], nrad[5];
  for (int n = 0; n < 5; ++n) {
    nu[n] = gnu[n];
    nrad[n] = static_cast<double>(gng[n]) * u.nH_avg * u.vol_rt;
  }
  double sigma_pht[3][5];
  photoion_sigma(nu, sigma_pht);
  {
    float* uu = c.u;
    float* ww = c.w;
    uu[0] = static_cast<float>(rrf::max(static_cast<double>(uu[0]), u.Dgas_min));
    const double rhog = static_cast<double>(uu[0]) * u.rho_unit;
    const double nH = rhog * XHy / mH_cgs;
    const double nHe = rhog * YHe / mHe_cgs;
    double nHI = rrf::min(static_cast<double>(ww[0]) * u.H_unit, nH);
    double nHeI = rrf::min(static_cast<double>(ww[1]) * u.He_unit, nHe);
    double nHeII = rrf::min(static_cast<double>(ww[2]) * u.He_unit, nHe - nHeI);
    double nHII = nH - nHI;
    double nHeIII = nHe - nHeI - nHeII;
    double ne = nHII + nHeII + 2 * nHeIII;
    double neff = nH + nHe + ne;
    double nff = nHII + nHeII + 4 * nHeIII;

    const float kin_f = (((0.0f + uu[1] * uu[1]) + uu[2] * uu[2]) + uu[3] * uu[3]) / uu[0] / 2.0f;
    const double kin = static_cast<double>(kin_f);
    double eth = (static_cast<double>(uu[4]) - kin) * u.eden_unit;
    double tmp = (g - 1) * eth / k_cgs / neff;
    if (tmp < u.temp_min) {
      tmp = u.temp_min;
      eth = neff * k_cgs * tmp / (g - 1);
    }
    int iT = 1 + static_cast<int>(std::round((std::log10(tmp) - lgT_min) / lgT_del));
    if (iT < 1) iT = 1;          // the serial port throws here (legacy only)
    if (iT > Ntemp) iT = Ntemp;
    const RateSetP rs = X.rates.legacy_physics ? rates_node(X.rates, iT) : rates_interp(X.rates, tmp);

    double ngam[5];
    for (int n = 0; n < 5; ++n) ngam[n] = static_cast<double>(gng[n]) * u.nH_avg;
    double gamma_pht[3][5], gamma_max[3][5];
    for (int inu = 0; inu < 5; ++inu)
      for (int sp = 0; sp < 3; ++sp)
        gamma_pht[sp][inu] = ngam[inu] * sigma_pht[sp][inu] * (u.len_rt / dtsec_rt);
    for (int sp = 0; sp < 3; ++sp)
      for (int inu = 0; inu < 5; ++inu) gamma_max[sp][inu] = gamma_pht[sp][inu];
    double gam_pht[3];
    auto sum_dim2 = [&](double out[3]) {
      for (int sp = 0; sp < 3; ++sp) {
        double t = 0;
        for (int inu = 0; inu < 5; ++inu) t += gamma_pht[sp][inu];
        out[sp] = t;
      }
    };
    sum_dim2(gam_pht);

    int n = 0;
    double ne_prev = static_cast<double>(3.40282346638528859812e+38F);
    double nH1 = 0, nHe1 = 0, nHe2 = 0;
    while (std::fabs(ne / ne_prev - 1) > static_cast<double>(1e-3f) && n < 5) {
      n = n + 1;
      ne_prev = ne;
      const double AH2 = rs.alpha[3] * ne * dtsec_rt;
      const double AHe2 = rs.alpha[4] * ne * dtsec_rt;
      const double AHe3 = rs.alpha[5] * ne * dtsec_rt;
      const double BH1 = (gam_pht[0] + rs.gcol[0] * ne) * dtsec_rt;
      const double BHe1 = (gam_pht[1] + rs.gcol[1] * ne) * dtsec_rt;
      const double BHe2 = (gam_pht[2] + rs.gcol[2] * ne) * dtsec_rt;
      const double sub1 = 1 + AHe3 + AHe2 + BHe2;
      const double sub2 = AHe2;
      const double sub3 = BHe1 - AHe3;
      const double sub4 = 1 + BHe1;
      const double sub5 = AHe3 * nHe;
      const double sub6 = sub1 * sub4 - sub2 * sub3;
      nH1 = (nHI + AH2 * nH) / (1 + AH2 + BH1);
      nHe1 = (nHeI * sub1 + (nHeII + sub5) * sub2) / sub6;
      nHe2 = (nHeI * sub3 + (nHeII + sub5) * sub4) / sub6;
      ne = (nH - nH1) + nHe2 + 2 * (nHe - nHe1 - nHe2);
      for (int inu = 0; inu < 5; ++inu) {
        const double num1 = (((0.0 + gamma_pht[0][inu] * nH1) + gamma_pht[1][inu] * nHe1) +
                             gamma_pht[2][inu] * nHe2) *
                            dtsec_rt * u.vol_gas;
        if (num1 > u.ngam_min) {
          const double x = nrad[inu] / num1;
          for (int sp = 0; sp < 3; ++sp)
            gamma_pht[sp][inu] = rrf::min(x * gamma_pht[sp][inu], gamma_max[sp][inu]);
        }
      }
      sum_dim2(gam_pht);
    }

    nHI = nH1;
    nHII = nH - nHI;
    nHeI = nHe1;
    nHeII = nHe2;
    nHeIII = nHe - nHeI - nHeII;
    ne = nHII + nHeII + 2 * nHeIII;
    neff = nH + nHe + ne;
    nff = nHII + nHeII + 4 * nHeIII;

    for (int inu = 0; inu < 5; ++inu) {
      const double d = rhog / u.rhob_avg;
      const double x = rrf::min(std::exp(-(d - rho_source) / rho_source), 1.0);
      const double num3[3] = {gamma_pht[0][inu] * nHI * dtsec_rt, gamma_pht[1][inu] * nHeI * dtsec_rt,
                              gamma_pht[2][inu] * nHeII * dtsec_rt};
      const double snum3 = ((0.0 + num3[0]) + num3[1]) + num3[2];
      nrad[inu] = rrf::max(nrad[inu] - x * snum3 * u.vol_gas, 0.0);
    }

    double lam_pht[3];
    for (int sp = 0; sp < 3; ++sp) {
      double t = 0;
      for (int inu = 0; inu < 5; ++inu) t += gamma_pht[sp][inu] * (nu[inu] - nu0[sp]) * eV2erg;
      lam_pht[sp] = t;
    }
    const double heat_rate = ((0.0 + lam_pht[0] * nHI) + lam_pht[1] * nHeI) + lam_pht[2] * nHeII;

    double sn_rate;
    if (tmp > 1e4 && rhog / u.rhob_avg > rho_source) {
      const double d = rhog / u.rhob_avg;
      const double x = rrf::max(1 - std::exp(-(d - rho_source) / rho_source), 0.0);
      const double tdyn = std::sqrt(3 * pi / (32 * G_cgs * rhog));
      sn_rate = x * eSN * eSF * rhog / tdyn * (c_cgs * c_cgs);
    } else {
      sn_rate = 0;
    }

    const double Tcmb = u.Tcmb;
    double cool_rate;
    if (X.rates.legacy_physics) {
      cool_rate = rs.cool[3] * ne * nHII + rs.cool[4] * ne * nHeII + rs.cool[5] * ne * nHeIII +
                  rs.cool[6] * ne * nHI + rs.cool[7] * ne * nHeI + rs.cool[8] * ne * nHeII +
                  rs.cool[9] * ne * nHI + rs.cool[10] * ne * nHeI + rs.cool[11] * ne * nHeII +
                  rs.cool[12] * ne * nff + rs.cool[13] * ne * (tmp - Tcmb) * rrf::pow4(Tcmb);
    } else {
      cool_rate = rs.cool[3] * ne * nHII + rs.cool[4] * ne * nHeII + rs.cool[5] * ne * nHeIII +
                  rs.cool[6] * ne * nHI + rs.cool[7] * ne * nHeI + rs.cool[8] * ne * nHeII +
                  rs.cool[9] * ne * nHI + rs.cool[10] * ne * ne * nHeII + rs.cool[11] * ne * nHeII +
                  rs.cool[12] * ne * nff + rs.cool[13] * ne * (tmp - Tcmb) * rrf::pow4(Tcmb) +
                  rs.cool[14] * ne * nHeII + rs.cool[15] * ne * nHeI;
    }
    {
      const double x = rhog / u.rhob_avg;
      if (x > 200) {
        cool_rate = cool_rate * rrf::pow2(200 / x);
        sn_rate = sn_rate * std::pow(200 / x, 1.5);
      }
    }
    const double net_rate = heat_rate + sn_rate - cool_rate;
    if (net_rate > 0) {
      eth = eth + net_rate * dtsec_rt;
    } else if (net_rate < 0) {
      eth = eth / (1 + std::log(1 - net_rate * dtsec_rt / eth));
    }
    tmp = (g - 1) * eth / k_cgs / neff;
    if (tmp < u.temp_min) {
      tmp = u.temp_min;
      eth = neff * k_cgs * tmp / (g - 1);
    }
    uu[4] = static_cast<float>(kin + eth / u.eden_unit);
    ww[0] = static_cast<float>(rrf::min(rrf::max(nHI, 0.0), nH) / u.H_unit);
    ww[1] = static_cast<float>(rrf::min(rrf::max(nHeI, 0.0), nHe) / u.He_unit);
    ww[2] = static_cast<float>(rrf::min(rrf::max(nHeII, 0.0), nHe) / u.He_unit);
    for (int m = 0; m < 5; ++m) us[m] = us[m] + static_cast<double>(uu[m]);
    for (int m = 0; m < 4; ++m) ws[m] = ws[m] + static_cast<double>(ww[m]);
    for (int sp = 0; sp < 3; ++sp) {
      c.gamma[sp] = static_cast<float>(static_cast<double>(c.gamma[sp]) + gam_pht[sp] / Nc);
      c.lambda[sp] = static_cast<float>(static_cast<double>(c.lambda[sp]) + lam_pht[sp] / Nc);
    }
  }
  for (int n = 0; n < 5; ++n)
    gfg[n] = static_cast<float>(static_cast<double>(gfg[n]) + nrad[n] / u.vol_rt / u.nH_avg);
  for (int m = 0; m < 5; ++m) us[m] = us[m] / Nc;
  for (int m = 0; m < 4; ++m) ws[m] = ws[m] / Nc;
  store_gas_state(c, gas_state_from_avg(us, ws, u), u);
  for (int m = 0; m < 5; ++m) {
    if (static_cast<double>(gng[m]) < u.ngam_min) {
      gfg[m] = 0;
      gng[m] = 0;
    } else {
      gfg[m] = gfg[m] / gng[m];
      gng[m] = gng[m] * gfg[m];
    }
  }
  return true;
}

}  // namespace rrk
