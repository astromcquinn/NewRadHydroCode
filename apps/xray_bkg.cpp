// The X-ray background of xray_heating = 1 without a simulation
// (docs/xray_heating.md): the SFR density of the source model, the X-ray
// emissivity, the rates per atom, and a one-zone thermal history of gas at
// the mean density outside the ionised regions (X-ray heating and
// ionisation, Compton coupling to the CMB, adiabatic cooling, case B
// recombination; helium stays neutral). What the preheating model predicts
// for a set of parameters.
//
// usage: xray_bkg [params.txt ...] [--set key=value ...] [--zmin 5]
//                 [--zout-step 0.5] [--xe0 2e-4] [--spectrum z]
//   params files: `key = value` lines; the xray_* keys and the source-model
//   keys (sfr_model = abundance, am_*, l_per_sfr, tvir_min, m_min) are used,
//   others ignored (so a run's params.txt or data/uv_lf/am_params_*.txt
//   can be passed as is).
//   --xe0: the ionised fraction at xray_zmax (the residual of
//   recombination; the simulations start fully neutral).
//   --spectrum z: print J_E(E) at z for a neutral IGM instead.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "radray/constants.hpp"
#include "radray/halo_sources.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/xray_background.hpp"

using namespace radray;

namespace {

std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const auto b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

void set_param(XrayParams& xp, HaloSourceParams& hp, const std::string& line) {
  const auto eq = line.find('=');
  if (eq == std::string::npos) return;
  const std::string key = trim(line.substr(0, eq)), val = trim(line.substr(eq + 1));
  if (key.empty()) return;
  if (!set_xray_param(xp, key, val)) set_halo_param(hp, key, val);
}

double hubble(double z) {
  using namespace cp;
  const double a = 1 / (1 + z);
  return pc::H0_cgs * hubble0 *
         std::sqrt(omegam / (a * a * a) + omegal + omegar / (a * a * a * a));
}

// Case B HII recombination, Hui & Gnedin (1997)
double alpha_B(double T) {
  const double l = 2 * 157807.0 / T;
  return 2.753e-14 * std::pow(l, 1.5) / std::pow(1 + std::pow(l / 2.740, 0.407), 2.242);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    XrayParams xp;
    xp.on = true;
    HaloSourceParams hp;
    double zmin = 5, zout = 0.5, xe0 = 2e-4, zspec = -1;
    for (int a = 1; a < argc; ++a) {
      const std::string s = argv[a];
      if (s == "--set" && a + 1 < argc) {
        set_param(xp, hp, argv[++a]);
      } else if (s == "--zmin" && a + 1 < argc) {
        zmin = std::stod(argv[++a]);
      } else if (s == "--zout-step" && a + 1 < argc) {
        zout = std::stod(argv[++a]);
      } else if (s == "--xe0" && a + 1 < argc) {
        xe0 = std::stod(argv[++a]);
      } else if (s == "--spectrum" && a + 1 < argc) {
        zspec = std::stod(argv[++a]);
      } else if (!s.empty() && s[0] == '-') {
        throw std::runtime_error("unknown option " + s);
      } else {
        std::ifstream f(s);
        if (!f) throw std::runtime_error("cannot read " + s);
        std::string line;
        while (std::getline(f, line)) {
          const auto h = line.find('#');
          set_param(xp, hp, h == std::string::npos ? line : line.substr(0, h));
        }
      }
    }
    // the one-zone tool supplies its own absorber fractions
    const XrayParams::Opacity opacity = xp.opacity;
    xp.opacity = XrayParams::Opacity::Box;
    const XrayBackground bg(xp, hp.model);

    std::printf("# xray_bkg: log L_X/SFR = %.3f in [%.3g, %.3g] keV, emitted %.3g-%.3g keV, alpha = %.3g\n",
                xp.log_lx_sfr, xp.e0_keV, xp.eband_keV, xp.e0_keV, xp.emax_keV, xp.alpha);
    std::printf("# mass function ST A = %g a = %g p = %g; tvir_min = %g K, l_per_sfr = %g\n",
                xp.st.A, xp.st.a, xp.st.p, xp.tvir_min >= 0 ? xp.tvir_min : hp.model.Tvir_min,
                hp.model.L_per_sfr);

    if (zspec >= 0) {
      const double xf[3] = {1, 1, 0};
      std::vector<double> E, J;
      bg.spectrum(zspec, xf, E, J);
      std::printf("# J_E at z = %g, neutral IGM\n# E[keV] J_E[erg/s/cm^2/sr/keV]\n", zspec);
      for (std::size_t i = 0; i < E.size(); i += 10) std::printf("%.5g %.5e\n", E[i], J[i]);
      return 0;
    }

    using namespace pc;
    using namespace cp;
    const double nH0 = rhoc_cgs * omegab * hubble0 * hubble0 * XHy / mH_cgs;
    const double nHe0 = rhoc_cgs * omegab * hubble0 * hubble0 * YHe / mHe_cgs;
    // Compton coefficient of rt_rates.cpp (coef_cool(14)), in double
    const double cC = 4 * sS_cgs / c_cgs * (8.0 / 3.0) * 3.36e-24 / (gamma_ideal - 1);
    // start: adiabatic since decoupling, the code's temp_min (calc_sim_var)
    double z = xp.z_max;
    double T = Tcmb0 * (1 + z) * (1 + z) / folded::temp_x;
    double x = xe0;
    double edep_tot = 0;  // heat per H atom, eV
    const double ds = 0.0005, ds_rates = 0.0025;
    XrayRates r{};
    double s_next_rates = -1, s_next_out = std::log(1 + z);
    std::printf("# one-zone gas at the mean density, opacity %s; heat = X-ray heat deposited per H atom\n",
                opacity == XrayParams::Opacity::None ? "none" : "of the one-zone neutral IGM");
    std::printf("# z M_min[Msun] sfrd[Msun/yr/Mpc^3] eps_X[erg/s/Mpc^3] gamma_HI[1/s] phi_HI[erg/s]"
                " Hx_per_H[erg/s] x_e T[K] T_cmb[K] T_adiab[K] heat[eV]\n");
    double s = std::log(1 + z);
    const double s_end = std::log(1 + zmin);
    while (s > s_end - 1e-12) {
      z = std::exp(s) - 1;
      if (s_next_rates < 0 || s <= s_next_rates) {
        double xf[3] = {1 - x, 1, 0};
        if (opacity == XrayParams::Opacity::None) xf[0] = xf[1] = xf[2] = 0;
        r = bg.rates(z, xf);
        s_next_rates = s - ds_rates;
      }
      const double nH = nH0 * std::pow(1 + z, 3), nHe = nHe0 * std::pow(1 + z, 3);
      const double nHI = (1 - x) * nH, ne = x * nH;
      double fheat, fHI, fHeI;
      rrk::xray_secondary(x, fheat, fHI, fHeI);
      const double edep = r.phi[0] * nHI + r.phi[1] * nHe;
      const double Hx = fheat * edep;
      if (s <= s_next_out + 1e-9) {
        const double Tad = Tcmb0 * (1 + z) * (1 + z) / folded::temp_x;
        std::printf("%.3f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", z, bg.m_min_msun(z),
                    r.sfrd, r.eps, r.gamma[0], r.phi[0], Hx / nH, x, T, Tcmb0 * (1 + z), Tad,
                    edep_tot);
        s_next_out = std::log(1 + std::max(zmin, z - zout));
      }
      const double H = hubble(z);
      const double dt = ds / H;
      const double Tg = Tcmb0 * (1 + z);
      const double ion = (1 - x) * (r.gamma[0] + (nHI > 0 ? fHI * edep / (HI_eV * eV_cgs * nHI) : 0));
      const double dxdt = ion - alpha_B(T) * ne * x;
      const double ntot = nH * (1 + x) + nHe;
      const double dTdt = -2 * H * T + 2.0 / 3.0 * (Hx - cC * ne * (T - Tg) * std::pow(Tg, 4)) / (k_cgs * ntot) -
                          T * nH * dxdt / ntot;
      T = std::max(T + dTdt * dt, 1.0);
      x = std::min(std::max(x + dxdt * dt, 0.0), 1.0);
      edep_tot += Hx / nH * dt / eV_cgs;
      s -= ds;
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "xray_bkg: " << e.what() << "\n";
    return 1;
  }
}
