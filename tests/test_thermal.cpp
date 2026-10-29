// thermal_init = recfast (thermal_history.hpp, docs/thermal_history.md): the
// RECFAST table and its interpolation, the residual-electron Compton step
// against the exact relaxation (with and without the dual-energy field, an
// ionised cell untouched), and calc_sim_var's floors (bit-identical by
// default, the recfast floors otherwise), and recfast_table = FILE (the
// table of another cosmology; last, since it replaces the compiled table).
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>

#include "check.hpp"
#include "radray/config.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/grid.hpp"
#include "radray/rt_rates.hpp"
#include "radray/thermal_history.hpp"

using namespace radray;

int main() {
  using namespace pc;
  using namespace cp;

  // --- the table (data/thermal/recfast_sim_cosmology.txt) ---
  CHECK_CLOSE(recfast_T(100.0), 1.681108e+02, 1e-6);
  CHECK_CLOSE(recfast_T(20.0), 9.351560e+00, 1e-6);
  CHECK_CLOSE(recfast_xe(100.0), 2.660750e-04, 1e-6);
  CHECK_CLOSE(recfast_T(3000.0), Tcmb0 * 3001, 1e-12);  // above the table: T_CMB
  double prev = 0;
  for (double z = 0; z < 1500; z += 0.7) {  // T rises with z, below T_CMB
    const double T = recfast_T(z);
    CHECK(T > prev && T <= Tcmb0 * (1 + z) * 1.0001);
    CHECK(recfast_xe(z) > 1e-4 && recfast_xe(z) < 1.2);
    prev = T;
  }
  // below z ~ 30 the adiabat of decoupled gas: T/(1+z)^2 ~ const
  CHECK_CLOSE(recfast_T(10) / 121, recfast_T(20) / 441, 3e-2);

  // --- calc_sim_var ---
  {
    Config cfg = Config::make(8, 4, 8, 4.0);
    SimClock clk;
    clk.redshift = 300;
    clk.scalefactor = 1 / 301.0;
    Units u0, u1;
    calc_sim_var(cfg, clk, u0);
    CHECK(u0.Thydro_min == 1 / u0.temp1_unit);  // the default: bit-identical
    CHECK(u0.temp_min == u0.Tcmb && u0.temp_max == u0.Tcmb);  // z > 216: T_CMB
    cfg.thermal_init = 1;
    cfg.thermal_floor = 0.01;
    calc_sim_var(cfg, clk, u1);
    CHECK(u1.temp_min == 1.0);  // legacy_physics = 1: the rate table's edge
    CHECK(u1.temp_max == 1e10);
    CHECK_CLOSE(u1.Thydro_min * u1.temp1_unit, 0.01, 1e-12);
    cfg.legacy_physics = false;
    calc_sim_var(cfg, clk, u1);
    CHECK(u1.temp_min == 0.01);
  }

  // --- the Compton step on one cell ---
  {
    const Config cfg = Config::make(8, 4, 8, 4.0);
    const double z = 50, a = 1 / (1 + z);
    SimClock clk;
    clk.redshift = z;
    clk.scalefactor = a;
    Units u;
    calc_sim_var(cfg, clk, u);
    AtomicRates ar;
    calc_atomic_rates(ar, true);
    const double C14 = ar.coef_cool(14, 1);
    // C14 = 4 a_rad sigma_T k / (m_e c) to the precision of its literals
    CHECK_CLOSE(C14, 4 * (4 * sS_cgs / c_cgs) * sT_cgs * k_cgs / (9.10938e-28 * c_cgs), 2e-3);
    const double T0 = 20, dt = 3e14;
    for (int with_e = 0; with_e < 2; ++with_e) {
      GasGrid g;
      g.allocate(1);
      const double D = u.Dgas_avg;
      const double nH = D * u.rho_unit * XHy / mH_cgs, nHe = D * u.rho_unit * YHe / mHe_cgs;
      const double eth = (nH + nHe) * k_cgs * T0 / (gamma_ideal - 1) / u.eden_unit;
      g.u[0] = static_cast<float>(D);
      g.u[1] = static_cast<float>(D * 1e-3);  // some momentum: kin kept
      g.u[2] = g.u[3] = 0;
      const float kin = g.u[1] * g.u[1] / g.u[0] / 2.0f;
      g.u[4] = static_cast<float>(static_cast<double>(kin) + eth);
      g.w[0] = static_cast<float>(D * XHy);
      g.w[1] = static_cast<float>(D * YHe);
      g.w[2] = g.w[3] = 0;
      if (with_e) g.e.assign(1, static_cast<float>(eth));
      CHECK(compton_residual_step(u, a, dt, C14, g) == 1);
      const double e1 = with_e ? g.e[0] : static_cast<double>(g.u[4]) - kin;
      const double T1 = e1 / eth * T0;
      const double Tg = Tcmb0 / a, xe = recfast_xe(z);
      const double rate = (gamma_ideal - 1) * C14 * xe * nH * std::pow(Tg, 4) / (k_cgs * (nH + nHe));
      const double Tx = Tg + (T0 - Tg) * std::exp(-rate * dt);
      std::printf("Compton step (e = %d): T %.3f -> %.4f K, exact %.4f (t_C = %.3g yr)\n", with_e, T0,
                  T1, Tx, 1 / rate / yr2sec);
      // (1e-5: the float neutral fraction leaves n_e ~ 1e-8 n_H of its own)
      CHECK_CLOSE(T1, Tx, 1e-4);
      CHECK(T1 > T0 && T1 < Tg);
    }
    // an ionised cell: its own electrons exceed the residual ones
    GasGrid g;
    g.allocate(1);
    g.u[0] = static_cast<float>(u.Dgas_avg);
    g.u[4] = 1e-3f;
    g.w[0] = g.w[1] = g.w[2] = 0;
    const float u4 = g.u[4];
    CHECK(compton_residual_step(u, a, dt, C14, g) == 0);
    CHECK(g.u[4] == u4);
  }

  // --- recfast_table = FILE (make_recfast_table.py --out) ---
  {
    Cosmology planck;
    planck.omegam = static_cast<float>(0.305);  // params.txt values are rounded to real(4)
    planck.omegab = static_cast<float>(0.048);
    planck.hubble0 = static_cast<float>(0.68);
    planck.derive();
    const Cosmology orig;
    CHECK(recfast_table_matches(orig) && !recfast_table_matches(planck));
    const char* path = "test_thermal_recfast.txt";
    auto write = [&](bool with_cosmo) {
      std::ofstream f(path);
      f << "# a test table\n";
      if (with_cosmo) f << "# cosmology omegam=0.305 omegab=0.048 hubble0=0.68\n";
      f << "# z  x_e  T_gas[K]  T_cmb[K]\n  200.00 2.0e-04 4.0e+02 5.4e+02\n  100.00 1.0e-04 1.0e+02 2.7e+02\n"
           " 2000.00 1.0e+00 5.0e+03 5.4e+03\n";
    };
    auto throws = [&](const Cosmology& c) {
      try {
        load_recfast_table(path, c);
      } catch (const std::runtime_error&) {
        return true;
      }
      return false;
    };
    write(false);
    CHECK(throws(planck));  // no cosmology line
    write(true);
    CHECK(throws(orig));  // another cosmology than the run's
    CHECK(recfast_table_matches(orig));  // the compiled table still in use
    CHECK(!throws(planck));
    CHECK(recfast_table_matches(planck) && !recfast_table_matches(orig));
    CHECK(recfast_T(100.0) == 1e2 && recfast_xe(200.0) == 2e-4);  // the rows, sorted
    // log-linear in ln(1+z) between the rows
    const double w = (std::log1p(150.0) - std::log1p(100.0)) / (std::log1p(200.0) - std::log1p(100.0));
    CHECK_CLOSE(recfast_T(150.0), std::exp(std::log(1e2) * (1 - w) + std::log(4e2) * w), 1e-14);
    CHECK(recfast_T(50.0) == 1e2);  // below the table: the lowest row
    CHECK(recfast_T(3000.0) == Tcmb0 * 3001);  // above it: T_CMB
    std::remove(path);
  }

  TEST_MAIN_END();
}
