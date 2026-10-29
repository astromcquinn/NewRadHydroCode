// thermal_init = recfast (thermal_history.hpp, docs/thermal_history.md).
#include "radray/thermal_history.hpp"

#include <cmath>

#include "radray/constants.hpp"
#include "radray/parallel.hpp"

namespace radray {

namespace {
#include "recfast_table.inc"

// log-linear interpolation of column col in ln(1+z)
double table(double z, int col) {
  const double s = std::log1p(z);
  if (z <= kRecfast[0][0]) return kRecfast[0][col];
  for (int i = 1; i < kRecfastRows; ++i) {
    if (z <= kRecfast[i][0]) {
      const double s0 = std::log1p(kRecfast[i - 1][0]), s1 = std::log1p(kRecfast[i][0]);
      const double w = (s - s0) / (s1 - s0);
      return std::exp(std::log(kRecfast[i - 1][col]) * (1 - w) + std::log(kRecfast[i][col]) * w);
    }
  }
  return -1;  // above the table
}
}  // namespace

double recfast_T(double z) {
  const double T = table(z, 2);
  return T > 0 ? T : cp::Tcmb0 * (1 + z);
}

double recfast_xe(double z) {
  const double x = table(z, 1);
  return x > 0 ? x : kRecfast[kRecfastRows - 1][1];
}

std::int64_t compton_residual_step(const Units& u, double a, double dt_sec, double C14,
                                   GasGrid& gas) {
  using namespace pc;
  using namespace cp;
  const double z = 1 / a - 1;
  const double Tg = Tcmb0 / a;
  const double xres = recfast_xe(z);
  const double g = gamma_ideal;
  const double Tg4 = Tg * Tg * Tg * Tg;
  const bool has_e = !gas.e.empty();
  const std::int64_t nc = gas.ncell;
  std::int64_t changed = 0;
#pragma omp parallel for schedule(static) reduction(+ : changed) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c) {
    float* uu = &gas.u[c * 5];
    const float* ww = &gas.w[c * 4];
    // the densities and the thermal energy as photo_gas_cell (rt_cell.hpp)
    const double rhog = static_cast<double>(uu[0]) * u.rho_unit;
    if (!(rhog > 0)) continue;
    const double nH = rhog * XHy / mH_cgs;
    const double nHe = rhog * YHe / mHe_cgs;
    const double nHI = std::fmin(static_cast<double>(ww[0]) * u.H_unit, nH);
    const double nHeI = std::fmin(static_cast<double>(ww[1]) * u.He_unit, nHe);
    const double nHeII = std::fmin(static_cast<double>(ww[2]) * u.He_unit, nHe - nHeI);
    const double ne = (nH - nHI) + nHeII + 2 * (nHe - nHeI - nHeII);
    const double ne_add = xres * nH - ne;
    if (!(ne_add > 0)) continue;
    const double neff = nH + nHe + ne;
    const float kin_f = (((0.0f + uu[1] * uu[1]) + uu[2] * uu[2]) + uu[3] * uu[3]) / uu[0] / 2.0f;
    const double kin = static_cast<double>(kin_f);
    const double eth = has_e ? static_cast<double>(gas.e[c]) * u.eden_unit
                             : (static_cast<double>(uu[4]) - kin) * u.eden_unit;
    if (!(eth > 0)) continue;
    const double rate = (g - 1) * C14 * ne_add * Tg4 / (k_cgs * neff);  // 1/t_C
    const double eeq = neff * k_cgs * Tg / (g - 1);
    const double enew = eeq + (eth - eeq) * std::exp(-rate * dt_sec);
    uu[4] = static_cast<float>(kin + enew / u.eden_unit);
    if (has_e) gas.e[c] = static_cast<float>(enew / u.eden_unit);
    changed += 1;
  }
  return changed;
}

}  // namespace radray
