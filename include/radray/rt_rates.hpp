// Atomic rates and photoionisation cross sections.
// Port of rt_rates.f90 (calc_atomic_rates, photoion_sigma); METAL off.
#pragma once

#include <array>
#include <vector>

#include "radray/constants.hpp"

namespace radray {

// Tabulated rates on the temperature grid T_i = 10**((i-1)*lgT_del+lgT_min),
// i = 1..Ntemp. Accessors take Fortran (1-based) indices.
struct AtomicRates {
  // coef_cool(16,Ntemp) (slot 16: P10, legacy_physics = 0 only), alpha_rec(6,Ntemp), gamma_col(3,Ntemp)
  std::vector<double> coef_cool_ = std::vector<double>(16 * cp::Ntemp, 0.0);
  std::vector<double> alpha_rec_ = std::vector<double>(6 * cp::Ntemp, 0.0);
  std::vector<double> gamma_col_ = std::vector<double>(3 * cp::Ntemp, 0.0);

  double& coef_cool(int n, int iT) { return coef_cool_[(n - 1) + 16 * (iT - 1)]; }
  double& alpha_rec(int n, int iT) { return alpha_rec_[(n - 1) + 6 * (iT - 1)]; }
  double& gamma_col(int n, int iT) { return gamma_col_[(n - 1) + 3 * (iT - 1)]; }
  double coef_cool(int n, int iT) const { return coef_cool_[(n - 1) + 16 * (iT - 1)]; }
  double alpha_rec(int n, int iT) const { return alpha_rec_[(n - 1) + 6 * (iT - 1)]; }
  double gamma_col(int n, int iT) const { return gamma_col_[(n - 1) + 3 * (iT - 1)]; }
};

// All rates at one temperature. Index n-1 holds Fortran slot n:
// alpha[n-1] = alpha_rec(n), gcol[n-1] = gamma_col(n), cool[n-1] = coef_cool(n).
struct RateSet {
  std::array<double, 6> alpha{};
  std::array<double, 3> gcol{};
  std::array<double, 16> cool{};
};

// rt_rates.f90: calc_atomic_rates. legacy_physics = false applies the
// physics fixes P1-P8 (docs/port_notes.md section 5b) to the tables.
void calc_atomic_rates(AtomicRates& r, bool legacy_physics = true);

// The formulas of calc_atomic_rates evaluated at any temperature T [K];
// at a table node it returns exactly the table entries.
RateSet atomic_rates_at(double T, bool legacy_physics = true);

// Table entries at node iT (1..Ntemp): the legacy nearest-node lookup.
RateSet rates_node(const AtomicRates& r, int iT);

// P9 (legacy_physics = 0; harness rt_rates.f90 rates_interp): rates at
// temperature tmp, interpolated linearly in log R versus log T between the
// two enclosing nodes; the edge node outside the table. See rt_rates.cpp.
RateSet rates_interp(const AtomicRates& r, double tmp);
// Same, between nodes i and i+1 (1 <= i < Ntemp) with weight w in [0,1).
RateSet rates_interp(const AtomicRates& r, int i, double w);
// One rate between node values a and b (harness rate_lerp).
double rate_lerp(double a, double b, double w);

// rt_rates.f90: photoion_sigma(nu) -> sigma(3,Nfreq) in cm^2.
// Returned as sigma[species 0..2][freq 0..4].
using SigmaTable = std::array<std::array<double, kNfreq>, 3>;
SigmaTable photoion_sigma(const std::array<double, kNfreq>& nu);

}  // namespace radray
