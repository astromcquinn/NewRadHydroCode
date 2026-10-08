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
  // coef_cool(15,Ntemp), alpha_rec(6,Ntemp), gamma_col(3,Ntemp)
  std::vector<double> coef_cool_ = std::vector<double>(15 * cp::Ntemp, 0.0);
  std::vector<double> alpha_rec_ = std::vector<double>(6 * cp::Ntemp, 0.0);
  std::vector<double> gamma_col_ = std::vector<double>(3 * cp::Ntemp, 0.0);

  double& coef_cool(int n, int iT) { return coef_cool_[(n - 1) + 15 * (iT - 1)]; }
  double& alpha_rec(int n, int iT) { return alpha_rec_[(n - 1) + 6 * (iT - 1)]; }
  double& gamma_col(int n, int iT) { return gamma_col_[(n - 1) + 3 * (iT - 1)]; }
  double coef_cool(int n, int iT) const { return coef_cool_[(n - 1) + 15 * (iT - 1)]; }
  double alpha_rec(int n, int iT) const { return alpha_rec_[(n - 1) + 6 * (iT - 1)]; }
  double gamma_col(int n, int iT) const { return gamma_col_[(n - 1) + 3 * (iT - 1)]; }
};

// rt_rates.f90: calc_atomic_rates
void calc_atomic_rates(AtomicRates& r);

// rt_rates.f90: photoion_sigma(nu) -> sigma(3,Nfreq) in cm^2.
// Returned as sigma[species 0..2][freq 0..4].
using SigmaTable = std::array<std::array<double, kNfreq>, 3>;
SigmaTable photoion_sigma(const std::array<double, kNfreq>& nu);

}  // namespace radray
