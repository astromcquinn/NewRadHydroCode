// The Ewald force of a point mass in the periodic box, for the PM force
// tests (test_pm_refine.cpp, test_amrex_pm.cpp).
#pragma once

#include <cmath>

#include "radray/constants.hpp"

namespace pm_ewald {

using namespace radray;

// Gradient of the periodic Green's function psi of the unit box,
// lap(psi) = delta - 1 (Ewald split, alpha = 2; converged to ~1e-10).
inline void ewald_grad(const double r[3], double g[3]) {
  const double alpha = 2.0, pi = pc::pi;
  g[0] = g[1] = g[2] = 0;
  for (int nx = -3; nx <= 3; ++nx)
    for (int ny = -3; ny <= 3; ++ny)
      for (int nz = -3; nz <= 3; ++nz) {
        const double d[3] = {r[0] - nx, r[1] - ny, r[2] - nz};
        const double dd = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (dd > 3.5) continue;
        const double w = (std::erfc(alpha * dd) + 2 * alpha * dd / std::sqrt(pi) *
                                                      std::exp(-alpha * alpha * dd * dd)) /
                         (4 * pi * dd * dd * dd);
        for (int m = 0; m < 3; ++m) g[m] += d[m] * w;
      }
  const int K = 8;
  for (int mx = -K; mx <= K; ++mx)
    for (int my = -K; my <= K; ++my)
      for (int mz = -K; mz <= K; ++mz) {
        if (mx == 0 && my == 0 && mz == 0) continue;
        const double k[3] = {2 * pi * mx, 2 * pi * my, 2 * pi * mz};
        const double k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
        const double w = std::exp(-k2 / (4 * alpha * alpha)) / k2 *
                         std::sin(k[0] * r[0] + k[1] * r[1] + k[2] * r[2]);
        for (int m = 0; m < 3; ++m) g[m] += w * k[m];
      }
}

inline double wrap(double x, double L) { return x - L * std::floor(x / L + 0.5); }

// Exact force (gas grid units) at separation dx (gas cells) from a point
// mass carrying all the dark matter, (1-fb) Ngas^3, scalefactor a:
// lap(phi) = (2a/3)(rho - mean), F = -grad(phi).
inline void ewald_force(const double dx[3], int Ngas, double a, double F[3]) {
  const double L = Ngas;
  const double ru[3] = {wrap(dx[0], L) / L, wrap(dx[1], L) / L, wrap(dx[2], L) / L};
  double g[3];
  ewald_grad(ru, g);
  const double M = (1 - cp::omegab / cp::omegam) * L * L * L;
  for (int m = 0; m < 3; ++m) F[m] = -(2 * a / 3) * (M / (L * L)) * g[m];
}

}  // namespace pm_ewald
