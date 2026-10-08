// hmf_table: analytic halo mass functions (mass_function.hpp) for comparison
// with simulated FoF catalogs (tools/analysis/hmf_compare.py).
//
//   hmf_table z [logM_min logM_max dlogM] [--tinker-zmax Z]
//
// Masses in Msun (not Msun/h), densities per comoving Mpc^3 (not h^-3), in
// the simulation cosmology (cosmo_parameters.f90). Columns:
//   logM  dn/dlog10M[PS]  dn/dlog10M[ST]  dn/dlog10M[Tinker08]
//   n(>M)[PS]  n(>M)[ST]  n(>M)[Tinker08]
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "radray/mass_function.hpp"

using namespace radray;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: hmf_table z [logM_min logM_max dlogM] [--tinker-zmax Z]\n");
    return 2;
  }
  double z = std::atof(argv[1]), lo = 8, hi = 13, d = 0.1, tzmax = 1e30;
  int pos = 0;
  for (int i = 2; i < argc; ++i) {
    if (!std::strcmp(argv[i], "--tinker-zmax") && i + 1 < argc) {
      tzmax = std::atof(argv[++i]);
    } else {
      const double v = std::atof(argv[i]);
      if (pos == 0) lo = v;
      else if (pos == 1) hi = v;
      else if (pos == 2) d = v;
      ++pos;
    }
  }
  MassFunctionCosmology c;
  MassFunction mf(c);
  mf.tinker_zmax_evol = tzmax;
  const double h = c.hubble0;
  using K = MassFunction::Kind;
  const K kinds[3] = {K::PressSchechter, K::ShethTormen, K::Tinker08};
  std::printf("# z = %g; logM [Msun]; dn/dlog10M [Mpc^-3]; n(>M) [Mpc^-3]; PS ST Tinker08\n", z);
  // dn/dlog10M on a fine grid from lo to 10^16.5 Msun, then the cumulative
  // n(>M) by the trapezoid rule from the top (one pass per mass function).
  const double fine = std::min(d, 0.02), top = 16.5;
  const int nf = static_cast<int>(std::ceil((top - lo) / fine));
  std::vector<std::array<double, 3>> dn(nf + 1), cum(nf + 1);
  for (int i = 0; i <= nf; ++i)
    for (int k = 0; k < 3; ++k)  // dn/dlog10M = ln10 dn/dlnM; (Mpc/h)^-3 -> Mpc^-3
      dn[i][k] = std::log(10.0) * mf.dndlnM(std::pow(10.0, lo + i * fine) * h, z, kinds[k]) * h * h * h;
  cum[nf] = {0, 0, 0};
  for (int i = nf - 1; i >= 0; --i)
    for (int k = 0; k < 3; ++k) cum[i][k] = cum[i + 1][k] + 0.5 * (dn[i][k] + dn[i + 1][k]) * fine;
  for (double lm = lo; lm <= hi + 1e-9; lm += d) {
    const int i = std::min(nf, static_cast<int>(std::lround((lm - lo) / fine)));
    std::printf("%.3f", lm);
    for (int k = 0; k < 3; ++k) std::printf(" %.6e", dn[i][k]);
    for (int k = 0; k < 3; ++k) std::printf(" %.6e", cum[i][k]);
    std::printf("\n");
  }
  return 0;
}
