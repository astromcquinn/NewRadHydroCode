// Analytic halo mass functions (Press-Schechter, Sheth-Tormen) for
// validating the halo finder. Linear power spectrum: BBKS transfer
// function with the Sugiyama (1995) shape parameter, normalised to sigma8;
// linear growth from the Heath (1977) integral. Accurate to ~10-20%, which
// is all the "rough" comparison of the halo finder needs.
#pragma once

#include "radray/constants.hpp"

namespace radray {

struct MassFunctionCosmology {
  // cosmo_parameters.f90 (sigma8, nsinit are not in constants.hpp)
  double omegam = cp::omegam;
  double omegab = cp::omegab;
  double hubble0 = cp::hubble0;
  double sigma8 = static_cast<double>(0.80f);
  double ns = static_cast<double>(0.96f);
};

class MassFunction {
 public:
  explicit MassFunction(const MassFunctionCosmology& c = {});

  // Linear growth factor normalised to D(0) = 1 (flat LCDM).
  double growth(double z) const;
  // rms linear overdensity at z=0 in a top hat of mass M (Msun/h)
  double sigma0(double M_msun_h) const;
  // Comoving mean matter density, Msun/h per (Mpc/h)^3
  double rho_mean() const;

  enum class Kind { PressSchechter, ShethTormen };
  // dn/dlnM in (Mpc/h)^-3 at mass M (Msun/h) and redshift z.
  double dndlnM(double M_msun_h, double z, Kind k) const;
  // Number of halos per (Mpc/h)^3 with M1 <= M < M2 (Msun/h).
  double n_between(double M1, double M2, double z, Kind k) const;

 private:
  MassFunctionCosmology c_;
  double norm_ = 1;
  double sigma_unnorm(double R) const;
};

}  // namespace radray
