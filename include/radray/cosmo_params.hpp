// Runtime cosmology and the preheating temperature floor
// (docs/production_run.md section 10; docs/zoom_turbulence.md sections 4.1
// item 5 and 4.5).
//
// The original fixes the cosmology at compile time (cosmo_parameters.f90:
// omegam = 0.27, omegab = 0.045, hubble0 = 0.70, sigma8 = 0.80, nsinit =
// 0.96, all real(4) literals stored in real(8) parameters). Here it is a
// run-time parameter set carried by Config::cosmo, with the original's
// values as the defaults: a run without the keys is bit-identical to the
// original, and so is a run that writes the original's values explicitly,
// because every value is rounded to real(4) as the original's literals are.
// The derived parameters follow the original's definitions: a flat universe
// omegal = 1 - omegam (omegar not subtracted, as in cosmo_parameters.f90),
// omegar = 4.48E-7*(1+0.69)*Tcmb0**4/hubble0**2 and the matter-temperature
// decoupling factor temp_x = 1E3*(omegab*hubble0**2)**0.4 of calc_sim_var.
// wde = -1 and Tcmb0 stay the original's constants.
#pragma once

#include <string>

#include "radray/constants.hpp"

namespace radray {

struct Cosmology {
  // params.txt keys omegam, omegab, hubble0, sigma8, ns
  double omegam = cp::omegam;
  double omegab = cp::omegab;
  double hubble0 = cp::hubble0;
  double sigma8 = cp::sigma8;  // the IC tool and the analytic mass functions only
  double ns = cp::nsinit;      // the IC tool and the analytic mass functions only
  // Derived (derive()); the defaults are the original's folded constants.
  double omegal = cp::omegal;
  double omegar = cp::omegar;
  double temp_x = folded::temp_x;

  // Sets omegal, omegar and temp_x from omegam, omegab and hubble0, and
  // checks the ranges (throws std::runtime_error). With the original's
  // inputs the derived values are the original's bit for bit.
  void derive();
  // omegam, omegab and hubble0 are the original's (sigma8 and ns do not
  // enter the simulation).
  bool is_original() const;
  // One line for the run logs, e.g. "omegam 0.3111 omegab 0.049 ...".
  std::string summary() const;
};

// The runtime cosmology's derived parameters as the original computes them,
// for any inputs (no shortcut to the folded constants): the unit test checks
// that these reproduce the folded constants for the original's inputs.
double cosmology_omegar(double hubble0);
double cosmology_temp_x(double omegab, double hubble0);

// Parses one params.txt entry of the runtime cosmology into c (omegam,
// omegab, hubble0, sigma8, ns; rounded to real(4)). Returns false for any
// other key. Call derive() after the last key.
bool set_cosmology_key(Cosmology& c, const std::string& key, const std::string& val);

}  // namespace radray
