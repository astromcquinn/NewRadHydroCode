#include "radray/cosmo_params.hpp"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace radray {

double cosmology_omegar(double hubble0) {
  // cosmo_parameters.f90: omegar = 4.48E-7*(1+0.69)*Tcmb0**4/hubble0**2,
  // evaluated left to right in real(8) after the real(4) product
  return cp::omegar_c1 * cp::Tcmb0_pow4 / (hubble0 * hubble0);
}

double cosmology_temp_x(double omegab, double hubble0) {
  // cosmo_tools.f90 calc_sim_var: x = 1E3*(ob*h0**2)**0.4 (0.4 a real(4)
  // literal; gfortran folds the power with MPFR, glibc's pow is correctly
  // rounded for the original's value, which the unit test checks)
  return 1e3 * std::pow(omegab * (hubble0 * hubble0), static_cast<double>(0.4f));
}

void Cosmology::derive() {
  if (!(omegam > 0 && omegam <= 1)) throw std::runtime_error("omegam must be in (0, 1]");
  if (!(omegab > 0 && omegab < omegam)) throw std::runtime_error("omegab must be in (0, omegam)");
  if (!(hubble0 > 0 && hubble0 < 2)) throw std::runtime_error("hubble0 must be in (0, 2)");
  if (!(sigma8 > 0)) throw std::runtime_error("sigma8 must be > 0");
  if (!(ns > 0)) throw std::runtime_error("ns must be > 0");
  omegal = 1 - omegam;
  if (is_original()) {
    // the folded constants of the original (equal to the runtime formulas
    // on glibc; kept so that the default never depends on the libm)
    omegar = cp::omegar;
    temp_x = folded::temp_x;
  } else {
    omegar = cosmology_omegar(hubble0);
    temp_x = cosmology_temp_x(omegab, hubble0);
  }
}

bool Cosmology::is_original() const {
  return omegam == cp::omegam && omegab == cp::omegab && hubble0 == cp::hubble0;
}

std::string Cosmology::summary() const {
  char buf[256];
  std::snprintf(buf, sizeof buf, "omegam %.7g omegab %.7g hubble0 %.7g sigma8 %.7g ns %.7g%s",
                omegam, omegab, hubble0, sigma8, ns, is_original() ? " (the original)" : "");
  return buf;
}

bool set_cosmology_key(Cosmology& c, const std::string& key, const std::string& val_in) {
  double* dst = nullptr;
  if (key == "omegam")
    dst = &c.omegam;
  else if (key == "omegab")
    dst = &c.omegab;
  else if (key == "hubble0")
    dst = &c.hubble0;
  else if (key == "sigma8")
    dst = &c.sigma8;
  else if (key == "ns")
    dst = &c.ns;
  else
    return false;
  std::string val = val_in;
  for (auto& ch : val)
    if (ch == 'd' || ch == 'D') ch = 'e';  // Fortran exponents
  // Rounded to real(4) directly from the decimal string (strtof), as
  // gfortran rounds the original's literals: omegam = 0.27 is the default.
  const char* s = val.c_str();
  char* end = nullptr;
  errno = 0;
  const float v = std::strtof(s, &end);
  while (end && (*end == ' ' || *end == '\t')) ++end;
  if (end == s || (end && *end) || errno == ERANGE)
    throw std::runtime_error(key + ": not a number: " + val_in);
  *dst = static_cast<double>(v);
  return true;
}

}  // namespace radray
