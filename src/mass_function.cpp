// Analytic mass functions (see include/radray/mass_function.hpp).
#include "radray/mass_function.hpp"

#include <algorithm>
#include <cmath>

namespace radray {

namespace {
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kRhoCrit = 2.775e11;  // Msun/h per (Mpc/h)^3
constexpr double kDeltaC = 1.686;
}  // namespace

MassFunction::MassFunction(const MassFunctionCosmology& c) : c_(c) {
  norm_ = c_.sigma8 / sigma_unnorm(8.0);
}

double MassFunction::rho_mean() const { return kRhoCrit * c_.omegam; }

double MassFunction::sigma_unnorm(double R) const {
  // BBKS T(q), q = k/Gamma (k in h/Mpc), Sugiyama (1995) shape parameter
  const double om = c_.omegam, ob = c_.omegab, h = c_.hubble0;
  const double gamma = om * h * std::exp(-ob * (1 + std::sqrt(2 * h) / om));
  const int n = 4000;
  const double lk0 = std::log(1e-5), lk1 = std::log(1e3);
  const double dl = (lk1 - lk0) / n;
  double s = 0;
  for (int i = 0; i <= n; ++i) {
    const double k = std::exp(lk0 + i * dl);
    const double q = k / gamma;
    const double T = std::log(1 + 2.34 * q) / (2.34 * q) *
                     std::pow(1 + 3.89 * q + std::pow(16.1 * q, 2) + std::pow(5.46 * q, 3) +
                                  std::pow(6.71 * q, 4),
                              -0.25);
    const double x = k * R;
    const double W = x < 1e-4 ? 1 : 3 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
    const double f = std::pow(k, c_.ns) * T * T * W * W * k * k * k;  // k^3 P(k) W^2
    s += (i == 0 || i == n) ? 0.5 * f : f;
  }
  return std::sqrt(s * dl / (2 * kPi * kPi));
}

double MassFunction::sigma0(double M) const {
  const double R = std::cbrt(3 * M / (4 * kPi * rho_mean()));
  return norm_ * sigma_unnorm(R);
}

double MassFunction::growth(double z) const {
  // D(a) proportional to H(a) * integral_0^a da' / (a' H(a'))^3
  const double om = c_.omegam, ol = 1 - om;
  auto D = [&](double a) {
    const int n = 4000;
    double s = 0;
    for (int i = 1; i <= n; ++i) {  // midpoint rule
      const double x = a * (i - 0.5) / n;
      const double E = std::sqrt(om / (x * x * x) + ol);
      s += 1 / std::pow(x * E, 3);
    }
    s *= a / n;
    return std::sqrt(om / (a * a * a) + ol) * s;
  };
  return D(1 / (1 + z)) / D(1.0);
}

double MassFunction::dndlnM(double M, double z, Kind kind) const {
  const double D = growth(z);
  const double s = sigma0(M) * D;
  const double eps = 0.02;
  const double s1 = sigma0(M * std::exp(-eps)), s2 = sigma0(M * std::exp(eps));
  const double dlnsinv_dlnM = -(std::log(s2) - std::log(s1)) / (2 * eps);
  return dndlnM_sigma(M, s, dlnsinv_dlnM, z, kind);
}

double MassFunction::dndlnM_sigma(double M, double s, double dlnsinv_dlnM, double z,
                                  Kind kind) const {
  const double nu = kDeltaC / s;
  double f;
  if (kind == Kind::PressSchechter) {
    f = std::sqrt(2 / kPi) * nu * std::exp(-0.5 * nu * nu);
  } else if (kind == Kind::Tinker08) {
    // f(sigma) = A [(sigma/b)^-a + 1] exp(-c/sigma^2), Delta = 200 (mean);
    // A = A0 (1+z)^-0.14, a = a0 (1+z)^-0.06, b = b0 (1+z)^-alpha,
    // log10 alpha = -(0.75/log10(Delta/75))^1.2. f(sigma) = nu f(nu) here,
    // so it replaces f directly.
    const double zz = std::min(z, tinker_zmax_evol);
    const double A0 = 0.186, a0 = 1.47, b0 = 2.57, c = 1.19, Delta = 200;
    const double alpha = std::pow(10.0, -std::pow(0.75 / std::log10(Delta / 75), 1.2));
    const double A = A0 * std::pow(1 + zz, -0.14);
    const double a = a0 * std::pow(1 + zz, -0.06);
    const double b = b0 * std::pow(1 + zz, -alpha);
    f = A * (std::pow(s / b, -a) + 1) * std::exp(-c / (s * s));
  } else {
    const double A = st.A, a = st.a, p = st.p;
    const double an2 = a * nu * nu;
    f = A * std::sqrt(2 * a / kPi) * (1 + std::pow(an2, -p)) * nu * std::exp(-0.5 * an2);
  }
  return rho_mean() / M * f * dlnsinv_dlnM;
}

double MassFunction::n_between(double M1, double M2, double z, Kind k) const {
  const int n = 32;
  const double l1 = std::log(M1), l2 = std::log(M2), dl = (l2 - l1) / n;
  double s = 0;
  for (int i = 0; i < n; ++i) s += dndlnM(std::exp(l1 + (i + 0.5) * dl), z, k);
  return s * dl;
}

}  // namespace radray
