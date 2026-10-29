// The uniform X-ray background of xray_heating = 1 (xray_background.hpp,
// docs/xray_heating.md).
#include "radray/xray_background.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "radray/constants.hpp"
#include "radray/halo_sources.hpp"
#include "radray/halos.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/parallel.hpp"
#include "radray/simulation.hpp"

namespace radray {

namespace {

using namespace pc;
using namespace cp;

constexpr double keV2erg = 1e3 * eV_cgs;
// ln M (Msun/h) range of the sigma table; the step is dndlnM's central
// difference (+-0.02 in ln M), so the derivative uses the neighbouring nodes
constexpr double kLnMlo = 13.8155105579642741;  // ln 1e6
constexpr double kLnMhi = 36.8413614879047;     // ln 1e16
constexpr double kDlnM = 0.02;
constexpr double kDlzTab = 0.01;  // ln(1+z) step of the rho_SFR table

// E(z)^2 = H^2/H0^2 of calc_sim_var (cosmology.cpp)
double hubble_e2(const Cosmology& co, double z) {
  const double a = 1 / (1 + z);
  return co.omegam / (a * a * a) + co.omegal + co.omegar / (a * a * a * a);
}
double hubble_cgs(const Cosmology& co, double z) { return H0_cgs * co.hubble0 * std::sqrt(hubble_e2(co, z)); }

double num(const std::string& key, const std::string& v) {
  try {
    std::size_t pos = 0;
    const double x = std::stod(v, &pos);
    if (v.find_first_not_of(" \t", pos) != std::string::npos) throw std::invalid_argument(v);
    return x;
  } catch (const std::exception&) {
    throw std::runtime_error(key + ": not a number: '" + v + "'");
  }
}

// int_0^1 e^{-x t} dt and int_0^1 t e^{-x t} dt
void exp_weights(double x, double& i0, double& i1) {
  if (x < 1e-3) {
    i0 = 1 - x / 2 + x * x / 6;
    i1 = 0.5 - x / 3 + x * x / 8;
  } else {
    const double e = std::exp(-x);
    i0 = (1 - e) / x;
    i1 = (1 - (1 + x) * e) / (x * x);
  }
}

}  // namespace

bool set_xray_param(XrayParams& p, const std::string& key, const std::string& val) {
  if (key == "xray_heating") {
    const double v = num(key, val);
    if (v != 0 && v != 1) throw std::runtime_error("xray_heating must be 0 or 1");
    p.on = v == 1;
  } else if (key == "xray_log_lx_sfr") {
    p.log_lx_sfr = num(key, val);
  } else if (key == "xray_e0") {
    p.e0_keV = num(key, val);
  } else if (key == "xray_eband") {
    p.eband_keV = num(key, val);
  } else if (key == "xray_emax") {
    p.emax_keV = num(key, val);
  } else if (key == "xray_alpha") {
    p.alpha = num(key, val);
  } else if (key == "xray_tvir_min") {
    p.tvir_min = num(key, val);
  } else if (key == "xray_m_min") {
    p.m_min = num(key, val);
  } else if (key == "xray_zmax") {
    p.z_max = num(key, val);
  } else if (key == "xray_dlog") {
    p.dlog = num(key, val);
    if (!(p.dlog > 0 && p.dlog <= 0.1)) throw std::runtime_error("xray_dlog must be in (0, 0.1]");
  } else if (key == "xray_st") {
    double v[3];
    std::stringstream ss(val);
    std::string tok;
    int n = 0;
    while (std::getline(ss, tok, ',')) {
      if (n >= 3) throw std::runtime_error("xray_st must be A,a,p");
      v[n++] = num(key, tok);
    }
    if (n != 3) throw std::runtime_error("xray_st must be A,a,p");
    p.st.A = v[0];
    p.st.a = v[1];
    p.st.p = v[2];
  } else if (key == "xray_opacity") {
    if (val == "box")
      p.opacity = XrayParams::Opacity::Box;
    else if (val == "neutral")
      p.opacity = XrayParams::Opacity::Neutral;
    else if (val == "none")
      p.opacity = XrayParams::Opacity::None;
    else
      throw std::runtime_error("xray_opacity must be box, neutral or none");
  } else {
    return false;
  }
  return true;
}

XrayBackground::XrayBackground(const XrayParams& p, const SourceModelParams& src, const Cosmology& co)
    : p_(p), co_(co), am_(src.am), l_per_sfr_(src.L_per_sfr), mf_(MassFunctionCosmology::from(co)) {
  if (src.model != SfrModel::AbundanceMatched)
    throw std::runtime_error(
        "xray_heating = 1 needs sfr_model = abundance (the X-ray emissivity follows its SFR)");
  if (!(p.e0_keV > 0 && p.e0_keV < p.eband_keV && p.eband_keV <= p.emax_keV))
    throw std::runtime_error("xray_* energies: need 0 < xray_e0 < xray_eband <= xray_emax");
  tvir_min_ = p.tvir_min >= 0 ? p.tvir_min : src.Tvir_min;
  m_min_ = p.m_min >= 0 ? p.m_min : src.M_min;
  mf_.st = p.st;

  // sigma0 on the ln M grid and dln(1/sigma)/dlnM by central differences
  // of the neighbouring nodes (+-kDlnM, as dndlnM)
  const int nM = static_cast<int>(std::lround((kLnMhi - kLnMlo) / kDlnM)) + 1;
  lnM0_ = kLnMlo;
  dlnM_ = kDlnM;
  std::vector<double> s(static_cast<std::size_t>(nM) + 2);
  for (int i = 0; i < nM + 2; ++i) s[i] = mf_.sigma0(std::exp(lnM0_ + (i - 1) * dlnM_));
  sig0_.resize(nM);
  dsig_.resize(nM);
  for (int i = 0; i < nM; ++i) {
    sig0_[i] = s[i + 1];
    dsig_[i] = -(std::log(s[i + 2]) - std::log(s[i])) / (2 * dlnM_);
  }

  // rho_SFR(z) on the ln(1+z) grid
  dlz_tab_ = kDlzTab;
  const int nz = static_cast<int>(std::ceil(std::log(1 + p_.z_max) / dlz_tab_)) + 1;
  lsfrd_.resize(nz);
  for (int k = 0; k < nz; ++k) {
    const double z = std::exp(k * dlz_tab_) - 1;
    lsfrd_[k] = std::log(std::max(sfrd_direct(z), 1e-300));
  }
}

double XrayBackground::m_min_msun(double z) const {
  // T = mu mH G M / (2 k R), M = 4 pi/3 Delta rho_crit(z) R^3 (halos.cpp)
  const double e2 = hubble_e2(co_, z);
  const double H = H0_cgs * co_.hubble0 * std::sqrt(e2);
  const double rho_crit = 3 * H * H / (8 * pi * G_cgs);
  const double om_z = co_.omegam * std::pow(1 + z, 3) / e2;
  const double dvir = delta_vir_bryan_norman(om_z);
  const double Mg = std::pow(2 * k_cgs * tvir_min_ / (p_.mu * mH_cgs * G_cgs), 1.5) /
                    std::sqrt(4 * pi * dvir * rho_crit / 3);
  return std::max(Mg / Msun_cgs, m_min_);
}

double XrayBackground::sfrd_direct(double z) const {
  const double h = co_.hubble0;
  const double D = mf_.growth(z);
  const int nM = static_cast<int>(sig0_.size());
  auto f = [&](int i) {  // dn/dlnM SFR per Mpc^3 at node i
    const double M = std::exp(lnM0_ + i * dlnM_);  // Msun/h
    const double dn = mf_.dndlnM_sigma(M, sig0_[i] * D, dsig_[i], z,
                                       MassFunction::Kind::ShethTormen);
    return dn * h * h * h * am_.luv(M / h, z) / l_per_sfr_;
  };
  // trapezoid from ln(M_min h) (linear in the first, partial interval)
  const double lmin = std::log(m_min_msun(z) * h);
  double x = (lmin - lnM0_) / dlnM_;
  if (x < 0) x = 0;
  const int i0 = static_cast<int>(std::floor(x));
  if (i0 >= nM - 1) return 0;
  const double w = x - i0;  // fraction of interval i0 below M_min
  const double fa = f(i0), fb = f(i0 + 1);
  const double fm = fa + w * (fb - fa);
  double sum = 0.5 * (fm + fb) * (1 - w) * dlnM_;
  double prev = fb;
  for (int i = i0 + 2; i < nM; ++i) {
    const double fi = f(i);
    sum += 0.5 * (prev + fi) * dlnM_;
    prev = fi;
  }
  return sum;
}

double XrayBackground::sfrd(double z) const {
  const double s = std::log(1 + z) / dlz_tab_;
  if (s < 0) return std::exp(lsfrd_.front());
  const int k = static_cast<int>(std::floor(s));
  if (k >= static_cast<int>(lsfrd_.size()) - 1) return z > p_.z_max ? 0 : std::exp(lsfrd_.back());
  const double w = s - k;
  return std::exp(lsfrd_[k] + w * (lsfrd_[k + 1] - lsfrd_[k]));
}

double XrayBackground::emissivity(double z) const {
  return std::pow(10.0, p_.log_lx_sfr) * sfrd(z);
}

double XrayBackground::shape(double E) const {
  if (E < p_.e0_keV || E > p_.emax_keV) return 0;
  const double a = p_.alpha;
  const double norm = std::fabs(a - 1) < 1e-12
                          ? std::log(p_.eband_keV / p_.e0_keV)
                          : (std::pow(p_.eband_keV, 1 - a) - std::pow(p_.e0_keV, 1 - a)) / (1 - a);
  return std::pow(E, -a) / norm;
}

double XrayBackground::shape_avg(double E1, double E2) const {
  const double a = std::max(E1, p_.e0_keV), b = std::min(E2, p_.emax_keV);
  if (!(b > a)) return 0;
  // (1/ln(E2/E1)) int_a^b shape(E) dE/E, shape = E^-alpha/norm
  const double al = p_.alpha;
  const double I = std::fabs(al) < 1e-12 ? std::log(b / a) : (std::pow(a, -al) - std::pow(b, -al)) / al;
  return I / std::log(E2 / E1) * shape(p_.e0_keV) * std::pow(p_.e0_keV, al);
}

void XrayBackground::spectrum(double z0, const double xfrac[3], std::vector<double>& E_keV,
                              std::vector<double>& J) const {
  const double d = p_.dlog;
  const double s0 = std::log(1 + z0), smax = std::log(1 + p_.z_max);
  E_keV.clear();
  J.clear();
  if (!(smax > s0)) return;
  // z_k: ln(1+z_k) = s0 + k d, k = 0..K; E_m = emax e^{(m - M) d}, m = 0..M
  // observed at z0, and up to M + K for the energies at earlier z_k (a
  // photon observed at E_i had E_{i+k} at z_k)
  const int K = static_cast<int>(std::ceil((smax - s0) / d));
  const double Elo = p_.e0_keV * (1 + z0) / (1 + p_.z_max);
  const int M = static_cast<int>(std::ceil(std::log(p_.emax_keV / Elo) / d));
  const int nE = M + K + 1;
  std::vector<double> E(nE), sig[3], shp(nE);
  for (int j = 0; j < 3; ++j) sig[j].resize(nE);
  for (int m = 0; m < nE; m += 5) {
    double nu[5];
    for (int q = 0; q < 5; ++q) nu[q] = 1e3 * p_.emax_keV * std::exp((m + q - M) * d);  // eV
    double sg[3][5];
    rrk::photoion_sigma(nu, sg);
    for (int q = 0; q < 5 && m + q < nE; ++q) {
      E[m + q] = nu[q] * 1e-3;
      for (int j = 0; j < 3; ++j) sig[j][m + q] = sg[j][q];
    }
  }
  // The spectral shape averaged over each bin in ln E (+-d/2), so that the
  // emission edges at e0 and emax enter with the fraction of the bin they
  // cover (point values would make an O(d) error at the edges)
  for (int m = 0; m < nE; ++m) shp[m] = shape_avg(E[m] * std::exp(-0.5 * d), E[m] * std::exp(0.5 * d));

  // Per redshift node: the emission factor c/(4 pi) (1+z0)^3 eps_c(z)/H(z)
  // per keV and cm^3 (comoving), and the absorption factor of species j,
  // dtau/dln(1+z) = c/H n_j,0 (1+z)^3 x_j sigma_j (dl = c dln(1+z)/H)
  const double nH0 = rhoc_cgs * co_.omegab * co_.hubble0 * co_.hubble0 * XHy / mH_cgs;
  const double nHe0 = rhoc_cgs * co_.omegab * co_.hubble0 * co_.hubble0 * YHe / mHe_cgs;
  const double n0[3] = {nH0 * xfrac[0], nHe0 * xfrac[1], nHe0 * xfrac[2]};
  const double Mpc3 = Mpc2cm * Mpc2cm * Mpc2cm;
  std::vector<double> gfac(K + 1), afac(K + 1);
  for (int k = 0; k <= K; ++k) {
    const double z = std::exp(s0 + k * d) - 1;
    const double H = hubble_cgs(co_, z);
    gfac[k] = z > p_.z_max ? 0 : c_cgs / (4 * pi) * std::pow(1 + z0, 3) * emissivity(z) / Mpc3 / H;
    afac[k] = c_cgs * (1 + z) * (1 + z) * (1 + z) / H;
  }
  auto kappa = [&](int k, int m) {
    return afac[k] * ((n0[0] * sig[0][m] + n0[1] * sig[1][m]) + n0[2] * sig[2][m]);
  };

  E_keV.resize(M + 1);
  J.resize(M + 1);
  for (int i = 0; i <= M; ++i) {
    double tau = 0, sum = 0;
    double ka = kappa(0, i);
    double ga = gfac[0] * shp[i];
    for (int k = 0; k < K; ++k) {
      const double kb = kappa(k + 1, i + k + 1);
      const double gb = gfac[k + 1] * shp[i + k + 1];
      const double dtau = 0.5 * (ka + kb) * d;
      double i0, i1;
      exp_weights(dtau, i0, i1);
      sum += std::exp(-tau) * (ga * (i0 - i1) + gb * i1) * d;
      tau += dtau;
      if (tau > 700) break;
      ka = kb;
      ga = gb;
    }
    E_keV[i] = E[i];
    J[i] = sum;
  }
}

XrayRates XrayBackground::rates(double z0, const double xfrac_box[3]) const {
  XrayRates r;
  r.z = z0;
  double xf[3] = {xfrac_box[0], xfrac_box[1], xfrac_box[2]};
  if (p_.opacity == XrayParams::Opacity::Neutral) {
    xf[0] = 1;
    xf[1] = 1;
    xf[2] = 0;
  } else if (p_.opacity == XrayParams::Opacity::None) {
    xf[0] = xf[1] = xf[2] = 0;
  }
  for (int j = 0; j < 3; ++j) r.xfrac[j] = xf[j];
  r.sfrd = sfrd(z0);
  r.eps = emissivity(z0);
  std::vector<double> E, J;
  spectrum(z0, xf, E, J);
  const int n = static_cast<int>(E.size());
  if (n < 2) return r;
  const double Ith[3] = {HI_eV * eV_cgs, HeI_eV * eV_cgs, HeII_eV * eV_cgs};
  // 4 pi int J/E sigma dE on the log grid: dE = E dlnE, trapezoid in ln E
  for (int i = 0; i < n; ++i) {
    double nu[5] = {1e3 * E[i], 1e3 * E[i], 1e3 * E[i], 1e3 * E[i], 1e3 * E[i]};
    double sg[3][5];
    rrk::photoion_sigma(nu, sg);
    const double w = (i == 0 || i == n - 1 ? 0.5 : 1.0) * 4 * pi * p_.dlog / keV2erg * J[i];
    for (int j = 0; j < 3; ++j) {
      r.gamma[j] += w * sg[j][0];
      r.phi[j] += w * sg[j][0] * (E[i] * keV2erg - Ith[j]);
    }
  }
  return r;
}

void box_ion_fractions(const Simulation& s, double x[3]) {
  const GasGrid& g = s.gas;
  const int N = g.N;
  const std::int64_t plane = static_cast<std::int64_t>(N) * N;
  std::vector<double> sums(static_cast<std::size_t>(N) * 4, 0.0);
#pragma omp parallel for schedule(static) if (par::worth(g.ncell))
  for (int k = 0; k < N; ++k) {
    double a[4] = {0, 0, 0, 0};
    for (std::int64_t c = k * plane; c < (k + 1) * plane; ++c) {
      a[0] += static_cast<double>(g.u[c * 5]);
      for (int q = 0; q < 3; ++q) a[q + 1] += static_cast<double>(g.w[c * 4 + q]);
    }
    for (int q = 0; q < 4; ++q) sums[static_cast<std::size_t>(k) * 4 + q] = a[q];
  }
  double t[4] = {0, 0, 0, 0};
  for (int k = 0; k < N; ++k)
    for (int q = 0; q < 4; ++q) t[q] += sums[static_cast<std::size_t>(k) * 4 + q];
  if (!(t[0] > 0)) {
    x[0] = x[1] = 1;
    x[2] = 0;
    return;
  }
  x[0] = t[1] / (XHy * t[0]);
  x[1] = t[2] / (YHe * t[0]);
  x[2] = t[3] / (YHe * t[0]);
}

}  // namespace radray
