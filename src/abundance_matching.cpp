// Abundance matching to observed UV luminosity functions
// (see include/radray/abundance_matching.hpp, docs/abundance_matching.md).
#include "radray/abundance_matching.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace radray {

namespace {
constexpr double kLn10 = 2.302585092994045684;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kPc = 3.0856775814913673e18;  // cm
constexpr double kM_lo = -30, kM_hi = -5;      // LF integration range
constexpr double kdM = 0.002;

// Linear interpolation of y(x) on an increasing grid; clamps outside.
double interp(const std::vector<double>& x, const std::vector<double>& y, double v) {
  if (v <= x.front()) return y.front();
  if (v >= x.back()) return y.back();
  const auto it = std::upper_bound(x.begin(), x.end(), v);
  const std::size_t i = static_cast<std::size_t>(it - x.begin());
  const double t = (v - x[i - 1]) / (x[i] - x[i - 1]);
  return y[i - 1] + t * (y[i] - y[i - 1]);
}

// Cumulative n(<M_UV) of an entry on a fine magnitude grid (log10 n).
struct CumulativeLF {
  std::vector<double> M, logn;
  explicit CumulativeLF(const LFEntry& e) {
    const int n = static_cast<int>(std::lround((kM_hi - kM_lo) / kdM));
    M.resize(n + 1);
    logn.resize(n + 1);
    double s = 0, prev = e.phi(kM_lo);
    for (int i = 0; i <= n; ++i) {
      M[i] = kM_lo + i * kdM;
      const double f = e.phi(M[i]);
      if (i > 0) s += 0.5 * (f + prev) * kdM;
      prev = f;
      logn[i] = std::log10(std::max(s, 1e-300));
    }
  }
  // Magnitude at which n(<M) = n (bisection on the monotone table).
  double M_at(double n) const {
    const double ln = std::log10(n);
    if (ln <= logn.front() || ln > logn.back()) return std::numeric_limits<double>::quiet_NaN();
    const auto it = std::lower_bound(logn.begin(), logn.end(), ln);
    const std::size_t i = static_cast<std::size_t>(it - logn.begin());
    if (i == 0) return M.front();
    const double t = (ln - logn[i - 1]) / (logn[i] - logn[i - 1]);
    return M[i - 1] + t * (M[i] - M[i - 1]);
  }
};

// Cumulative halo density n(>M) on a log grid (Msun), Mpc^-3
struct CumulativeHMF {
  std::vector<double> logM, logn;
  CumulativeHMF(const MassFunction& mf, double z, MassFunction::Kind k, double lo, double h) {
    const double hi = 16.5, d = 0.02;
    const int n = static_cast<int>(std::ceil((hi - lo) / d));
    logM.resize(n + 1);
    logn.resize(n + 1);
    std::vector<double> f(n + 1);
    for (int i = 0; i <= n; ++i) {
      logM[i] = lo + i * d;
      f[i] = mf.dndlnM(std::pow(10.0, logM[i]) * h, z, k);  // (Mpc/h)^-3
    }
    double s = 0;
    logn[n] = -300;
    for (int i = n - 1; i >= 0; --i) {
      s += 0.5 * (f[i] + f[i + 1]) * d * kLn10;
      logn[i] = std::log10(std::max(s * h * h * h, 1e-300));  // Mpc^-3
    }
  }
  double logn_at(double lm) const { return interp(logM, logn, lm); }
};

}  // namespace

// ---------------------------------------------------------------------------

double luv_from_muv(double M_uv) {
  const double d = 10 * kPc;
  return 4 * kPi * d * d * std::pow(10.0, -0.4 * (M_uv + 48.60));
}

double muv_from_luv(double L_uv) {
  const double d = 10 * kPc;
  return -2.5 * std::log10(L_uv / (4 * kPi * d * d)) - 48.60;
}

double LFEntry::phi(double M) const {
  const double dM = M - Mstar;
  if (form == Form::Schechter) {
    const double x = std::pow(10.0, -0.4 * dM);
    return 0.4 * kLn10 * phistar * std::pow(x, alpha + 1) * std::exp(-x);
  }
  return phistar / (std::pow(10.0, 0.4 * (alpha + 1) * dM) + std::pow(10.0, 0.4 * (beta + 1) * dM));
}

double LFEntry::n_brighter(double M_uv) const {
  const int n = std::max(1, static_cast<int>(std::ceil((M_uv - kM_lo) / kdM)));
  const double h = (M_uv - kM_lo) / n;
  double s = 0;
  for (int i = 0; i <= n; ++i) s += (i == 0 || i == n ? 0.5 : 1.0) * phi(kM_lo + i * h);
  return s * h;
}

std::vector<LFEntry> read_lf_table(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<LFEntry> out;
  std::string line;
  int ln = 0;
  while (std::getline(in, line)) {
    ++ln;
    const auto c = line.find('#');
    if (c != std::string::npos) line.resize(c);
    std::istringstream ss(line);
    LFEntry e;
    std::string form, status;
    if (!(ss >> e.z)) continue;  // blank or comment line
    if (!(ss >> form >> e.Mstar >> e.phistar >> e.alpha >> e.beta >> e.M_bright >> e.M_faint >>
          e.weight >> e.ref >> status))
      throw std::runtime_error(path + ":" + std::to_string(ln) + ": expected 11 columns");
    if (form == "schechter")
      e.form = LFEntry::Form::Schechter;
    else if (form == "dpl")
      e.form = LFEntry::Form::DoublePowerLaw;
    else
      throw std::runtime_error(path + ":" + std::to_string(ln) + ": form must be schechter|dpl");
    if (status != "confirmed" && status != "unconfirmed")
      throw std::runtime_error(path + ":" + std::to_string(ln) +
                               ": status must be confirmed|unconfirmed");
    e.confirmed = status == "confirmed";
    out.push_back(e);
  }
  return out;
}

double DustModel::A_uv(double M_obs) const {
  if (!enabled) return 0;
  const double beta = beta0 + dbeta_dM * (M_obs + 19.5);
  return std::max(0.0, irx_a + irx_b * beta);
}

double halo_n_above(const MassFunction& mf, double h, double M_msun, double z,
                    MassFunction::Kind k) {
  CumulativeHMF c(mf, z, k, std::log10(M_msun), h);
  return std::pow(10.0, c.logn.front());
}

std::vector<AMPoint> abundance_match(const LFEntry& e, const AMOptions& opt) {
  MassFunction mf(opt.cosmo);
  mf.tinker_zmax_evol = opt.tinker_zmax_evol;
  const double h = opt.cosmo.hubble0;
  const CumulativeHMF hmf(mf, e.z, opt.hmf, opt.logM_min - 0.1, h);
  const CumulativeLF clf(e);
  const double faint = opt.extend_faint_to != 0 ? std::max(e.M_faint, opt.extend_faint_to)
                                                : e.M_faint;
  std::vector<AMPoint> pts;
  const int n = static_cast<int>(std::lround((opt.logM_max - opt.logM_min) / opt.dlogM));
  for (int i = 0; i <= n; ++i) {
    AMPoint p;
    p.z = e.z;
    p.logM = opt.logM_min + i * opt.dlogM;
    p.n_halo = std::pow(10.0, hmf.logn_at(p.logM));
    p.M_obs = clf.M_at(p.n_halo);
    if (!std::isfinite(p.M_obs)) continue;
    p.M_int = p.M_obs - opt.dust.A_uv(p.M_obs);
    p.logL = std::log10(luv_from_muv(p.M_int));
    p.weight = e.weight;
    p.in_range = p.M_obs >= e.M_bright && p.M_obs <= faint;
    pts.push_back(p);
  }
  return pts;
}

// ---------------------------------------------------------------------------
// Parametric model
// ---------------------------------------------------------------------------

double LuvModel::dlogL_norm(double z) const {
  if (norm_z.empty() || !(z > z_hi)) return 0;
  double za = z_hi, da = 0;
  for (std::size_t k = 0; k < norm_z.size(); ++k) {
    if (z <= norm_z[k]) return da + (norm_dlogL[k] - da) * (z - za) / (norm_z[k] - za);
    za = norm_z[k];
    da = norm_dlogL[k];
  }
  return da;  // held beyond the last knot
}

double LuvModel::log_luv(double logM, double z) const {
  return log_luv_shape(logM, z) + dlogL_norm(z);
}

double LuvModel::log_luv_shape(double logM, double z) const {
  const double x = std::clamp(z, z_lo, z_hi) - z0;
  const double logLp = p[0] + p[1] * x + p[2] * x * x;
  const double a = p[3] + p[4] * x;
  const double logMp = p[5] + p[6] * x;
  const double b = p[7];
  const double u = logM - logMp;
  // log10(2 / (10^(-a u) + 10^(-b u))), evaluated without overflow
  const double ea = -a * u, eb = -b * u;
  const double mx = std::max(ea, eb);
  return logLp + std::log10(2.0) - mx -
         std::log10(std::pow(10.0, ea - mx) + std::pow(10.0, eb - mx));
}

double LuvModel::luv(double M, double z) const {
  if (!(M > 0)) return 0;
  return std::pow(10.0, log_luv(std::log10(M), z));
}

const char* LuvModel::name(int i) {
  static const char* n[kNpar] = {"am_l0", "am_l1", "am_l2", "am_a0",
                                 "am_a1", "am_m0", "am_m1", "am_b"};
  return n[i];
}

std::string LuvModel::to_params() const {
  std::ostringstream s;
  s.precision(10);
  for (int i = 0; i < kNpar; ++i) s << name(i) << " = " << p[i] << "\n";
  s << "am_z0 = " << z0 << "\n";
  if (z_lo > 0) s << "am_zlo = " << z_lo << "\n";
  if (z_hi < 1e29) s << "am_zhi = " << z_hi << "\n";
  if (!norm_z.empty()) {
    s << "am_norm_z = ";
    for (std::size_t k = 0; k < norm_z.size(); ++k) s << (k ? "," : "") << norm_z[k];
    s << "\nam_norm_dl = ";
    for (std::size_t k = 0; k < norm_dlogL.size(); ++k) s << (k ? "," : "") << norm_dlogL[k];
    s << "\n";
  }
  return s.str();
}

FitResult fit_luv_model(const std::vector<AMPoint>& pts, const LuvModel& start,
                        const std::array<bool, LuvModel::kNpar>& free) {
  std::vector<const AMPoint*> use;
  for (const auto& p : pts)
    if (p.in_range && p.weight > 0 && p.z <= start.z_hi) use.push_back(&p);
  FitResult r;
  r.model = start;
  r.npoints = static_cast<int>(use.size());
  if (use.empty()) return r;
  std::vector<int> idx;
  for (int i = 0; i < LuvModel::kNpar; ++i)
    if (free[i]) idx.push_back(i);
  const int nd = static_cast<int>(idx.size());
  double wsum = 0;
  for (const auto* p : use) wsum += p->weight;

  auto cost = [&](const std::vector<double>& v) {
    LuvModel m = start;
    for (int k = 0; k < nd; ++k) m.p[idx[k]] = v[k];
    double s = 0;
    for (const auto* p : use) {
      const double d = m.log_luv_shape(p->logM, p->z) - p->logL;
      s += p->weight * d * d;
    }
    return s / wsum;
  };

  std::vector<double> x0(nd);
  for (int k = 0; k < nd; ++k) x0[k] = start.p[idx[k]];

  // Nelder-Mead with restarts from the current best
  int iters = 0;
  for (int restart = 0; restart < 4 && nd > 0; ++restart) {
    std::vector<std::vector<double>> s(nd + 1, x0);
    for (int k = 0; k < nd; ++k) s[k + 1][k] += (idx[k] == 5 ? 0.5 : 0.2);
    std::vector<double> f(nd + 1);
    for (int i = 0; i <= nd; ++i) f[i] = cost(s[i]);
    for (int it = 0; it < 20000; ++it, ++iters) {
      std::vector<int> o(nd + 1);
      for (int i = 0; i <= nd; ++i) o[i] = i;
      std::sort(o.begin(), o.end(), [&](int a, int b) { return f[a] < f[b]; });
      const int best = o[0], worst = o[nd], second = o[nd - 1];
      if (std::fabs(f[worst] - f[best]) < 1e-14 * (1 + f[best])) break;
      std::vector<double> c(nd, 0);
      for (int i = 0; i <= nd; ++i)
        if (i != worst)
          for (int k = 0; k < nd; ++k) c[k] += s[i][k] / nd;
      auto along = [&](double t) {
        std::vector<double> y(nd);
        for (int k = 0; k < nd; ++k) y[k] = c[k] + t * (s[worst][k] - c[k]);
        return y;
      };
      const auto xr = along(-1);
      const double fr = cost(xr);
      if (fr < f[best]) {
        const auto xe = along(-2);
        const double fe = cost(xe);
        if (fe < fr) {
          s[worst] = xe;
          f[worst] = fe;
        } else {
          s[worst] = xr;
          f[worst] = fr;
        }
      } else if (fr < f[second]) {
        s[worst] = xr;
        f[worst] = fr;
      } else {
        const auto xc = along(fr < f[worst] ? -0.5 : 0.5);
        const double fc = cost(xc);
        if (fc < std::min(fr, f[worst])) {
          s[worst] = xc;
          f[worst] = fc;
        } else {
          for (int i = 0; i <= nd; ++i) {
            if (i == best) continue;
            for (int k = 0; k < nd; ++k) s[i][k] = s[best][k] + 0.5 * (s[i][k] - s[best][k]);
            f[i] = cost(s[i]);
          }
        }
      }
    }
    int b = 0;
    for (int i = 1; i <= nd; ++i)
      if (f[i] < f[b]) b = i;
    x0 = s[b];
  }
  for (int k = 0; k < nd; ++k) r.model.p[idx[k]] = x0[k];
  r.rms_dex = std::sqrt(cost(x0));
  r.iterations = iters;
  return r;
}

LuvModel fit_highz_norm(const std::vector<AMPoint>& pts, const LuvModel& shape) {
  LuvModel m = shape;
  m.norm_z.clear();
  m.norm_dlogL.clear();
  std::vector<double> zs;
  for (const auto& p : pts)
    if (p.in_range && p.weight > 0 && p.z > shape.z_hi) zs.push_back(p.z);
  std::sort(zs.begin(), zs.end());
  zs.erase(std::unique(zs.begin(), zs.end()), zs.end());
  for (double z : zs) {
    double sw = 0, sd = 0;
    for (const auto& p : pts) {
      if (!p.in_range || !(p.weight > 0) || p.z != z) continue;
      sw += p.weight;
      sd += p.weight * (p.logL - shape.log_luv_shape(p.logM, p.z));
    }
    m.norm_z.push_back(z);
    m.norm_dlogL.push_back(sd / sw);
  }
  return m;
}

std::vector<double> model_lf(const LuvModel& m, double z, const std::vector<double>& M_uv,
                             const AMOptions& opt) {
  MassFunction mf(opt.cosmo);
  mf.tinker_zmax_evol = opt.tinker_zmax_evol;
  const double h = opt.cosmo.hubble0;
  const CumulativeHMF hmf(mf, z, opt.hmf, 6.0, h);
  // Halo mass whose intrinsic L matches an observed magnitude (L monotone
  // in M for a, b > 0): bisection in log M.
  auto logM_of = [&](double M_obs) {
    const double logL = std::log10(luv_from_muv(M_obs - opt.dust.A_uv(M_obs)));
    double lo = 6, hi = 16.5;
    if (m.log_luv(lo, z) >= logL) return lo;
    if (m.log_luv(hi, z) <= logL) return hi;
    for (int it = 0; it < 100; ++it) {
      const double mid = 0.5 * (lo + hi);
      (m.log_luv(mid, z) < logL ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
  };
  std::vector<double> out;
  const double d = 0.01;
  for (double M : M_uv) {
    const double n1 = std::pow(10.0, hmf.logn_at(logM_of(M - d)));  // brighter
    const double n2 = std::pow(10.0, hmf.logn_at(logM_of(M + d)));
    out.push_back((n2 - n1) / (2 * d));
  }
  return out;
}

}  // namespace radray
