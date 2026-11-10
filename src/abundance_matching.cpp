// Abundance matching to observed UV luminosity functions
// (see include/radray/abundance_matching.hpp, docs/abundance_matching.md).
//
// The offline half of sfr_model = abundance: apps/am_fit reads the LF table
// (data/uv_lf/lf_params.txt, read_lf_table), matches n_halo(>M) = n_gal(<M_UV)
// at every entry's redshift (abundance_match), fits the smooth L_UV(M, z) of
// LuvModel to the matched points (fit_luv_model, fit_highz_norm) and writes
// its parameters as params.txt lines (LuvModel::to_params: am_l0 .. am_b,
// am_z0, am_zlo, am_zhi, am_norm_z, am_norm_dl). At run time only LuvModel is
// used: halo_sources.cpp assign_sources evaluates luv() per halo (the am_*
// keys are parsed by set_halo_param) and xray_background.cpp integrates it
// over the mass function. New code, not a port of anything in the Fortran.
//
// Units: halo masses in Msun (the mass function works in Msun/h and
// (Mpc/h)^-3, converted with h = cosmo.hubble0), number densities in
// comoving Mpc^-3 (not h^3 Mpc^-3), L_UV in erg/s/Hz, magnitudes AB.
#include "radray/abundance_matching.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace radray {

namespace {
constexpr double kLn10 = 2.302585092994045684;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kPc = 3.0856775814913673e18;  // cm
constexpr double kM_lo = -30, kM_hi = -5;      // LF integration range
// step of that grid (12500 points): the trapezoid of CumulativeLF and the
// upper bound of the step of LFEntry::n_brighter
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
// The mass function is evaluated in Msun/h and (Mpc/h)^-3 and converted
// with h; with an enabled SimHMF the analytic dn/dlnM is multiplied by its
// FoF correction R(N, z) (docs/abundance_matching.md section 1). The grid
// runs from 10^lo to 10^16.5 Msun in 0.02 dex and is integrated from the
// top down, so logn[n] = -300 stands for n(>10^16.5 Msun) = 0.
struct CumulativeHMF {
  std::vector<double> logM, logn;
  CumulativeHMF(const MassFunction& mf, double z, MassFunction::Kind k, double lo, double h,
                const SimHMF* sim = nullptr) {
    const double hi = 16.5, d = 0.02;
    const int n = static_cast<int>(std::ceil((hi - lo) / d));
    logM.resize(n + 1);
    logn.resize(n + 1);
    std::vector<double> f(n + 1);
    for (int i = 0; i <= n; ++i) {
      logM[i] = lo + i * d;
      f[i] = mf.dndlnM(std::pow(10.0, logM[i]) * h, z, k);  // (Mpc/h)^-3
      if (sim && sim->enabled) f[i] *= sim->factor(std::pow(10.0, logM[i]), z);
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

// M_UV <-> L_nu (erg/s/Hz): the AB definition at 10 pc (header)
double luv_from_muv(double M_uv) {
  const double d = 10 * kPc;
  return 4 * kPi * d * d * std::pow(10.0, -0.4 * (M_uv + 48.60));
}

double muv_from_luv(double L_uv) {
  const double d = 10 * kPc;
  return -2.5 * std::log10(L_uv / (4 * kPi * d * d)) - 48.60;
}

// phi(M_UV) in Mpc^-3 mag^-1 for the three forms of the table (header); the
// SchechterCurved factor applies faintward of M_UV = -16 only, with the
// curvature delta stored in beta.
double LFEntry::phi(double M) const {
  const double dM = M - Mstar;
  if (form == Form::Schechter || form == Form::SchechterCurved) {
    const double x = std::pow(10.0, -0.4 * dM);
    double f = 0.4 * kLn10 * phistar * std::pow(x, alpha + 1) * std::exp(-x);
    if (form == Form::SchechterCurved && M > -16)
      f *= std::pow(10.0, -0.4 * beta * (M + 16) * (M + 16));
    return f;
  }
  return phistar / (std::pow(10.0, 0.4 * (alpha + 1) * dM) + std::pow(10.0, 0.4 * (beta + 1) * dM));
}

// n(<M_UV) by the trapezoid rule over [kM_lo, M_uv], with the step shrunk
// from kdM so that whole intervals cover the range.
double LFEntry::n_brighter(double M_uv) const {
  const int n = std::max(1, static_cast<int>(std::ceil((M_uv - kM_lo) / kdM)));
  const double h = (M_uv - kM_lo) / n;
  double s = 0;
  for (int i = 0; i <= n; ++i) s += (i == 0 || i == n ? 0.5 : 1.0) * phi(kM_lo + i * h);
  return s * h;
}

// One LFEntry per non-comment line (the column order is in the header); a
// parse error names the file and the line.
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
    else if (form == "schechter_curv")
      e.form = LFEntry::Form::SchechterCurved;
    else
      throw std::runtime_error(path + ":" + std::to_string(ln) +
                               ": form must be schechter|dpl|schechter_curv");
    if (status != "confirmed" && status != "unconfirmed")
      throw std::runtime_error(path + ":" + std::to_string(ln) +
                               ": status must be confirmed|unconfirmed");
    e.confirmed = status == "confirmed";
    out.push_back(e);
  }
  return out;
}

// The FoF correction of the header: R = r_inf / (1 + (N/N_c)^-p) with
// N = M / m_part and N_c = Nc0 ((1+z)/7)^q; zero below nmin particles.
double SimHMF::factor(double M, double z) const {
  if (!enabled) return 1;
  const double N = M / mpart_msun;
  if (N < nmin) return 0;
  const double Nc = Nc0 * std::pow((1 + z) / 7.0, q);
  return r_inf / (1 + std::pow(N / Nc, -p));
}

// key = value lines of fit_fof_hmf.py's output ('#' starts a comment): keys
// mpart_msun (required), fof_r_inf, fof_p, fof_Nc0, fof_q, nmin.
SimHMF SimHMF::read(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  SimHMF s;
  s.enabled = true;
  std::string line;
  while (std::getline(in, line)) {
    const auto c = line.find('#');
    if (c != std::string::npos) line.resize(c);
    const auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    std::istringstream ks(line.substr(0, eq));
    std::string k;
    ks >> k;
    const double v = std::stod(line.substr(eq + 1));
    if (k == "mpart_msun") s.mpart_msun = v;
    else if (k == "fof_r_inf") s.r_inf = v;
    else if (k == "fof_p") s.p = v;
    else if (k == "fof_Nc0") s.Nc0 = v;
    else if (k == "fof_q") s.q = v;
    else if (k == "nmin") s.nmin = v;
    else throw std::runtime_error(path + ": unknown key " + k);
  }
  if (!(s.mpart_msun > 0)) throw std::runtime_error(path + ": mpart_msun missing");
  return s;
}

// Meurer et al. (1999) IRX-beta: A_UV = irx_a + irx_b beta_UV, with beta_UV
// linear in M_obs around -19.5 and no negative extinction.
double DustModel::A_uv(double M_obs) const {
  if (!enabled) return 0;
  const double beta = beta0 + dbeta_dM * (M_obs + 19.5);
  return std::max(0.0, irx_a + irx_b * beta);
}

// n(>M) is the bottom node of a CumulativeHMF starting at M (no SimHMF
// correction here).
double halo_n_above(const MassFunction& mf, double h, double M_msun, double z,
                    MassFunction::Kind k) {
  CumulativeHMF c(mf, z, k, std::log10(M_msun), h);
  return std::pow(10.0, c.logn.front());
}

// One LF entry: for each halo mass of the opt grid, n_halo(>M) from the
// cumulative mass function (whose grid starts 0.1 dex below logM_min, so the
// interpolation never clamps at the bottom) and the observed magnitude with
// the same cumulative density, M_obs = n_gal^-1(n_halo). Masses outside the
// LF table's density range (M_at gives NaN) produce no point. The intrinsic
// magnitude removes the dust; in_range / in_shape mark the points the high-z
// normalisation and the shape fit use (header).
std::vector<AMPoint> abundance_match(const LFEntry& e, const AMOptions& opt) {
  MassFunction mf(opt.cosmo);
  mf.tinker_zmax_evol = opt.tinker_zmax_evol;
  mf.st = opt.st;
  const double h = opt.cosmo.hubble0;
  const CumulativeHMF hmf(mf, e.z, opt.hmf, opt.logM_min - 0.1, h, &opt.sim_hmf);
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
    p.in_range = p.M_obs >= e.M_bright && p.M_obs <= e.M_faint;
    p.in_shape = p.M_obs >= e.M_bright && p.M_obs <= faint;
    pts.push_back(p);
  }
  return pts;
}

// ---------------------------------------------------------------------------
// Parametric model
// ---------------------------------------------------------------------------

// The high-z normalisation offset of the header: 0 up to z_hi, then piecewise
// linear through (z_hi, 0) and the knots, held at the last knot beyond it.
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

// The double power law of the header in log10, with L_p, a and M_p evaluated
// at the redshift clamped to [z_lo, z_hi].
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

// params.txt keys of p[0..7] (parsed by set_halo_param, halo_sources.cpp)
const char* LuvModel::name(int i) {
  static const char* n[kNpar] = {"am_l0", "am_l1", "am_l2", "am_a0",
                                 "am_a1", "am_m0", "am_m1", "am_b"};
  return n[i];
}

// params.txt lines; am_zlo / am_zhi only when they clamp, the knots only when
// set, so a default model round-trips to the same text
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
  // the points: in_shape, positive weight, at or below the shape's z_hi (the
  // entries above it only set the normalisation, fit_highz_norm)
  std::vector<const AMPoint*> use;
  for (const auto& p : pts)
    if (p.in_shape && p.weight > 0 && p.z <= start.z_hi) use.push_back(&p);
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

  // the weighted mean square of the log L residuals and the monotonicity
  // penalty; the fit minimises their sum (cost)
  auto cost_parts = [&](const std::vector<double>& v) {
    LuvModel m = start;
    for (int k = 0; k < nd; ++k) m.p[idx[k]] = v[k];
    double s = 0;
    for (const auto* p : use) {
      const double d = m.log_luv_shape(p->logM, p->z) - p->logL;
      s += p->weight * d * d;
    }
    // L(M) must increase with M at every redshift (abundance matching is
    // monotone): penalise b < 0 and a(z) < 0 over the clamped z range.
    double pen = std::min(0.0, m.p[7]) * std::min(0.0, m.p[7]);
    for (double z : {m.z_lo, m.z_hi}) {
      const double a = m.p[3] + m.p[4] * (std::clamp(z, m.z_lo, m.z_hi) - m.z0);
      pen += std::min(0.0, a) * std::min(0.0, a);
    }
    return std::make_pair(s / wsum, pen);
  };
  auto cost = [&](const std::vector<double>& v) {
    const auto [msq, pen] = cost_parts(v);
    return msq + 1e3 * pen;
  };

  std::vector<double> x0(nd);
  for (int k = 0; k < nd; ++k) x0[k] = start.p[idx[k]];

  // Nelder-Mead with restarts from the current best
  // Nelder-Mead: 4 restarts of up to 20000 iterations each, the initial
  // simplex displaced by 0.2 in every free parameter (0.5 dex in m0, the pivot
  // mass); the standard reflection (t = -1), expansion (-2), outside / inside
  // contraction (-0.5 / +0.5) and shrink towards the best vertex.
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
  // the weighted rms of the residuals, without the monotonicity penalty
  r.rms_dex = std::sqrt(cost_parts(x0).first);
  r.iterations = iters;
  return r;
}

// One knot per distinct redshift above z_hi among the in_range points: the
// weighted mean of logL_AM - log_luv_shape at that redshift (header).
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

// phi(M_UV) = dn/dM_UV as a central difference of the cumulative halo
// density over +-0.01 mag, with M_halo(M_UV) the inversion of L(M) after
// the dust is re-applied; the halo grid starts at 10^6 Msun (the lower
// bound of the inversion).
std::vector<double> model_lf(const LuvModel& m, double z, const std::vector<double>& M_uv,
                             const AMOptions& opt) {
  MassFunction mf(opt.cosmo);
  mf.tinker_zmax_evol = opt.tinker_zmax_evol;
  mf.st = opt.st;
  const double h = opt.cosmo.hubble0;
  const CumulativeHMF hmf(mf, z, opt.hmf, 6.0, h, &opt.sim_hmf);
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
