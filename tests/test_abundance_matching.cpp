// Tests of abundance matching and the AbundanceMatched source model
// (docs/abundance_matching.md):
//  * AB magnitude <-> L_nu conversion;
//  * LF integration (d n_brighter / dM = phi) for Schechter and DPL;
//  * abundance_match: n_halo(>M) == n_gal(<M_UV) at every matched point;
//  * model_lf: its integral equals the halo abundance it came from;
//  * fit_luv_model recovers the parameters of synthetic points;
//  * Tinker et al. (2008) mass function: sane relative to Sheth-Tormen;
//  * the source model: L, sfr, the xi_ion/f_esc scaling, the z clamp and
//    the params.txt keys;
//  * the high-z normalisation with the shape fixed (fit_highz_norm).
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "check.hpp"
#include "radray/abundance_matching.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/halo_sources.hpp"
#include "radray/halos.hpp"
#include "radray/mass_function.hpp"

using namespace radray;

namespace {

LFEntry schechter(double z) {
  LFEntry e;
  e.z = z;
  e.form = LFEntry::Form::Schechter;
  e.Mstar = -21.0;
  e.phistar = 1e-4;
  e.alpha = -2.1;
  e.M_bright = -22.5;
  e.M_faint = -16.5;
  return e;
}

void test_magnitudes() {
  // M_UV = 51.595 - 2.5 log10 L (AB zero point 48.60, 10 pc)
  CHECK_CLOSE(std::log10(luv_from_muv(-20.0)), (51.595 + 20.0) / 2.5, 2e-5);
  CHECK_CLOSE(muv_from_luv(luv_from_muv(-17.3)), -17.3, 1e-12);
}

void test_lf_integral() {
  for (int form = 0; form < 2; ++form) {
    LFEntry e = schechter(8);
    if (form == 1) {
      e.form = LFEntry::Form::DoublePowerLaw;
      e.beta = -4.0;
    }
    for (double M : {-22.0, -20.0, -17.0, -14.0}) {
      const double d = 0.01;
      const double dn = (e.n_brighter(M + d) - e.n_brighter(M - d)) / (2 * d);
      CHECK_CLOSE(dn, e.phi(M), 1e-3);
    }
  }
  // Curved Schechter (Bouwens+22): unchanged brightward of -16, suppressed
  // by 10^(-0.4 delta (M+16)^2) faintward
  {
    LFEntry c = schechter(6);
    c.form = LFEntry::Form::SchechterCurved;
    c.beta = 0.3;  // delta
    LFEntry s0 = schechter(6);
    CHECK(c.phi(-18) == s0.phi(-18));
    CHECK_CLOSE(c.phi(-13), s0.phi(-13) * std::pow(10.0, -0.4 * 0.3 * 9), 1e-12);
  }
  // Schechter phi at M* is 0.4 ln10 phi* e^-1
  const LFEntry e = schechter(8);
  CHECK_CLOSE(e.phi(e.Mstar), 0.4 * std::log(10.0) * e.phistar / std::exp(1.0), 1e-12);
}

void test_abundance_match() {
  AMOptions opt;
  opt.logM_min = 9;
  opt.logM_max = 12;
  opt.dlogM = 0.5;
  const LFEntry e = schechter(7);
  const auto pts = abundance_match(e, opt);
  CHECK(pts.size() == 7);
  MassFunction mf(opt.cosmo);
  double prevL = 0;
  for (const auto& p : pts) {
    // matched abundance; the cumulative tables are piecewise linear in log
    CHECK_CLOSE(e.n_brighter(p.M_obs), p.n_halo, 2e-3);
    CHECK_CLOSE(halo_n_above(mf, opt.cosmo.hubble0, std::pow(10.0, p.logM), e.z, opt.hmf),
                p.n_halo, 1e-2);
    CHECK(p.logL > prevL);  // monotone
    prevL = p.logL;
    CHECK(p.M_int == p.M_obs);  // dust off
  }
  // Dust makes the intrinsic magnitudes brighter where A_UV > 0
  opt.dust.enabled = true;
  opt.dust.beta0 = -1.8;
  const auto pd = abundance_match(e, opt);
  CHECK(pd.back().M_int < pd.back().M_obs);
}

void test_model_lf() {
  AMOptions opt;
  LuvModel m;  // defaults: a smooth, monotone relation
  const double z = 8;
  std::vector<double> M;
  for (double x = -26; x <= -15 + 1e-9; x += 0.05) M.push_back(x);
  const auto phi = model_lf(m, z, M, opt);
  // integral of phi from -26 to -15 == n_halo(>M_h(-15)) - n_halo(>M_h(-26))
  double s = 0;
  for (std::size_t i = 1; i < M.size(); ++i) s += 0.5 * (phi[i] + phi[i - 1]) * 0.05;
  auto logM_of = [&](double Muv) {
    const double lL = std::log10(luv_from_muv(Muv));
    double lo = 6, hi = 16.5;
    for (int it = 0; it < 100; ++it) {
      const double mid = 0.5 * (lo + hi);
      (m.log_luv(mid, z) < lL ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
  };
  MassFunction mf(opt.cosmo);
  const double h = opt.cosmo.hubble0;
  const double n15 = halo_n_above(mf, h, std::pow(10.0, logM_of(-15)), z, opt.hmf);
  const double n26 = halo_n_above(mf, h, std::pow(10.0, logM_of(-26)), z, opt.hmf);
  CHECK_CLOSE(s, n15 - n26, 1e-2);
}

void test_fit_recovery() {
  LuvModel truth;
  truth.p = {29.3, 0.05, -0.03, 1.5, 0.02, 11.4, -0.1, 0.4};
  truth.z_lo = 0;
  truth.z_hi = 100;
  std::vector<AMPoint> pts;
  for (double z : {6.0, 7.0, 8.0, 9.0, 10.0, 11.0}) {
    for (double lm = 9; lm <= 12.5 + 1e-9; lm += 0.1) {
      AMPoint p;
      p.z = z;
      p.logM = lm;
      p.logL = truth.log_luv(lm, z);
      p.in_shape = true;
      pts.push_back(p);
    }
  }
  LuvModel start = truth;
  start.p = {29.0, 0.0, 0.0, 1.3, 0.0, 11.6, 0.0, 0.3};
  std::array<bool, LuvModel::kNpar> free;
  free.fill(true);
  const FitResult r = fit_luv_model(pts, start, free);
  CHECK(r.rms_dex < 1e-4);
  for (int i = 0; i < LuvModel::kNpar; ++i) CHECK(std::fabs(r.model.p[i] - truth.p[i]) < 5e-3);
  // Fixed parameters stay fixed
  free.fill(false);
  free[0] = true;
  const FitResult r2 = fit_luv_model(pts, start, free);
  CHECK(r2.model.p[3] == start.p[3]);
}

void test_tinker() {
  MassFunction mf;
  using K = MassFunction::Kind;
  // Tinker (Delta=200m) and Sheth-Tormen agree to ~10-30% at cluster masses at z=0
  const double t = mf.dndlnM(1e14, 0, K::Tinker08), s = mf.dndlnM(1e14, 0, K::ShethTormen);
  CHECK(t > 0 && s > 0);
  CHECK(std::fabs(std::log(t / s)) < std::log(1.5));
  // At high z, finite and below ST for rare halos (the known ST excess)
  const double t10 = mf.dndlnM(1e11, 10, K::Tinker08), s10 = mf.dndlnM(1e11, 10, K::ShethTormen);
  CHECK(t10 > 0 && std::isfinite(t10));
  CHECK(t10 < s10);
  // Freezing the evolution changes the result above the freeze redshift only
  MassFunction mf2;
  mf2.tinker_zmax_evol = 3;
  CHECK(mf2.dndlnM(1e11, 2, K::Tinker08) == mf.dndlnM(1e11, 2, K::Tinker08));
  CHECK(mf2.dndlnM(1e11, 8, K::Tinker08) != mf.dndlnM(1e11, 8, K::Tinker08));
}

void test_highz_norm() {
  LuvModel shape;
  shape.p = {29.3, 0.05, -0.03, 1.5, 0.02, 11.4, -0.1, 0.4};
  shape.z_lo = 6;
  shape.z_hi = 10;
  // No knots: L(M) held at its z_hi form above z_hi
  CHECK(shape.log_luv(10, 13) == shape.log_luv(10, 10));
  // Synthetic JWST points: shape + 0.3 dex at z = 11, -0.2 dex at z = 14,
  // plus an out-of-range point that must be ignored
  std::vector<AMPoint> pts;
  for (double z : {11.0, 14.0})
    for (double lm = 9.5; lm <= 11.5 + 1e-9; lm += 0.25) {
      AMPoint p;
      p.z = z;
      p.logM = lm;
      p.logL = shape.log_luv_shape(lm, z) + (z == 11 ? 0.3 : -0.2);
      p.in_range = true;
      pts.push_back(p);
    }
  AMPoint bad;
  bad.z = 11;
  bad.logM = 8;
  bad.logL = 0;
  bad.in_range = false;
  pts.push_back(bad);
  const LuvModel m = fit_highz_norm(pts, shape);
  CHECK(m.norm_z.size() == 2);
  CHECK(std::fabs(m.norm_dlogL[0] - 0.3) < 1e-12 && std::fabs(m.norm_dlogL[1] + 0.2) < 1e-12);
  // The shape is unchanged; zero offset up to z_hi, linear between knots,
  // held beyond the last
  CHECK(m.p == shape.p);
  CHECK(m.dlogL_norm(9) == 0 && m.dlogL_norm(10) == 0);
  CHECK(std::fabs(m.dlogL_norm(10.5) - 0.15) < 1e-12);
  CHECK(std::fabs(m.dlogL_norm(12.5) - 0.05) < 1e-12);
  CHECK(m.dlogL_norm(20) == m.dlogL_norm(14));
  CHECK(std::fabs(m.log_luv(10, 11) - (shape.log_luv(10, 10) + 0.3)) < 1e-12);
  // params round trip
  HaloSourceParams hp;
  std::istringstream in(m.to_params());
  std::string line;
  while (std::getline(in, line)) {
    const auto eq = line.find(" = ");
    CHECK(set_halo_param(hp, line.substr(0, eq), line.substr(eq + 3)));
  }
  CHECK(hp.model.am.norm_z == m.norm_z);
  for (std::size_t k = 0; k < m.norm_z.size(); ++k)
    CHECK(std::fabs(hp.model.am.norm_dlogL[k] - m.norm_dlogL[k]) < 1e-9);
  CHECK(std::fabs(hp.model.am.log_luv(10.3, 12) - m.log_luv(10.3, 12)) < 1e-8);
}

void test_sim_hmf() {
  SimHMF sh;
  sh.enabled = true;
  sh.mpart_msun = 5e7;
  sh.r_inf = 2.0;
  sh.p = 0.7;
  sh.Nc0 = 200;
  sh.q = 2.5;
  // zero below nmin particles, r_inf/2 at N = N_c, -> r_inf at large N
  CHECK(sh.factor(19 * 5e7, 6) == 0);
  CHECK_CLOSE(sh.factor(200 * 5e7, 6), 1.0, 1e-12);
  CHECK_CLOSE(sh.factor(200 * 5e7 * std::pow(10.0 / 7.0, 2.5), 9), 1.0, 1e-12);
  CHECK(sh.factor(1e16, 6) > 1.99);
  // Matching against the corrected HMF: n(>M) is the integral of ST x factor
  AMOptions a, b;
  a.hmf = b.hmf = MassFunction::Kind::ShethTormen;
  b.sim_hmf = sh;
  a.logM_min = b.logM_min = 9.5;
  a.logM_max = b.logM_max = 10.5;
  a.dlogM = b.dlogM = 0.5;
  const LFEntry e = schechter(7);
  const auto pa = abundance_match(e, a), pb = abundance_match(e, b);
  CHECK(pa.size() == pb.size());
  for (std::size_t i = 0; i < pa.size(); ++i) {
    // n(>M) of the corrected HMF, by direct integration of ST x factor
    MassFunction mf(a.cosmo);
    const double h = a.cosmo.hubble0;
    double n = 0;
    const double d = 0.002;
    for (double lm = pb[i].logM; lm < 16.5; lm += d) {
      const double m1 = lm + 0.5 * d;
      n += mf.dndlnM(std::pow(10.0, m1) * h, e.z, a.hmf) * h * h * h *
           sh.factor(std::pow(10.0, m1), e.z) * d * std::log(10.0);
    }
    CHECK_CLOSE(pb[i].n_halo, n, 2e-2);
  }
}

HaloUnits units_at(double z) {
  const Config cfg = Config::make(64);
  SimClock clk;
  init_variables(Cosmology{}, z, clk);
  Units u;
  calc_sim_var(cfg, clk, u);
  return make_halo_units(cfg, clk, u);
}

void test_source_model() {
  const HaloUnits hu = units_at(8);
  std::vector<Halo> halos(3);
  const double M[3] = {3e11, 2e10, 1e9};
  for (int i = 0; i < 3; ++i) {
    halos[i].id = i;
    halos[i].N = 100 - i;
    halos[i].M = M[i];
    halos[i].T = 1e5;
    halos[i].x[0] = 0.1 * (i + 1);
  }
  halos[2].T = 5e3;  // below Tvir_min: not a source

  HaloSourceParams hp;
  CHECK(set_halo_param(hp, "sfr_model", "abundance"));
  CHECK(set_halo_param(hp, "am_l0", "29.3"));
  CHECK(set_halo_param(hp, "am_a0", "1.5"));
  CHECK(set_halo_param(hp, "am_zhi", "10"));
  CHECK(set_halo_param(hp, "xi_ion", "25.2"));
  CHECK(hp.model.f_esc == cp::fResc);  // default: the original's escape fraction
  CHECK(!set_halo_param(hp, "am_bogus", "1"));
  CHECK(hp.model.model == SfrModel::AbundanceMatched);
  CHECK(hp.model.am.p[0] == 29.3 && hp.model.am.z_hi == 10);

  double xi_pop2 = 0;
  for (double v : cp::pop2_init) xi_pop2 += v;
  const double scale = std::pow(10.0, 25.2) / xi_pop2;  // f_esc == fResc
  CHECK(std::fabs(scale - 1) < 1e-3);

  const SourceCatalog cat = assign_sources(halos, hu, hp.model);
  CHECK(cat.Nsource == 2);
  for (int i = 0; i < 2; ++i) {
    const double Luv = hp.model.am.luv(M[i], 8);
    CHECK(cat.halo[i].L == static_cast<float>(Luv * scale));
    CHECK(cat.halo[i].sfr == static_cast<float>(Luv / hp.model.L_per_sfr));
  }
  CHECK(cat.halo[2].L == 0 && cat.halo[2].sfr == 0);

  // N_ion,esc = f_esc xi_ion L_UV: doubling f_esc doubles L, +0.3 dex in xi
  // gives 10^0.3
  HaloSourceParams h2 = hp;
  h2.model.f_esc = 0.4;
  const SourceCatalog c2 = assign_sources(halos, hu, h2.model);
  CHECK_CLOSE(c2.halo[0].L, 2.0 * cat.halo[0].L, 1e-6);
  CHECK(c2.halo[0].sfr == cat.halo[0].sfr);
  h2 = hp;
  h2.model.log_xi_ion = 25.5;
  const SourceCatalog c3 = assign_sources(halos, hu, h2.model);
  CHECK_CLOSE(c3.halo[0].L, std::pow(10.0, 0.3) * cat.halo[0].L, 1e-6);

  // f_esc(z) = min(f_esc_max, f_esc ((1+z)/(1+f_esc_z0))^beta); beta = 0
  // (default) is the constant f_esc exactly
  CHECK(f_esc_at(hp.model, 8) == hp.model.f_esc);
  h2 = hp;
  CHECK(set_halo_param(h2, "f_esc", "0.1"));
  CHECK(set_halo_param(h2, "f_esc_beta", "2"));
  CHECK(set_halo_param(h2, "f_esc_z0", "6"));
  CHECK_CLOSE(f_esc_at(h2.model, 6), 0.1, 1e-12);
  CHECK_CLOSE(f_esc_at(h2.model, 13), 0.4, 1e-12);
  CHECK(f_esc_at(h2.model, 40) == 1.0);  // capped at f_esc_max = 1
  CHECK(set_halo_param(h2, "f_esc_max", "0.3"));
  CHECK(f_esc_at(h2.model, 13) == 0.3);
  const SourceCatalog c4 = assign_sources(halos, hu, h2.model);  // z = 8
  CHECK_CLOSE(c4.halo[0].L, cat.halo[0].L * 0.1 * std::pow(9.0 / 7.0, 2) / cp::fResc, 1e-6);

  // The relation is held fixed above am_zhi
  CHECK(hp.model.am.luv(1e10, 12) == hp.model.am.luv(1e10, 10));
  CHECK(hp.model.am.luv(1e10, 9) != hp.model.am.luv(1e10, 10));
}

}  // namespace

int main() {
  test_magnitudes();
  test_lf_integral();
  test_abundance_match();
  test_model_lf();
  test_fit_recovery();
  test_tinker();
  test_highz_norm();
  test_sim_hmf();
  test_source_model();
  TEST_MAIN_END();
}
