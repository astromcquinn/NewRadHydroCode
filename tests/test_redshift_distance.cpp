// Unit tests of redshift_per_distance (Config::redshift_per_distance,
// redshift_fscale in cosmology.hpp; docs/production_run.md): the radiation is
// redshifted once per hydro step by fscale = (a - da)/a (the original, flag
// off), or by fscale^s with s = min(nstep_rt*dl_rt/dl_c, 1) the distance it
// actually moved over the light distance (flag on, c_reduced < 1).
//   * rt_path_fraction matches the calc_timestep geometry, is c_reduced to
//     within the nint rounding, and the rt_substep_carry total is exact;
//   * redshift_fscale: flag off and c_reduced = 1 give (a - da)/a exactly;
//     flag on gives fscale^s; s is capped at 1;
//   * one ray (ray_redshifting) and one cell spectrum (moment_redshift, M1):
//     the spectrum is exactly redshift_spectrum with fscale^s (flag on) or
//     fscale (flag off); with no photons in the lowest group the photon
//     number is conserved and the energy falls by that factor (to the
//     group-range clamp of redshift_spectrum, < 1% of ln), so
//     ln(E1/E0)/ln(fscale) = s ~ c_reduced with the flag and 1 without;
//     with c_reduced = 1 the flag changes nothing (bit-identical).
#include <cmath>
#include <cstring>
#include <memory>

#include "check.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/rt.hpp"
#include "radray/rt_moments.hpp"

using namespace radray;

namespace {

// Independent copy of the calc_timestep geometry: comoving light distance
// of the step [a1, a1 + da] and the RT cell size, in cm.
double light_dist(double a1, double da) {
  return 2 * pc::c_cgs / (pc::H0_cgs * cp::hubble0 * std::sqrt(cp::omegam)) *
         (std::sqrt(a1 + da) - std::sqrt(a1));
}
double rt_cell(const Config& cfg) { return cfg.box / cp::hubble0 * pc::Mpc2cm / cfg.Nrt; }

// One hydro step of the RT clock: calc_timestep (fixed da) from the current
// scalefactor, then the clock at the end of the step, where post_rt
// redshifts.
void step(Simulation& s, double da) {
  calc_timestep(s.cfg, da, s.clk);
  s.clk.scalefactor += s.clk.da_sim;
  s.clk.redshift = 1 / s.clk.scalefactor - 1;
  calc_sim_var(s.cfg, s.clk, s.units);
}

std::unique_ptr<Simulation> make_sim(int N, double z, double c_reduced, bool flag) {
  auto s = std::make_unique<Simulation>();
  init_simulation(*s, N, 4, z);
  s->cfg.c_reduced = c_reduced;
  s->cfg.redshift_per_distance = flag;
  return s;
}

// A da (fraction of a) at z = 6.5, N = 32 for which nint rounds the light
// cells at c_reduced down (s < 1) with a few substeps.
double da_rounding_down(double c_reduced) {
  for (int k = 1; k < 2000; ++k) {
    auto s = make_sim(32, 6.5, c_reduced, true);
    const double da = s->clk.scalefactor * 1e-3 * k;
    step(*s, da);
    if (s->clk.nstep_rt >= 3 && rt_path_fraction(s->cfg, s->clk) < 0.99 * c_reduced) return da;
  }
  return 0;
}

void test_path_fraction() {
  for (double c : {0.1, 0.03, 0.3, 1.0}) {
    auto sp = make_sim(32, 6.5, c, true);
    Simulation& s = *sp;
    const double a1 = s.clk.scalefactor, da = 0.02 * a1;
    step(s, da);
    const double dl_c = light_dist(a1, da), dl_rt = rt_cell(s.cfg);
    const double sp_ = rt_path_fraction(s.cfg, s.clk);
    CHECK(s.clk.nstep_rt == static_cast<int>(std::lround(dl_c / dl_rt * c)));
    CHECK_CLOSE(sp_, s.clk.nstep_rt * dl_rt / dl_c, 1e-12);
    // c_reduced to within the nint rounding (half a cell)
    CHECK(std::fabs(sp_ - c) <= 0.5 * dl_rt / dl_c * (1 + 1e-9));
  }
  // rt_substep_carry: per step s fluctuates, the total path is c_reduced
  // times the total light distance to within half a cell
  {
    auto sp = make_sim(32, 8.0, 0.1, true);
    Simulation& s = *sp;
    s.cfg.rt_substep_carry = true;
    double path = 0, light = 0;
    for (int n = 0; n < 40; ++n) {
      const double a1 = s.clk.scalefactor, da = 0.0037 * a1;
      step(s, da);
      path += s.clk.nstep_rt * rt_cell(s.cfg);
      light += light_dist(a1, da);
      CHECK_CLOSE(rt_path_fraction(s.cfg, s.clk) * light_dist(a1, da),
                  s.clk.nstep_rt * rt_cell(s.cfg), 1e-12);
    }
    CHECK(std::fabs(path - 0.1 * light) <= 0.5 * rt_cell(s.cfg) * (1 + 1e-9));
  }
}

void test_fscale() {
  const double c = 0.1;
  const double da = da_rounding_down(c);
  CHECK(da > 0);
  for (int flag = 0; flag <= 1; ++flag) {
    for (double cr : {c, 1.0}) {
      auto sp = make_sim(32, 6.5, cr, flag == 1);
      Simulation& s = *sp;
      step(s, da);
      const double fscale = (s.clk.scalefactor - s.clk.da_sim) / s.clk.scalefactor;
      const double f = redshift_fscale(s.cfg, s.clk);
      if (flag == 0 || cr >= 1) {
        CHECK(f == fscale);  // the original, bit for bit
      } else {
        const double sp_ = rt_path_fraction(s.cfg, s.clk);
        CHECK(sp_ < 1);
        CHECK_CLOSE(f, std::pow(fscale, sp_), 1e-15);
        // ln f / ln fscale = s ~ c_reduced: the radiation loses c/c_hat
        // times less energy per hydro step than with the original
        CHECK_CLOSE(std::log(f) / std::log(fscale), sp_, 1e-12);
        CHECK(std::fabs(sp_ - c) < 0.5 * c);
      }
    }
  }
  // Cap: c_reduced close to 1 and a step of less than half a cell of light:
  // the nstep_rt >= 1 guard makes the rays cover more than the light
  // distance (s > 1); the factor stays the per-time one
  {
    auto sp = make_sim(32, 6.5, 0.9, true);
    Simulation& s = *sp;
    const double a1 = s.clk.scalefactor;
    double da = 1e-3 * a1;
    while (light_dist(a1, da) / rt_cell(s.cfg) * 0.9 > 0.4) da *= 0.5;
    step(s, da);
    CHECK(s.clk.nstep_rt == 1);
    CHECK(rt_path_fraction(s.cfg, s.clk) > 1);
    CHECK(redshift_fscale(s.cfg, s.clk) ==
          (s.clk.scalefactor - s.clk.da_sim) / s.clk.scalefactor);
  }
}

// Photon number and energy of a spectrum
void totals(const float* nu, const float* ng, double& n, double& e) {
  n = e = 0;
  for (int g = 0; g < kNfreq; ++g) {
    n += ng[g];
    e += static_cast<double>(nu[g]) * ng[g];
  }
}

const float kNg[kNfreq] = {0.0f, 1.0f, 0.6f, 0.4f, 0.2f};  // none in group 1: no loss

// One ray through ray_redshifting; returns its final nu, ngam
void redshift_one_ray(double cr, bool flag, double da, float* nu, float* ng, double& fs_used,
                      double& fscale, double& ngam_min) {
  auto sp = make_sim(16, 6.5, cr, flag);
  Simulation& s = *sp;
  s.rays.Nrbuf = s.rays.Nfree_max / 4;
  ray_buffer_allocate(s);
  step(s, da);
  std::int64_t ir;
  ray_create(s, 1, &ir);
  Ray& r = s.rays.at(ir);
  r.nside = 4;
  for (int g = 0; g < kNfreq; ++g) {
    r.nu[g] = static_cast<float>(cp::freq_init[g]);
    r.ngam[g] = kNg[g];
  }
  ray_redshifting(s);
  std::memcpy(nu, r.nu, sizeof r.nu);
  std::memcpy(ng, r.ngam, sizeof r.ngam);
  fs_used = redshift_fscale(s.cfg, s.clk);
  fscale = (s.clk.scalefactor - s.clk.da_sim) / s.clk.scalefactor;
  ngam_min = s.units.ngam_min;
}

// One cell spectrum through moment_redshift (M1)
void redshift_one_cell(double cr, bool flag, double da, float* nu, float* ng, double& fs_used,
                       double& fscale, double& ngam_min) {
  auto sp = make_sim(8, 6.5, cr, flag);
  Simulation& s = *sp;
  s.cfg.rt_method = RTMethod::M1;
  moment_init(s);
  step(s, da);
  MomentState& M = *s.moments;
  const std::int64_t c = s.rt.cell(3, 4, 5);
  for (int g = 0; g < kNfreq; ++g) {
    M.n[g * M.ncell + c] = kNg[g];
    M.e[g * M.ncell + c] = static_cast<float>(cp::freq_init[g] * kNg[g]);
  }
  moment_redshift(s);
  for (int g = 0; g < kNfreq; ++g) {
    ng[g] = M.n[g * M.ncell + c];
    nu[g] = ng[g] > 0 ? static_cast<float>(static_cast<double>(M.e[g * M.ncell + c]) / ng[g])
                      : 0.0f;
  }
  fs_used = redshift_fscale(s.cfg, s.clk);
  fscale = (s.clk.scalefactor - s.clk.da_sim) / s.clk.scalefactor;
  ngam_min = s.units.ngam_min;
}

// redshift_spectrum of the initial spectrum with factor f
void expected(double f, double ngam_min, float* nu, float* ng) {
  const RedshiftTable T = make_redshift_table();
  double dnu[kNfreq], dng[kNfreq], keep[kNfreq], down[kNfreq];
  for (int g = 0; g < kNfreq; ++g) {
    dnu[g] = static_cast<float>(cp::freq_init[g]);
    dng[g] = kNg[g];
  }
  redshift_spectrum(T, f, ngam_min, dnu, dng, keep, down);
  for (int g = 0; g < kNfreq; ++g) {
    nu[g] = static_cast<float>(dnu[g]);
    ng[g] = static_cast<float>(dng[g]);
  }
}

template <class F>
void test_spectrum(F redshift_one) {
  double n0, e0;
  float nu0[kNfreq];
  for (int g = 0; g < kNfreq; ++g) nu0[g] = static_cast<float>(cp::freq_init[g]);
  totals(nu0, kNg, n0, e0);
  const double da = da_rounding_down(0.1);

  float nu_off[kNfreq], ng_off[kNfreq], nu_on[kNfreq], ng_on[kNfreq], nu_x[kNfreq], ng_x[kNfreq];
  double fs_off, fs_on, fscale, ngam_min;
  // flag off, c_reduced = 0.1: redshifted by fscale (per time)
  redshift_one(0.1, false, da, nu_off, ng_off, fs_off, fscale, ngam_min);
  double n1, e1;
  totals(nu_off, ng_off, n1, e1);
  CHECK(fs_off == fscale);
  expected(fscale, ngam_min, nu_x, ng_x);
  for (int g = 0; g < kNfreq; ++g) {
    CHECK_CLOSE(ng_off[g], ng_x[g], 1e-6);
    CHECK_CLOSE(nu_off[g], nu_x[g], 1e-6);  // moments: nu = e/n of floats
  }
  CHECK_CLOSE(n1, n0, 1e-6);
  CHECK_CLOSE(std::log(e1 / e0) / std::log(fscale), 1.0, 0.02);

  // flag on: by fscale^s, s ~ 0.1 (per distance)
  redshift_one(0.1, true, da, nu_on, ng_on, fs_on, fscale, ngam_min);
  totals(nu_on, ng_on, n1, e1);
  CHECK(fs_on != fscale);
  const double s_path = std::log(fs_on) / std::log(fscale);
  CHECK(s_path > 0.05 && s_path < 0.1);
  expected(fs_on, ngam_min, nu_x, ng_x);
  for (int g = 0; g < kNfreq; ++g) {
    CHECK_CLOSE(ng_on[g], ng_x[g], 1e-6);
    CHECK_CLOSE(nu_on[g], nu_x[g], 1e-6);
  }
  CHECK_CLOSE(n1, n0, 1e-6);
  // the per-step energy loss is the per-distance one, s ~ 0.1 times the
  // per-time one
  CHECK_CLOSE(std::log(e1 / e0) / std::log(fscale), s_path, 0.02);

  // c_reduced = 1: the flag changes nothing, bit for bit
  redshift_one(1.0, false, da, nu_off, ng_off, fs_off, fscale, ngam_min);
  redshift_one(1.0, true, da, nu_on, ng_on, fs_on, fscale, ngam_min);
  CHECK(fs_on == fs_off);
  CHECK(std::memcmp(nu_on, nu_off, sizeof nu_on) == 0);
  CHECK(std::memcmp(ng_on, ng_off, sizeof ng_on) == 0);
}

}  // namespace

int main() {
  test_path_fraction();
  test_fscale();
  test_spectrum(redshift_one_ray);
  test_spectrum(redshift_one_cell);
  TEST_MAIN_END();
}
