// Test of the sub-cycled gravity of the refined levels (pm_level_subcycle,
// docs/stage_g.md section 6): a light particle on a circular orbit around a
// clump of 511 particles at one point (they move as one body), inside the
// finest of two refined levels (the clump at the centre of level 2, the
// orbit 6 level-2 cells from it), in a static universe (a = 1, da = 0) with
// a cold, light gas, through radhydro's base step (two hydro calls around
// the PM step) at a fixed base step of 1/64 of the orbital period, for one
// orbit. The orbiter's energy relative to the clump, E = v^2/2 + Phi(r)
// with Phi from the measured radial force (direction-averaged), is compared
// step by step with a sub-cycled run at half the base step (the energy's
// swing around the orbit from the mesh's angular anisotropy is the same in
// every run and cancels). Without sub-cycling level 2's particles are
// kicked once per base step; with it 2 x 2^2 = 8 times (two hydro calls of
// four level-2 substeps), so the leapfrog's energy error should fall by
// about 8^2 = 64 (the request's (2^l)^2 = 16 counted one hydro call per
// base step). Checks a factor of at least 16 and prints a checksum of the
// sub-cycled energies (the same on any rank count).
//
// usage: test_amrex_pm_subcycle [--mgs n] [--orbits n] [--steps-per-orbit n]
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>

#include "../src/amrex/gas_state.hpp"
#include "check.hpp"
#include "radray/amrex_state.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/hydro.hpp"
#include "radray/hydro_unsplit.hpp"
#include "radray/pm.hpp"
#include "radray/pm_backend.hpp"
#include "radray/rt_backend.hpp"
#include "radray/simulation.hpp"

using namespace radray;

namespace {
namespace ax = amrex;

int g_mgs = 8;
int g_orbits = 1;
int g_spo = 64;            // base steps per orbit (at 16 the cold gas's infall sub-cycles the hydro)
constexpr int kN = 32;     // the gas grid and level-0 mesh
constexpr int kNdm = 8;    // 512 particles: 511 in the clump, 1 orbiter
const double kXc[3] = {4.013, 4.027, 3.991};  // the clump (DM-grid units: level-0 cell 16 = 4.0)
constexpr double kR = 0.375;                  // the orbit radius (DM units: 1.5 level-0 cells, 6 level-2 cells)

struct Run {
  std::vector<double> E;  // the orbiter's energy after every base step
  double E0 = 0, ekin = 0;
  std::int64_t substeps = 0;
};

Run run(bool sub, int spo) {
  Simulation sim;
  Config& c = sim.cfg;
  c = Config::make(kN, 4, kNdm, 4.0);
  c.rt_backend = RTBackend::AMReX;
  c.hydro_backend = HydroBackend::AMReX;
  c.pm_backend = PMBackend::AMReX;
  c.hydro_scheme = HydroScheme::Unsplit;
  c.amrex_state = true;
  c.amrex_max_grid_size = g_mgs;
  c.amrex_hydro_max_grid_size = 0;
  c.amrex_max_level = 2;
  c.amrex_refine_mass = 1e9;
  c.amrex_n_error_buf = 0;
  const int rb[12] = {8, 8, 8, 23, 23, 23, 12, 12, 12, 19, 19, 19};
  for (int m = 0; m < 12; ++m) c.amrex_refine_box[m] = rb[m];
  c.amrex_refine_nbox = 2;
  c.pm_level_solve = 1;
  c.pm_level_subcycle = sub ? 1 : 0;
  c.pm_mlmg_iters = 30;
  c.pm_mlmg_tol = 1e-11;
  c.dual_energy = 0;
  c.N_rt_gas = 1;
  amrex_state_init(sim);
  gs::GasState* st = gs::active();
  CHECK(st != nullptr);
  // a cold gas at rest, 1e-3 of the matter density: its gravity is negligible
  for (ax::MFIter mfi(st->G); mfi.isValid(); ++mfi) {
    auto g = st->G.array(mfi);
    ax::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
      g(i, j, k, gs::kU) = 1e-3f;
      for (int m = 1; m < 4; ++m) g(i, j, k, gs::kU + m) = 0.0f;
      g(i, j, k, gs::kU + 4) = 1e-9f;
      g(i, j, k, gs::kE) = 1e-9f;
      for (int m = 0; m < 4; ++m) g(i, j, k, gs::kW + m) = 1e-4f;
    });
  }
  sim.gas.N = kN;
  sim.gas.e_on = true;
  Units un;
  un.Dgas_avg = 1e-3;
  un.Dgas_min = 1e-10;
  un.Thydro_min = 1e-12;
  un.Tgas_min = 0;
  un.Tgas_max = 1e30;
  HydroState hs;

  // the hierarchy: hydro calls of negligible length regrid to the boxes
  {
    SimClock c0;
    c0.scalefactor = 1.0;
    c0.dt_sim = 1e-12;
    for (int it = 0; it < 3 && st->finest < 2; ++it)  // (a regrid adds one level to an empty hierarchy)
      hydro_half_step(c, un, c0, HydroHalf::First, sim.gas, hs);
    CHECK(st->finest == 2);
  }

  Domains d;
  init_domain_gas(c, d);
  const pm::PMConfig pc = pm::PMConfig::make(c, true);
  // the radial force of the clump, measured on a PM state of its own (the
  // clump and the orbiter together at the clump: the orbiter's share is
  // 1/512), averaged over directions in the orbit's plane: a_r(r) in DM-grid
  // units, and Phi(r) = -int a_r
  const int nr = 64, nphi = 16;
  const double r0 = 0.15, r1 = 0.6;
  std::vector<double> rr(nr), ar(nr, 0.0);
  {
    pm::PMState sm;
    pm::init_pm(sm, pc, d);
    sm.xv.assign(static_cast<std::size_t>(6 * sm.pc.Npart), 0.0f);
    for (std::int64_t p = 0; p < sm.pc.Npart; ++p)
      for (int m = 0; m < 3; ++m) sm.xv[static_cast<std::size_t>(6 * p + m)] = static_cast<float>(kXc[m]);
    amrex_pm_load_particles(sm);
    amrex_pm_fields(sm, nullptr, 1.0);
    std::vector<float> xs;
    for (int q = 0; q < nr; ++q) {
      rr[q] = r0 + (r1 - r0) * q / (nr - 1);
      for (int k = 0; k < nphi; ++k) {
        const double ph = 2 * M_PI * (k + 0.5) / nphi;
        xs.push_back(static_cast<float>(kXc[0] + rr[q] * std::cos(ph)));
        xs.push_back(static_cast<float>(kXc[1] + rr[q] * std::sin(ph)));
        xs.push_back(static_cast<float>(kXc[2]));
      }
    }
    std::vector<double> F;
    std::vector<int> lev;
    amrex_pm_force_at(sm, xs, F, &lev);
    for (int q = 0; q < nr; ++q)
      for (int k = 0; k < nphi; ++k) {
        const int n = q * nphi + k;
        CHECK(lev[static_cast<std::size_t>(n)] == 2);
        const double ph = 2 * M_PI * (k + 0.5) / nphi;
        ar[q] += (F[3 * n] * std::cos(ph) + F[3 * n + 1] * std::sin(ph)) / sm.pc.x_dm_m / nphi;
      }
  }
  auto phi = [&](double r) {  // -int_{r0}^{r} a_r dr (trapezoid, linear in between)
    double p = 0;
    for (int q = 0; q + 1 < nr; ++q) {
      const double a = rr[q], b = rr[q + 1];
      if (r <= a) break;
      const double hi = std::fmin(r, b);
      const double fa = ar[q], fb = ar[q] + (ar[q + 1] - ar[q]) * (hi - a) / (b - a);
      p -= 0.5 * (fa + fb) * (hi - a);
    }
    return p;
  };
  auto ar_at = [&](double r) {
    for (int q = 0; q + 1 < nr; ++q)
      if (r <= rr[q + 1]) return ar[q] + (ar[q + 1] - ar[q]) * (r - rr[q]) / (rr[q + 1] - rr[q]);
    return ar[nr - 1];
  };
  // the orbiter (the last particle) on the circular orbit in the x-y plane
  const double vc = std::sqrt(kR * std::fabs(ar_at(kR)));
  const double period = 2 * M_PI * kR / vc;
  pm::PMState s;
  pm::init_pm(s, pc, d);
  const std::int64_t np = s.pc.Npart;
  s.xv.assign(static_cast<std::size_t>(6 * np), 0.0f);
  for (std::int64_t p = 0; p < np; ++p)
    for (int m = 0; m < 3; ++m) s.xv[static_cast<std::size_t>(6 * p + m)] = static_cast<float>(kXc[m]);
  {
    float* p = s.xv_p(np);
    p[0] = static_cast<float>(kXc[0] + kR);
    p[4] = static_cast<float>(vc);
  }
  amrex_pm_load_particles(s);

  Run out;
  auto energy = [&]() {
    amrex_pm_gather_particles(s);
    double e = 0;
    if (io_rank()) {
      const float* cl = s.xv_p(1);
      const float* o = s.xv_p(np);
      double r2 = 0, v2 = 0;
      for (int m = 0; m < 3; ++m) {
        double dx = static_cast<double>(o[m]) - cl[m];
        dx -= kNdm * std::floor(dx / kNdm + 0.5);
        r2 += dx * dx;
        const double dv = static_cast<double>(o[3 + m]) - cl[3 + m];
        v2 += dv * dv;
      }
      e = 0.5 * v2 + phi(std::sqrt(r2));
    }
    ax::ParallelDescriptor::Bcast(&e, 1, ax::ParallelDescriptor::IOProcessorNumber());
    return e;
  };
  out.E0 = energy();
  out.ekin = 0.5 * vc * vc;
  SimClock clk;
  clk.scalefactor = 1.0;
  clk.redshift = 0.0;
  clk.da_sim = 0.0;
  clk.dt_sim = period / spo;
  const bool pm_sub = amrex_pm_subcycled(s);
  CHECK(pm_sub == sub);
  for (int n = 0; n < g_orbits * spo; ++n) {
    if (pm_sub) amrex_pm_subcycle_call(s, clk);
    hydro_half_step(c, un, clk, HydroHalf::First, sim.gas, hs);
    pm::particle_mesh_unsplit(s, nullptr, nullptr, 0.0, clk, un);
    if (pm_sub) amrex_pm_subcycle_call(s, clk);
    hydro_half_step(c, un, clk, HydroHalf::Second, sim.gas, hs);
    if (pm_sub) amrex_pm_subcycle_finish(s, clk);
    out.E.push_back(energy());
    if (io_rank() && std::getenv("RADRAY_PM_LEVEL_DEBUG"))
      std::printf("  %s step %d: E %.6e (E0 %.6e), substeps %lld\n", sub ? "sub" : "base", n, out.E.back(), out.E0,
                  static_cast<long long>(hs.ucounters.substeps));
  }
  out.substeps = hs.ucounters.substeps;
  if (io_rank() && std::getenv("RADRAY_PM_LEVEL_DEBUG")) std::printf("%s\n", amrex_pm_timers(s).c_str());
  amrex_hydro_release();
  amrex_state_release();
  return out;
}

// The energy error against the reference run (ref_per steps per step):
// max |E(n) - E_ref(ref_per (n + 1) - 1)| / (v^2/2) over the base steps.
// The energy itself also swings by a few per cent around the orbit (the
// mesh force's angular anisotropy against the direction-averaged Phi),
// the same in every run: the difference removes it.
double err_vs(const Run& r, const Run& ref, int ref_per) {
  double m = 0;
  for (std::size_t n = 0; n < r.E.size(); ++n)
    m = std::fmax(m, std::fabs(r.E[n] - ref.E[static_cast<std::size_t>(ref_per) * (n + 1) - 1]));
  return m / r.ekin;
}

}  // namespace

int main(int argc, char** argv) {
  backend_init(argc, argv);
  struct Fin {
    ~Fin() { backend_finalize(); }
  } fin;
  for (int a = 1; a + 1 < argc; ++a) {
    if (std::strcmp(argv[a], "--mgs") == 0) g_mgs = std::atoi(argv[a + 1]);
    if (std::strcmp(argv[a], "--orbits") == 0) g_orbits = std::atoi(argv[a + 1]);
    if (std::strcmp(argv[a], "--steps-per-orbit") == 0) g_spo = std::atoi(argv[a + 1]);
  }
  if (io_rank())
    std::printf("test_amrex_pm_subcycle: ranks %d, max_grid_size %d, %d orbits of %d base steps\n", backend_nranks(),
                g_mgs, g_orbits, g_spo);
  const Run r0 = run(false, g_spo);
  const Run r1 = run(true, g_spo);
  const Run rr = run(true, 2 * g_spo);  // the reference: sub-cycled at half the base step
  const double d0 = err_vs(r0, rr, 2), d1 = err_vs(r1, rr, 2);
  if (io_rank()) {
    std::printf("orbit energy error against the reference (sub-cycled, half the base step), max |E - E_ref| / "
                "(v^2/2): pm_level_subcycle 0: %.3e, 1: %.3e (ratio %.1f)\n",
                d0, d1, d0 / d1);
    std::printf("hydro substeps: %lld, %lld and %lld (reference)\n", static_cast<long long>(r0.substeps),
                static_cast<long long>(r1.substeps), static_cast<long long>(rr.substeps));
    unsigned long long h = 1469598103934665603ULL;  // FNV-1a of the energies (rank independence)
    for (double e : r1.E) {
      unsigned char b[8];
      std::memcpy(b, &e, 8);
      for (unsigned char x : b) h = (h ^ x) * 1099511628211ULL;
    }
    std::printf("energies checksum (pm_level_subcycle 1): %016llx\n", h);
  }
  CHECK(d1 < d0 / 16);
  TEST_MAIN_END();
}
