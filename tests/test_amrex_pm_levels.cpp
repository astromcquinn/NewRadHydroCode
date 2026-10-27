// Tests of the level Poisson solve (pm_level_solve = 1, Stage G,
// docs/stage_g.md) on a gas state with refined levels built here (no hydro):
//   * the force of a single particle against the Ewald sum with level 1
//     over the whole box (the periodic level solve must match the uniform
//     pm_refine = 2 mesh), with level 1 over a region around the particle
//     and with two nested levels (the resolution of the finest level
//     covering a point, the error growing only at the level boundaries);
//   * no self-force on a refined level;
//   * the level-1-everywhere forces against the uniform pm_refine = 2
//     backend at the same points (the same discrete problem: the MLMG
//     solution against the FFT's);
//   * total momentum of DM + gas over one step with level 1 everywhere
//     (conserved to the solver tolerance) and with a region (reported, and
//     bounded).
// Runs on any number of ranks (every rank holds the gathered results).
//
// usage: test_amrex_pm_levels [--mgs n]
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>

#include "../src/amrex/gas_state.hpp"
#include "check.hpp"
#include "pm_ewald.hpp"
#include "radray/amrex_state.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/hydro.hpp"
#include "radray/pm.hpp"
#include "radray/pm_backend.hpp"
#include "radray/rt_backend.hpp"
#include "radray/simulation.hpp"

using namespace radray;

namespace {
using namespace pm_ewald;
namespace ax = amrex;

int g_mgs = 8;
constexpr int kN = 32;  // the gas grid (and the level-0 mesh, pm_refine = 1)

// A region of level-0 cells [lo, hi] (inclusive) per level, refined to
// level l: the test's hierarchy
struct Region {
  int lo[3], hi[3];
};

// The simulation with the gas state (amrex_state = 1) and the refined
// levels of `regions` (regions[l - 1] is level l), zero gas
struct World {
  std::unique_ptr<Simulation> sim;
  std::vector<Region> regions;

  World(int Ndm, const std::vector<Region>& regs) : regions(regs) {
    sim = std::make_unique<Simulation>();
    Config& c = sim->cfg;
    c = Config::make(kN, 4, Ndm, 4.0);
    c.rt_backend = RTBackend::AMReX;
    c.hydro_backend = HydroBackend::AMReX;
    c.pm_backend = PMBackend::AMReX;
    c.hydro_scheme = HydroScheme::Unsplit;
    c.amrex_state = true;
    c.amrex_max_grid_size = g_mgs;
    c.amrex_hydro_max_grid_size = 0;
    c.amrex_max_level = static_cast<int>(regs.size());
    c.pm_level_solve = 1;
    c.pm_mlmg_iters = 30;
    c.pm_mlmg_tol = 1e-11;
    c.N_rt_gas = 1;
    amrex_state_init(*sim);
    gs::GasState* st = gs::active();
    CHECK(st != nullptr);
    st->fine.clear();
    for (std::size_t l = 1; l <= regs.size(); ++l) {
      const Region& R = regs[l - 1];
      gs::GasState::FineLevel F;
      const ax::Box b(ax::IntVect(R.lo[0] << l, R.lo[1] << l, R.lo[2] << l),
                      ax::IntVect(((R.hi[0] + 1) << l) - 1, ((R.hi[1] + 1) << l) - 1, ((R.hi[2] + 1) << l) - 1));
      F.ba = ax::BoxArray(b);
      F.ba.maxSize(2 * g_mgs);
      F.dm = ax::DistributionMapping(F.ba);
      F.G.define(F.ba, F.dm, gs::kNC, 0);
      F.G.setVal(0.0f);
      st->fine.push_back(std::move(F));
    }
    st->finest = static_cast<int>(regs.size());
  }
  ~World() { amrex_state_release(); }
};

// The PM state of the world: host arrays on every rank, the particles of
// xv0 (DM-grid units) loaded onto the level
void make_pm(const World& w, pm::PMState& s, const std::vector<float>& xv0) {
  Domains d;
  init_domain_gas(w.sim->cfg, d);
  pm::PMConfig pc = pm::PMConfig::make(w.sim->cfg, true);
  pm::init_pm(s, pc, d);
  CHECK(static_cast<std::int64_t>(xv0.size()) == 6 * s.pc.Npart);
  s.xv = xv0;
  amrex_pm_load_particles(s);
}

// The force at random points at distance d (gas cells) from x0, against
// Ewald: the rms relative error per distance bin
struct ForceStats {
  std::vector<double> d, rms;
};

ForceStats force_errors(pm::PMState& s, const float x0[3], double a) {
  ForceStats st;
  std::mt19937 rng(11);
  std::normal_distribution<double> G(0, 1);
  const int nbin = 28, npt = 400;
  for (int b = 0; b < nbin; ++b) {
    const double d = 0.125 * std::pow(2.0, b / 4.0);
    if (d > kN / 4.0) break;
    std::vector<float> xs(static_cast<std::size_t>(3) * npt);
    std::vector<double> dx(static_cast<std::size_t>(3) * npt);
    for (int n = 0; n < npt; ++n) {
      double e[3] = {G(rng), G(rng), G(rng)};
      const double en = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
      for (int c = 0; c < 3; ++c) {
        double x = x0[c] + d * e[c] / en / kN;
        x -= std::floor(x);
        xs[3 * n + c] = static_cast<float>(x);
        dx[3 * n + c] = (static_cast<double>(xs[3 * n + c]) - x0[c]) * kN;
      }
    }
    std::vector<double> F;
    std::vector<int> lev;
    amrex_pm_force_at(s, xs, F, &lev);
    double se = 0, sl[4] = {0, 0, 0, 0};
    int nl[4] = {0, 0, 0, 0};
    for (int n = 0; n < npt; ++n) {
      double Fe[3];
      ewald_force(&dx[3 * n], kN, a, Fe);
      const double fe = std::sqrt(Fe[0] * Fe[0] + Fe[1] * Fe[1] + Fe[2] * Fe[2]);
      double err2 = 0;
      for (int c = 0; c < 3; ++c) err2 += (F[3 * n + c] - Fe[c]) * (F[3 * n + c] - Fe[c]);
      se += err2 / (fe * fe);
      sl[lev[n]] += err2 / (fe * fe);
      nl[lev[n]] += 1;
    }
    st.d.push_back(d);
    st.rms.push_back(std::sqrt(se / npt));
    if (io_rank() && std::getenv("RADRAY_PM_LEVEL_DEBUG"))
      std::printf("   d %7.4f: per level (n, rms): L0 (%d, %.3e) L1 (%d, %.3e) L2 (%d, %.3e)\n", d, nl[0],
                  nl[0] ? std::sqrt(sl[0] / nl[0]) : 0.0, nl[1], nl[1] ? std::sqrt(sl[1] / nl[1]) : 0.0, nl[2],
                  nl[2] ? std::sqrt(sl[2] / nl[2]) : 0.0);
  }
  return st;
}

double rms_at(const ForceStats& st, double d) {
  for (std::size_t b = 0; b < st.d.size(); ++b)
    if (std::fabs(st.d[b] - d) < 1e-9) return st.rms[b];
  CHECK(false);
  return 1;
}

const float kX0[3] = {0.31372f, 0.65213f, 0.18944f};  // gas cells (10.0, 20.9, 6.1)

// The single particle's fields on a hierarchy; the self-force checked
ForceStats single_particle(const std::vector<Region>& regs, const char* label) {
  const double a = 1.0;
  World w(1, regs);
  pm::PMState s;
  std::vector<float> xv(6, 0.0f);
  for (int c = 0; c < 3; ++c) xv[c] = kX0[c];
  make_pm(w, s, xv);
  amrex_pm_fields(s, nullptr, a);
  std::vector<float> xs(kX0, kX0 + 3);
  std::vector<double> F;
  amrex_pm_force_at(s, xs, F);
  double Fe[3];
  const double dref[3] = {1, 0, 0};
  ewald_force(dref, kN, a, Fe);
  const double fself = std::sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
  if (io_rank()) std::printf("%s: self-force / F(1 gas cell) = %.2e\n", label, fself / std::fabs(Fe[0]));
  CHECK(fself < 1e-5 * std::fabs(Fe[0]));
  return force_errors(s, kX0, a);
}

// The uniform backend with pm_refine = r (no state), the same particle
ForceStats uniform(int r, pm::PMState& s) {
  const double a = 1.0;
  Domains d;
  init_domain_gas(Config::make(kN), d);
  pm::PMConfig pc = pm::PMConfig::make(kN, 1, true, r);
  pc.backend = PMBackend::AMReX;
  pc.amrex_max_grid_size = g_mgs * r;
  pm::init_pm(s, pc, d);
  for (int c = 0; c < 3; ++c) s.xv_p(1)[c] = kX0[c];
  std::vector<float> u(5 * s.pc.Ncell, 0.0f);
  amrex_pm_fields(s, u.data(), a);
  // the same points and bins as force_errors, through the host's force_at
  ForceStats st;
  std::mt19937 rng(11);
  std::normal_distribution<double> G(0, 1);
  const int nbin = 28, npt = 400;
  for (int b = 0; b < nbin; ++b) {
    const double dd = 0.125 * std::pow(2.0, b / 4.0);
    if (dd > kN / 4.0) break;
    double se = 0;
    for (int n = 0; n < npt; ++n) {
      double e[3] = {G(rng), G(rng), G(rng)};
      const double en = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
      float xt[3];
      double dx[3];
      for (int c = 0; c < 3; ++c) {
        double x = kX0[c] + dd * e[c] / en / kN;
        x -= std::floor(x);
        xt[c] = static_cast<float>(x);
        dx[c] = (static_cast<double>(xt[c]) - kX0[c]) * kN;
      }
      double F[3], Fe[3];
      pm::force_at(s, xt, F);
      for (int c = 0; c < 3; ++c) F[c] /= r;
      ewald_force(dx, kN, a, Fe);
      const double fe = std::sqrt(Fe[0] * Fe[0] + Fe[1] * Fe[1] + Fe[2] * Fe[2]);
      double err2 = 0;
      for (int c = 0; c < 3; ++c) err2 += (F[c] - Fe[c]) * (F[c] - Fe[c]);
      se += err2 / (fe * fe);
    }
    st.d.push_back(dd);
    st.rms.push_back(std::sqrt(se / npt));
  }
  return st;
}

void test_single_particle() {
  pm::PMState u1, u2;
  const ForceStats s1 = uniform(1, u1);
  const ForceStats s2 = uniform(2, u2);
  // level 1 over the whole box; over 16^3 level-0 cells around the particle;
  // two levels, the second over the central 8^3
  const Region whole{{0, 0, 0}, {kN - 1, kN - 1, kN - 1}};
  const Region r16{{2, 13, 0}, {17, 28, 15}};
  const Region r8{{6, 17, 2}, {13, 24, 9}};
  const ForceStats w1 = single_particle({whole}, "level 1 everywhere");
  const ForceStats p1 = single_particle({r16}, "level 1 over 16^3 cells");
  const ForceStats p2 = single_particle({r16, r8}, "levels 1 and 2 over 16^3 and 8^3 cells");
  if (io_rank()) {
    std::printf(" d [gas cells]   rms |F-F_ewald|/|F_ewald|   uniform r=1  uniform r=2  L1 everywhere  L1 region   L1+L2 region\n");
    for (std::size_t b = 0; b < s1.d.size(); ++b)
      std::printf("   %8.4f                               %9.3e    %9.3e    %9.3e      %9.3e   %9.3e\n", s1.d[b],
                  s1.rms[b], s2.rms[b], w1.rms[b], p1.rms[b], p2.rms[b]);
  }
  // level 1 everywhere: the uniform r = 2 mesh's accuracy at every distance
  for (std::size_t b = 0; b < s2.d.size(); ++b) CHECK(w1.rms[b] < 1.05 * s2.rms[b] + 1e-4);
  // a region: inside it (d <= 4) the r = 2 accuracy, well away from the
  // boundary (the particle is >= 6 cells from every face); at 8 cells the
  // points straddle the boundary: no worse than the r = 1 mesh
  for (double d : {1.0, 2.0, 4.0}) CHECK(rms_at(p1, d) < 1.25 * rms_at(s2, d) + 2e-3);
  CHECK(rms_at(p1, 8.0) < 1.25 * rms_at(s1, 8.0) + 2e-3);
  // two levels: inside the inner region (d <= 2) a 4x finer mesh than the
  // base: the criteria of the uniform r = 4 mesh (test_amrex_pm.cpp)
  CHECK(rms_at(p2, 1.0) < 0.04);
  CHECK(rms_at(p2, 2.0) < 0.012);
  CHECK(rms_at(p2, 1.0) < 1.5 * rms_at(s1, 4.0) + 2e-3);
  for (double d : {1.0, 2.0}) CHECK(rms_at(p2, d) < 0.7 * rms_at(p1, d) + 1e-3);
}

// Level 1 everywhere against the uniform pm_refine = 2 backend: the same
// particles and the same discrete problem (the CIC on the 64^3 mesh, the
// 7-point Laplacian, the 4-point gradient); only the solver differs (MLMG
// to 1e-11 against the FFT, both in float potentials)
void test_vs_uniform() {
  const double a = 0.5;
  const int Ndm = 16;
  std::vector<float> xv(6 * static_cast<std::size_t>(Ndm) * Ndm * Ndm, 0.0f);
  std::mt19937 rng(3);
  std::uniform_real_distribution<float> U(0, static_cast<float>(Ndm));
  for (std::size_t p = 0; p < xv.size() / 6; ++p)
    for (int c = 0; c < 3; ++c) xv[6 * p + c] = U(rng);
  // the uniform r = 2 backend
  pm::PMState u;
  {
    Domains d;
    init_domain_gas(Config::make(kN), d);
    pm::PMConfig pc = pm::PMConfig::make(kN, Ndm, true, 2);
    pc.backend = PMBackend::AMReX;
    pc.amrex_max_grid_size = 2 * g_mgs;
    pm::init_pm(u, pc, d);
    u.xv = xv;
    std::vector<float> ug(5 * u.pc.Ncell, 0.0f);
    amrex_pm_fields(u, ug.data(), a);
  }
  // the level solve with level 1 everywhere
  const Region whole{{0, 0, 0}, {kN - 1, kN - 1, kN - 1}};
  World w(Ndm, {whole});
  pm::PMState s;
  make_pm(w, s, xv);
  amrex_pm_fields(s, nullptr, a);
  // the forces at random points
  const int npt = 2000;
  std::vector<float> xs(static_cast<std::size_t>(3) * npt);
  std::uniform_real_distribution<float> V(0, 1);
  for (int n = 0; n < 3 * npt; ++n) xs[n] = V(rng);
  std::vector<double> F;
  amrex_pm_force_at(s, xs, F);
  double worst = 0, fmax = 0;
  for (int n = 0; n < npt; ++n) {
    double Fu[3];
    pm::force_at(u, &xs[3 * n], Fu);
    for (int m = 0; m < 3; ++m) {
      worst = std::fmax(worst, std::fabs(Fu[m] / 2 - F[3 * n + m]));
      fmax = std::fmax(fmax, std::fabs(Fu[m] / 2));
    }
  }
  if (io_rank()) std::printf("level 1 everywhere vs uniform pm_refine = 2: max |dF| / max |F| = %.2e\n", worst / fmax);
  CHECK(worst < 1e-4 * fmax);
}

// One step from rest: random particles and gas; the total momentum of DM
// and gas (level 0's composite) must stay zero: to the solver tolerance
// with level 1 everywhere, bounded with a region
void test_momentum() {
  const int Ndm = 16;
  const Region whole{{0, 0, 0}, {kN - 1, kN - 1, kN - 1}};
  const Region r16{{2, 13, 0}, {17, 28, 15}};
  const Region r8{{6, 17, 2}, {13, 24, 9}};
  struct Case {
    std::vector<Region> regs;
    const char* label;
    double tol;
  };
  const Case cases[] = {{{whole}, "level 1 everywhere", 1e-6},
                        {{whole, whole}, "levels 1 and 2 everywhere", 1e-6},
                        {{r16}, "level 1 over 16^3 cells", 1e-2},
                        {{r16, r8}, "levels 1 and 2 over 16^3 and 8^3 cells", 1e-2}};
  for (const Case& cs : cases) {
    World w(Ndm, cs.regs);
    std::vector<float> xv(6 * static_cast<std::size_t>(Ndm) * Ndm * Ndm, 0.0f);
    std::mt19937 rng(5);
    std::uniform_real_distribution<float> U(0, static_cast<float>(Ndm));
    for (std::size_t p = 0; p < xv.size() / 6; ++p)
      for (int c = 0; c < 3; ++c) xv[6 * p + c] = U(rng);
    pm::PMState s;
    make_pm(w, s, xv);
    // the gas: a random density per level-0 cell (the same on every rank),
    // at rest, E = e = 1; the fine levels take their parents' values
    gs::GasState* st = gs::active();
    std::vector<float> rho(static_cast<std::size_t>(kN) * kN * kN);
    std::uniform_real_distribution<double> R(0.05, 0.3);
    for (float& x : rho) x = static_cast<float>(R(rng));
    auto fill = [&](gs::FFab& G, int l) {
      for (ax::MFIter mfi(G); mfi.isValid(); ++mfi) {
        auto g = G.array(mfi);
        const ax::Box bx = mfi.validbox();
        ax::LoopOnCpu(bx, [&](int i, int j, int k) {
          const std::size_t c = static_cast<std::size_t>(i >> l) +
                                static_cast<std::size_t>(kN) * ((j >> l) + static_cast<std::size_t>(kN) * (k >> l));
          g(i, j, k, gs::kU) = rho[c];
          for (int m = 1; m < 4; ++m) g(i, j, k, gs::kU + m) = 0.0f;
          g(i, j, k, gs::kU + 4) = 1.0f;
          g(i, j, k, gs::kE) = 1.0f;
        });
      }
    };
    fill(st->G, 0);
    for (int l = 1; l <= st->finest; ++l) fill(st->fine[static_cast<std::size_t>(l) - 1].G, l);
    SimClock clk;
    clk.dt_sim = 0.01;
    clk.scalefactor = 0.5;
    Units un;
    un.Dgas_min = 1e-10;
    un.Tgas_min = 0;
    un.Tgas_max = 1e30;
    amrex_pm_step_unsplit(s, nullptr, nullptr, 1e-3, clk, un);
    amrex_pm_gather_particles(s);
    // rank 0's mirror: the particles; every rank: its gas cells
    const double mdm = (1 - cp::omegab / cp::omegam) * s.pc.dm_mass_ratio;
    double P[3] = {0, 0, 0}, sabs = 0;
    if (io_rank())
      for (std::int64_t ip = 1; ip <= s.pc.Npart; ++ip)
        for (int m = 0; m < 3; ++m) {
          const double p = mdm * s.xv_p(ip)[3 + m] * s.pc.x_dm_m;
          P[m] += p;
          sabs += std::fabs(p);
        }
    for (ax::MFIter mfi(st->G); mfi.isValid(); ++mfi) {
      auto g = st->G.const_array(mfi);
      ax::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
        for (int m = 0; m < 3; ++m) {
          const double p = g(i, j, k, gs::kU + 1 + m);
          P[m] += p;
          sabs += std::fabs(p);
        }
      });
    }
    ax::ParallelDescriptor::ReduceRealSum(P, 3);
    ax::ParallelDescriptor::ReduceRealSum(sabs);
    const double err = std::fmax(std::fabs(P[0]), std::fmax(std::fabs(P[1]), std::fabs(P[2]))) / sabs;
    if (io_rank()) std::printf("momentum, %s: |P| / sum|p| = %.2e (amax %.3g arms %.3g)\n", cs.label, err, clk.amax_dm, clk.arms_dm);
    CHECK(err < cs.tol);
    CHECK(clk.amax_dm > 0 && clk.arms_dm > 0);
  }
}

}  // namespace

int main(int argc, char** argv) {
  backend_init(argc, argv);
  struct Fin {
    ~Fin() { backend_finalize(); }
  } fin;
  for (int a = 1; a + 1 < argc; ++a)
    if (std::strcmp(argv[a], "--mgs") == 0) g_mgs = std::atoi(argv[a + 1]);
  if (io_rank()) std::printf("test_amrex_pm_levels: ranks %d, max_grid_size %d\n", backend_nranks(), g_mgs);
  test_single_particle();
  test_vs_uniform();
  test_momentum();
  TEST_MAIN_END();
}
