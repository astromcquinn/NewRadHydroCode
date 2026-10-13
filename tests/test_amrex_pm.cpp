// Tests of the AMReX PM backend (pm_backend = amrex, docs/stage_e_pm.md),
// the counterparts of test_pm_refine.cpp on the level:
//   * the force of a single particle against the Ewald sum of a point mass
//     in the periodic box, for pm_refine = 1, 2, 4 at the same gas grid
//     (the backend's mesh force gathered onto the host, interpolated with
//     the host's force_at);
//   * no self-force;
//   * the gas kick (the backend's kick_gas) against Ewald at the gas cell
//     centres;
//   * the backend's mesh force against the serial code's on the same state
//     (round-off: the deposit order and the FFT);
//   * total momentum of DM + gas conserved by one backend step, for every
//     kick variant (MF, non-MF, unsplit) and pm_refine = 1, 2, 4.
// Runs on any number of ranks (every rank holds the gathered state and
// checks it); the output is rank 0's.
//
// usage: test_amrex_pm [--fft pencil|amrex] [--mgs n]
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "check.hpp"
#include "pm_ewald.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/hydro.hpp"
#include "radray/pm.hpp"
#include "radray/pm_backend.hpp"
#include "radray/rt_backend.hpp"

using namespace radray;

namespace {
using namespace pm_ewald;

int g_fft = 0, g_mgs = 16;

void make_state(pm::PMState& s, int Ngas, int Ndm, int refine, bool mf, PMBackend be) {
  Domains d;
  init_domain_gas(Config::make(Ngas), d);
  pm::PMConfig pc = pm::PMConfig::make(Ngas, Ndm, mf, refine);
  pc.backend = be;
  pc.amrex_max_grid_size = g_mgs;
  pc.amrex_fft = g_fft;
  pm::init_pm(s, pc, d);
}

// The backend's fields of a single particle: self-force and the force at
// random points at distance d (gas cells) against Ewald
struct ForceStats {
  int refine;
  std::vector<double> d, rms;
};

ForceStats single_particle(int Ngas, int refine) {
  const double a = 1.0;
  pm::PMState s;
  make_state(s, Ngas, 1, refine, true, PMBackend::AMReX);
  const float x0[3] = {0.31372f, 0.65213f, 0.18944f};
  for (int c = 0; c < 3; ++c) s.xv_p(1)[c] = x0[c];
  std::vector<float> u(5 * s.pc.Ncell, 0.0f);
  amrex_pm_fields(s, u.data(), a);
  CHECK(s.fmesh_valid);

  double F[3], Fe[3];
  pm::force_at(s, s.xv_p(1), F);
  const double dref[3] = {1, 0, 0};
  ewald_force(dref, Ngas, a, Fe);
  const double fself = std::sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) / refine;
  if (io_rank()) std::printf("pm_refine = %d: self-force / F(1 gas cell) = %.2e\n", refine, fself / std::fabs(Fe[0]));
  CHECK(fself < 1e-5 * std::fabs(Fe[0]));

  ForceStats st;
  st.refine = refine;
  std::mt19937 rng(11);
  std::normal_distribution<double> G(0, 1);
  const int nbin = 28, npt = 400;
  for (int b = 0; b < nbin; ++b) {
    const double d = 0.125 * std::pow(2.0, b / 4.0);
    if (d > Ngas / 4.0) break;
    double se = 0;
    for (int n = 0; n < npt; ++n) {
      double e[3] = {G(rng), G(rng), G(rng)};
      const double en = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
      float xt[3];
      for (int c = 0; c < 3; ++c) {
        double x = x0[c] + d * e[c] / en / Ngas;
        x -= std::floor(x);
        xt[c] = static_cast<float>(x);
      }
      double dx[3];
      for (int c = 0; c < 3; ++c) dx[c] = (static_cast<double>(xt[c]) - x0[c]) * Ngas;
      pm::force_at(s, xt, F);
      for (int c = 0; c < 3; ++c) F[c] /= refine;
      ewald_force(dx, Ngas, a, Fe);
      const double fe = std::sqrt(Fe[0] * Fe[0] + Fe[1] * Fe[1] + Fe[2] * Fe[2]);
      double err2 = 0;
      for (int c = 0; c < 3; ++c) err2 += (F[c] - Fe[c]) * (F[c] - Fe[c]);
      se += err2 / (fe * fe);
    }
    st.d.push_back(d);
    st.rms.push_back(std::sqrt(se / npt));
  }
  return st;
}

double rms_at(const ForceStats& st, double d) {
  for (std::size_t b = 0; b < st.d.size(); ++b)
    if (std::fabs(st.d[b] - d) < 1e-9) return st.rms[b];
  CHECK(false);
  return 1;
}

void test_single_particle() {
  const int Ngas = 32;
  const ForceStats s1 = single_particle(Ngas, 1);
  const ForceStats s2 = single_particle(Ngas, 2);
  const ForceStats s4 = single_particle(Ngas, 4);
  if (io_rank()) {
    std::printf(" d [gas cells]   rms |F-F_ewald|/|F_ewald|   r=1        r=2        r=4\n");
    for (std::size_t b = 0; b < s1.d.size(); ++b)
      std::printf("   %8.4f                               %9.3e  %9.3e  %9.3e\n", s1.d[b], s1.rms[b],
                  s2.rms[b], s4.rms[b]);
  }
  // the criteria of test_pm_refine.cpp
  for (double d : {4.0, 8.0}) {
    CHECK(rms_at(s2, d / 2) < 1.25 * rms_at(s1, d) + 2e-3);
    CHECK(rms_at(s4, d / 4) < 1.5 * rms_at(s1, d) + 2e-3);
  }
  CHECK(rms_at(s1, 4.0) < 0.04);
  CHECK(rms_at(s2, 2.0) < 0.04);
  CHECK(rms_at(s4, 1.0) < 0.04);
  CHECK(rms_at(s1, 8.0) < 0.012);
  CHECK(rms_at(s2, 4.0) < 0.012);
  CHECK(rms_at(s4, 2.0) < 0.012);
  CHECK(rms_at(s2, 4.0) < 0.5 * rms_at(s1, 4.0));
}

// The backend's potential and mesh force against the serial code's on the
// same particles and gas (round-off)
void test_vs_serial() {
  for (int r : {1, 2}) {
    const int Ngas = 32, Ndm = 16 * r;
    pm::PMState a, b;
    make_state(a, Ngas, Ndm, r, true, PMBackend::Serial);
    make_state(b, Ngas, Ndm, r, true, PMBackend::AMReX);
    std::mt19937 rng(3);
    std::uniform_real_distribution<float> U(0, static_cast<float>(Ndm));
    std::uniform_real_distribution<double> R(0.05, 0.3);
    for (std::int64_t ip = 1; ip <= a.pc.Npart; ++ip)
      for (int c = 0; c < 3; ++c) a.xv_p(ip)[c] = b.xv_p(ip)[c] = U(rng);
    std::vector<float> u(5 * a.pc.Ncell, 0.0f);
    for (std::int64_t c = 0; c < a.pc.Ncell; ++c) u[5 * c] = static_cast<float>(R(rng));
    pm::link_dm(a);
    pm::calc_density_field(a, u.data());
    pm::calc_gravity_field(a, 0.5);
    amrex_pm_fields(b, u.data(), 0.5);
    double worst_phi = 0, worst_f = 0, fmax = 0, pmax = 0;
    for (int k = 1; k <= a.pc.Nm; ++k)
      for (int j = 1; j <= a.pc.Nm; ++j)
        for (int i = 1; i <= a.pc.Nm; ++i) {
          const double pa = a.density[a.dcell(i, j, k)], pb = b.density[b.dcell(i, j, k)];
          worst_phi = std::fmax(worst_phi, std::fabs(pa - pb));
          pmax = std::fmax(pmax, std::fabs(pa));
        }
    // the force at the particles
    for (std::int64_t ip = 1; ip <= a.pc.Npart; ++ip) {
      double Fa[3], Fb[3];
      pm::force_at(a, a.xv_p(ip), Fa);
      pm::force_at(b, b.xv_p(ip), Fb);
      for (int m = 0; m < 3; ++m) {
        worst_f = std::fmax(worst_f, std::fabs(Fa[m] - Fb[m]));
        fmax = std::fmax(fmax, std::fabs(Fa[m]));
      }
    }
    if (io_rank())
      std::printf("backend vs serial, pm_refine = %d: potential max |dphi| / max |phi| = %.2e, force "
                  "max |dF| / max |F| = %.2e\n",
                  r, worst_phi / pmax, worst_f / fmax);
    CHECK(worst_phi < 1e-5 * pmax);
    CHECK(worst_f < 1e-5 * fmax);
  }
}

// The gas kick of a single particle's field (MF): v_gas += F dt, against
// Ewald at the gas cell centres 3-8 gas cells away
void test_gas_force() {
  const int Ngas = 32;
  const double a = 1.0;
  for (int r : {1, 2}) {
    pm::PMState s;
    make_state(s, Ngas, 1, r, true, PMBackend::AMReX);
    const float x0[3] = {0.31372f, 0.65213f, 0.18944f};
    for (int c = 0; c < 3; ++c) s.xv_p(1)[c] = x0[c];
    std::vector<float> u(5 * s.pc.Ncell, 0.0f), v(3 * s.pc.Ncell, 0.0f);
    SimClock clk;
    clk.dt_sim = 1;
    clk.scalefactor = a;
    Units un;
    amrex_pm_step(s, u.data(), v.data(), clk, un);
    double se = 0, n = 0;
    for (int k = 1; k <= Ngas; ++k)
      for (int j = 1; j <= Ngas; ++j)
        for (int i = 1; i <= Ngas; ++i) {
          const double dx[3] = {wrap(i - 0.5 - x0[0] * Ngas, Ngas), wrap(j - 0.5 - x0[1] * Ngas, Ngas),
                                wrap(k - 0.5 - x0[2] * Ngas, Ngas)};
          const double d = std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
          if (d < 3 || d > 8) continue;
          double Fe[3];
          ewald_force(dx, Ngas, a, Fe);
          const float* vv = &v[3 * s.gcell(i, j, k)];
          double e2 = 0, f2 = 0;
          for (int m = 0; m < 3; ++m) {
            e2 += (vv[m] - Fe[m]) * (vv[m] - Fe[m]);
            f2 += Fe[m] * Fe[m];
          }
          se += e2 / f2;
          n += 1;
        }
    if (io_rank()) std::printf("gas kick, pm_refine = %d: rms error 3-8 gas cells %.3e\n", r, std::sqrt(se / n));
    CHECK(std::sqrt(se / n) < (r == 1 ? 0.02 : 0.01));
  }
}

// One backend step of random particles and gas from rest: the total
// momentum, in mesh units, must stay zero to round-off
void test_momentum() {
  struct Case {
    int Ngas, Ndm, r;
    int kind;  // 0 MF, 1 non-MF, 2 unsplit
  };
  const Case cases[] = {{16, 8, 1, 0}, {16, 16, 1, 0},  {16, 32, 2, 0}, {16, 16, 2, 0},
                        {16, 32, 2, 1}, {16, 32, 2, 2}, {16, 64, 4, 0}, {16, 16, 1, 2}, {32, 32, 1, 1}};
  for (const Case& cs : cases) {
    pm::PMState s;
    make_state(s, cs.Ngas, cs.Ndm, cs.r, cs.kind != 1, PMBackend::AMReX);
    std::mt19937 rng(5);
    std::uniform_real_distribution<float> U(0, static_cast<float>(cs.Ndm));
    std::uniform_real_distribution<double> R(0.05, 0.3);
    for (std::int64_t ip = 1; ip <= s.pc.Npart; ++ip)
      for (int c = 0; c < 6; ++c) s.xv_p(ip)[c] = c < 3 ? U(rng) : 0.0f;
    std::vector<float> u(5 * s.pc.Ncell, 0.0f), v(3 * s.pc.Ncell, 0.0f);
    for (std::int64_t c = 0; c < s.pc.Ncell; ++c) {
      u[5 * c] = static_cast<float>(R(rng));
      u[5 * c + 4] = 1.0f;
    }
    const std::vector<float> u0 = u;
    SimClock clk;
    clk.dt_sim = 0.01;
    clk.scalefactor = 0.5;
    Units un;
    un.Dgas_min = 1e-10;
    un.Tgas_min = 0;
    un.Tgas_max = 1e30;
    if (cs.kind == 2)
      amrex_pm_step_unsplit(s, u.data(), nullptr, 1e-3, clk, un);
    else
      amrex_pm_step(s, u.data(), v.data(), clk, un);
    const double mdm = (1 - cp::omegab / cp::omegam) * s.pc.N_dm_m_c;
    const int r = cs.r;
    double P[3] = {0, 0, 0}, sabs = 0;
    for (std::int64_t ip = 1; ip <= s.pc.Npart; ++ip)
      for (int m = 0; m < 3; ++m) {
        const double p = mdm * s.xv_p(ip)[3 + m] * s.pc.x_dm_m;
        P[m] += p;
        sabs += std::fabs(p);
      }
    for (std::int64_t c = 0; c < s.pc.Ncell; ++c)
      for (int m = 0; m < 3; ++m) {
        const double dpg = cs.kind == 0 ? static_cast<double>(u[5 * c]) * v[3 * c + m]
                                        : static_cast<double>(u[5 * c + 1 + m]) - u0[5 * c + 1 + m];
        const double p = dpg * r * r * r * r;
        P[m] += p;
        sabs += std::fabs(p);
      }
    const double err = std::fmax(std::fabs(P[0]), std::fmax(std::fabs(P[1]), std::fabs(P[2]))) / sabs;
    if (io_rank())
      std::printf("momentum: Ngas %d, Ndm %d, pm_refine %d, %s: |P|/sum|p| = %.2e (amax %.3g arms %.3g)\n",
                  cs.Ngas, cs.Ndm, cs.r, cs.kind == 0 ? "MF" : cs.kind == 1 ? "non-MF" : "unsplit", err,
                  clk.amax_dm, clk.arms_dm);
    CHECK(err < 2e-6);
    CHECK(clk.amax_dm > 0 && clk.arms_dm > 0 && clk.arms_dm <= clk.amax_dm * std::sqrt(3.0));
  }
}

}  // namespace

int main(int argc, char** argv) {
  backend_init(argc, argv);
  struct Fin {
    ~Fin() { backend_finalize(); }
  } fin;
  for (int a = 1; a + 1 < argc; ++a) {
    if (std::strcmp(argv[a], "--fft") == 0) g_fft = std::strcmp(argv[a + 1], "amrex") == 0 ? 1 : 0;
    if (std::strcmp(argv[a], "--mgs") == 0) g_mgs = std::atoi(argv[a + 1]);
  }
  if (io_rank()) std::printf("test_amrex_pm: ranks %d, fft %s, max_grid_size %d\n", backend_nranks(), g_fft ? "amrex" : "pencil", g_mgs);
  test_single_particle();
  test_vs_serial();
  test_gas_force();
  test_momentum();
  TEST_MAIN_END();
}
