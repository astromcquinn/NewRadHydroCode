// Tests of the refined PM mesh, pm_refine > 1 (docs/pm_refine.md):
//   * the force of a single particle against the Ewald sum of a point mass
//     in the periodic box (minus the mean density, as the PM solves), for
//     pm_refine = 1, 2, 4 at the same gas grid: the force is accurate down
//     to a fixed number of MESH cells, i.e. r times closer in gas cells;
//   * no self-force;
//   * the gas kick (average of the mesh forces over the gas cell / r)
//     against the same Ewald force at the gas cell centres;
//   * total momentum of DM + gas is conserved by one kick (spread and
//     average are adjoint), for every kick variant (MF, non-MF, unsplit);
//   * pm_refine = 1 through the new code paths equals the original layout.
//
// usage: test_pm_refine [--table FILE]   (FILE: the force table of the
// single-particle test, for tools/analysis/pm_refine_figures.py)
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

using namespace radray;

namespace {
using namespace pm_ewald;

void make_state(pm::PMState& s, int Ngas, int Ndm, int refine, bool mf) {
  Domains d;
  init_domain_gas(Config::make(Ngas), d);
  pm::init_pm(s, pm::PMConfig::make(Ngas, Ndm, mf, refine), d);
}

// ewald_grad, wrap, ewald_force: tests/pm_ewald.hpp (shared with test_amrex_pm)

struct ForceStats {
  int refine;
  std::vector<double> d, rms, bias;  // per distance bin (gas cells)
};

// One particle (Ndm = 1) at x0, no gas. The force at random points at
// distance d from it, in gas cells, against Ewald.
ForceStats single_particle(int Ngas, int refine, std::FILE* table) {
  const double a = 1.0;
  pm::PMState s;
  make_state(s, Ngas, 1, refine, true);
  const float x0[3] = {0.31372f, 0.65213f, 0.18944f};  // DM units: the box is [0,1)
  for (int c = 0; c < 3; ++c) s.xv_p(1)[c] = x0[c];
  std::vector<float> u(5 * s.pc.Ncell, 0.0f), v(3 * s.pc.Ncell, 0.0f);
  pm::link_dm(s);
  pm::calc_density_field(s, u.data());
  pm::calc_gravity_field(s, a);

  // no self-force
  double F[3], Fe[3];
  pm::force_at(s, s.xv_p(1), F);
  const double dref[3] = {1, 0, 0};
  ewald_force(dref, Ngas, a, Fe);
  const double fself = std::sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) / refine;
  std::printf("pm_refine = %d: self-force / F(1 gas cell) = %.2e\n", refine, fself / std::fabs(Fe[0]));
  CHECK(fself < 1e-5 * std::fabs(Fe[0]));

  ForceStats st;
  st.refine = refine;
  std::mt19937 rng(11);
  std::normal_distribution<double> G(0, 1);
  std::uniform_real_distribution<double> U01(0, 1);
  const int nbin = 28, npt = 400;
  for (int b = 0; b < nbin; ++b) {
    const double d = 0.125 * std::pow(2.0, b / 4.0);  // 0.125 .. 16 gas cells
    if (d > Ngas / 4.0) break;
    double se = 0, sb = 0;
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
      for (int c = 0; c < 3; ++c) F[c] /= refine;  // mesh -> gas cells
      ewald_force(dx, Ngas, a, Fe);
      const double fe = std::sqrt(Fe[0] * Fe[0] + Fe[1] * Fe[1] + Fe[2] * Fe[2]);
      double err2 = 0, rad = 0;
      for (int c = 0; c < 3; ++c) {
        err2 += (F[c] - Fe[c]) * (F[c] - Fe[c]);
        rad += F[c] * Fe[c];
      }
      se += err2 / (fe * fe);
      sb += rad / (fe * fe) - 1;
    }
    st.d.push_back(d);
    st.rms.push_back(std::sqrt(se / npt));
    st.bias.push_back(sb / npt);
    if (table)
      std::fprintf(table, "%d %d %.6f %.6e %.6e\n", Ngas, refine, d, st.rms.back(), st.bias.back());
  }
  (void)U01;
  return st;
}

double rms_at(const ForceStats& st, double d) {
  for (std::size_t b = 0; b < st.d.size(); ++b)
    if (std::fabs(st.d[b] - d) < 1e-9) return st.rms[b];
  CHECK(false);
  return 1;
}

void test_single_particle(std::FILE* table) {
  const int Ngas = 32;
  const ForceStats s1 = single_particle(Ngas, 1, table);
  const ForceStats s2 = single_particle(Ngas, 2, table);
  const ForceStats s4 = single_particle(Ngas, 4, table);
  std::printf(" d [gas cells]   rms |F-F_ewald|/|F_ewald|   r=1        r=2        r=4\n");
  for (std::size_t b = 0; b < s1.d.size(); ++b)
    std::printf("   %8.4f                               %9.3e  %9.3e  %9.3e\n", s1.d[b], s1.rms[b],
                s2.rms[b], s4.rms[b]);
  // The error is a function of the distance in mesh cells: pm_refine = r
  // at d/r gas cells is as accurate as pm_refine = 1 at d (to within the
  // sampling noise and the far-field periodic terms).
  // Measured (Ngas = 32): r = 1 at 4 / 8 gas cells 2.6e-2 / 9.4e-3,
  // r = 2 at 2 / 4 gas cells 2.8e-2 / 9.1e-3, r = 4 at 1 / 2 gas cells
  // 3.6e-2 / 1.05e-2.
  for (double d : {4.0, 8.0}) {
    CHECK(rms_at(s2, d / 2) < 1.25 * rms_at(s1, d) + 2e-3);
    CHECK(rms_at(s4, d / 4) < 1.5 * rms_at(s1, d) + 2e-3);
  }
  // and the absolute accuracy of the 4-point PM: < 4% at 4 mesh cells,
  // < 1.2% at 8 mesh cells
  CHECK(rms_at(s1, 4.0) < 0.04);
  CHECK(rms_at(s2, 2.0) < 0.04);
  CHECK(rms_at(s4, 1.0) < 0.04);
  CHECK(rms_at(s1, 8.0) < 0.012);
  CHECK(rms_at(s2, 4.0) < 0.012);
  CHECK(rms_at(s4, 2.0) < 0.012);
  // pm_refine = 2 is ~3x more accurate than pm_refine = 1 at fixed distance
  CHECK(rms_at(s2, 4.0) < 0.5 * rms_at(s1, 4.0));
}

// The gas kick of a single particle's field: v_gas += F dt (MF), compared
// with Ewald at the gas cell centres beyond 3 gas cells.
void test_gas_force() {
  const int Ngas = 32;
  const double a = 1.0;
  for (int r : {1, 2}) {
    pm::PMState s;
    make_state(s, Ngas, 1, r, true);
    const float x0[3] = {0.31372f, 0.65213f, 0.18944f};
    for (int c = 0; c < 3; ++c) s.xv_p(1)[c] = x0[c];
    std::vector<float> u(5 * s.pc.Ncell, 0.0f), v(3 * s.pc.Ncell, 0.0f);
    SimClock clk;
    clk.dt_sim = 1;
    Units un;
    pm::link_dm(s);
    pm::calc_density_field(s, u.data());
    pm::calc_gravity_field(s, a);
    pm::calc_kick(s, u.data(), v.data(), clk, un);
    double se = 0, n = 0, worst = 0;
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
          worst = std::fmax(worst, std::sqrt(e2 / f2));
        }
    std::printf("gas kick, pm_refine = %d: rms error 3-8 gas cells %.3e, max %.3e\n", r,
                std::sqrt(se / n), worst);
    // r = 1 at 3-8 cells is 3-8 mesh cells; r = 2: 6-16 mesh cells, plus
    // the averaging over the gas cell (a 2-cell top-hat)
    CHECK(std::sqrt(se / n) < (r == 1 ? 0.02 : 0.01));
  }
}

// One kick of random particles and gas from rest: the total momentum,
// in mesh units, must stay zero to round-off.
void test_momentum() {
  struct Case {
    int Ngas, Ndm, r;
    int kind;  // 0 MF, 1 non-MF, 2 unsplit
  };
  const Case cases[] = {{16, 8, 1, 0}, {16, 16, 1, 0},  {16, 32, 2, 0}, {16, 16, 2, 0},
                        {16, 32, 2, 1}, {16, 32, 2, 2}, {16, 64, 4, 0}, {16, 16, 1, 2}};
  for (const Case& cs : cases) {
    pm::PMState s;
    make_state(s, cs.Ngas, cs.Ndm, cs.r, cs.kind != 1);
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
    Units un;
    un.Dgas_min = 1e-10;
    un.Tgas_min = 0;
    un.Tgas_max = 1e30;
    pm::link_dm(s);
    pm::calc_density_field(s, u.data());
    pm::calc_gravity_field(s, 0.5);
    if (cs.kind == 2)
      pm::calc_kick_unsplit(s, u.data(), nullptr, 1e-3, clk, un);
    else
      pm::calc_kick(s, u.data(), v.data(), clk, un);
    const double mdm = (1 - cp::omegab / cp::omegam) * s.pc.dm_mass_ratio;  // mesh units
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
        // gas momentum change in gas units: rho dv (MF) or d(rho v)
        const double dpg = cs.kind == 0 ? static_cast<double>(u[5 * c]) * v[3 * c + m]
                                        : static_cast<double>(u[5 * c + 1 + m]) - u0[5 * c + 1 + m];
        const double p = dpg * r * r * r * r;  // mass r^3, velocity *r: mesh units
        P[m] += p;
        sabs += std::fabs(p);
      }
    const double err = std::fmax(std::fabs(P[0]), std::fmax(std::fabs(P[1]), std::fabs(P[2]))) / sabs;
    std::printf("momentum: Ngas %d, Ndm %d, pm_refine %d, %s: |P|/sum|p| = %.2e\n", cs.Ngas, cs.Ndm,
                cs.r, cs.kind == 0 ? "MF" : cs.kind == 1 ? "non-MF" : "unsplit", err);
    CHECK(err < 2e-6);
  }
}

// The parallel link_dm of pm_refine > 1 builds exactly the serial
// (original) linked lists: ll, hoc, hoc_cpu identical
void test_link_parallel() {
  pm::PMState s;
  make_state(s, 32, 64, 2, true);  // 64^3 particles on the 64^3 mesh
  std::mt19937 rng(9);
  std::uniform_real_distribution<float> U(0, 64);
  for (std::int64_t ip = 1; ip <= s.pc.Npart; ++ip)
    for (int c = 0; c < 3; ++c) s.xv_p(ip)[c] = U(rng);
  // clustered particles too: many in one row and one cell
  for (std::int64_t ip = 1; ip <= s.pc.Npart; ip += 7) {
    s.xv_p(ip)[1] = 3.25f;
    s.xv_p(ip)[2] = 17.5f;
    if (ip % 3 == 0) s.xv_p(ip)[0] = 40.1f;
  }
  pm::link_dm(s);
  const int N = s.pc.Nm;
  std::vector<std::int64_t> ll(s.pc.Npart + 1, 0), hoc(s.pc.Nmcell, 0), hc(N * N, 0);
  for (std::int64_t ip = 1; ip <= s.pc.Npart; ++ip) {
    const int j = int(double(s.xv_p(ip)[1]) * s.pc.x_dm_m) % N;
    const int k = int(double(s.xv_p(ip)[2]) * s.pc.x_dm_m) % N;
    ll[ip] = hc[j + N * k];
    hc[j + N * k] = ip;
  }
  for (int k = 0; k < N; ++k)
    for (int j = 0; j < N; ++j)
      for (std::int64_t ip = hc[j + N * k]; ip > 0;) {
        const int i = int(double(s.xv_p(ip)[0]) * s.pc.x_dm_m) % N;
        const std::int64_t next = ll[ip];
        ll[ip] = hoc[i + N * (j + std::int64_t(N) * k)];
        hoc[i + N * (j + std::int64_t(N) * k)] = ip;
        ip = next;
      }
  CHECK(ll == s.ll);
  CHECK(hoc == s.hoc);
  CHECK(hc == s.hoc_cpu);
  std::printf("link_dm (pm_refine 2, %lld particles): lists identical to the serial build: %s\n",
              static_cast<long long>(s.pc.Npart), (ll == s.ll && hoc == s.hoc) ? "yes" : "NO");
}

// pm_refine = 1 builds exactly the original configuration
void test_refine1_config() {
  const pm::PMConfig a = pm::PMConfig::make(32, 16, true, 1);
  CHECK(a.Nm == 32 && a.N_dm_m == 2 && a.N_dm_m_c == 8 && a.x_dm_m == 2.0 && a.x_dm_gas == 2.0);
  CHECK(a.nmr == a.ngasr && a.Nmcell == a.Ncell);
  const pm::PMConfig b = pm::PMConfig::make(32, 64, true, 2);
  CHECK(b.Nm == 64 && b.N_dm_m == 1 && b.x_dm_m == 1.0 && b.x_dm_gas == 0.5);
  // more particles than mesh cells (Stage G's zoom tests, docs/stage_g.md):
  // the particle mass a double below 1
  const pm::PMConfig c = pm::PMConfig::make(32, 64, true, 1);
  CHECK(c.Nm == 32 && c.N_dm_m == 0 && c.x_dm_m == 0.5 && c.dm_mass_ratio == 0.125 && c.x_dm_gas == 0.5);
  CHECK(a.dm_mass_ratio == 8.0 && b.dm_mass_ratio == 1.0);
  // neither a multiple nor a divisor: refused
  bool threw = false;
  try {
    pm::PMConfig::make(32, 48, true, 1);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);
  // domains scaled to the mesh, same order
  pm::PMState s1, s2;
  make_state(s1, 16, 8, 1, true);
  make_state(s2, 16, 32, 2, true);
  CHECK(s1.dom.box.size() == s2.dom.box.size());
  for (std::size_t n = 0; n < s1.dom.box.size(); ++n)
    for (int e = 0; e < 6; e += 2) {
      CHECK(s2.dom.box[n][e] == 2 * s1.dom.box[n][e] - 1);
      CHECK(s2.dom.box[n][e + 1] == 2 * s1.dom.box[n][e + 1]);
    }
}

}  // namespace

int main(int argc, char** argv) {
  std::FILE* table = nullptr;
  for (int a = 1; a + 1 < argc; ++a)
    if (std::strcmp(argv[a], "--table") == 0) {
      table = std::fopen(argv[a + 1], "w");
      if (table) std::fprintf(table, "# Ngas pm_refine d_gas_cells rms_rel_err mean_radial_bias\n");
    }
  test_refine1_config();
  test_link_parallel();
  test_single_particle(table);
  test_gas_force();
  test_momentum();
  if (table) std::fclose(table);
  TEST_MAIN_END();
}
