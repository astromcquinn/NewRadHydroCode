// Tests of the on-the-fly halo finder and source model:
//  * FoF partition == brute-force O(N^2) FoF on random, Gaussian-blob
//    (including blobs across the periodic boundaries) and Soneira-Peebles
//    particle sets, for several cell-grid sizes;
//  * determinism (repeat runs; 1/2/4 OpenMP threads when built with OpenMP);
//  * halo properties (periodic centre of mass, mean velocity, virial R/T,
//    SO refinement, potential-minimum centre);
//  * the mass function of a synthetic halo population drawn from
//    Sheth-Tormen, recovered by FoF;
//  * the source models.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include "check.hpp"
#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/halo_sources.hpp"
#include "radray/halos.hpp"
#include "radray/mass_function.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace radray;
using i64 = std::int64_t;

namespace {

struct Particles {
  std::vector<float> xv;  // 6 per particle
  i64 n() const { return static_cast<i64>(xv.size() / 6); }
  void add(double x, double y, double z, double L, double vx = 0, double vy = 0,
           double vz = 0) {
    auto w = [L](double a) {
      a = std::fmod(a, L);
      if (a < 0) a += L;
      float f = static_cast<float>(a);
      if (f >= static_cast<float>(L)) f = 0;
      return f;
    };
    xv.insert(xv.end(), {w(x), w(y), w(z), static_cast<float>(vx), static_cast<float>(vy),
                         static_cast<float>(vz)});
  }
};

void add_uniform(Particles& P, i64 n, double L, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> U(0, L);
  for (i64 i = 0; i < n; ++i) P.add(U(rng), U(rng), U(rng), L);
}

void add_blob(Particles& P, i64 n, const double c[3], double sigma, double L,
              std::mt19937_64& rng, double v0 = 0) {
  std::normal_distribution<double> G(0, sigma);
  std::normal_distribution<double> V(0, 1);
  for (i64 i = 0; i < n; ++i)
    P.add(c[0] + G(rng), c[1] + G(rng), c[2] + G(rng), L, v0 + V(rng), -v0 + V(rng), V(rng));
}

// Soneira & Peebles (1978) hierarchical point process.
void add_soneira_peebles(Particles& P, int levels, int eta, double lambda, double R0,
                         double L, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> U(0, 1);
  auto in_sphere = [&](double R, double out[3]) {
    for (;;) {
      double x = 2 * U(rng) - 1, y = 2 * U(rng) - 1, z = 2 * U(rng) - 1;
      if (x * x + y * y + z * z <= 1) {
        out[0] = R * x;
        out[1] = R * y;
        out[2] = R * z;
        return;
      }
    }
  };
  std::vector<std::array<double, 3>> cur = {{U(rng) * L, U(rng) * L, U(rng) * L}};
  double R = R0;
  for (int l = 0; l < levels; ++l) {
    std::vector<std::array<double, 3>> next;
    for (const auto& c : cur)
      for (int e = 0; e < eta; ++e) {
        double d[3];
        in_sphere(R, d);
        next.push_back({c[0] + d[0], c[1] + d[1], c[2] + d[2]});
      }
    cur.swap(next);
    R /= lambda;
  }
  for (const auto& c : cur) P.add(c[0], c[1], c[2], L);
}

bool same_partition(const FofGroups& a, const FofGroups& b) {
  return a.label == b.label && a.group_label == b.group_label && a.offset == b.offset &&
         a.members == b.members;
}

// ---------------------------------------------------------------------------

void test_fof_vs_bruteforce() {
  std::mt19937_64 rng(12345);
  const double L = 16;
  struct Case {
    const char* name;
    Particles P;
    double ll;
  };
  std::vector<Case> cases;
  {
    Case c{"uniform", {}, 0};
    add_uniform(c.P, 3000, L, rng);
    c.ll = 0.6 * L / std::cbrt(3000.0);  // near percolation: big groups
    cases.push_back(std::move(c));
  }
  {
    Case c{"blobs+boundary", {}, 0};
    add_uniform(c.P, 1500, L, rng);
    const double centres[][3] = {{0.05, 8, 8},    {15.95, 15.9, 0.1}, {0, 0, 0},
                                 {8, 15.99, 7},   {3, 4, 5},          {10, 2, 15.8},
                                 {12, 12, 12},    {6, 0.2, 15.9}};
    for (const auto& ce : centres) add_blob(c.P, 150 + 20 * (cases.size()), ce, 0.25, L, rng);
    c.ll = 0.2 * L / std::cbrt(static_cast<double>(c.P.n()));
    cases.push_back(std::move(c));
  }
  {
    Case c{"soneira-peebles", {}, 0};
    for (int t = 0; t < 4; ++t) add_soneira_peebles(c.P, 5, 4, 1.9, 4.0, L, rng);
    add_uniform(c.P, 500, L, rng);
    c.ll = 0.2 * L / std::cbrt(static_cast<double>(c.P.n()));
    cases.push_back(std::move(c));
  }
  for (auto& c : cases) {
    const i64 n = c.P.n();
    for (i64 nmin : {i64{1}, i64{5}, i64{20}}) {
      const FofGroups ref = fof_partition_bruteforce(c.P.xv.data(), n, 6, L, c.ll, nmin);
      for (int ncmax : {0, 1, 2, 3, 7}) {
        const FofGroups got = fof_partition(c.P.xv.data(), n, 6, L, c.ll, nmin, ncmax);
        const bool ok = same_partition(ref, got);
        if (!ok)
          std::fprintf(stderr, "FoF mismatch: %s nmin=%lld ncmax=%d\n", c.name,
                       static_cast<long long>(nmin), ncmax);
        CHECK(ok);
      }
      if (nmin == 20)
        std::printf("  %-16s N=%lld ll=%.3f groups(>=20)=%lld\n", c.name,
                    static_cast<long long>(n), c.ll, static_cast<long long>(ref.ngroups()));
    }
  }
  // Groups across the boundary are really joined: blob at the corner
  {
    Particles P;
    const double cc[3] = {0, 0, 0};
    add_blob(P, 400, cc, 0.1, L, rng);
    const FofGroups g = fof_partition(P.xv.data(), P.n(), 6, L, 0.2, 1);
    CHECK(g.ngroups() >= 1);
    i64 big = 0;
    for (i64 k = 0; k < g.ngroups(); ++k) big = std::max(big, g.size(k));
    CHECK(big >= 390);
  }
}

HaloUnits test_units(int Ndm, double z) {
  // Ngas = 2*Ndm, as in radhydro_global.f90
  const Config cfg = Config::make(2 * Ndm);
  SimClock clk;
  init_variables(z, clk);
  Units u;
  calc_sim_var(cfg, clk, u);
  return make_halo_units(cfg, clk, u);
}

std::vector<char> serialise(const std::vector<Halo>& h) {
  std::vector<char> b(h.size() * sizeof(Halo));
  if (!h.empty()) std::memcpy(b.data(), h.data(), b.size());
  return b;
}

void test_determinism() {
  std::mt19937_64 rng(777);
  const int Ndm = 32;
  const double L = Ndm;
  Particles P;
  add_uniform(P, 20000, L, rng);
  for (int t = 0; t < 30; ++t) {
    std::uniform_real_distribution<double> U(0, L);
    const double c[3] = {U(rng), U(rng), U(rng)};
    add_blob(P, 50 + 37 * t, c, 0.3 + 0.02 * t, L, rng, 3);
  }
  for (int t = 0; t < 3; ++t) add_soneira_peebles(P, 5, 4, 1.9, 4, L, rng);
  const HaloUnits hu = test_units(Ndm, 8);
  for (int variant = 0; variant < 3; ++variant) {
    HaloFinderParams p;
    if (variant == 1) p.so_refine = true;
    if (variant == 2) p.centre = HaloCentre::PotentialMinimum;
    const auto ref = serialise(find_halos(P.xv.data(), P.n(), hu, p));
    CHECK(serialise(find_halos(P.xv.data(), P.n(), hu, p)) == ref);
#ifdef _OPENMP
    const int nt0 = omp_get_max_threads();
    for (int nt : {1, 2, 3, 4, 8}) {
      omp_set_num_threads(nt);
      const bool same = serialise(find_halos(P.xv.data(), P.n(), hu, p)) == ref;
      if (!same) std::fprintf(stderr, "non-deterministic with %d threads (variant %d)\n", nt, variant);
      CHECK(same);
    }
    omp_set_num_threads(nt0);
#endif
  }
#ifdef _OPENMP
  std::printf("  determinism: checked 1,2,3,4,8 threads\n");
#else
  std::printf("  determinism: serial build (repeat runs only)\n");
#endif
}

void test_properties() {
  std::mt19937_64 rng(4242);
  const int Ndm = 32;
  const double L = Ndm;
  const HaloUnits hu = test_units(Ndm, 6);
  Particles P;
  add_uniform(P, 4000, L, rng);  // sparse background, never linked
  // Halo across the corner, mean velocity (5,-5,0)
  const double c[3] = {0.1, 31.95, 0.0};
  add_blob(P, 3000, c, 0.15, L, rng, 5);
  HaloFinderParams p;
  const auto h = find_halos(P.xv.data(), P.n(), hu, p);
  CHECK(!h.empty());
  if (h.empty()) return;
  const Halo& H = h[0];
  CHECK(H.N >= 2990 && H.N <= 3010);
  CHECK_CLOSE(H.M, H.N * hu.mpart_msun, 1e-12);
  for (int d = 0; d < 3; ++d) {
    double dx = H.x[d] / hu.pos_to_mpch - c[d];
    dx -= L * std::nearbyint(dx / L);
    CHECK(std::fabs(dx) < 0.02);
    CHECK(H.x[d] >= 0 && H.x[d] < L * hu.pos_to_mpch);
  }
  CHECK(std::fabs(H.v[0] / hu.vel_to_kms - 5) < 0.1);
  CHECK(std::fabs(H.v[1] / hu.vel_to_kms + 5) < 0.1);
  // Virial definitions
  const double dvir = delta_vir_bryan_norman(hu.omegam_z);
  CHECK(dvir > 170 && dvir < 18 * pc::pi * pc::pi);
  const double Rphys = H.R / hu.hubble0 * pc::Mpc2cm / (1 + hu.redshift);
  const double rho = H.M * pc::Msun_cgs / (4.0 / 3 * pc::pi * Rphys * Rphys * Rphys);
  CHECK_CLOSE(rho, dvir * hu.rho_crit, 1e-10);
  const double T = 0.59 * pc::mH_cgs * pc::G_cgs * H.M * pc::Msun_cgs / (2 * pc::k_cgs * Rphys);
  CHECK_CLOSE(H.T, T, 1e-12);
  std::printf("  halo: N=%lld M=%.3e Msun R=%.4f Mpc/h T=%.3e K  (z=6, mp=%.3e)\n",
              static_cast<long long>(H.N), H.M, H.R, H.T, hu.mpart_msun);
  // Potential-minimum centre close to the blob centre
  HaloFinderParams pp;
  pp.centre = HaloCentre::PotentialMinimum;
  const auto hp = find_halos(P.xv.data(), P.n(), hu, pp);
  for (int d = 0; d < 3; ++d) {
    double dx = hp[0].x[d] / hu.pos_to_mpch - c[d];
    dx -= L * std::nearbyint(dx / L);
    CHECK(std::fabs(dx) < 0.1);
  }
  // SO: mean density inside R_SO >= Delta*rho_crit, and R_SO ~ the radius
  // where the Gaussian blob's enclosed density drops to the threshold
  HaloFinderParams ps;
  ps.so_refine = true;
  const auto hs = find_halos(P.xv.data(), P.n(), hu, ps);
  CHECK(!hs.empty());
  if (!hs.empty()) {
    const Halo& S = hs[0];
    const double Rp = S.R / hu.hubble0 * pc::Mpc2cm / (1 + hu.redshift);
    const double rho_so = S.M * pc::Msun_cgs / (4.0 / 3 * pc::pi * Rp * Rp * Rp);
    CHECK(rho_so >= dvir * hu.rho_crit * (1 - 1e-9));
    CHECK(S.N > 0 && S.N <= 3100);
    std::printf("  SO:   N=%lld R=%.4f Mpc/h (FoF N=%lld)\n", static_cast<long long>(S.N),
                S.R, static_cast<long long>(S.N_fof));
  }
}

// Synthetic population: halos drawn from Sheth-Tormen, each a singular
// isothermal sphere truncated at its Bryan & Norman virial radius, on a
// Poisson background; FoF must recover the input mass function.
void test_mass_function() {
  const int Ndm = 64;
  const double L = Ndm;
  const double z = 3;
  const HaloUnits hu = test_units(Ndm, z);
  const double h = hu.hubble0;
  const double Vbox = std::pow(L * hu.pos_to_mpch, 3);  // (Mpc/h)^3
  MassFunction mf;
  std::mt19937_64 rng(2026);
  const double mp_h = hu.mpart_msun * h;  // Msun/h
  const double M1 = 30 * mp_h, M2 = 3e12;
  const int nb = 24;
  const double dl = std::log(M2 / M1) / nb;
  const double dvir = delta_vir_bryan_norman(hu.omegam_z);
  Particles P;
  i64 nhalo_in = 0, nhp = 0;
  std::uniform_real_distribution<double> U(0, 1);
  for (int b = 0; b < nb; ++b) {
    const double ma = M1 * std::exp(b * dl), mb = ma * std::exp(dl);
    const double mean = mf.n_between(ma, mb, z, MassFunction::Kind::ShethTormen) * Vbox;
    std::poisson_distribution<int> pois(mean);
    const int nh = pois(rng);
    for (int k = 0; k < nh; ++k) {
      const double M = ma * std::exp(U(rng) * dl);  // Msun/h
      const i64 np = static_cast<i64>(std::llround(M / mp_h));
      const double R = virial_radius_mpch(np * hu.mpart_msun, dvir, hu) / hu.pos_to_mpch;
      const double c[3] = {U(rng) * L, U(rng) * L, U(rng) * L};
      for (i64 q = 0; q < np; ++q) {
        const double r = R * U(rng);  // SIS: M(<r) ~ r
        const double mu = 2 * U(rng) - 1, ph = 2 * pc::pi * U(rng);
        const double s = std::sqrt(1 - mu * mu);
        P.add(c[0] + r * s * std::cos(ph), c[1] + r * s * std::sin(ph), c[2] + r * mu, L);
      }
      ++nhalo_in;
      nhp += np;
    }
  }
  const i64 ntot = static_cast<i64>(Ndm) * Ndm * Ndm;
  std::mt19937_64 rng2(99);
  add_uniform(P, std::max<i64>(ntot - nhp, 0), L, rng2);

  HaloFinderParams p;
  p.min_npart = 20;
  const auto halos = find_halos(P.xv.data(), P.n(), hu, p);
  std::printf("  mass function (z=%g, %lld particles, %lld input halos, %zu FoF halos >= 20 p)\n",
              z, static_cast<long long>(P.n()), static_cast<long long>(nhalo_in), halos.size());
  std::printf("    %10s %10s %8s %8s %8s\n", "M1[Msun/h]", "M2", "FoF", "ST", "PS");
  const int nbin = 6;
  const double dlb = std::log(M2 / M1) / nbin;
  for (int b = 0; b < nbin; ++b) {
    const double ma = M1 * std::exp(b * dlb), mb = ma * std::exp(dlb);
    int nf = 0;
    for (const auto& H : halos)
      if (H.M * h >= ma && H.M * h < mb) ++nf;
    const double st = mf.n_between(ma, mb, z, MassFunction::Kind::ShethTormen) * Vbox;
    const double ps = mf.n_between(ma, mb, z, MassFunction::Kind::PressSchechter) * Vbox;
    std::printf("    %10.3e %10.3e %8d %8.1f %8.1f\n", ma, mb, nf, st, ps);
    if (st >= 20) {
      const double ratio = nf / st;
      const double tol = 0.3 + 3 / std::sqrt(st);  // 3 sigma Poisson + 30% systematics
      if (std::fabs(ratio - 1) > tol)
        std::fprintf(stderr, "mass function bin %d: FoF/ST = %.3f\n", b, ratio);
      CHECK(std::fabs(ratio - 1) <= tol);
    }
  }
  // Mass function sanity: ST > PS at high mass, ST < PS... (just finite)
  CHECK(mf.dndlnM(1e10, 6, MassFunction::Kind::ShethTormen) > 0);
  CHECK_CLOSE(mf.growth(0), 1.0, 1e-12);
  CHECK(mf.growth(3) < 0.4 && mf.growth(3) > 0.3);  // ~0.32 for Om=0.27
}

void test_source_models() {
  const int Ndm = 32;
  const HaloUnits hu = test_units(Ndm, 8);
  std::vector<Halo> halos;
  for (int i = 0; i < 50; ++i) {
    Halo H;
    H.id = i;
    H.N = 2000 - 37 * i;
    H.M = H.N * hu.mpart_msun;
    H.x[0] = 0.1 * i;
    H.R = virial_radius_mpch(H.M, delta_vir_bryan_norman(hu.omegam_z), hu);
    H.T = virial_temperature(H.M, H.R, 0.59, hu);
    halos.push_back(H);
  }
  SourceModelParams p;
  p.Tvir_min = halos[30].T;  // first 31 halos are sources
  const SourceCatalog c = assign_sources(halos, hu, p);
  CHECK(c.Nsource == 31);
  CHECK(c.halo.size() == 50);
  double sfr = 0;
  for (i64 i = 0; i < c.Nsource; ++i) {
    sfr += c.halo[i].sfr;
    CHECK_CLOSE(c.halo[i].L, 8.0e27 * c.halo[i].sfr, 1e-6);
    CHECK(c.halo[i].mstar > 0 && c.halo[i].fgas > 0);
    if (i > 0) CHECK(c.halo[i].sfr <= c.halo[i - 1].sfr);
  }
  for (std::size_t i = c.Nsource; i < c.halo.size(); ++i) CHECK(c.halo[i].L == 0);
  const double vbox = std::pow(cp::box_default / cp::hubble0, 3);
  CHECK_CLOSE(sfr, global_sfr_density(8, p, nullptr) * vbox, 1e-5);
  // Hernquist & Springel fit: 0.013 at z=0 by construction, peak at z~5
  CHECK_CLOSE(global_sfr_density(0, p, nullptr), 0.013, 1e-3);  // omegar: chi(0)!=1
  double zpk = 0, best = 0;
  for (double zz = 0; zz < 15; zz += 0.05)
    if (global_sfr_density(zz, p, nullptr) > best) {
      best = global_sfr_density(zz, p, nullptr);
      zpk = zz;
    }
  std::printf("  HS03 fit: peak %.3f Msun/yr/Mpc^3 at z=%.2f; rho_sfr(z=8)=%.3e\n", best, zpk,
              global_sfr_density(8, p, nullptr));
  CHECK(zpk > 4 && zpk < 7);
  // Table: linear interpolation
  std::vector<SfrHistoryRow> t(3);
  t[0].z = 10;
  t[0].sfr = 0.01;
  t[1].z = 8;
  t[1].sfr = 0.03;
  t[2].z = 6;
  t[2].sfr = 0.05;
  CHECK_CLOSE(global_sfr_density(9, p, &t), 0.02, 1e-12);
  CHECK(global_sfr_density(11, p, &t) == 0);
  // Dynamical model
  SourceModelParams pd;
  pd.model = SfrModel::Dynamical;
  pd.Tvir_min = 0;
  const SourceCatalog d = assign_sources(halos, hu, pd);
  CHECK(d.Nsource == 50);
  const double fb = cp::omegab / cp::omegam;
  const double rho_vir = delta_vir_bryan_norman(hu.omegam_z) * hu.rho_crit;
  const double tdyn = std::sqrt(3 * pc::pi / (32 * pc::G_cgs * rho_vir)) * pc::sec2yr;
  CHECK_CLOSE(d.halo[0].sfr, 0.05 * fb * halos[0].M / tdyn, 1e-6);
  std::printf("  dynamical: M=%.3e Msun -> sfr=%.3e Msun/yr (t_dyn=%.3e yr)\n", halos[0].M,
              d.halo[0].sfr, tdyn);
}

}  // namespace

int main() {
  std::printf("fof vs brute force\n");
  test_fof_vs_bruteforce();
  std::printf("determinism\n");
  test_determinism();
  std::printf("properties\n");
  test_properties();
  std::printf("source models\n");
  test_source_models();
  std::printf("mass function\n");
  test_mass_function();
  TEST_MAIN_END();
}
