// The uniform mean free path test of tools/analysis/moment_transfer.cpp
// (docs/rt_moments.md section 7) for the rays: the transfer function between
// the sources and the photons a cell holds during a substep, T(k) = sqrt(P_N /
// P_L), to compare with arctan(k lambda) / (k lambda).
//
// The ray tracing is the solver's own (rt_step_1: thresholds, pool guard,
// emission, cast / split / merge and the deposit; rt_step_3: attenuation by
// fgam and the move). Only the physics is replaced: instead of rt_step_2
// (photo_cell) every cell transmits fgam = exp(-1/lambda) of the photons
// deposited in it, in every group, so a ray loses exp(-1/lambda) per cell of
// path. The gas of the RT grid is made fully ionised (tau = 0), which makes
// the deposit the path share of the original everywhere.
//
//   ray_transfer <problem_dir> <out_prefix> --lambda L --nsrc n [--substeps m]
//                [--navg a] [--seed s] [--ngam-mean g] [--set key=value ...]
//
// problem_dir: params.txt and a uniform gas.bin of N^3 cells (sources.txt is
// not read). lambda in cells. nsrc random cells get L += 1 each (the
// generator and seed of moment_transfer.cpp: the same cells). The source
// luminosity is scaled so that the mean photon count per cell in steady state
// is --ngam-mean (default 0.1 per H atom, about that of the production runs
// after overlap; the ray deletion threshold ngam_min = 1e-6 is absolute).
// substeps: default 20 lambda (>= 400). Writes out_prefix.L, out_prefix.N (the
// photons deposited in the last substep, summed over the groups: what
// radiative_physics sees) and out_prefix.Navg (the mean over the last navg
// substeps, default 16), float32, N^3, x fastest, and prints the photon
// budget and the live rays per cell.
//
// Build (after the library): g++ -O2 -fopenmp -ffp-contract=off -Iinclude
//   tools/analysis/ray_transfer.cpp build/libradray.a -lfftw3f_omp -lfftw3f -o ray_transfer
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "radray/cosmology.hpp"
#include "radray/io.hpp"
#include "radray/kernels/ray_kernels.hpp"
#include "radray/parallel.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/rt.hpp"
#include "radray/rt_backend.hpp"
#include "radray/rt_moments.hpp"
#include "radray/simulation.hpp"

int main(int argc, char** argv) {
  using namespace radray;
  backend_init(argc, argv);
  struct Fin {
    ~Fin() { backend_finalize(); }
  } fin;
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <problem_dir> <out_prefix> --lambda L --nsrc n [--substeps m] [--navg a] "
                         "[--seed s] [--ngam-mean g] [--set key=value ...]\n", argv[0]);
    return 2;
  }
  const std::string pdir = argv[1], out = argv[2];
  double lambda = 0, ngam_mean = 0.1;
  long nsrc = 0;
  int substeps = 0, navg = 16;
  unsigned seed = 1;
  std::vector<std::string> sets;
  for (int a = 3; a < argc; ++a) {
    const std::string k = argv[a];
    if (a + 1 >= argc) {
      std::fprintf(stderr, "ray_transfer: %s needs a value\n", k.c_str());
      return 2;
    }
    const char* v = argv[++a];
    if (k == "--lambda") lambda = std::atof(v);
    else if (k == "--nsrc") nsrc = std::atol(v);
    else if (k == "--substeps") substeps = std::atoi(v);
    else if (k == "--navg") navg = std::atoi(v);
    else if (k == "--seed") seed = static_cast<unsigned>(std::atoi(v));
    else if (k == "--ngam-mean") ngam_mean = std::atof(v);
    else if (k == "--set") sets.emplace_back(v);
    else {
      std::fprintf(stderr, "ray_transfer: unknown argument '%s'\n", k.c_str());
      return 2;
    }
  }
  if (!(lambda > 0) || nsrc <= 0) {
    std::fprintf(stderr, "ray_transfer: --lambda and --nsrc are required\n");
    return 2;
  }
  if (substeps <= 0) substeps = std::max(400, static_cast<int>(20 * lambda));
  navg = std::max(1, std::min(navg, substeps));

  try {
    ProblemParams P = read_params(pdir + "/params.txt");
    for (const std::string& kv : sets) {
      const auto eq = kv.find('=');
      if (eq == std::string::npos) throw std::runtime_error("--set needs key=value: " + kv);
      set_param(P, kv.substr(0, eq), kv.substr(eq + 1));
    }
    // the set-up of apps/rt_driver.cpp
    auto sim = std::make_unique<Simulation>();
    Simulation& s = *sim;
    init_simulation(s, P.N, P.nray_max_factor, P.z_init, P.box);
    s.cfg.legacy_rt = (P.legacy_rt == 1);
    s.cfg.legacy_physics = (P.legacy_physics == 1);
    s.cfg.c_reduced = P.c_reduced;
    s.cfg.det_parallel = (P.deterministic_parallel == 1);
    s.cfg.ray_rotate = P.ray_rotate;
    s.cfg.ray_rotate_seed = static_cast<std::uint64_t>(P.ray_rotate_seed);
    s.cfg.rt_cells_per_substep = P.rt_cells_per_substep;
    s.cfg.photon_conserve = Config::photon_conserve_effective(P.photon_conserve, s.cfg.legacy_physics);
    s.cfg.photon_consume_dense =
        Config::photon_consume_dense_effective(P.photon_consume_dense, s.cfg.legacy_physics);
    s.cfg.ray_deposit = P.ray_deposit;
    copy_rt_method(P.rtm, s.cfg);
    if (use_moments(s.cfg)) throw std::runtime_error("ray_transfer is for rt_method = rays");
    copy_cosmology_param(P.rtm, s.cfg);
    init_variables(s.cfg.cosmo, P.z_init, s.clk);
    calc_sim_var(s.cfg, s.clk, s.units);
    copy_rt_backend_param(P.rtm, s.cfg);
    check_rt_backend_param(s.cfg);
    check_rt_cells_param(s.cfg);
    if (use_amrex(s.cfg)) throw std::runtime_error("ray_transfer is for rt_backend = serial");
    par::set_threads(P.threads);
    if (!s.cfg.legacy_physics) calc_atomic_rates(s.rates, s.cfg.legacy_physics);
    s.raylog = nullptr;
    s.out = nullptr;
    read_gas(pdir + "/gas.bin", s.gas);

    const int N = P.N;
    const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
    // the sources of moment_transfer.cpp (nsrc > 0)
    std::vector<float> L(nc, 0.0f);
    {
      std::mt19937_64 rng(seed);
      std::uniform_int_distribution<std::int64_t> u(0, nc - 1);
      for (long n = 0; n < nsrc; ++n) L[u(rng)] += 1.0f;
    }

    // one hydro step's clock (the substep length), as rt_driver; the clock is
    // then held: every substep is the same
    calc_timestep(s.cfg, P.da_sim, s.clk);
    advance_clock_hydro_step(s.cfg, s.clk, s.units);
    s.clk.istep_rt = 1;

    // luminosity per unit L for the wanted steady-state mean: per substep the
    // sources emit sum(L) q photons and the cells absorb f of what they hold
    const double fabs_ = 1 - std::exp(-1.0 / lambda);
    double q1 = 0;
    {
      double ng[rrk::NF];
      rrk::source_photons(1.0f, rt_params(s), ng);
      for (int g = 0; g < rrk::NF; ++g) q1 += ng[g];
    }
    const double luv = ngam_mean * fabs_ * static_cast<double>(nc) / (static_cast<double>(nsrc) * q1);
    std::vector<Source> sources;
    for (int k = 0; k < N; ++k)
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
          const float l = L[i + static_cast<std::int64_t>(N) * (j + static_cast<std::int64_t>(N) * k)];
          if (l > 0) {
            Source so;
            so.i = i + 1;
            so.j = j + 1;
            so.k = k + 1;
            so.Luv = static_cast<float>(l * luv);
            so.sfr = 1.0f;
            sources.push_back(so);
          }
        }
    calc_halo(s, sources);
    prep_rt(s);
    RTGrid& G = s.rt;
    // fully ionised gas: tau = 0 in calc_rt_1, the deposit is the path share
    std::fill(G.fHI.begin(), G.fHI.end(), 0.0f);
    std::fill(G.fHeI.begin(), G.fHeI.end(), 0.0f);
    std::fill(G.fHeII.begin(), G.fHeII.end(), 0.0f);
    const float surv = static_cast<float>(std::exp(-1.0 / lambda));

    std::printf("# N = %d, lambda = %g cells, %zu source cells (%ld draws), nray_max_factor = %d, "
                "Nray_max = %lld, threads = %d\n# substeps = %d, navg = %d, Luv per unit L = %.6e, "
                "expected mean photons per cell = %g, det_parallel = %d, ray_rotate = %d\n",
                N, lambda, sources.size(), nsrc, P.nray_max_factor, static_cast<long long>(s.cfg.Nray_max),
                par::max_threads(), substeps, navg, luv, ngam_mean, int(s.cfg.det_parallel), s.cfg.ray_rotate);
    std::fflush(stdout);

    std::vector<double> avg(nc, 0.0);
    std::vector<float> last(nc, 0.0f);
    for (int it = 1; it <= substeps; ++it) {
      rt_step_1(s);  // emission, cast / split / merge, the deposit
      const bool acc = it > substeps - navg;
      double tot = 0;
#pragma omp parallel for schedule(static) reduction(+ : tot)
      for (int k = 0; k < N; ++k)
        for (int j = 0; j < N; ++j)
          for (int i = 0; i < N; ++i) {
            const std::int64_t c = G.cell(i + 1, j + 1, k + 1);
            const float* ng = G.ngam_p(c);
            double sum = 0;
            for (int g = 0; g < kNfreq; ++g) sum += static_cast<double>(ng[g]);
            tot += sum;
            const std::int64_t o = i + static_cast<std::int64_t>(N) * (j + static_cast<std::int64_t>(N) * k);
            if (acc) avg[o] += sum / navg;
            if (it == substeps) last[o] = static_cast<float>(sum);
            float* fg = G.fgam_p(c);
            for (int g = 0; g < kNfreq; ++g) fg[g] = surv;  // in place of radiative_physics
          }
      if (it % 20 == 0 || it == substeps || it <= 3) {
        const std::int64_t nr = s.rays.Nray - s.rays.save.nfree;
        std::printf("substep %5d  photons per cell %.6f (steady state %.6f)  live rays per cell %.2f  "
                    "pool limited %lld forced %lld\n",
                    it, tot / static_cast<double>(nc), ngam_mean, static_cast<double>(nr) / static_cast<double>(nc),
                    static_cast<long long>(s.rays.npool_limited), static_cast<long long>(s.rays.npool_forced));
        std::fflush(stdout);
      }
      rt_step_3(s);  // attenuation by fgam, the move
      s.clk.istep_rt = it + 1;
    }

    auto dump = [&](const std::string& f, const float* v) {
      FILE* fp = std::fopen(f.c_str(), "wb");
      if (!fp) throw std::runtime_error("cannot write " + f);
      std::fwrite(v, sizeof(float), static_cast<std::size_t>(nc), fp);
      std::fclose(fp);
    };
    std::vector<float> avgf(avg.begin(), avg.end());
    dump(out + ".L", L.data());
    dump(out + ".N", last.data());
    dump(out + ".Navg", avgf.data());
    std::printf("# %s\n", format_rt_timers(s.timers).c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ray_transfer: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
