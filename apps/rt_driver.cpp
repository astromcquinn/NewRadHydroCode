// RT-only driver implementing docs/validation_spec.md.
//
// usage: rt_driver <problem_dir> <outdir> [--nsteps n] [--set key=value ...]
//
// --set overrides (or adds) one params.txt entry, e.g. --set legacy_rt=0 or
// --set legacy_physics=0.
//
// For each "hydro step" this does what radhydro.f90's main loop does with
// hydro and PM switched off:
//   calc_sim_var; calc_timestep (fixed da_sim); radhydro clock update
//   (two half steps, each followed by calc_sim_var); calc_halo (sources
//   from sources.txt, as add_halo_to_rt); radiative_transfer; summarize_rt
// and then writes <outdir>/step_%04d/{<field>.bin,stats.txt}.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "radray/cosmology.hpp"
#include "radray/io.hpp"
#include "radray/rt.hpp"
#include "radray/simulation.hpp"

int main(int argc, char** argv) {
  using namespace radray;
  if (argc < 3) {
    std::fprintf(stderr,
                 "usage: %s <problem_dir> <outdir> [--nsteps n] [--set key=value ...]\n",
                 argv[0]);
    return 2;
  }
  const std::string pdir = argv[1];
  const std::string odir = argv[2];
  int nsteps_override = -1;
  std::vector<std::string> sets;
  for (int a = 3; a < argc; ++a) {
    const std::string s = argv[a];
    if (s == "--nsteps" && a + 1 < argc) {
      nsteps_override = std::atoi(argv[++a]);
    } else if (s == "--set" && a + 1 < argc) {
      sets.emplace_back(argv[++a]);
    } else {
      std::fprintf(stderr, "rt_driver: unknown argument '%s'\n", s.c_str());
      return 2;
    }
  }

  try {
    const auto t0 = std::chrono::steady_clock::now();
    ProblemParams P = read_params(pdir + "/params.txt");
    for (const std::string& kv : sets) {
      const auto eq = kv.find('=');
      if (eq == std::string::npos) throw std::runtime_error("--set needs key=value: " + kv);
      set_param(P, kv.substr(0, eq), kv.substr(eq + 1));
    }
    if (nsteps_override >= 0) P.nsteps = nsteps_override;

    std::filesystem::create_directories(odir);
    std::ofstream raylog(odir + "/ray.log");

    auto sim = std::make_unique<Simulation>();
    Simulation& s = *sim;
    init_simulation(s, P.N, P.nray_max_factor, P.z_init);
    s.cfg.legacy_rt = (P.legacy_rt == 1);
    s.cfg.legacy_physics = (P.legacy_physics == 1);
    // The rate tables were built by init_rt with the default (legacy) physics.
    if (!s.cfg.legacy_physics) calc_atomic_rates(s.rates, s.cfg.legacy_physics);
    s.raylog = &raylog;
    s.out = &std::cout;

    read_gas(pdir + "/gas.bin", s.gas);
    const std::vector<Source> sources = read_sources(pdir + "/sources.txt");

    std::cout << " Problem    : " << pdir << "\n"
              << " N          : " << P.N << "\n"
              << " Nray_max   : " << s.cfg.Nray_max << "\n"
              << " z_init     : " << P.z_init << "\n"
              << " da_sim     : " << P.da_sim << "\n"
              << " nsteps     : " << P.nsteps << "\n"
              << " legacy_rt  : " << P.legacy_rt << "\n"
              << " legacy_phys: " << P.legacy_physics << "\n"
              << " sources    : " << sources.size() << "\n";

    calc_sim_var(s.cfg, s.clk, s.units);
    for (int step = 1; step <= P.nsteps; ++step) {
      calc_sim_var(s.cfg, s.clk, s.units);
      const double z1 = s.clk.redshift;
      calc_timestep(s.cfg, P.da_sim, s.clk);
      advance_clock_hydro_step(s.cfg, s.clk, s.units);
      calc_halo(s, sources);
      radiative_transfer(s);

      std::cout << "\n Step number  : " << s.clk.istep << "\n"
                << " Redshift     : " << z1 << " -> " << s.clk.redshift << "\n"
                << " RT           : " << s.clk.nstep_rt << "  dtyr_rt = "
                << static_cast<float>(s.clk.dtyr_rt) << "\n";
      summarize_rt(s);

      char name[32];
      std::snprintf(name, sizeof name, "step_%04d", s.clk.istep);
      const std::string sdir = odir + "/" + name;
      std::filesystem::create_directories(sdir);
      write_step_output(sdir, s);
    }
    const double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "\n Simulation completed in " << secs << " s\n";
  } catch (const std::exception& e) {
    std::fprintf(stderr, "rt_driver: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
