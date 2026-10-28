// RT-only driver implementing docs/validation_spec.md.
//
// usage: rt_driver <problem_dir> <outdir> [--nsteps n] [--set key=value ...]
//
// --set overrides (or adds) one params.txt entry, e.g. --set legacy_rt=0 or
// --set legacy_physics=0. Ray snapshots for the compression study:
// --set ray_dump_dir=DIR --set ray_dump_every=K (or ray_dump_at="k1 k2").
// Lossy-precision experiment (changes results; research only):
// --set ray_quant_bits=m [--set ray_quant_dither=0] [--set ray_quant_struct=1].
// Moment RT instead of the rays (docs/rt_moments.md): --set rt_method=m1|otvet.
// Rays crossing k RT cells per substep (docs/rt_multicell.md):
// --set rt_cells_per_substep=k (legacy_rt = 0); --set rt_substep_carry=1;
// --set photon_conserve=1 (the B17 fix of photo_cell's photon limit);
// --set photon_consume_dense=0|1 (P11; both follow legacy_physics by default);
// --set ray_deposit=exposure.
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
#include "radray/ray_dump.hpp"
#include "radray/parallel.hpp"
#include "radray/plane_parallel.hpp"
#include "radray/rt.hpp"
#include "radray/rt_backend.hpp"
#include "radray/recomb_rad.hpp"
#include "radray/rt_moments.hpp"
#include "radray/simulation.hpp"

int main(int argc, char** argv) {
  using namespace radray;
  backend_init(argc, argv);  // rt_backend = amrex: amrex::Initialize (MPI)
  struct Fin {
    ~Fin() { backend_finalize(); }
  } fin;
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

    const bool io = io_rank();  // with MPI only rank 0 writes files
    if (io) std::filesystem::create_directories(odir);
    backend_barrier();
    std::ofstream raylog;
    if (io) raylog.open(odir + "/ray.log");

    auto sim = std::make_unique<Simulation>();
    Simulation& s = *sim;
    init_simulation(s, P.N, P.nray_max_factor, P.z_init, P.box);
    s.cfg.legacy_rt = (P.legacy_rt == 1);
    s.cfg.legacy_physics = (P.legacy_physics == 1);
    s.cfg.c_reduced = P.c_reduced;
    s.cfg.det_parallel = (P.deterministic_parallel == 1);
    s.cfg.ray_rotate = P.ray_rotate;
    s.cfg.ray_rotate_seed = static_cast<std::uint64_t>(P.ray_rotate_seed);
    s.cfg.redshift_per_distance = (P.redshift_per_distance == 1);
    s.cfg.rt_cells_per_substep = P.rt_cells_per_substep;  // docs/rt_multicell.md
    s.cfg.photon_conserve = Config::photon_conserve_effective(P.photon_conserve, s.cfg.legacy_physics);  // B17
    s.cfg.photon_consume_dense =
        Config::photon_consume_dense_effective(P.photon_consume_dense, s.cfg.legacy_physics);  // P11
    s.cfg.ray_deposit = P.ray_deposit;                     // share | exposure
    s.cfg.cool_suppress_overdensity = P.cool_suppress_overdensity;
    if (P.rt_substep_carry) {
      if (s.cfg.legacy_rt) throw std::runtime_error("rt_substep_carry requires legacy_rt = 0");
      s.cfg.rt_substep_carry = true;
    }
    copy_rt_method(P.rtm, s.cfg);  // rt_method = rays | m1 | otvet
    copy_rt_backend_param(P.rtm, s.cfg);  // rt_backend = serial | amrex
    copy_plane_param(P.rtm, s.cfg);       // plane_parallel (plane_parallel.hpp)
    check_recomb_param(s.cfg);     // recomb_case, recomb_rad (recomb_rad.hpp)
    check_rt_backend_param(s.cfg);
    check_rt_cells_param(s.cfg);
    check_plane_param(s.cfg);
    plane_init(s);  // plane_parallel = 1: the face sources (the RT runs every step here)
    par::set_threads(P.threads);
    // The rate tables were built by init_rt with the default (legacy) physics.
    if (!s.cfg.legacy_physics) calc_atomic_rates(s.rates, s.cfg.legacy_physics);
    s.raylog = io ? &raylog : nullptr;
    RayDumpConfig dump;
    dump.dir = P.ray_dump_dir;
    dump.every = P.ray_dump_every;
    dump.at = parse_substep_list(P.ray_dump_at);
    ray_dump_enable(s, dump);
    s.ray_quant.bits = P.ray_quant_bits;
    s.ray_quant.dither = (P.ray_quant_dither != 0);
    s.ray_quant.structural = (P.ray_quant_struct != 0);
    s.ray_quant.offbits = P.ray_quant_offbits;
    if (s.ray_quant.enabled())
      std::cout << " RESEARCH: ray_quant bits=" << P.ray_quant_bits << " dither="
                << P.ray_quant_dither << " struct=" << P.ray_quant_struct
                << " (results differ from the original on purpose)\n";
    s.out = io ? &std::cout : nullptr;

    read_gas(pdir + "/gas.bin", s.gas);
    const std::vector<Source> sources = read_sources(pdir + "/sources.txt");

    if (io) std::cout << " Problem    : " << pdir << "\n"
              << " N          : " << P.N << "\n"
              << " Nray_max   : " << s.cfg.Nray_max << "\n"
              << " z_init     : " << P.z_init << "\n"
              << " da_sim     : " << P.da_sim << "\n"
              << " nsteps     : " << P.nsteps << "\n"
              << " legacy_rt  : " << P.legacy_rt << "\n"
              << " legacy_phys: " << P.legacy_physics << "\n"
              << " c_reduced  : " << P.c_reduced << "\n"
              << " det_par    : " << P.deterministic_parallel << "\n"
              << " ray_rotate : " << P.ray_rotate << "\n"
              << " z_per_dist : " << P.redshift_per_distance << "\n"
              << (P.rt_cells_per_substep != 1 || P.rt_substep_carry
                      ? " rt_cells   : " + std::to_string(P.rt_cells_per_substep) +
                            " per substep, substep carry " + std::to_string(P.rt_substep_carry) + "\n"
                      : std::string())
              << (s.cfg.photon_conserve != !s.cfg.legacy_physics
                      ? " photon_cons: " + std::to_string(int(s.cfg.photon_conserve)) + " (B17 fix; set explicitly)\n"
                      : std::string())
              << (s.cfg.photon_consume_dense != !s.cfg.legacy_physics
                      ? " consume_den: " + std::to_string(int(s.cfg.photon_consume_dense)) +
                            " (P11 fix; set explicitly)\n"
                      : std::string())
              << (P.ray_deposit ? " ray_deposit: exposure\n" : "")
              << " rt_method  : " << rt_method_summary(s.cfg) << "\n"
              << (plane_summary(s.cfg).empty() ? std::string() : " plane      : " + plane_summary(s.cfg) + "\n")
              << (P.box != 5 ? " box        : " + std::to_string(P.box) + " Mpc/h\n" : std::string())
              << " backend    : " << rt_backend_summary(s.cfg) << ", ranks " << backend_nranks() << "\n"
              << " threads    : " << par::max_threads() << "\n"
              << " sources    : " << sources.size() << "\n";

    calc_sim_var(s.cfg, s.clk, s.units);
    if (io) plane_write_spectrum(odir, s);  // plane_parallel = 1: the bins and fluxes
    for (int step = 1; step <= P.nsteps; ++step) {
      calc_sim_var(s.cfg, s.clk, s.units);
      const double z1 = s.clk.redshift;
      calc_timestep(s.cfg, P.da_sim, s.clk);
      advance_clock_hydro_step(s.cfg, s.clk, s.units);
      calc_halo(s, sources);
      radiative_transfer(s);

      if (io)
        std::cout << "\n Step number  : " << s.clk.istep << "\n"
                  << " Redshift     : " << z1 << " -> " << s.clk.redshift << "\n"
                  << " RT           : " << s.clk.nstep_rt << "  dtyr_rt = "
                  << static_cast<float>(s.clk.dtyr_rt) << "\n";
      summarize_rt(s);

      if (io) {
        char name[32];
        std::snprintf(name, sizeof name, "step_%04d", s.clk.istep);
        const std::string sdir = odir + "/" + name;
        std::filesystem::create_directories(sdir);
        write_step_output(sdir, s);
      }
    }
    backend_barrier();
    const double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (io)
      std::cout << "\n Simulation completed in " << secs << " s\n " << format_rt_timers(s.timers)
                << "\n ray pool guard: split limited in " << s.rays.npool_limited
                << " substeps, merging forced in " << s.rays.npool_forced << "\n"
                << (use_amrex(s.cfg) ? " " + amrex_rt_timers(s) + "\n" : std::string());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "rt_driver: error: %s\n", e.what());
    backend_abort(1);  // more than one rank: end the job (the others would wait)
    return 1;
  }
  return 0;
}
