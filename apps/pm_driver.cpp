// Gravity-only (PM) driver: the C++ counterpart of
// reference/harness_pm/src/pm_harness.f90.
//
// usage: pm_driver <problem_dir> <outdir> [--nsteps n] [--quiet] [--set key=value ...]
//
// --set overrides a params.txt key (pm_backend = serial | amrex and the
// amrex_pm_* parameters of docs/stage_e_pm.md among them; with
// pm_backend = amrex the driver runs under MPI and rank 0 writes).
//
// Reads <problem_dir>/params.txt and the initial conditions
// <problem_dir>/{xv_dm.init,u_gas.init}. For each step it does what
// radhydro.f90's main loop does with hydro and RT switched off:
//   calc_sim_var; calc_timestep (complete: PM/hydro constraints, ra_max, or
//   a fixed da_sim); half clock step + calc_sim_var; particle_mesh; half clock
//   step + calc_sim_var; summarize_dm; summarize_gas (the full routine of
//   the hydro port); calc_power_spectrum every pk_every steps
// and writes <outdir>/step_%04d/{xv_dm,u_gas,v_gas,phi}.bin + stats.txt
// (step_0000: the initial state, plus w_gas.bin) and <outdir>/Pk/pk.<ztag>.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "radray/cosmology.hpp"
#include "radray/fft.hpp"
#include "radray/hydro.hpp"
#include "radray/io.hpp"
#include "radray/pm.hpp"
#include "radray/pm_backend.hpp"
#include "radray/rt_backend.hpp"

namespace {

struct PMParams {
  int Ngas = 32;
  int Ndm = -1;
  int pm_refine = 1;
  int mf = 1;
  double z_init = 300;
  double da_sim = 0;
  int nsteps = 10;
  int pk_every = 0;
  radray::Config backend;  // pm_backend, amrex_pm_max_grid_size, amrex_pm_fft, amrex_max_grid_size
};

void set_pm_param(PMParams& p, const std::string& k, std::istringstream& vs) {
  if (k == "Ngas") vs >> p.Ngas;
  else if (k == "Ndm") vs >> p.Ndm;
  else if (k == "pm_refine") vs >> p.pm_refine;
  else if (k == "mf") vs >> p.mf;
  else if (k == "z_init") vs >> p.z_init;
  else if (k == "da_sim") vs >> p.da_sim;
  else if (k == "nsteps") vs >> p.nsteps;
  else if (k == "pk_every") vs >> p.pk_every;
  else {
    std::string v;
    vs >> v;
    if (radray::set_pm_backend_param(p.backend, k, v)) return;
    if (radray::set_rt_backend_param(p.backend, k, v)) return;  // amrex_max_grid_size
    std::cerr << "pm_driver: ignoring parameter " << k << "\n";
  }
}

PMParams read_pm_params(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot open " + path);
  PMParams p;
  std::string line;
  while (std::getline(f, line)) {
    const auto h = line.find('#');
    if (h != std::string::npos) line.resize(h);
    const auto e = line.find('=');
    if (e == std::string::npos) continue;
    std::istringstream ks(line.substr(0, e)), vs(line.substr(e + 1));
    std::string k;
    ks >> k;
    set_pm_param(p, k, vs);
  }
  if (p.Ndm < 0) p.Ndm = p.Ngas / 2;
  return p;
}

}  // namespace

int main(int argc, char** argv) {
  using namespace radray;
  backend_init(argc, argv);  // pm_backend = amrex: amrex::Initialize (MPI)
  struct Fin {
    ~Fin() { backend_finalize(); }
  } fin;
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <problem_dir> <outdir> [--nsteps n] [--quiet] [--set key=value ...]\n",
                 argv[0]);
    return 2;
  }
  const std::string pdir = argv[1];
  const std::string odir = argv[2];
  int nsteps_override = -1;
  bool quiet = false;
  std::vector<std::string> sets;
  for (int a = 3; a < argc; ++a) {
    const std::string s = argv[a];
    if (s == "--nsteps" && a + 1 < argc) nsteps_override = std::atoi(argv[++a]);
    else if (s == "--quiet") quiet = true;
    else if (s == "--set" && a + 1 < argc) sets.push_back(argv[++a]);
  }

  try {
    const auto t0 = std::chrono::steady_clock::now();
    PMParams P = read_pm_params(pdir + "/params.txt");
    for (const std::string& kv : sets) {
      const auto e = kv.find('=');
      if (e == std::string::npos) throw std::runtime_error("--set expects key=value: " + kv);
      std::istringstream vs(kv.substr(e + 1));
      set_pm_param(P, kv.substr(0, e), vs);
    }
    if (nsteps_override >= 0) P.nsteps = nsteps_override;
    const bool io = io_rank();  // rank 0 writes (pm_backend = amrex)
    std::ostream* out = quiet || !io ? nullptr : &std::cout;

    if (io) std::filesystem::create_directories(odir + "/Pk");

    Config cfg = Config::make(P.Ngas, 64, P.Ndm);
    cfg.pm_refine = P.pm_refine;
    copy_pm_backend_param(P.backend, cfg);
    cfg.amrex_max_grid_size = P.backend.amrex_max_grid_size;
    check_pm_backend_param(cfg);
    if (out) *out << " pm_driver: " << pm_backend_summary(cfg) << ", ranks = " << backend_nranks() << "\n";
    pm::PMState pm;
    SimClock clk;
    Units units;
    Domains dom;
    GasGrid gas;
    HydroState hs;

    // Initialize (order as in radhydro.f90)
    init_variables(P.z_init, clk);
    init_domain_gas(cfg, dom);
    init_hydro(cfg, gas, hs);
    pm::init_pm(pm, pm::PMConfig::make(cfg, P.mf != 0), dom);
    pm.gather_potential = true;  // phi.bin below
    fft::fft_3d('f', pm.pc.Nm, pm.pc.Nm, pm.pc.Nm, pm.density.data());  // init_fft
    calc_sim_var(cfg, clk, units);

    pm::read_init_conditions(pdir, pm, gas.u, gas.w, units);

    calc_sim_var(cfg, clk, units);
    pm::summarize_dm(pm, units, out);
    summarize_gas(cfg, units, clk, gas, hs, nullptr);
    std::vector<float>& u_gas = gas.u;
    std::vector<float>& v_gas = hs.v;
    std::vector<float>& w_gas = gas.w;

    auto write_step = [&](bool initial) {
      if (!io) return;
      char name[32];
      std::snprintf(name, sizeof name, "step_%04d", clk.istep);
      const std::string d = odir + "/" + name;
      std::filesystem::create_directories(d);
      write_raw_f32(d + "/xv_dm.bin", pm.xv);
      write_raw_f32(d + "/u_gas.bin", u_gas);
      write_raw_f32(d + "/v_gas.bin", v_gas);
      // the potential on the PM mesh (Nm^3; the gas grid with pm_refine = 1)
      std::vector<float> phi(pm.pc.Nmcell);
      for (int k = 1; k <= pm.pc.Nm; ++k)
        for (int j = 1; j <= pm.pc.Nm; ++j)
          for (int i = 1; i <= pm.pc.Nm; ++i) phi[pm.mcell(i, j, k)] = pm.density[pm.dcell(i, j, k)];
      write_raw_f32(d + "/phi.bin", phi);
      if (initial) write_raw_f32(d + "/w_gas.bin", w_gas);

      const pm::DMSummary ds = pm::summarize_dm(pm, units, nullptr);
      std::FILE* fp = std::fopen((d + "/stats.txt").c_str(), "w");
      if (!fp) throw std::runtime_error("cannot write stats.txt in " + d);
      auto d1 = [&](const std::string& k, double x) { std::fprintf(fp, "%s = %.17g\n", k.c_str(), x); };
      std::fprintf(fp, "istep = %d\n", clk.istep);
      std::fprintf(fp, "nstep_rt = %d\n", clk.nstep_rt);
      d1("redshift", clk.redshift);
      d1("scalefactor", clk.scalefactor);
      d1("year", clk.year);
      d1("da_sim", clk.da_sim);
      d1("dt_sim", clk.dt_sim);
      d1("dt_yr", clk.dt_yr);
      d1("amax_dm", clk.amax_dm);
      d1("arms_dm", clk.arms_dm);
      d1("vmax_dm", clk.vmax_dm);
      d1("vrms_dm", clk.vrms_dm);
      d1("cmax_gas", clk.cmax_gas);
      d1("vmax_gas", clk.vmax_gas);
      d1("vrms_gas", clk.vrms_gas);
      const char* ax = "xyz";
      for (int m = 0; m < 3; ++m) {
        d1(std::string("dm_") + ax[m] + "_avg", ds.xavg[m]);
        d1(std::string("dm_") + ax[m] + "_min", ds.xmin[m]);
        d1(std::string("dm_") + ax[m] + "_max", ds.xmax[m]);
      }
      for (int m = 0; m < 3; ++m) {
        d1(std::string("dm_v") + ax[m] + "_avg", ds.vavg[m]);
        d1(std::string("dm_v") + ax[m] + "_min", ds.vmin[m]);
        d1(std::string("dm_v") + ax[m] + "_max", ds.vmax[m]);
      }
      std::fclose(fp);
    };
    write_step(true);

    double t_pm = 0;
    std::vector<double> pk_dm, pk_gas;
    for (int step = 1; step <= P.nsteps; ++step) {
      calc_sim_var(cfg, clk, units);
      calc_timestep(cfg, units, P.da_sim > 0 ? TimestepMode::Fixed : TimestepMode::Original,
                    P.da_sim, clk, out);

      // radhydro: forward hydro skipped
      clk.year = clk.year + clk.dt_yr / 2;
      clk.scalefactor = clk.scalefactor + clk.da_sim / 2;
      clk.redshift = 1 / clk.scalefactor - 1;
      calc_sim_var(cfg, clk, units);

      const auto tp0 = std::chrono::steady_clock::now();
      pm::particle_mesh(pm, u_gas.data(), v_gas.data(), clk, units);
      t_pm += std::chrono::duration<double>(std::chrono::steady_clock::now() - tp0).count();

      // reverse hydro skipped
      clk.year = clk.year + clk.dt_yr / 2;
      clk.scalefactor = clk.scalefactor + clk.da_sim / 2;
      clk.redshift = 1 / clk.scalefactor - 1;
      calc_sim_var(cfg, clk, units);

      pm::summarize_dm(pm, units, out);
      summarize_gas(cfg, units, clk, gas, hs, nullptr);

      if (P.pk_every > 0 && clk.istep % P.pk_every == 0) {
        pm::calc_power_spectrum(pm, u_gas.data(), cfg.box, pk_dm, pk_gas);
        if (io) pm::write_power_spectrum(odir + "/Pk/pk." + pm::ztag(clk.redshift), pm, pk_dm, pk_gas);
      }
      write_step(false);
    }

    const double total =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (io) {
      std::printf(" pm_driver: particle_mesh wall time %10.4f s\n", t_pm);
      std::printf(" pm_driver: total wall time %10.4f s\n", total);
      if (use_amrex_pm(cfg)) std::printf(" %s\n", amrex_pm_timers(pm).c_str());
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "pm_driver: %s\n", e.what());
    backend_abort(1);  // more than one rank: end the job (the others would wait)
    return 1;
  }
  return 0;
}
