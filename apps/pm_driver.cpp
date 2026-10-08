// Gravity-only (PM) driver: the C++ counterpart of
// reference/harness_pm/src/pm_harness.f90.
//
// usage: pm_driver <problem_dir> <outdir> [--nsteps n] [--quiet]
//
// Reads <problem_dir>/params.txt and the initial conditions
// <problem_dir>/{xv_dm.init,u_gas.init}. For each step it does what
// radhydro.f90's main loop does with hydro and RT switched off:
//   calc_sim_var; calc_timestep (complete: PM/hydro constraints, ra_max, or
//   a fixed da_sim); half clock step + calc_sim_var; particle_mesh; half clock
//   step + calc_sim_var; summarize_dm; summarize_gas (time-step part);
//   calc_power_spectrum every pk_every steps
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

#include "radray/cosmology.hpp"
#include "radray/fft.hpp"
#include "radray/io.hpp"
#include "radray/pm.hpp"

namespace {

struct PMParams {
  int Ngas = 32;
  int Ndm = -1;
  int mf = 1;
  double z_init = 300;
  double da_sim = 0;
  int nsteps = 10;
  int pk_every = 0;
};

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
    if (k == "Ngas") vs >> p.Ngas;
    else if (k == "Ndm") vs >> p.Ndm;
    else if (k == "mf") vs >> p.mf;
    else if (k == "z_init") vs >> p.z_init;
    else if (k == "da_sim") vs >> p.da_sim;
    else if (k == "nsteps") vs >> p.nsteps;
    else if (k == "pk_every") vs >> p.pk_every;
    else std::cerr << "pm_driver: ignoring parameter " << k << "\n";
  }
  if (p.Ndm < 0) p.Ndm = p.Ngas / 2;
  return p;
}

}  // namespace

int main(int argc, char** argv) {
  using namespace radray;
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <problem_dir> <outdir> [--nsteps n] [--quiet]\n", argv[0]);
    return 2;
  }
  const std::string pdir = argv[1];
  const std::string odir = argv[2];
  int nsteps_override = -1;
  bool quiet = false;
  for (int a = 3; a < argc; ++a) {
    const std::string s = argv[a];
    if (s == "--nsteps" && a + 1 < argc) nsteps_override = std::atoi(argv[++a]);
    else if (s == "--quiet") quiet = true;
  }

  try {
    const auto t0 = std::chrono::steady_clock::now();
    PMParams P = read_pm_params(pdir + "/params.txt");
    if (nsteps_override >= 0) P.nsteps = nsteps_override;
    std::ostream* out = quiet ? nullptr : &std::cout;

    std::filesystem::create_directories(odir + "/Pk");

    const Config cfg = Config::make(P.Ngas);
    pm::PMState pm;
    SimClock clk;
    Units units;
    std::vector<float> u_gas(5 * cfg.Ncell, 0.0f), v_gas(3 * cfg.Ncell, 0.0f),
        w_gas(4 * cfg.Ncell, 0.0f);

    // Initialize (order as in radhydro.f90)
    init_variables(P.z_init, clk);
    pm::init_pm(pm, pm::PMConfig::make(P.Ngas, P.Ndm, P.mf != 0));
    fft::fft_3d('f', P.Ngas, P.Ngas, P.Ngas, pm.density.data());  // init_fft
    calc_sim_var(cfg, clk, units);

    pm::read_init_conditions(pdir, pm, u_gas, w_gas, units);

    calc_sim_var(cfg, clk, units);
    pm::summarize_dm(pm, units, out);
    pm::GasSpeeds gs = pm::summarize_gas_speeds(P.Ngas, u_gas.data(), v_gas.data(),
                                                 w_gas.data(), units);

    auto write_step = [&](bool initial) {
      char name[32];
      std::snprintf(name, sizeof name, "step_%04d", clk.istep);
      const std::string d = odir + "/" + name;
      std::filesystem::create_directories(d);
      write_raw_f32(d + "/xv_dm.bin", pm.xv);
      write_raw_f32(d + "/u_gas.bin", u_gas);
      write_raw_f32(d + "/v_gas.bin", v_gas);
      std::vector<float> phi(cfg.Ncell);
      for (int k = 1; k <= P.Ngas; ++k)
        for (int j = 1; j <= P.Ngas; ++j)
          for (int i = 1; i <= P.Ngas; ++i) phi[pm.gcell(i, j, k)] = pm.density[pm.dcell(i, j, k)];
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
      d1("amax_dm", pm.amax_dm);
      d1("arms_dm", pm.arms_dm);
      d1("vmax_dm", pm.vmax_dm);
      d1("vrms_dm", pm.vrms_dm);
      d1("cmax_gas", gs.cmax_gas);
      d1("vmax_gas", gs.vmax_gas);
      d1("vrms_gas", gs.vrms_gas);
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
      pm::calc_timestep_full(cfg, pm, gs, clk, units, P.da_sim, out);

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
      gs = pm::summarize_gas_speeds(P.Ngas, u_gas.data(), v_gas.data(), w_gas.data(), units);

      if (P.pk_every > 0 && clk.istep % P.pk_every == 0) {
        pm::calc_power_spectrum(pm, u_gas.data(), pk_dm, pk_gas);
        pm::write_power_spectrum(odir + "/Pk/pk." + pm::ztag(clk.redshift), pm, pk_dm, pk_gas);
      }
      write_step(false);
    }

    const double total =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf(" pm_driver: particle_mesh wall time %10.4f s\n", t_pm);
    std::printf(" pm_driver: total wall time %10.4f s\n", total);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "pm_driver: %s\n", e.what());
    return 1;
  }
  return 0;
}
