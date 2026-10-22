// Driver of the hydro accuracy test suite (docs/hydro_tests.md).
//
// usage: hydro_suite_driver <problem_dir> <outdir> [--set key=value ...] [--quiet]
//
// It calls the unmodified solver (src/hydro.cpp, src/pm.cpp) in the same
// order as apps/hydro_driver (mode = hydro) or apps/radhydro without RT
// (mode = cosmo). Only the *driver-level* choices differ from those apps, so
// that physics tests can be run cleanly:
//
//   * output times. mode = hydro: the step is the code's Courant step
//     dt = 2*dt_gas_factor/cmax_gas (= calc_timestep, TimestepMode::Courant)
//     times `courant` (default 1, i.e. the code's own step), shortened so
//     that the run lands exactly on the requested code times `t_out`. The
//     clock (a, z) is advanced consistently (calc_da_from_dt). mode = cosmo:
//     the original step dt = min(2*dt_gas_factor/cmax_gas,
//     dt_pm_factor/sqrt(amax_dm)), da <= a*ra_max/(1-ra_max) (as
//     calc_timestep, TimestepMode::Original, for istep > 0), additionally
//     capped at da/a <= dlna_max (test-only switch, default = ra_max/(1-ra_max),
//     i.e. off) and shortened to land on the output redshifts `z_out`.
//   * the initial gas state is read complete (u_gas.bin with the thermal
//     energy, optional w_gas.bin, v_gas.bin), not through
//     read_init_conditions (which sets T = Tgas_min).
//   * no radiative transfer, no sources (adiabatic gas).
//
// The hydro equations in the code's supercomoving variables contain no
// expansion terms for gamma = 5/3 (docs/port_notes_hydro.md §4.4), so in
// mode = hydro the code time t = sum(dt_sim) is the Newtonian time of the
// test problem; only the floors (T >= 1 K, Dgas_min) and c_limit depend on a.
//
// params.txt keys:
//   N (required), mode = hydro|cosmo, z_init (10), courant (1),
//   t_out = t1 t2 ...  (hydro), z_out = z1 z2 ... (cosmo), Ndm (N/2),
//   box (5, Mpc/h), dlna_max (ra_max/(1-ra_max)), max_steps (100000),
//   write_fields (1), pk (1: P(k) at every cosmo output)
// Inputs: u_gas.bin (5,N^3), optional w_gas.bin (4,N^3) (default neutral),
//   optional v_gas.bin (3,N^3) (default 0); cosmo: xv_dm.init (6,Ndm^3).
// Outputs: <outdir>/log.txt (one line per step), <outdir>/out_%03d/
//   {info.txt, u_gas.bin, v_gas.bin[, xv_dm.bin, pk.txt]}, out_000 = initial.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "radray/constants.hpp"
#include "radray/cosmology.hpp"
#include "radray/fft.hpp"
#include "radray/hydro.hpp"
#include "radray/hydro_backend.hpp"
#include "radray/hydro_unsplit.hpp"
#include "radray/rt_backend.hpp"
#include "radray/io.hpp"
#include "radray/pm.hpp"

namespace {

using namespace radray;
using Clock = std::chrono::steady_clock;

std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const auto b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

std::vector<double> parse_list(std::string v) {
  for (auto& c : v)
    if (c == ',') c = ' ';
  std::istringstream in(v);
  std::vector<double> out;
  double x;
  while (in >> x) out.push_back(x);
  return out;
}

struct Params {
  int N = -1;
  int Ndm = -1;
  int pm_refine = 1;  // PM mesh = pm_refine*N (docs/pm_refine.md)
  bool cosmo = false;
  double z_init = 10;
  double courant = 1;
  double box = 5;
  double dlna_max = stepp::ra_max / (1 - stepp::ra_max);
  std::vector<double> t_out, z_out;
  int max_steps = 100000;
  int write_fields = 1;
  int pk = 1;
  Config hydro;  // hydro_scheme and the unsplit parameters (docs/hydro_unsplit.md)
};

void set_param(Params& p, const std::string& k, const std::string& v) {
  if (k == "N") p.N = std::stoi(v);
  else if (k == "Ndm") p.Ndm = std::stoi(v);
  else if (k == "pm_refine") p.pm_refine = std::stoi(v);
  else if (k == "mode") {
    if (v == "hydro") p.cosmo = false;
    else if (v == "cosmo") p.cosmo = true;
    else throw std::runtime_error("mode must be hydro or cosmo");
  } else if (k == "z_init") p.z_init = std::stod(v);
  else if (k == "courant") p.courant = std::stod(v);
  else if (k == "box") p.box = std::stod(v);
  else if (k == "dlna_max") p.dlna_max = std::stod(v);
  else if (k == "t_out") p.t_out = parse_list(v);
  else if (k == "z_out") p.z_out = parse_list(v);
  else if (k == "max_steps") p.max_steps = std::stoi(v);
  else if (k == "write_fields") p.write_fields = std::stoi(v);
  else if (k == "pk") p.pk = std::stoi(v);
  else if (set_hydro_scheme_param(p.hydro, k, v)) {
  } else if (set_hydro_backend_param(p.hydro, k, v)) {  // hydro_backend = amrex, the hierarchy (docs/stage_f.md)
  } else std::cerr << "hydro_suite_driver: ignoring parameter " << k << "\n";
}

void read_f32(const std::string& path, std::vector<float>& a) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  const auto size = std::filesystem::file_size(path);
  if (size != a.size() * sizeof(float))
    throw std::runtime_error(path + ": wrong size " + std::to_string(size) + ", expected " +
                             std::to_string(a.size() * sizeof(float)));
  in.read(reinterpret_cast<char*>(a.data()), static_cast<std::streamsize>(size));
}

struct Totals {
  double mass = 0, mom[3] = {0, 0, 0}, energy = 0, thermal = 0, thermal_e = 0;
};

// Lab-frame totals (as hydro_driver's stats.txt)
Totals totals(const GasGrid& gas, const HydroState& hs) {
  Totals t;
  for (std::int64_t c = 0; c < gas.ncell; ++c) {
    const float* u = &gas.u[c * 5];
    const float* v = &hs.v[c * 3];
    const double d = u[0];
    double p[3], ekg = 0, ekl = 0;
    for (int m = 0; m < 3; ++m) {
      p[m] = static_cast<double>(u[1 + m]) + d * v[m];
      ekg += static_cast<double>(u[1 + m]) * u[1 + m];
      ekl += p[m] * p[m];
      t.mom[m] += p[m];
    }
    t.mass += d;
    t.thermal += u[4] - ekg / d / 2;
    t.energy += u[4] - ekg / d / 2 + ekl / d / 2;
    if (!gas.e.empty()) t.thermal_e += gas.e[c];
  }
  return t;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <problem_dir> <outdir> [--set key=value ...] [--quiet]\n",
                 argv[0]);
    return 2;
  }
  const std::string pdir = argv[1], odir = argv[2];
  bool quiet = false;
  std::vector<std::string> sets;
  for (int a = 3; a < argc; ++a) {
    const std::string s = argv[a];
    if (s == "--set" && a + 1 < argc) sets.emplace_back(argv[++a]);
    else if (s == "--quiet") quiet = true;
    else {
      std::fprintf(stderr, "unknown argument %s\n", s.c_str());
      return 2;
    }
  }
  backend_init(argc, argv);  // hydro_backend = amrex: amrex::Initialize (MPI)
  struct Fin {
    ~Fin() { backend_finalize(); }
  } fin;
  const bool io = io_rank();
  try {
    const auto t_start = Clock::now();
    Params P;
    {
      std::ifstream in(pdir + "/params.txt");
      if (!in) throw std::runtime_error("cannot open " + pdir + "/params.txt");
      std::string line;
      while (std::getline(in, line)) {
        const auto h = line.find('#');
        if (h != std::string::npos) line.resize(h);
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        set_param(P, trim(line.substr(0, eq)), trim(line.substr(eq + 1)));
      }
    }
    for (const auto& kv : sets) {
      const auto eq = kv.find('=');
      if (eq == std::string::npos) throw std::runtime_error("--set needs key=value");
      set_param(P, trim(kv.substr(0, eq)), trim(kv.substr(eq + 1)));
    }
    if (P.N <= 0) throw std::runtime_error("N missing");
    if (!P.cosmo && P.t_out.empty()) throw std::runtime_error("t_out missing");
    if (P.cosmo && P.z_out.empty()) throw std::runtime_error("z_out missing");
    std::ostream* out = quiet || !io ? nullptr : &std::cout;
    if (io) std::filesystem::create_directories(odir);

    Config cfg = Config::make(P.N, 1, P.Ndm, P.box);
    cfg.pm_refine = P.pm_refine;
    copy_hydro_scheme(P.hydro, cfg);
    copy_hydro_backend_param(P.hydro, cfg);
    check_hydro_backend_param(cfg);
    const bool unsplit = cfg.hydro_scheme == HydroScheme::Unsplit;
    SimClock clk;
    Units un;
    GasGrid gas;
    HydroState hs;
    Domains dom;
    pm::PMState pms;

    init_variables(cfg.cosmo, P.z_init, clk);
    init_domain_gas(cfg, dom);
    init_hydro(cfg, gas, hs);
    if (P.cosmo) {
      pm::init_pm(pms, pm::PMConfig::make(cfg), dom);
      fft::fft_3d('f', pms.pc.Nm, pms.pc.Nm, pms.pc.Nm, pms.density.data());  // init_fft
    }
    calc_sim_var(cfg, clk, un);

    read_f32(pdir + "/u_gas.bin", gas.u);
    if (std::filesystem::exists(pdir + "/w_gas.bin")) {
      read_f32(pdir + "/w_gas.bin", gas.w);
    } else {
      for (std::int64_t c = 0; c < gas.ncell; ++c) {
        const double rho = gas.u[c * 5];
        gas.w[c * 4 + 0] = static_cast<float>(cp::XHy * rho);
        gas.w[c * 4 + 1] = static_cast<float>(cp::YHe * rho);
        gas.w[c * 4 + 2] = 0;
        gas.w[c * 4 + 3] = 0;
      }
    }
    if (std::filesystem::exists(pdir + "/v_gas.bin")) read_f32(pdir + "/v_gas.bin", hs.v);
    if (P.cosmo) read_f32(pdir + "/xv_dm.init", pms.xv);

    std::ofstream log(odir + "/log.txt");
    log << "# istep t a z dt_sim cmax_gas amax_dm mass mom_x mom_y mom_z energy thermal "
           "wall_step[s]\n";
    log.precision(10);

    double t = 0;
    int iout = 0;
    double wall_hydro = 0, wall_pm = 0;
    auto write_output = [&](double wall) {
      char name[32];
      std::snprintf(name, sizeof name, "out_%03d", iout);
      const std::string d = odir + "/" + name;
      if (!io) {  // every rank holds the host state (Stage E1 model); rank 0 writes
        ++iout;
        return;
      }
      std::filesystem::create_directories(d);
      if (P.write_fields) {
        write_raw_f32(d + "/u_gas.bin", gas.u);
        write_raw_f32(d + "/v_gas.bin", hs.v);
        if (!gas.e.empty()) write_raw_f32(d + "/e_gas.bin", gas.e);
        if (P.cosmo) write_raw_f32(d + "/xv_dm.bin", pms.xv);
      }
      if (P.cosmo && P.pk) {
        std::vector<double> pk_dm, pk_gas;
        pm::calc_power_spectrum(pms, gas.u.data(), cfg.box, pk_dm, pk_gas);
        pm::write_power_spectrum(d + "/pk.txt", pms, pk_dm, pk_gas);
      }
      const Totals T = totals(gas, hs);
      std::FILE* fp = std::fopen((d + "/info.txt").c_str(), "w");
      if (!fp) throw std::runtime_error("cannot write " + d);
      auto w = [&](const char* k, double x) { std::fprintf(fp, "%s = %.17g\n", k, x); };
      std::fprintf(fp, "istep = %d\n", clk.istep);
      w("t", t);
      w("scalefactor", clk.scalefactor);
      w("redshift", clk.redshift);
      w("temp1_unit", un.temp1_unit);
      w("vel_unit", un.vel_unit);
      w("x_unit", un.x_unit);
      w("time_unit", un.time_unit);
      w("cmax_gas", clk.cmax_gas);
      w("mass", T.mass);
      std::fprintf(fp, "momentum = %.17g %.17g %.17g\n", T.mom[0], T.mom[1], T.mom[2]);
      w("energy", T.energy);
      w("thermal", T.thermal);
      if (unsplit) {
        w("thermal_e", T.thermal_e);
        const UnsplitCounters& U = hs.ucounters;
        std::fprintf(fp, "unsplit_substeps = %lld\nunsplit_max_substeps = %lld\n"
                         "unsplit_cells_from_e = %lld\nunsplit_cells_reset = %lld\n"
                         "unsplit_cells_synced = %lld\n"
                         "unsplit_floor_T = %lld\nunsplit_floor_rho = %lld\n",
                     static_cast<long long>(U.substeps), static_cast<long long>(U.max_substeps),
                     static_cast<long long>(U.cells_from_e), static_cast<long long>(U.cells_reset),
                     static_cast<long long>(U.cells_synced),
                     static_cast<long long>(U.floor_T), static_cast<long long>(U.floor_rho));
      }
      w("wall", wall);
      w("wall_hydro", wall_hydro);
      w("wall_pm", wall_pm);
      std::fclose(fp);
      ++iout;
    };

    // Initial state: summaries as the apps do (summarize_gas clamps the state)
    calc_sim_var(cfg, clk, un);
    init_unsplit_state(cfg, un, gas, hs);  // hydro_scheme = unsplit: lab frame, gas.e
    if (P.cosmo) pm::summarize_dm(pms, un, nullptr);
    summarize_gas(cfg, un, clk, gas, hs, nullptr);
    write_output(0);

    const double inf = std::numeric_limits<double>::infinity();
    std::size_t next = 0;
    const std::size_t nout = P.cosmo ? P.z_out.size() : P.t_out.size();
    while (next < nout) {
      if (clk.istep >= P.max_steps) throw std::runtime_error("max_steps reached");
      const auto ts = Clock::now();
      calc_sim_var(cfg, clk, un);
      const double a = clk.scalefactor;
      bool hit = false;
      if (!P.cosmo) {
        // the code's Courant step (calc_timestep: dt_vg = 2*dt_gas_factor/cmax_gas)
        double dt = P.courant * 2 * stepp::dt_gas_factor / clk.cmax_gas;
        const double remain = P.t_out[next] - t;
        if (dt >= remain * (1 - 1e-12)) {
          dt = remain;
          hit = true;
        } else if (dt > 0.75 * remain) {
          dt = 0.5 * remain;  // avoid a tiny last step
        }
        clk.dt_sim = dt;
        calc_da_from_dt(cfg.cosmo, a, clk.da_sim, dt, clk.dt_yr);
      } else {
        const double dt_vg = 2 * stepp::dt_gas_factor / clk.cmax_gas;
        const double dt_ap = clk.amax_dm > 0 ? stepp::dt_pm_factor / std::sqrt(clk.amax_dm) : inf;
        double da;
        calc_da_from_dt(cfg.cosmo, a, da, P.courant * std::min(dt_vg, dt_ap), clk.dt_yr);
        da = std::min(da, a * stepp::ra_max / (1 - stepp::ra_max));
        da = std::min(da, a * P.dlna_max);
        const double a_out = 1 / (1 + P.z_out[next]);
        const double remain = a_out - a;
        if (da >= remain * (1 - 1e-12)) {
          da = remain;
          hit = true;
        } else if (da > 0.75 * remain) {
          da = 0.5 * remain;
        }
        clk.da_sim = da;
        calc_dt_from_da(cfg.cosmo, a, da, clk.dt_sim, clk.dt_yr);
      }
      clk.dt_sec = clk.dt_yr * pc::yr2sec;
      clk.istep = clk.istep + 1;

      const auto th = Clock::now();
      if (!P.cosmo) {
        radhydro_hydro_step(cfg, clk, un, gas, hs);
        wall_hydro += std::chrono::duration<double>(Clock::now() - th).count();
      } else {
        // radhydro.f90 radhydro without RT (apps/radhydro: hydro_half_step
        // routes to the legacy sweeps or to the unsplit step)
        hydro_half_step(cfg, un, clk, HydroHalf::First, gas, hs);
        auto t1 = Clock::now();
        wall_hydro += std::chrono::duration<double>(t1 - th).count();
        clk.year = clk.year + clk.dt_yr / 2;
        clk.scalefactor = clk.scalefactor + clk.da_sim / 2;
        clk.redshift = 1 / clk.scalefactor - 1;
        calc_sim_var(cfg, clk, un);
        if (unsplit)
          pm::particle_mesh_unsplit(pms, gas.u.data(), gas.e.data(),
                                    cfg.dual_energy ? cfg.dual_eta1 : 0.0, clk, un);
        else
          pm::particle_mesh(pms, gas.u.data(), hs.v.data(), clk, un);
        auto t2 = Clock::now();
        wall_pm += std::chrono::duration<double>(t2 - t1).count();
        hydro_half_step(cfg, un, clk, HydroHalf::Second, gas, hs);
        wall_hydro += std::chrono::duration<double>(Clock::now() - t2).count();
        clk.year = clk.year + clk.dt_yr / 2;
        clk.scalefactor = clk.scalefactor + clk.da_sim / 2;
        clk.redshift = 1 / clk.scalefactor - 1;
        calc_sim_var(cfg, clk, un);
      }
      t += clk.dt_sim;
      if (P.cosmo) pm::summarize_dm(pms, un, nullptr);
      summarize_gas(cfg, un, clk, gas, hs, nullptr);
      if (P.cosmo && hit) {
        clk.scalefactor = 1 / (1 + P.z_out[next]);  // remove round-off of the a sum
        clk.redshift = P.z_out[next];
        calc_sim_var(cfg, clk, un);
      }
      const double wall = std::chrono::duration<double>(Clock::now() - ts).count();
      const Totals T = totals(gas, hs);
      log << clk.istep << " " << t << " " << clk.scalefactor << " " << clk.redshift << " "
          << clk.dt_sim << " " << clk.cmax_gas << " " << clk.amax_dm << " " << T.mass << " "
          << T.mom[0] << " " << T.mom[1] << " " << T.mom[2] << " " << T.energy << " "
          << T.thermal << " " << wall << "\n";
      if (out && clk.istep % 20 == 0) {
        std::printf(" step %6d t = %.5f z = %.5f dt = %.4g wall %.3f s\n", clk.istep, t,
                    clk.redshift, clk.dt_sim, wall);
        std::fflush(stdout);
      }
      if (hit) {
        write_output(std::chrono::duration<double>(Clock::now() - t_start).count());
        ++next;
      }
    }
    const double total = std::chrono::duration<double>(Clock::now() - t_start).count();
    if (!io) return 0;
    std::printf(" hydro_suite_driver: %s: %d steps, t = %.6g, z = %.6g; wall %.2f s "
                "(hydro %.2f, pm %.2f)\n",
                pdir.c_str(), clk.istep, t, clk.redshift, total, wall_hydro, wall_pm);
    if (use_amrex_hydro(cfg)) std::printf(" %s\n", amrex_hydro_timers().c_str());
    std::ofstream(odir + "/timing.txt") << "steps = " << clk.istep << "\nwall = " << total
                                        << "\nwall_hydro = " << wall_hydro
                                        << "\nwall_pm = " << wall_pm << "\n";
  } catch (const std::exception& e) {
    std::fprintf(stderr, "hydro_suite_driver: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
