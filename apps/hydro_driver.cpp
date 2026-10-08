// Hydro-only driver: the C++ counterpart of
// reference/harness_hydro/src/hydro_harness.f90 (see its header and
// docs/port_notes_hydro.md).
//
// usage: hydro_driver <problem_dir> <outdir> [--nsteps n] [--quiet]
//
// Reads  <problem_dir>/params.txt (N, z_init, nsteps, timestep = fixed|courant,
//        da_sim), u_gas.bin (5,N^3), w_gas.bin (4,N^3), optional v_gas.bin
//        (3,N^3).
// Writes <outdir>/step_%04d/{u_gas,w_gas,v_gas}.bin, stats.txt,
//        gas_summary.txt for the initial state (step_0000) and every step,
//        and <outdir>/domain_gas.txt.
//
// Per step, as radhydro.f90 without PM and RT:
//   calc_sim_var; calc_timestep (hydro part); radhydro (transposes and the
//   six sweeps, two half clock steps); summarize_gas.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "radray/cosmology.hpp"
#include "radray/hydro.hpp"
#include "radray/io.hpp"

namespace {

using namespace radray;

struct HydroParams {
  int N = -1;
  double z_init = 10;
  int nsteps = 10;
  TimestepMode mode = TimestepMode::Fixed;
  double da_sim = 1e-4;
};

std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const auto b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

HydroParams read_hydro_params(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  HydroParams p;
  std::string line;
  while (std::getline(in, line)) {
    const auto h = line.find('#');
    if (h != std::string::npos) line = line.substr(0, h);
    const auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    const std::string key = trim(line.substr(0, eq));
    std::string val = trim(line.substr(eq + 1));
    if (key == "timestep") {
      if (val == "fixed")
        p.mode = TimestepMode::Fixed;
      else if (val == "courant")
        p.mode = TimestepMode::Courant;
      else
        throw std::runtime_error("timestep must be fixed or courant: " + val);
      continue;
    }
    for (auto& ch : val)
      if (ch == 'd' || ch == 'D') ch = 'e';
    if (key == "N")
      p.N = std::stoi(val);
    else if (key == "z_init")
      p.z_init = std::stod(val);
    else if (key == "da_sim")
      p.da_sim = std::stod(val);
    else if (key == "nsteps")
      p.nsteps = std::stoi(val);
    else
      std::cerr << "hydro_driver: ignoring parameter " << key << "\n";
  }
  if (p.N <= 0) throw std::runtime_error("params.txt: N missing");
  return p;
}

void read_field(const std::string& path, std::vector<float>& a) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  const auto size = std::filesystem::file_size(path);
  if (size != a.size() * sizeof(float))
    throw std::runtime_error(path + ": wrong size " + std::to_string(size) + ", expected " +
                             std::to_string(a.size() * sizeof(float)));
  in.read(reinterpret_cast<char*>(a.data()), static_cast<std::streamsize>(size));
  if (!in) throw std::runtime_error("read error on " + path);
}

void write_domain_gas_test(const std::string& fn, Domains dom) {
  std::FILE* fp = std::fopen(fn.c_str(), "w");
  if (!fp) throw std::runtime_error("cannot write " + fn);
  std::fprintf(fp, "Ndomain_gas = %d\n", dom.Ndomain);
  std::fprintf(fp, "N_dom_gas =  %d %d %d\n", dom.Nx_dom, dom.Nx_dom, dom.Nx_dom);
  for (int n = 1; n <= dom.Ndomain; ++n) {
    std::fprintf(fp, "indx_dom_gas_%d =", n);
    for (int x : dom.box[n - 1]) std::fprintf(fp, " %d", x);
    std::fprintf(fp, "\ndomain_gas_%d =", n);
    for (int x : dom.domain[n - 1]) std::fprintf(fp, " %d", x);
    std::fprintf(fp, "\n");
  }
  for (int pass = 1; pass <= 3; ++pass) {
    std::fprintf(fp, "order_pass_%d =", pass);
    for (int n = 1; n <= dom.Ndomain; ++n) {
      const int idom = set_domain_gas(dom);
      std::fprintf(fp, " %d", idom + 1);
      end_domain_gas(dom, idom);
    }
    std::fprintf(fp, "\n");
  }
  std::fprintf(fp, "status_after =");
  for (int n = 1; n <= dom.Ndomain; ++n) std::fprintf(fp, " %d", dom.domain[n - 1][0]);
  std::fprintf(fp, "\n");
  std::fclose(fp);
}

void write_step_output(const std::string& dir, const SimClock& clk, const Units& un,
                       const GasGrid& gas, const HydroState& hs) {
  write_raw_f32(dir + "/u_gas.bin", gas.u);
  write_raw_f32(dir + "/w_gas.bin", gas.w);
  write_raw_f32(dir + "/v_gas.bin", hs.v);

  // Conserved totals, exactly as hydro_harness.f90: write_step_output
  double mass = 0, mom[3] = {0, 0, 0}, ener = 0, ener_grid = 0, eth_tot = 0;
  double wsum[4] = {0, 0, 0, 0}, vsum[3] = {0, 0, 0};
  for (std::int64_t c = 0; c < gas.ncell; ++c) {
    const float* u = &gas.u[c * 5];
    const float* w = &gas.w[c * 4];
    const float* v = &hs.v[c * 3];
    const double d = u[0];
    const double px = static_cast<double>(u[1]) + d * static_cast<double>(v[0]);
    const double py = static_cast<double>(u[2]) + d * static_cast<double>(v[1]);
    const double pz = static_cast<double>(u[3]) + d * static_cast<double>(v[2]);
    const double ekg = (static_cast<double>(u[1]) * u[1] + static_cast<double>(u[2]) * u[2] +
                        static_cast<double>(u[3]) * u[3]) /
                       d / 2;
    const double ekl = (px * px + py * py + pz * pz) / d / 2;
    mass = mass + d;
    mom[0] = mom[0] + px;
    mom[1] = mom[1] + py;
    mom[2] = mom[2] + pz;
    ener_grid = ener_grid + static_cast<double>(u[4]);
    eth_tot = eth_tot + (static_cast<double>(u[4]) - ekg);
    ener = ener + ((static_cast<double>(u[4]) - ekg) + ekl);
    for (int m = 0; m < 4; ++m) wsum[m] = wsum[m] + static_cast<double>(w[m]);
    for (int m = 0; m < 3; ++m) vsum[m] = vsum[m] + static_cast<double>(v[m]);
  }

  std::FILE* fp = std::fopen((dir + "/stats.txt").c_str(), "w");
  if (!fp) throw std::runtime_error("cannot write stats.txt in " + dir);
  auto d1 = [&](const char* k, double x) { std::fprintf(fp, "%s = %.17g\n", k, x); };
  std::fprintf(fp, "istep = %d\n", clk.istep);
  d1("redshift", clk.redshift);
  d1("scalefactor", clk.scalefactor);
  d1("year", clk.year);
  d1("da_sim", clk.da_sim);
  d1("dt_sim", clk.dt_sim);
  d1("dt_yr", clk.dt_yr);
  d1("dt_hy", hs.dt_hy);
  d1("c_limit", hs.c_limit);
  d1("cmax_gas", hs.cmax_gas);
  d1("vmax_gas", hs.vmax_gas);
  d1("vrms_gas", hs.vrms_gas);
  d1("vel_unit", un.vel_unit);
  d1("temp1_unit", un.temp1_unit);
  d1("mass", mass);
  std::fprintf(fp, "momentum = %.17g %.17g %.17g\n", mom[0], mom[1], mom[2]);
  d1("energy", ener);
  d1("energy_grid", ener_grid);
  d1("thermal", eth_tot);
  std::fprintf(fp, "w_total = %.17g %.17g %.17g %.17g\n", wsum[0], wsum[1], wsum[2], wsum[3]);
  std::fprintf(fp, "v_total = %.17g %.17g %.17g\n", vsum[0], vsum[1], vsum[2]);
  std::fclose(fp);
}

void write_summary(const std::string& dir, const GasSummary& s) {
  std::ofstream f(dir + "/gas_summary.txt");
  f << "summary = " << format_gas_record(s.record) << "\n";
}

std::string step_dir(const std::string& odir, int istep) {
  char name[32];
  std::snprintf(name, sizeof name, "step_%04d", istep);
  const std::string d = odir + "/" + name;
  std::filesystem::create_directories(d);
  return d;
}

}  // namespace

int main(int argc, char** argv) {
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
    if (s == "--nsteps" && a + 1 < argc)
      nsteps_override = std::atoi(argv[++a]);
    else if (s == "--quiet")
      quiet = true;
  }
  std::ostream* out = quiet ? nullptr : &std::cout;

  try {
    const auto t0 = std::chrono::steady_clock::now();
    HydroParams P = read_hydro_params(pdir + "/params.txt");
    if (nsteps_override >= 0) P.nsteps = nsteps_override;
    std::filesystem::create_directories(odir);

    const Config cfg = Config::make(P.N);
    SimClock clk;
    Units un;
    GasGrid gas;
    HydroState hs;
    Domains dom;

    init_variables(P.z_init, clk);
    init_domain_gas(cfg, dom, out);
    init_hydro(cfg, gas, hs);
    write_domain_gas_test(odir + "/domain_gas.txt", dom);
    calc_sim_var(cfg, clk, un);

    read_field(pdir + "/u_gas.bin", gas.u);
    read_field(pdir + "/w_gas.bin", gas.w);
    if (std::filesystem::exists(pdir + "/v_gas.bin")) read_field(pdir + "/v_gas.bin", hs.v);

    if (out)
      *out << " Problem    : " << pdir << "\n N          : " << P.N
           << "\n z_init     : " << P.z_init << "\n timestep   : "
           << (P.mode == TimestepMode::Fixed ? "fixed" : "courant")
           << "\n da_sim     : " << P.da_sim << "\n nsteps     : " << P.nsteps << "\n";

    std::string sdir = step_dir(odir, clk.istep);
    calc_sim_var(cfg, clk, un);
    write_summary(sdir, summarize_gas(cfg, un, clk, gas, hs, out));
    write_step_output(sdir, clk, un, gas, hs);

    double thy = 0;
    for (int step = 1; step <= P.nsteps; ++step) {
      const auto t1 = std::chrono::steady_clock::now();
      calc_sim_var(cfg, clk, un);
      calc_timestep_hydro(un, P.mode, P.da_sim, /*amax_dm=*/0.0, clk, hs, out);
      radhydro_hydro_step(cfg, clk, un, gas, hs);
      const GasSummary sum = summarize_gas(cfg, un, clk, gas, hs, out);
      const double dt =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
      thy += dt;
      sdir = step_dir(odir, clk.istep);
      write_summary(sdir, sum);
      write_step_output(sdir, clk, un, gas, hs);
      if (out) {
        char buf[128];
        std::snprintf(buf, sizeof buf, " hydro_driver: step %5d wall %10.4f s\n", clk.istep, dt);
        *out << buf;
      }
    }
    const HydroCounters& c = hs.counters;
    const double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf(
        " hydro_driver: %s: %d steps, hydro %.3f s, total %.3f s\n"
        "   counters: smooth_cells_1 %lld, smooth_cells_2 %lld, c_limit branches %lld, "
        "freeze %lld, check_temp %lld\n"
        "   euler pencils %lld (subcycled %lld, max nt %lld), advect pencils %lld "
        "(subcycled %lld, max nt %lld)\n",
        pdir.c_str(), P.nsteps, thy, secs, static_cast<long long>(c.smooth_cells_1),
        static_cast<long long>(c.smooth_cells_2), static_cast<long long>(c.climit_branch),
        static_cast<long long>(c.freeze_smooth), static_cast<long long>(c.check_temp_fix),
        static_cast<long long>(c.euler_pencils), static_cast<long long>(c.euler_subcycled),
        static_cast<long long>(c.euler_max_nt), static_cast<long long>(c.advect_pencils),
        static_cast<long long>(c.advect_subcycled), static_cast<long long>(c.advect_max_nt));
  } catch (const std::exception& e) {
    std::fprintf(stderr, "hydro_driver: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
