// Full radiation-hydrodynamics driver: the C++ counterpart of the original
// radhydro.f90 main program (and of reference/harness_full, the gfortran
// build of it). See docs/full_run.md.
//
// usage: radhydro <problem_dir> <outdir> [--set key=value ...] [--quiet]
//
// Inputs (problem_dir):
//   params.txt  `key = value` lines (see RunParams below); --set overrides
//   xv_dm.init, u_gas.init   initial conditions (tools/ics/make_ics.py;
//                            directory `ics`, default problem_dir)
//   zdata1.txt, zdata2.txt   output redshifts, read as init_io does (the
//                            params zdata1/zdata2 override them)
//   Halo/       legacy catalogs + sfr.txt (source_mode = legacy)
//
// Sequence (radhydro.f90):
//   init_variables, init_domain_gas, init_domain_rt, init_hydro, init_pm,
//   init_rt, init_ray_tracing, init_io (zdata), init_fft, calc_sim_var;
//   read_init_conditions; calc_sim_var; summarize_dm/gas/rt;
//   do: calc_sim_var; calc_timestep (+ data flags); radhydro (transposes and
//       forward sweeps, half clock step, particle_mesh, reverse sweeps and
//       transposes, half clock step, [idata1 > 1: sources (read_halo or
//       on-the-fly FoF), calc_halo, radiative_transfer], summaries);
//       write_data_2; calc_power_spectrum + write_data_1;
//       exit when redshift <= zdata1(Ndata1).
//
// Outputs (outdir): step_%04d/stats.txt every step (step_0000 = initial
// state); the fields u_gas, w_gas, v_gas, xv_dm and the rt_grid fields at
// the steps selected by fields_every / fields_at_data1 and at the last step;
// Log/hydro.log (the summarize_gas unit-100 records, byte-identical to the
// Fortran); Pk/pk.<ztag>; write_data_1/2 with the F1 fix (DM/xv_dm,
// Gas/{u,v,w}_gas, RT/rt_grid, Ray/{ray,ray_index,htoc_ray} .<ztag>, see
// write_data_files and docs/full_run.md F1); timing.txt.
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
#include <vector>

#include "radray/checkpoint.hpp"
#include "radray/cosmology.hpp"
#include "radray/fft.hpp"
#include "radray/halo_sources.hpp"
#include "radray/hydro.hpp"
#include "radray/io.hpp"
#include "radray/parallel.hpp"
#include "radray/pm.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/rt.hpp"
#include "radray/simulation.hpp"

namespace {

using namespace radray;
using Clock = std::chrono::steady_clock;

std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const auto b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

double to_double(std::string v) {
  for (auto& ch : v)
    if (ch == 'd' || ch == 'D') ch = 'e';
  return std::stod(v);
}

struct RunParams {
  // radhydro_global.f90 sizes (compile-time in the original and in the
  // Fortran harness: -DNGRID, -DNDM, -DNRAY_MAX_FACTOR, -DBOX)
  int N = 32;               // Ngas = Nrt
  int Ndm = -1;             // default N/2
  double box = 5;           // Mpc/h
  int nray_max_factor = 64;
  // cosmo_parameters.f90 zinit (300 in the original)
  double z_init = 300;
  // Output schedule: zdata1 / zdata2 lists (else <problem>/zdata1.txt ...)
  std::string zdata1, zdata2;
  // Switches
  int legacy_rt = 1;
  int legacy_physics = 1;
  double c_reduced = 1.0;   // reduced speed of light (calc_timestep cscale)
  int deterministic_parallel = 0;  // parallel ray tracing (docs/parallel.md)
  int threads = 0;          // OpenMP threads (0: OpenMP default)
  double ra_max = -1;       // < 0: the original (real(4) 0.1); legacy_rt = 0 only
  int rt_substep_carry = 0;  // 1: exact total RT substeps (Config::rt_substep_carry)
  // Restart (docs/production_run.md): checkpoint every n steps (0: never),
  // keep only the latest one (1) or all (0)
  int restart_every = 0;
  int restart_keep_latest = 1;
  // Movie frames independent of the hydro step: every movie_dz in redshift,
  // starting at movie_z0 (default: the first zdata1 redshift, where RT
  // starts), written from inside the RT substeps (0: off)
  double movie_dz = 0;
  double movie_z0 = -1;
  std::string movie_fields = "fHI,temp";
  // Output volume: the step_* field files to write (comma list; empty: all
  // 16), and the write_data_1/2 files: all (F1 fix, default), norays
  // (without Ray/, which holds every ray: 96 B each) or dm (only xv_dm, as
  // the original with one thread)
  std::string fields;
  std::string data_files = "all";
  // Sources
  HaloSourceParams halo;
  std::string ics;          // directory of xv_dm.init/u_gas.init
  // Output and control
  int fields_every = 0;     // write fields every n steps (0: never)
  int fields_at_data1 = 1;  // write fields at the zdata1 output steps
  int max_steps = 1000;     // the original's timing(1000) array
};

void set_run_param(RunParams& P, const std::string& key_in, const std::string& val_in) {
  const std::string key = trim(key_in), val = trim(val_in);
  if (key == "N")
    P.N = std::stoi(val);
  else if (key == "Ndm")
    P.Ndm = std::stoi(val);
  else if (key == "box")
    P.box = to_double(val);
  else if (key == "nray_max_factor")
    P.nray_max_factor = std::stoi(val);
  else if (key == "z_init")
    P.z_init = to_double(val);
  else if (key == "zdata1")
    P.zdata1 = val;
  else if (key == "zdata2")
    P.zdata2 = val;
  else if (key == "legacy_rt")
    P.legacy_rt = std::stoi(val);
  else if (key == "legacy_physics")
    P.legacy_physics = std::stoi(val);
  else if (key == "ra_max")
    P.ra_max = to_double(val);
  else if (key == "rt_substep_carry")
    P.rt_substep_carry = std::stoi(val);
  else if (key == "restart_every")
    P.restart_every = std::stoi(val);
  else if (key == "restart_keep_latest")
    P.restart_keep_latest = std::stoi(val);
  else if (key == "movie_dz")
    P.movie_dz = to_double(val);
  else if (key == "movie_z0")
    P.movie_z0 = to_double(val);
  else if (key == "movie_fields")
    P.movie_fields = val;
  else if (key == "fields")
    P.fields = val;
  else if (key == "data_files") {
    if (val != "all" && val != "norays" && val != "dm")
      throw std::runtime_error("data_files must be all, norays or dm");
    P.data_files = val;
  }
  else if (key == "deterministic_parallel") {
    P.deterministic_parallel = std::stoi(val);
    if (P.deterministic_parallel != 0 && P.deterministic_parallel != 1)
      throw std::runtime_error("deterministic_parallel must be 0 or 1");
  } else if (key == "threads")
    P.threads = std::stoi(val);
  else if (key == "c_reduced") {
    P.c_reduced = to_double(val);
    if (!(P.c_reduced > 0)) throw std::runtime_error("c_reduced must be > 0");
  } else if (key == "ics")
    P.ics = val;
  else if (key == "fields_every")
    P.fields_every = std::stoi(val);
  else if (key == "fields_at_data1")
    P.fields_at_data1 = std::stoi(val);
  else if (key == "max_steps")
    P.max_steps = std::stoi(val);
  else if (!set_halo_param(P.halo, key, val))
    std::cerr << "radhydro: ignoring parameter " << key << "\n";
  if ((key == "legacy_rt" && P.legacy_rt != 0 && P.legacy_rt != 1) ||
      (key == "legacy_physics" && P.legacy_physics != 0 && P.legacy_physics != 1))
    throw std::runtime_error(key + " must be 0 or 1");
}

// io_tools.f90 init_io: `read(un,*,end=1) z`, one value per record (the
// rest of a record is skipped). A params list "z1 z2 ..." (or commas) is
// equivalent to a file with one value per line.
std::vector<double> read_zdata(const std::string& list, const std::string& path,
                               bool required) {
  std::vector<double> z;
  if (!list.empty()) {
    std::string s = list;
    for (auto& ch : s)
      if (ch == ',') ch = ' ';
    std::istringstream in(s);
    std::string tok;
    while (in >> tok) z.push_back(to_double(tok));
    return z;
  }
  std::ifstream in(path);
  if (!in) {
    if (required) throw std::runtime_error("cannot open " + path);
    return z;
  }
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::string tok;
    if (ls >> tok) z.push_back(to_double(tok));
  }
  return z;
}

void write_f32(const std::string& path, const float* p, std::size_t n) {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot write " + path);
  out.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n * sizeof(float)));
}

// io_tools.f90 write_data_1 / write_data_2 with the F1 fix: all the files
// the original writes when Ncpu >= 9 (with Ncpu = 1 it writes only xv_dm,
// docs/full_run.md F1), in the original formats (stream, Fortran layout):
//   data 2: DM/xv_dm, Gas/u_gas, Gas/v_gas, Gas/w_gas
//   data 1: the same plus RT/rt_grid (rt_grid_type, 36 real(4) per cell),
//           and, if Nray > 0, Ray/ray (ray_type, 96 B), Ray/ray_index
//           (int8), Ray/htoc_ray (int8 (4,Nrt^3)).
// Differences: Ray/ray and Ray/ray_index stop at the highest slot ever used
// (the original writes all Nray_max entries; the rest are zeros and the
// identity); D3D/zreion is not written (calc_gas_3 is not ported). The Ray
// files are created empty when Nray = 0, as in the original.
void write_data_files(const std::string& odir, const std::string& tag, const Simulation& s,
                      const std::vector<float>& xv, const std::vector<float>& vgas, bool data1,
                      const std::string& which) {
  write_f32(odir + "/DM/xv_dm." + tag, xv.data(), xv.size());
  if (which == "dm") return;
  std::filesystem::create_directories(odir + "/Gas");
  write_f32(odir + "/Gas/u_gas." + tag, s.gas.u.data(), s.gas.u.size());
  write_f32(odir + "/Gas/v_gas." + tag, vgas.data(), vgas.size());
  write_f32(odir + "/Gas/w_gas." + tag, s.gas.w.data(), s.gas.w.size());
  if (!data1) return;
  const RTGrid& G = s.rt;
  std::filesystem::create_directories(odir + "/RT");
  {
    std::vector<float> buf(static_cast<std::size_t>(36 * G.ncell));
    for (std::int64_t c = 0; c < G.ncell; ++c) {
      float* b = &buf[static_cast<std::size_t>(36 * c)];
      b[0] = G.rhog[c];
      b[1] = G.fHI[c];
      b[2] = G.fHeI[c];
      b[3] = G.fHeII[c];
      b[4] = G.temp[c];
      b[5] = G.Luv[c];
      b[6] = G.sfr[c];
      for (int m = 0; m < 3; ++m) {
        b[7 + m] = G.dx[c * 3 + m];
        b[10 + m] = G.gamma[c * 3 + m];
        b[13 + m] = G.lambda[c * 3 + m];
      }
      for (int m = 0; m < kNfreq; ++m) {
        b[16 + m] = G.nu[c * kNfreq + m];
        b[21 + m] = G.ngam[c * kNfreq + m];
        b[26 + m] = G.fgam[c * kNfreq + m];
        b[31 + m] = G.tau[c * kNfreq + m];
      }
    }
    write_f32(odir + "/RT/rt_grid." + tag, buf.data(), buf.size());
  }
  if (which == "norays") return;
  // the original opens (creates) the Ray files and writes them if Nray > 0
  std::filesystem::create_directories(odir + "/Ray");
  {
    struct FRay {  // radhydro_global.f90 ray_type
      std::int32_t nside, ipix;
      std::int64_t next, prev;
      float r, dr, x[3], dx[3], nu[kNfreq], ngam[kNfreq];
    };
    static_assert(sizeof(FRay) == 96, "Fortran ray_type layout");
    std::ofstream out(odir + "/Ray/ray." + tag, std::ios::binary);
    std::ofstream oix(odir + "/Ray/ray_index." + tag, std::ios::binary);
    if (!out || !oix) throw std::runtime_error("cannot write Ray/ files");
    // deterministic_parallel: storage slots depend on thread timing, so the
    // rays are written renumbered in the canonical order (cells k,j,i, each
    // chain head to tail): new index n for the n-th ray, ray_index = identity
    std::vector<std::int64_t> order, newidx;
    if (s.cfg.det_parallel && s.rays.Nray > 0) {
      newidx.assign(static_cast<std::size_t>(s.rays.hwm()) + 1, 0);
      for (std::int64_t c = 0; c < G.ncell; ++c)
        for (std::int64_t ir = G.htoc_ray(3, c); ir > 0; ir = s.rays.at(ir).next) {
          order.push_back(ir);
          newidx[ir] = static_cast<std::int64_t>(order.size());
        }
    }
    auto remap = [&](std::int64_t ir) { return newidx.empty() || ir <= 0 ? ir : newidx[ir]; };
    const std::int64_t top =
        s.rays.Nray > 0 ? (newidx.empty() ? s.rays.hwm() : static_cast<std::int64_t>(order.size()))
                        : 0;
    for (std::int64_t n = 1; n <= top; ++n) {
      const std::int64_t ir = newidx.empty() ? n : order[n - 1];
      const Ray& r = s.rays.at(ir);
      FRay f{};
      f.nside = r.nside;
      f.ipix = r.ipix;
      f.next = remap(r.next);
      f.prev = remap(r.prev);
      f.r = r.r;
      f.dr = r.dr;
      for (int m = 0; m < 3; ++m) {
        f.x[m] = r.x[m];
        f.dx[m] = r.dx[m];
      }
      for (int m = 0; m < kNfreq; ++m) {
        f.nu[m] = r.nu[m];
        f.ngam[m] = r.ngam[m];
      }
      out.write(reinterpret_cast<const char*>(&f), sizeof f);
      const std::int64_t idx =
          !newidx.empty() ? n
                          : (ir < static_cast<std::int64_t>(s.rays.ray_index_.size())
                                 ? static_cast<std::int64_t>(s.rays.ray_index_[ir])
                                 : ir);
      oix.write(reinterpret_cast<const char*>(&idx), sizeof idx);
    }
    std::ofstream oh(odir + "/Ray/htoc_ray." + tag, std::ios::binary);
    if (s.rays.Nray > 0) {
      std::vector<std::int64_t> h(G.htoc.size());
      for (std::size_t m = 0; m < h.size(); ++m) h[m] = remap(G.htoc[m]);
      oh.write(reinterpret_cast<const char*>(h.data()),
               static_cast<std::streamsize>(h.size() * sizeof(std::int64_t)));
    }
  }
}

// Peak resident set size (VmHWM of /proc/self/status) in MB, -1 if unknown
double peak_rss_mb() {
  std::ifstream in("/proc/self/status");
  std::string line;
  while (std::getline(in, line))
    if (line.rfind("VmHWM:", 0) == 0) return std::stod(line.substr(6)) / 1024.0;
  return -1;
}

struct Timers {
  double hydro = 0, pm = 0, halos = 0, rt = 0, summaries = 0, io = 0, pk = 0;
};

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <problem_dir> <outdir> [--set key=value ...] [--quiet] [--restart <checkpoint>]\n",
                 argv[0]);
    return 2;
  }
  const std::string pdir = argv[1];
  const std::string odir = argv[2];
  std::vector<std::string> sets;
  bool quiet = false;
  std::string restart_from;
  for (int a = 3; a < argc; ++a) {
    const std::string s = argv[a];
    if (s == "--set" && a + 1 < argc)
      sets.emplace_back(argv[++a]);
    else if (s == "--restart" && a + 1 < argc)
      restart_from = argv[++a];
    else if (s == "--quiet")
      quiet = true;
    else {
      std::fprintf(stderr, "radhydro: unknown argument '%s'\n", s.c_str());
      return 2;
    }
  }

  try {
    const auto t0 = Clock::now();
    RunParams P;
    P.halo.dir_halo = pdir + "/Halo/";
    {
      std::ifstream in(pdir + "/params.txt");
      if (!in) throw std::runtime_error("cannot open " + pdir + "/params.txt");
      std::string line;
      while (std::getline(in, line)) {
        const auto h = line.find('#');
        if (h != std::string::npos) line.resize(h);
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        set_run_param(P, line.substr(0, eq), line.substr(eq + 1));
      }
    }
    for (const std::string& kv : sets) {
      const auto eq = kv.find('=');
      if (eq == std::string::npos) throw std::runtime_error("--set needs key=value: " + kv);
      set_run_param(P, kv.substr(0, eq), kv.substr(eq + 1));
    }
    // relative paths of the problem are relative to the problem directory
    auto in_problem = [&](const std::string& p) {
      if (p.empty() || p[0] == '/') return p;
      return pdir + "/" + p;
    };
    const std::string ics_dir = P.ics.empty() ? pdir : in_problem(P.ics);
    if (!P.halo.dir_halo.empty() && P.halo.dir_halo[0] != '/' &&
        P.halo.dir_halo.rfind(pdir, 0) != 0)
      P.halo.dir_halo = in_problem(P.halo.dir_halo);
    if (!P.halo.model.sfr_history_path.empty())
      P.halo.model.sfr_history_path = in_problem(P.halo.model.sfr_history_path);

    std::filesystem::create_directories(odir + "/Log");
    std::filesystem::create_directories(odir + "/Pk");
    std::filesystem::create_directories(odir + "/DM");
    // a restart appends to the logs of the run directory
    const auto logmode = restart_from.empty() ? std::ios::out : std::ios::app;
    std::ofstream raylog(odir + "/Log/ray.log", logmode);
    std::ofstream hydrolog(odir + "/Log/hydro.log", logmode);
    std::ofstream nullout;
    std::ostream* out = quiet ? nullptr : &std::cout;

    // ------------------------------------------------------------------
    // Initialize (radhydro.f90 order)
    // ------------------------------------------------------------------
    auto sim = std::make_unique<Simulation>();
    Simulation& s = *sim;
    s.cfg = Config::make(P.N, P.nray_max_factor, P.Ndm, P.box);
    s.cfg.legacy_rt = (P.legacy_rt == 1);
    s.cfg.legacy_physics = (P.legacy_physics == 1);
    s.cfg.c_reduced = P.c_reduced;
    s.cfg.det_parallel = (P.deterministic_parallel == 1);
    if (s.cfg.det_parallel && s.cfg.legacy_rt)
      throw std::runtime_error("deterministic_parallel = 1 requires legacy_rt = 0");
    par::set_threads(P.threads);
    if (P.ra_max > 0) {
      if (s.cfg.legacy_rt)
        throw std::runtime_error("ra_max can only be changed with legacy_rt = 0");
      if (!(P.ra_max < 1)) throw std::runtime_error("ra_max must be in (0,1)");
      s.cfg.ra_max = P.ra_max;
    }
    if (P.rt_substep_carry) {
      if (s.cfg.legacy_rt) throw std::runtime_error("rt_substep_carry requires legacy_rt = 0");
      s.cfg.rt_substep_carry = true;
    }
    s.raylog = &raylog;
    s.out = out;
    const Config& cfg = s.cfg;

    Domains gdom;
    HydroState hs;
    pm::PMState pm;

    init_variables(P.z_init, s.clk);                  // init_variables
    init_domain_gas(cfg, gdom);                        // init_domain_gas
    init_domain_rt(cfg, s.dom);                        // init_domain_rt
    init_hydro(cfg, s.gas, hs);                        // init_hydro
    pm::init_pm(pm, pm::PMConfig::make(cfg), gdom);    // init_pm
    init_rt(s);                                        // init_rt (rates)
    init_ray_tracing(s);                               // init_ray_tracing
    // init_io: output redshifts
    const std::vector<double> zdata1 = read_zdata(P.zdata1, pdir + "/zdata1.txt", true);
    const std::vector<double> zdata2 = read_zdata(P.zdata2, pdir + "/zdata2.txt", false);
    if (zdata1.empty()) throw std::runtime_error("no zdata1 output redshifts");
    std::vector<double> adata1, adata2;
    for (double z : zdata1) adata1.push_back(1 / (1 + z));
    for (double z : zdata2) adata2.push_back(1 / (1 + z));
    const int Ndata1 = static_cast<int>(zdata1.size());
    const int Ndata2 = static_cast<int>(zdata2.size());
    DriverState dst;
    int& idata1 = dst.idata1;
    int& idata2 = dst.idata2;
    bool& data1_out = dst.data1_out;
    bool& data2_out = dst.data2_out;
    fft::fft_3d('f', cfg.Ngas, cfg.Ngas, cfg.Ngas, pm.density.data());  // init_fft
    calc_sim_var(cfg, s.clk, s.units);

    if (restart_from.empty()) pm::read_init_conditions(ics_dir, pm, s.gas.u, s.gas.w, s.units);

    if (out)
      *out << " radhydro: N = " << cfg.Ngas << ", Ndm = " << cfg.Ndm << ", box = " << cfg.box
           << " Mpc/h, z_init = " << P.z_init << ", Nray_max = " << cfg.Nray_max
           << "\n legacy_rt = " << P.legacy_rt << ", legacy_physics = " << P.legacy_physics << ", c_reduced = " << P.c_reduced
           << ", deterministic_parallel = " << P.deterministic_parallel
           << ", threads = " << par::max_threads() << ", sources = "
           << (P.halo.mode == HaloSourceParams::Mode::Legacy ? "legacy" : "onthefly")
           << ", zdata1: " << Ndata1 << " redshifts down to " << zdata1.back() << "\n";

    Timers T;
    std::int64_t& Nhalo = dst.Nhalo;
    auto summaries = [&]() {
      const auto ts = Clock::now();
      pm::summarize_dm(pm, s.units, out);
      const GasSummary gsum = summarize_gas(cfg, s.units, s.clk, s.gas, hs, out);
      hydrolog << format_gas_record(gsum.record) << "\n";
      summarize_rt(s);
      T.summaries += std::chrono::duration<double>(Clock::now() - ts).count();
    };

    auto write_stats = [&](const std::string& d) {
      // global diagnostics over the RT grid, k,j,i order, real(8) sums
      double xv = 0, xm = 0, mass = 0, tv = 0, tm = 0, xhe3 = 0;
      const RTGrid& G = s.rt;
      for (std::int64_t c = 0; c < G.ncell; ++c) {
        const double fHI = G.fHI[c], rho = G.rhog[c], T_ = G.temp[c];
        xv = xv + (1 - fHI);
        xm = xm + rho * (1 - fHI);
        mass = mass + rho;
        tv = tv + T_;
        tm = tm + rho * T_;
        xhe3 = xhe3 + (1 - static_cast<double>(G.fHeI[c]) - static_cast<double>(G.fHeII[c]));
      }
      const double nc = cfg.nrtcr;
      std::FILE* fp = std::fopen((d + "/stats.txt").c_str(), "w");
      if (!fp) throw std::runtime_error("cannot write stats.txt in " + d);
      auto i1 = [&](const char* k, long long x) { std::fprintf(fp, "%s = %lld\n", k, x); };
      auto d1 = [&](const char* k, double x) { std::fprintf(fp, "%s = %.17g\n", k, x); };
      const RayStat& rv = s.rays.save;
      i1("istep", s.clk.istep);
      i1("idata1", idata1);
      i1("idata2", idata2);
      i1("nstep_rt", s.clk.nstep_rt);
      i1("Nhalo", Nhalo);
      i1("Nsource", s.Nsource);
      i1("nray_active", s.rays.Nray - rv.nfree);
      i1("ncst1", rv.ncst1);
      i1("nmrg1", rv.nmrg1);
      i1("nspl1", rv.nspl1);
      i1("ndel1", rv.ndel1);
      d1("redshift", s.clk.redshift);
      d1("scalefactor", s.clk.scalefactor);
      d1("year", s.clk.year);
      d1("da_sim", s.clk.da_sim);
      d1("dt_sim", s.clk.dt_sim);
      d1("dt_yr", s.clk.dt_yr);
      d1("cmax_gas", s.clk.cmax_gas);
      d1("vmax_gas", s.clk.vmax_gas);
      d1("vrms_gas", s.clk.vrms_gas);
      d1("c_limit", hs.c_limit);
      d1("amax_dm", s.clk.amax_dm);
      d1("arms_dm", s.clk.arms_dm);
      d1("vmax_dm", s.clk.vmax_dm);
      d1("vrms_dm", s.clk.vrms_dm);
      d1("ngam_rays_total", ngam_rays_total(s));
      d1("Luv", s.rt_save.Luv);
      d1("sfr", s.rt_save.sfr);
      d1("xHII_V", xv / nc);
      d1("xHII_M", xm / mass);
      d1("xHeIII_V", xhe3 / nc);
      d1("T_V", tv / nc);
      d1("T_M", tm / mass);
      std::fclose(fp);
    };

    auto write_fields = [&](const std::string& d) {
      const RTGrid& G = s.rt;
      auto W = [&](const char* name, const std::vector<float>& v) {
        if (!P.fields.empty() && ("," + P.fields + ",").find(std::string(",") + name + ",") ==
                                     std::string::npos)
          return;
        write_f32(d + "/" + name + ".bin", v.data(), v.size());
      };
      W("fHI", G.fHI);
      W("fHeI", G.fHeI);
      W("fHeII", G.fHeII);
      W("temp", G.temp);
      W("rhog", G.rhog);
      W("ngam", G.ngam);
      W("nu", G.nu);
      W("fgam", G.fgam);
      W("tau", G.tau);
      W("gamma", G.gamma);
      W("lambda", G.lambda);
      W("dx", G.dx);
      W("u_gas", s.gas.u);
      W("w_gas", s.gas.w);
      W("v_gas", hs.v);
      W("xv_dm", pm.xv);
    };

    auto step_dir = [&]() {
      char name[32];
      std::snprintf(name, sizeof name, "step_%04d", s.clk.istep);
      const std::string d = odir + "/" + name;
      std::filesystem::create_directories(d);
      return d;
    };

    // Movie frames (movie_dz > 0): selected rt_grid fields every movie_dz
    // in redshift, from inside the RT substeps (Simulation::substep_hook),
    // or at the end of a hydro step without RT.
    std::vector<std::string> mfields;
    {
      std::string f = P.movie_fields;
      for (auto& ch : f)
        if (ch == ',') ch = ' ';
      std::istringstream in(f);
      std::string tok;
      while (in >> tok) mfields.push_back(tok);
    }
    std::ofstream movielog;
    auto write_frame = [&](const Simulation& ss) {
      const RTGrid& G = ss.rt;
      char name[64];
      for (const std::string& fld : mfields) {
        const std::vector<float>* v = nullptr;
        if (fld == "fHI") v = &G.fHI;
        else if (fld == "fHeI") v = &G.fHeI;
        else if (fld == "fHeII") v = &G.fHeII;
        else if (fld == "temp") v = &G.temp;
        else if (fld == "rhog") v = &G.rhog;
        else if (fld == "ngam") v = &G.ngam;
        else if (fld == "gamma") v = &G.gamma;
        else throw std::runtime_error("movie_fields: unknown field " + fld);
        std::snprintf(name, sizeof name, "/movie/frame_%05d_%s.bin", dst.movie_frame, fld.c_str());
        write_f32(odir + name, v->data(), v->size());
      }
      double xv = 0, xm = 0, mass = 0;
      for (std::int64_t c = 0; c < G.ncell; ++c) {
        xv = xv + (1 - static_cast<double>(G.fHI[c]));
        xm = xm + static_cast<double>(G.rhog[c]) * (1 - static_cast<double>(G.fHI[c]));
        mass = mass + static_cast<double>(G.rhog[c]);
      }
      char line[256];
      std::snprintf(line, sizeof line, "%6d %6d %6d %.8f %.8e %.8e\n", dst.movie_frame,
                    ss.clk.istep, ss.clk.istep_rt, ss.clk.redshift, xv / cfg.nrtcr, xm / mass);
      movielog << line << std::flush;
      dst.movie_frame += 1;
      dst.movie_next -= P.movie_dz;
    };
    if (P.movie_dz > 0) {
      std::filesystem::create_directories(odir + "/movie");
      movielog.open(odir + "/movie/frames.txt", logmode);
      if (restart_from.empty()) {
        movielog << "# frame istep istep_rt redshift xHII_V xHII_M (frame_NNNNN_<field>.bin: "
                    "float32 (N,N,N) rt_grid fields)\n";
        dst.movie_next = P.movie_z0 > 0 ? P.movie_z0 : zdata1.front();
      }
      s.substep_hook = [&](Simulation& ss) {
        while (ss.clk.redshift <= dst.movie_next) write_frame(ss);
      };
    }

    // Restart: the complete state of the end of a step (checkpoint.hpp)
    if (!restart_from.empty()) {
      read_checkpoint(restart_from, s, hs, pm, dst);
      if (out) *out << " radhydro: restarted from " << restart_from << " at step " << s.clk.istep << "\n";
    } else {
      // Summarize the initial state
      calc_sim_var(cfg, s.clk, s.units);
      summaries();
      const std::string d = step_dir();
      write_stats(d);
      if (P.fields_every > 0) write_fields(d);
    }

    std::ofstream timing(odir + "/timing.txt", logmode);
    if (restart_from.empty())
      timing << "# istep redshift nstep_rt Nsource nray_active t_hydro t_pm t_halos t_rt "
                "t_summaries t_io t_total [s, this step]\n";

    // ------------------------------------------------------------------
    // Time stepping
    // ------------------------------------------------------------------
    std::vector<double> pk_dm, pk_gas;
    for (;;) {
      if (s.clk.istep >= P.max_steps)
        throw std::runtime_error("max_steps reached before the last zdata1 redshift");
      const Timers T0 = T;
      const auto ts0 = Clock::now();

      calc_sim_var(cfg, s.clk, s.units);
      const TimestepInfo ti =
          calc_timestep(cfg, s.units, TimestepMode::Original, 0.0, s.clk, out);
      // Set IO flags (calc_timestep)
      if (ti.a2 >= adata1[idata1 - 1]) data1_out = true;
      if (idata2 <= Ndata2)
        if (ti.a2 >= adata2[idata2 - 1]) data2_out = true;

      // radhydro: forward hydro step with transposes in between
      auto th = Clock::now();
      transpose_ij(s.gas, hs);  // x,y,z -> y,x,z
      transpose_ik(s.gas, hs);  // y,x,z -> z,x,y
      hydro_step(cfg, s.units, s.clk, SweepOrder::Forward, s.gas, hs);
      transpose_ik(s.gas, hs);  // z,x,y -> y,x,z
      hydro_step(cfg, s.units, s.clk, SweepOrder::Forward, s.gas, hs);
      transpose_ij(s.gas, hs);  // y,x,z -> x,y,z
      hydro_step(cfg, s.units, s.clk, SweepOrder::Forward, s.gas, hs);
      T.hydro += std::chrono::duration<double>(Clock::now() - th).count();

      s.clk.year = s.clk.year + s.clk.dt_yr / 2;
      s.clk.scalefactor = s.clk.scalefactor + s.clk.da_sim / 2;
      s.clk.redshift = 1 / s.clk.scalefactor - 1;
      calc_sim_var(cfg, s.clk, s.units);

      // Gravity for gas and dark matter (gas ordered x,y,z, same as DM)
      th = Clock::now();
      pm::particle_mesh(pm, s.gas.u.data(), hs.v.data(), s.clk, s.units);
      T.pm += std::chrono::duration<double>(Clock::now() - th).count();

      // Reverse hydro step with transposes in between
      th = Clock::now();
      hydro_step(cfg, s.units, s.clk, SweepOrder::Reverse, s.gas, hs);
      transpose_ij(s.gas, hs);  // x,y,z -> y,x,z
      hydro_step(cfg, s.units, s.clk, SweepOrder::Reverse, s.gas, hs);
      transpose_ik(s.gas, hs);  // y,x,z -> z,x,y
      hydro_step(cfg, s.units, s.clk, SweepOrder::Reverse, s.gas, hs);
      transpose_ik(s.gas, hs);  // z,x,y -> y,x,z
      transpose_ij(s.gas, hs);  // y,x,z -> x,y,z
      T.hydro += std::chrono::duration<double>(Clock::now() - th).count();

      s.clk.year = s.clk.year + s.clk.dt_yr / 2;
      s.clk.scalefactor = s.clk.scalefactor + s.clk.da_sim / 2;
      s.clk.redshift = 1 / s.clk.scalefactor - 1;
      calc_sim_var(cfg, s.clk, s.units);

      // Ray tracing radiative transfer (after the first zdata1 output)
      if (idata1 > 1) {
        th = Clock::now();
        const SourceCatalog cat = update_sources(s, pm.xv.data(), pm.pc.Npart, P.halo, out);
        Nhalo = static_cast<std::int64_t>(cat.halo.size());
        T.halos += std::chrono::duration<double>(Clock::now() - th).count();
        if (cat.Nsource > 0) {
          th = Clock::now();
          radiative_transfer(s);
          T.rt += std::chrono::duration<double>(Clock::now() - th).count();
        }
      }

      // Summarize variables
      summaries();

      // IO
      th = Clock::now();
      const std::string d = step_dir();
      bool fields = (P.fields_every > 0 && s.clk.istep % P.fields_every == 0);
      if (data2_out) {  // write_data_2 (F1 fix: all files)
        write_data_files(odir, pm::ztag(s.clk.redshift), s, pm.xv, hs.v, false, P.data_files);
        idata2 = idata2 + 1;
        data2_out = false;
      }
      if (data1_out) {
        const auto tp = Clock::now();
        pm::calc_power_spectrum(pm, s.gas.u.data(), cfg.box, pk_dm, pk_gas);
        pm::write_power_spectrum(odir + "/Pk/pk." + pm::ztag(s.clk.redshift), pm, pk_dm,
                                 pk_gas);
        T.pk += std::chrono::duration<double>(Clock::now() - tp).count();
        // write_data_1 (F1 fix: all files)
        write_data_files(odir, pm::ztag(s.clk.redshift), s, pm.xv, hs.v, true, P.data_files);
        idata1 = idata1 + 1;
        data1_out = false;
        if (P.fields_at_data1) fields = true;
      }
      const bool last = s.clk.redshift <= zdata1[Ndata1 - 1];
      write_stats(d);
      if (fields || last) write_fields(d);
      T.io += std::chrono::duration<double>(Clock::now() - th).count();

      const double tstep = std::chrono::duration<double>(Clock::now() - ts0).count();
      {
        char buf[512];
        std::snprintf(buf, sizeof buf,
                      "%5d %10.6f %6d %6lld %9lld %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",
                      s.clk.istep, s.clk.redshift, s.clk.nstep_rt,
                      static_cast<long long>(s.Nsource),
                      static_cast<long long>(s.rays.Nray - s.rays.save.nfree),
                      T.hydro - T0.hydro, T.pm - T0.pm, T.halos - T0.halos, T.rt - T0.rt,
                      T.summaries - T0.summaries, T.io - T0.io, tstep);
        timing << buf << std::flush;
        if (quiet) {
          std::printf(" step %4d  z = %9.5f  nstep_rt = %5d  Nsource = %5lld  %.3f s\n",
                      s.clk.istep, s.clk.redshift, s.clk.nstep_rt,
                      static_cast<long long>(s.Nsource), tstep);
          std::fflush(stdout);
        }
      }

      // Frames not reached inside RT substeps (steps without RT)
      if (P.movie_dz > 0)
        while (s.clk.redshift <= dst.movie_next) write_frame(s);

      // Checkpoint (complete state at the end of this step)
      if (P.restart_every > 0 && (s.clk.istep % P.restart_every == 0 || last)) {
        std::filesystem::create_directories(odir + "/restart");
        char name[64];
        if (P.restart_keep_latest)
          std::snprintf(name, sizeof name, "/restart/checkpoint.bin");
        else
          std::snprintf(name, sizeof name, "/restart/checkpoint_%04d.bin", s.clk.istep);
        write_checkpoint(odir + name, s, hs, pm, dst);
      }

      // Stop
      if (last) break;
    }

    const double total = std::chrono::duration<double>(Clock::now() - t0).count();
    std::printf(
        " radhydro: %d steps, total %.2f s: hydro %.2f, pm %.2f (+ P(k) %.2f), halos %.2f, "
        "rt %.2f, summaries %.2f, io %.2f; peak RSS %.1f MB\n",
        s.clk.istep, total, T.hydro, T.pm, T.pk, T.halos, T.rt, T.summaries, T.io,
        peak_rss_mb());
    std::printf(" %s\n ray pool guard: split limited in %lld substeps, merging forced in %lld\n",
                format_rt_timers(s.timers).c_str(), static_cast<long long>(s.rays.npool_limited),
                static_cast<long long>(s.rays.npool_forced));
    timing << "# " << format_rt_timers(s.timers) << "\n";
    timing << "# total " << total << " hydro " << T.hydro << " pm " << T.pm << " pk " << T.pk
           << " halos " << T.halos << " rt " << T.rt << " summaries " << T.summaries << " io "
           << T.io << " peak_rss_MB " << peak_rss_mb() << "\n";
  } catch (const std::exception& e) {
    std::fprintf(stderr, "radhydro: error: %s (peak RSS %.1f MB)\n", e.what(), peak_rss_mb());
    return 1;
  }
  return 0;
}
