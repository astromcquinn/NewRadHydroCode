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
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
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
#include "radray/hydro_unsplit.hpp"
#include "radray/io.hpp"
#include "radray/parallel.hpp"
#include "radray/plane_parallel.hpp"
#include "radray/pm.hpp"
#include "radray/pm_backend.hpp"
#include "radray/ray_dump.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/rt.hpp"
#include "radray/amrex_state.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/hydro_backend.hpp"
#include "radray/rt_backend.hpp"
#include "radray/recomb_rad.hpp"
#include "radray/thermal_history.hpp"
#include "radray/xray_background.hpp"
#include "radray/rt_moments.hpp"
#include "radray/rt_tensor_diag.hpp"
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
  int pm_refine = 1;        // PM mesh = pm_refine*N (docs/pm_refine.md); 1 = original
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
  // c_reduced_schedule = c0 x1 c1 x2 c2 ... (new, default empty = off,
  // docs/production_run.md section 1): the reduced speed of light by
  // ionisation state. Each hydro step uses c_k for xk <= x < x(k+1), with
  // x = x_HII,V of the box at the end of the previous step (c0 below x1,
  // and before the RT has run: x = 0 while the RT grid is unwritten).
  // Overrides c_reduced. Stateless, so a restart reproduces it.
  // c_reduced_schedule_interp = 1 (default 0, the steps above): c is
  // continuous in x, linear between the midpoints of the plateaus
  // ((x_k + x_{k+1})/2 with x_0 = 0 and x_{n+1} = 1, value c_k there) and
  // constant below the first and above the last midpoint, so the mean of
  // the step schedule is kept but the photon reservoir n_gamma ~ ndot
  // lambda / c is never emptied in one step (the Gamma spikes of
  // docs/rt_suite_512.md section 4).
  std::vector<double> c_sched_c, c_sched_x;
  int c_sched_interp = 0;
  int deterministic_parallel = 0;  // parallel ray tracing (docs/parallel.md)
  int threads = 0;          // OpenMP threads (0: OpenMP default)
  double ra_max = -1;       // < 0: the original (real(4) 0.1); legacy_rt = 0 only
  int rt_substep_carry = 0;  // 1: exact total RT substeps (Config::rt_substep_carry)
  int ray_rotate = 0;        // B16 fix: Config::ray_rotate (0 = original)
  long long ray_rotate_seed = 0;
  int redshift_per_distance = 0;  // Config::redshift_per_distance (0 = original)
  int rt_cells_per_substep = 1;   // Config::rt_cells_per_substep (1 = original)
  int photon_conserve = -1;       // Config::photon_conserve (B17 fix): -1 = !legacy_physics, else 0 | 1
  int ray_deposit = 0;            // Config::ray_deposit (0 = share, the original; 1 = exposure)
  double cool_suppress_overdensity = 200.0;  // Config::cool_suppress_overdensity (200 = original; 0 = off)
  // thermal_init = floor | recfast (Config::thermal_init, docs/thermal_history.md)
  // and thermal_floor in K (< 0: 1 K with floor, 0.1 K with recfast)
  int thermal_init = 0;
  double thermal_floor = -1;
  // Restart (docs/production_run.md): checkpoint every n steps (0: never),
  // keep only the latest one (1) or all (0)
  int restart_every = 0;
  int restart_keep_latest = 1;
  // > 0: numbered checkpoints (checkpoint_NNNN.bin, as restart_keep_latest = 0),
  // of which only the newest restart_keep are kept (0: as restart_keep_latest)
  int restart_keep = 0;
  // Movie frames independent of the hydro step: every movie_dz in redshift,
  // starting at movie_z0 (default: the first zdata1 redshift, where RT
  // starts), written from inside the RT substeps (0: off)
  double movie_dz = 0;
  double movie_z0 = -1;
  std::string movie_fields = "fHI,temp";
  // Slab-only movie frames (docs/movies.md): 0 (default) writes the full
  // cube of each field; t > 0 writes only the slabs of t cells centred on
  // the box centre (first cell Nrt/2 - t/2) normal to each axis in
  // movie_slab_axes, to movie/frame_NNNNN_<field>.slab, geometry in
  // movie/slab.txt
  int movie_slab = 0;
  std::string movie_slab_axes = "xyz";
  // movie_every_substep = 1: a frame after every RT substep (movie_dz > 0
  // still turns the movie on; for the tests of the substep-averaged fields)
  int movie_every_substep = 0;
  // Output volume: the step_* field files to write (comma list; empty: all
  // 16), and the write_data_1/2 files: all (F1 fix, default), norays
  // (without Ray/, which holds every ray: 96 B each), nodm (norays without
  // DM/xv_dm: the particles stay only in the checkpoints), dm (only xv_dm,
  // as the original with one thread) or gas (only Gas/u_gas)
  std::string fields;
  std::string data_files = "all";
  // Sources
  HaloSourceParams halo;
  // X-ray background heating (xray_*; docs/xray_heating.md); off by default
  XrayParams xray;
  std::string ics;          // directory of xv_dm.init/u_gas.init
  int ics_check = 1;        // 1: z_init and the cosmology must match <ics>/ics_info.txt
  // Output and control
  int fields_every = 0;     // write fields every n steps (0: never)
  // fields_levels = 1: with every field output, one plane at the finest
  // resolution of the RT levels (step_NNNN/levels_slice.bin and .txt; the
  // refined runs' field files are the level-0 composite), normal to
  // fields_levels_axis (x|y|z) through base cell fields_levels_index
  // (-1: Nrt/2), the middle fine plane of that cell
  int fields_levels = 0;
  int fields_levels_axis = 2;
  int fields_levels_index = -1;
  int fields_at_data1 = 1;  // write fields at the zdata1 output steps
  // tensor_dump = 1: the Eddington-tensor diagnostic at the zdata1 output
  // steps (rt_tensor_diag.hpp; tensor.bin, otvet_h.bin, sources.txt)
  int tensor_dump = 0;
  int max_steps = 1000;     // the original's timing(1000) array
  int run_steps = 0;        // > 0: stop after this many steps of this invocation
  int step_timers = 0;      // 1: print the backends' cumulative phase timers after every step (rank 0)
                            // (timing runs from a checkpoint; 0: run to the end)
  // hydro_scheme = legacy | unsplit and the unsplit parameters
  // (docs/hydro_unsplit.md), parsed by set_hydro_scheme_param
  Config hydro;
  // rt_method = rays | m1 | otvet and the moment-solver parameters
  // (docs/rt_moments.md), parsed by set_rt_method_param
  Config rtm;
  // Ray snapshots for the compression study (radray/ray_dump.hpp); off by
  // default, no effect on results.
  RayDumpConfig ray_dump;
  // Lossy ray-precision experiment (radray/ray_quant.hpp); research only.
  RayQuantConfig ray_quant;
};

void set_run_param(RunParams& P, const std::string& key_in, const std::string& val_in) {
  const std::string key = trim(key_in), val = trim(val_in);
  if (key == "N")
    P.N = std::stoi(val);
  else if (key == "Ndm")
    P.Ndm = std::stoi(val);
  else if (key == "pm_refine")
    P.pm_refine = std::stoi(val);
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
  else if (key == "ray_rotate") {
    P.ray_rotate = std::stoi(val);
    if (P.ray_rotate < 0 || P.ray_rotate > 3) throw std::runtime_error("ray_rotate must be 0..3");
  } else if (key == "ray_rotate_seed")
    P.ray_rotate_seed = std::stoll(val);
  else if (key == "rt_cells_per_substep") {
    P.rt_cells_per_substep = std::stoi(val);
    if (P.rt_cells_per_substep < 1) throw std::runtime_error("rt_cells_per_substep must be >= 1");
  }
  else if (key == "ray_deposit") {
    if (val == "share" || val == "0")
      P.ray_deposit = 0;
    else if (val == "exposure" || val == "1")
      P.ray_deposit = 1;
    else
      throw std::runtime_error("ray_deposit must be share or exposure");
  }
  else if (key == "cool_suppress_overdensity") {
    P.cool_suppress_overdensity = std::stod(val);
    if (!(P.cool_suppress_overdensity >= 0)) throw std::runtime_error("cool_suppress_overdensity must be >= 0 (0 = off)");
  }
  else if (key == "thermal_init") {
    if (val == "floor")
      P.thermal_init = 0;
    else if (val == "recfast")
      P.thermal_init = 1;
    else
      throw std::runtime_error("thermal_init must be floor or recfast");
  } else if (key == "thermal_floor") {
    P.thermal_floor = to_double(val);
    if (!(P.thermal_floor > 0)) throw std::runtime_error("thermal_floor must be > 0 (K)");
  } else if (key == "photon_conserve") {
    P.photon_conserve = std::stoi(val);  // -1: follow legacy_physics (the default)
    if (P.photon_conserve < -1 || P.photon_conserve > 1)
      throw std::runtime_error("photon_conserve must be 0 or 1");
  }
  else if (key == "redshift_per_distance") {
    P.redshift_per_distance = std::stoi(val);
    if (P.redshift_per_distance != 0 && P.redshift_per_distance != 1)
      throw std::runtime_error("redshift_per_distance must be 0 or 1");
  }
  else if (key == "restart_every")
    P.restart_every = std::stoi(val);
  else if (key == "restart_keep_latest")
    P.restart_keep_latest = std::stoi(val);
  else if (key == "restart_keep") {
    P.restart_keep = std::stoi(val);
    if (P.restart_keep < 0) throw std::runtime_error("restart_keep must be >= 0");
  }
  else if (key == "movie_dz")
    P.movie_dz = to_double(val);
  else if (key == "movie_z0")
    P.movie_z0 = to_double(val);
  else if (key == "movie_fields")
    P.movie_fields = val;
  else if (key == "movie_every_substep")
    P.movie_every_substep = std::stoi(val);
  else if (key == "movie_slab") {
    P.movie_slab = std::stoi(val);
    if (P.movie_slab < 0) throw std::runtime_error("movie_slab must be >= 0");
  }
  else if (key == "movie_slab_axes") {
    if (val.empty()) throw std::runtime_error("movie_slab_axes: empty");
    for (std::size_t i = 0; i < val.size(); ++i)
      if ((val[i] != 'x' && val[i] != 'y' && val[i] != 'z') || val.find(val[i]) != i)
        throw std::runtime_error("movie_slab_axes must be distinct letters of xyz, e.g. xyz or z");
    P.movie_slab_axes = val;
  }
  else if (key == "fields")
    P.fields = val;
  else if (key == "data_files") {
    if (val != "all" && val != "norays" && val != "nodm" && val != "dm" && val != "gas")
      throw std::runtime_error("data_files must be all, norays, nodm, dm or gas");
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
  } else if (key == "c_reduced_schedule") {
    // c0 x1 c1 x2 c2 ... (spaces or commas); "" or "none" turns it off
    std::string t = val;
    for (char& ch : t)
      if (ch == ',') ch = ' ';
    std::istringstream in(t);
    std::vector<double> v;
    std::string w;
    while (in >> w) {
      if (w == "none") continue;
      v.push_back(to_double(w));
    }
    P.c_sched_c.clear();
    P.c_sched_x.clear();
    if (!v.empty()) {
      if (v.size() % 2 != 1) throw std::runtime_error("c_reduced_schedule: give c0 x1 c1 x2 c2 ...");
      for (std::size_t i = 0; i < v.size(); ++i) (i % 2 ? P.c_sched_x : P.c_sched_c).push_back(v[i]);
      for (double c : P.c_sched_c)
        if (!(c > 0)) throw std::runtime_error("c_reduced_schedule: speeds must be > 0");
      for (std::size_t i = 0; i < P.c_sched_x.size(); ++i)
        if (!(P.c_sched_x[i] > 0 && P.c_sched_x[i] <= 1) || (i > 0 && !(P.c_sched_x[i] > P.c_sched_x[i - 1])))
          throw std::runtime_error("c_reduced_schedule: thresholds must increase within (0, 1]");
    }
  } else if (key == "c_reduced_schedule_interp")
    P.c_sched_interp = std::stoi(val);
  else if (key == "ics")
    P.ics = val;
  else if (key == "fields_every")
    P.fields_every = std::stoi(val);
  else if (key == "fields_levels") {
    P.fields_levels = std::stoi(val);
    if (P.fields_levels != 0 && P.fields_levels != 1) throw std::runtime_error("fields_levels must be 0 or 1");
  } else if (key == "fields_levels_axis") {
    if (val != "x" && val != "y" && val != "z") throw std::runtime_error("fields_levels_axis must be x, y or z");
    P.fields_levels_axis = val[0] - 'x';
  } else if (key == "fields_levels_index")
    P.fields_levels_index = std::stoi(val);
  else if (key == "fields_at_data1")
    P.fields_at_data1 = std::stoi(val);
  else if (key == "tensor_dump")
    P.tensor_dump = std::stoi(val);
  else if (key == "max_steps")
    P.max_steps = std::stoi(val);
  else if (key == "run_steps")
    P.run_steps = std::stoi(val);
  else if (key == "step_timers")
    P.step_timers = std::stoi(val);
  else if (set_rt_method_param(P.rtm, key, val)) {
    // moment RT (docs/rt_moments.md)
  } else if (set_plane_param(P.rtm, key, val)) {
    // plane-parallel RT-domain mode and the step caps (docs/plane_parallel.md)
  } else if (set_rt_backend_param(P.rtm, key, val)) {
    // rt_backend = serial | amrex (docs/stage_d.md)
  } else if (set_pm_backend_param(P.rtm, key, val)) {
  } else if (set_amrex_state_param(P.rtm, key, val)) {
    // amrex_state = 0 | 1 (docs/stage_e_state.md)
    // pm_backend = serial | amrex (docs/stage_e_pm.md)
  } else if (set_hydro_scheme_param(P.hydro, key, val)) {
    // hydro scheme (docs/hydro_unsplit.md)
  } else if (set_hydro_backend_param(P.hydro, key, val)) {
    // hydro_backend = serial | amrex (docs/stage_e_hydro.md)
  } else if (key == "ray_dump_dir")
    P.ray_dump.dir = val;
  else if (key == "ray_dump_every")
    P.ray_dump.every = std::stoll(val);
  else if (key == "ray_dump_at")
    P.ray_dump.at = parse_substep_list(val);
  else if (key == "ray_quant_bits")
    P.ray_quant.bits = std::stoi(val);
  else if (key == "ray_quant_dither")
    P.ray_quant.dither = std::stoi(val) != 0;
  else if (key == "ray_quant_struct")
    P.ray_quant.structural = std::stoi(val) != 0;
  else if (key == "ray_quant_offbits")
    P.ray_quant.offbits = std::stoi(val);
  else if (key == "ics_check")
    P.ics_check = std::stoi(val);
  else if (set_cosmology_param(P.rtm, key, val)) {
    // omegam, omegab, hubble0, sigma8, ns, T_min, z_floor (cosmology.hpp)
  } else if (!set_xray_param(P.xray, key, val) && !set_halo_param(P.halo, key, val))
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

// A raw float32 file written in pieces, e.g. from the slabs streamed to
// rank 0 (docs/rank0_streaming.md): the same bytes as write_f32 of the
// whole array. Only the rank that opens it writes; the others pass through.
struct F32Out {
  std::ofstream out;
  std::string path;
  F32Out(const std::string& p, bool open) : path(p) {
    if (!open) return;
    out.open(p, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + p);
  }
  void write(const float* p, std::size_t n) {
    if (out.is_open()) out.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n * sizeof(float)));
  }
  void zeros(std::size_t n) {
    if (!out.is_open()) return;
    const std::vector<float> z(std::min<std::size_t>(n, std::size_t(1) << 24), 0.0f);
    for (std::size_t done = 0; done < n;) {
      const std::size_t m = std::min(z.size(), n - done);
      write(z.data(), m);
      done += m;
    }
  }
  void close() {
    if (!out.is_open()) return;
    out.close();
    if (!out) throw std::runtime_error("error writing " + path);
  }
};

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
// The star-formation sink's particles (docs/stage_f.md section 2): DM/stars.<ztag>
// holds int64 n, then n records {int64 id; float x[3], v[3], mass, a_form}
// sorted by formation time and position (independent of the rank count);
// sink.txt the count and the total mass
void write_stars(const std::string& path, std::vector<pm::Star> st) {
  std::sort(st.begin(), st.end(), [](const pm::Star& a, const pm::Star& b) {
    if (a.a_form != b.a_form) return a.a_form < b.a_form;
    for (int c = 0; c < 3; ++c)
      if (a.x[c] != b.x[c]) return a.x[c] < b.x[c];
    return a.id < b.id;
  });
  std::filesystem::create_directories(std::filesystem::path(path).parent_path());
  std::ofstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write " + path);
  const std::int64_t n = static_cast<std::int64_t>(st.size());
  f.write(reinterpret_cast<const char*>(&n), sizeof n);
  static_assert(sizeof(pm::Star) == 40, "Star must be 40 packed bytes");
  f.write(reinterpret_cast<const char*>(st.data()), static_cast<std::streamsize>(st.size() * sizeof(pm::Star)));
}
void write_sink_summary(const std::string& path, const pm::PMState& pm) {
  double mass = 0;
  for (const pm::Star& st : pm.stars) mass += static_cast<double>(st.mass);
  std::ofstream f(path);
  f << "nstar = " << pm.stars.size() << "\nnstar_ever = " << pm.nstar << "\nmass = " << std::setprecision(17) << mass
    << "\n";
}

// The Ray files of write_data_1 (part of write_data_files; the host rays)
void write_ray_files(const std::string& odir, const std::string& tag, const Simulation& s) {
  const RTGrid& G = s.rt;
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
    if (s.cfg.det_parallel && s.rays.Nray > 0 && !G.htoc.empty()) {
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
    if (s.rays.Nray > 0 && G.htoc.empty()) {
      // amrex_state = 1: no host grid; the host chains are empty (the rays
      // live on the level), so htoc is zero (docs/rank0_streaming.md)
      const std::vector<std::int64_t> z(static_cast<std::size_t>(4 * G.N) * G.N, 0);
      for (int k = 0; k < G.N; ++k)
        oh.write(reinterpret_cast<const char*>(z.data()), static_cast<std::streamsize>(z.size() * sizeof(std::int64_t)));
    } else if (s.rays.Nray > 0) {
      std::vector<std::int64_t> h(G.htoc.size());
      for (std::size_t m = 0; m < h.size(); ++m) h[m] = remap(G.htoc[m]);
      oh.write(reinterpret_cast<const char*>(h.data()),
               static_cast<std::streamsize>(h.size() * sizeof(std::int64_t)));
    }
  }
}

void write_data_files(const std::string& odir, const std::string& tag, const Simulation& s,
                      const std::vector<float>& xv, const std::vector<float>& vgas, bool data1,
                      const std::string& which) {
  if (which == "gas") {  // only Gas/u_gas (the Lya forest velocities of igm_stats.py)
    std::filesystem::create_directories(odir + "/Gas");
    write_f32(odir + "/Gas/u_gas." + tag, s.gas.u.data(), s.gas.u.size());
    return;
  }
  if (which != "nodm") write_f32(odir + "/DM/xv_dm." + tag, xv.data(), xv.size());
  if (which == "dm") return;
  std::filesystem::create_directories(odir + "/Gas");
  write_f32(odir + "/Gas/u_gas." + tag, s.gas.u.data(), s.gas.u.size());
  write_f32(odir + "/Gas/v_gas." + tag, vgas.data(), vgas.size());
  write_f32(odir + "/Gas/w_gas." + tag, s.gas.w.data(), s.gas.w.size());
  if (!s.gas.e.empty())  // hydro_scheme = unsplit: the internal energy density
    write_f32(odir + "/Gas/e_gas." + tag, s.gas.e.data(), s.gas.e.size());
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
  if (which == "norays" || which == "nodm") return;
  write_ray_files(odir, tag, s);

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
  backend_init(argc, argv);  // rt_backend = amrex: amrex::Initialize (MPI)
  struct Fin {
    ~Fin() { backend_finalize(); }
  } fin;
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

    // With MPI (rt_backend = amrex) every rank runs this driver on the full
    // serial state and only rank 0 writes files; an unopened ofstream
    // discards what the other ranks write to it.
    const bool io = io_rank();
    if (io) {
      std::filesystem::create_directories(odir + "/Log");
      std::filesystem::create_directories(odir + "/Pk");
      std::filesystem::create_directories(odir + "/DM");
    }
    backend_barrier();
    // a restart appends to the logs of the run directory
    const auto logmode = restart_from.empty() ? std::ios::out : std::ios::app;
    std::ofstream raylog, hydrolog;
    if (io) {
      raylog.open(odir + "/Log/ray.log", logmode);
      hydrolog.open(odir + "/Log/hydro.log", logmode);
    }
    std::ofstream nullout;
    std::ostream* out = (quiet || !io) ? nullptr : &std::cout;

    // ------------------------------------------------------------------
    // Initialize (radhydro.f90 order)
    // ------------------------------------------------------------------
    auto sim = std::make_unique<Simulation>();
    Simulation& s = *sim;
    s.cfg = Config::make(P.N, P.nray_max_factor, P.Ndm, P.box);
    s.cfg.pm_refine = P.pm_refine;
    s.cfg.legacy_rt = (P.legacy_rt == 1);
    s.cfg.legacy_physics = (P.legacy_physics == 1);
    s.cfg.c_reduced = P.c_sched_c.empty() ? P.c_reduced : P.c_sched_c[0];
    s.cfg.det_parallel = (P.deterministic_parallel == 1);
    s.cfg.ray_rotate = P.ray_rotate;
    s.cfg.ray_rotate_seed = static_cast<std::uint64_t>(P.ray_rotate_seed);
    s.cfg.redshift_per_distance = (P.redshift_per_distance == 1);
    s.cfg.rt_cells_per_substep = P.rt_cells_per_substep;  // docs/rt_multicell.md
    s.cfg.photon_conserve = Config::photon_conserve_effective(P.photon_conserve, s.cfg.legacy_physics);  // B17
    s.cfg.ray_deposit = P.ray_deposit;                     // share | exposure
    s.cfg.cool_suppress_overdensity = P.cool_suppress_overdensity;
    s.cfg.thermal_init = P.thermal_init;  // docs/thermal_history.md
    s.cfg.thermal_floor = P.thermal_floor > 0 ? P.thermal_floor : (P.thermal_init == 1 ? 0.1 : 1.0);
    copy_hydro_scheme(P.hydro, s.cfg);
    // dual_e_self defaults to 1 with thermal_init = recfast (docs/hydro_unsplit.md
    // section 7, U7; user decision 2026-10-27); an explicit key wins
    if (!s.cfg.dual_e_self_set && s.cfg.thermal_init == 1) s.cfg.dual_e_self = true;
    copy_rt_method(P.rtm, s.cfg);
    copy_cosmology_param(P.rtm, s.cfg);  // omegam, ..., T_min, z_floor (cosmology.hpp)
    copy_rt_backend_param(P.rtm, s.cfg);  // rt_backend = serial | amrex
    copy_pm_backend_param(P.rtm, s.cfg);  // pm_backend = serial | amrex
    copy_hydro_backend_param(P.hydro, s.cfg);  // hydro_backend = serial | amrex
    copy_amrex_state_param(P.rtm, s.cfg);  // amrex_state = 0 | 1
    copy_plane_param(P.rtm, s.cfg);        // plane_parallel, dt_max_myr, rt_substeps_max
    check_recomb_param(s.cfg);  // recomb_rad (recomb_rad.hpp)
    check_rt_backend_param(s.cfg);
    check_pm_backend_param(s.cfg);
    check_hydro_backend_param(s.cfg);
    check_amrex_state_param(s.cfg);
    check_rt_cells_param(s.cfg);
    check_plane_param(s.cfg);
    if (s.cfg.plane_parallel && s.cfg.plane_z_re < 0)
      throw std::runtime_error("plane_parallel = 1 needs plane_z_re (the turn-on redshift) in radhydro");
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
    s.raylog = io ? &raylog : nullptr;
    s.out = out;
    const Config& cfg = s.cfg;
    // gammaHI_avg in fields or movie_fields: the substep-averaged Gamma_HI
    // (Simulation::gamma_avg, rt.hpp; docs/igm_diagnostics.md). The serial
    // RT, or the moment RT on rt_backend = amrex, whose level holds the sums
    // (amrex_gamma_avg_field gathers them; not the ray levels).
    const bool gavg_fields = ("," + P.fields + ",").find(",gammaHI_avg,") != std::string::npos;
    const bool gavg_movie = ("," + P.movie_fields + ",").find(",gammaHI_avg,") != std::string::npos;
    if (gavg_fields || gavg_movie) {
      if (use_amrex(cfg) && !use_moments(cfg))
        throw std::runtime_error(
            "gammaHI_avg (the substep-averaged Gamma_HI) is implemented for rt_backend = serial and for "
            "rt_method = m1 | otvet on rt_backend = amrex");
      s.gamma_avg.on = true;
    }
    // xray_heating = 1 (docs/xray_heating.md): a uniform X-ray background
    // from the SFR of the abundance-matched source model. Its rates are
    // computed once per hydro step with RT, at the middle of the step, with
    // the box's mass-weighted absorber fractions; Log/xray.txt lists them.
    // the residual-electron Compton step acts on the host gas grid only: not
    // on the levels of amrex_state = 1, nor on refined hydro levels
    if (cfg.thermal_init == 1 && (cfg.amrex_state || (use_amrex_hydro(cfg) && cfg.amrex_max_level > 0)))
      throw std::runtime_error(
          "thermal_init = recfast is implemented for the host gas grid only (amrex_state = 0, no refined "
          "hydro levels)");
    // The RECFAST table (recfast_table.inc) was computed for the original
    // cosmology; another cosmology needs its own table
    if (cfg.thermal_init == 1 && !cfg.cosmo.is_original())
      throw std::runtime_error(
          "thermal_init = recfast: the RECFAST table is for the original cosmology (omegam 0.27, omegab "
          "0.045, hubble0 0.70), not " + cfg.cosmo.summary());
    std::unique_ptr<XrayBackground> xbg;
    std::ofstream xraylog;
    if (P.xray.on) {
      if (use_amrex(cfg) || cfg.amrex_state)
        throw std::runtime_error("xray_heating is implemented for rt_backend = serial only");
      xbg = std::make_unique<XrayBackground>(P.xray, P.halo.model, cfg.cosmo);
      if (io) {
        xraylog.open(odir + "/Log/xray.txt", logmode);
        if (restart_from.empty())
          xraylog << "# step z sfrd[Msun/yr/Mpc^3] eps_X[erg/s/Mpc^3] xHI xHeI xHeII"
                     " gamma_HI gamma_HeI gamma_HeII[1/s] phi_HI phi_HeI phi_HeII[erg/s]\n";
      }
      s.xray.on = true;
      s.xray_hook = [&xbg, &xraylog](Simulation& sim) {
        // radiative_transfer has reset the clock to the start of the step
        const double a = sim.clk.scalefactor + sim.clk.da_sim / 2;
        double xf[3];
        box_ion_fractions(sim, xf);
        const XrayRates r = xbg->rates(1 / a - 1, xf);
        for (int j = 0; j < 3; ++j) {
          sim.xray.gamma[j] = r.gamma[j];
          sim.xray.phi[j] = r.phi[j];
        }
        if (xraylog.is_open()) {
          xraylog << sim.clk.istep << std::setprecision(8) << " " << r.z << " " << r.sfrd << " "
                  << r.eps;
          for (int j = 0; j < 3; ++j) xraylog << " " << r.xfrac[j];
          for (int j = 0; j < 3; ++j) xraylog << " " << r.gamma[j];
          for (int j = 0; j < 3; ++j) xraylog << " " << r.phi[j];
          xraylog << "\n" << std::flush;
        }
      };
    }

    Domains gdom;
    HydroState hs;
    pm::PMState pm;

    init_variables(cfg.cosmo, P.z_init, s.clk);       // init_variables
    init_domain_gas(cfg, gdom);                        // init_domain_gas
    init_domain_rt(cfg, s.dom);                        // init_domain_rt
    // amrex_state = 1 (amrex_state.hpp): the gas state on the levels, built
    // before any level; no rank keeps the host grids (the outputs stream
    // the levels to rank 0, docs/rank0_streaming.md), and the other ranks
    // never allocate the host PM arrays (at pm_refine = 2
    // with 1024^3 particles they are 45 GB per rank)
    const bool on_level = use_amrex_state(cfg);
    init_hydro(cfg, s.gas, hs);                        // init_hydro
    pm::init_pm(pm, pm::PMConfig::make(cfg), gdom, !(on_level && !io));  // init_pm
    init_rt(s);                                        // init_rt (rates)
    init_ray_tracing(s);                               // init_ray_tracing
    plane_init(s);                                     // plane_parallel = 1: the face sources
    if (cfg.sf_overdensity > 0 && !on_level)
      throw std::runtime_error("sf_overdensity > 0 (the star-formation sink) needs amrex_state = 1 in radhydro");
    if (on_level) {
      amrex_state_init(s);
      // the gas grid lives on the level; no rank keeps a host copy (rank 0
      // streams it for the outputs, docs/rank0_streaming.md)
      for (auto* v : {&s.gas.u, &s.gas.w, &hs.v}) {
        v->clear();
        v->shrink_to_fit();
      }
    }
    // the particles' mirror on rank 0 (pm.xv): current after the load and
    // after a gather, stale after a PM step
    bool xv_fresh = true;
    auto sync_xv = [&]() {
      if (!on_level || xv_fresh) return;
      amrex_pm_gather_particles(pm);
      xv_fresh = true;
    };
    // the particles into rank 0's xv and the sink's stars into rank 0's
    // list (collective), before an output or a checkpoint. The grids are
    // not copied: every output streams them from the levels
    // (docs/rank0_streaming.md).
    auto sync_mirror = [&]() {
      if (!on_level) return;
      sync_xv();
      if (cfg.sf_overdensity > 0) amrex_pm_gather_stars(pm);
    };
    // the star-formation sink: the hydro's events of the last call become
    // star particles on the PM level (collective)
    std::int64_t nstar_step = 0;
    auto sink_to_stars = [&]() {
      if (!(cfg.sf_overdensity > 0)) return;
      std::vector<SinkEvent> ev;
      amrex_hydro_take_sink_events(ev);
      nstar_step += amrex_pm_add_stars(pm, ev, s.clk.scalefactor);
    };
    ray_dump_enable(s, P.ray_dump);                    // (study only; off by default)
    s.ray_quant = P.ray_quant;                         // (research only; off by default)
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
    if (!(on_level && !io)) fft::fft_3d('f', pm.pc.Nm, pm.pc.Nm, pm.pc.Nm, pm.density.data());  // init_fft
    if (on_level && !io) {
      // amrex_state = 1: the particles live on the level; rank 0 keeps the
      // host PM arrays as the mirror (the halo finder, P(k), the outputs)
      for (auto* v : {&pm.ll, &pm.hoc, &pm.hoc_cpu}) {
        v->clear();
        v->shrink_to_fit();
      }
      for (auto* v : {&pm.xv, &pm.density, &pm.fmesh}) {
        v->clear();
        v->shrink_to_fit();
      }
    }
    calc_sim_var(cfg, s.clk, s.units);
    if (P.tensor_dump && use_amrex(cfg))
      throw std::runtime_error("tensor_dump is implemented for the serial RT backend only");

    // thermal_init = recfast: the ICs at the RECFAST temperature of z_init
    // amrex_state = 1: the thermal energy of the initial conditions, for
    // the streamed copy of u_gas.init onto the level (ics_to_level)
    double ics_Tinit = 0;
    std::vector<float> empty_grid;  // read_init_conditions: no host gas grid
    if (restart_from.empty()) {
      if (P.ics_check) check_ics_info(ics_dir, P.z_init, cfg.cosmo);  // z_init, cosmology
      const double T_init_K = cfg.thermal_init == 1 ? recfast_T(s.clk.redshift) : -1.0;
      ics_Tinit = pm::ics_tinit(s.units, T_init_K);
      pm::read_init_conditions(ics_dir, pm, on_level ? empty_grid : s.gas.u, s.gas.w, s.units, T_init_K);
    }
    // amrex_state = 1 (docs/rank0_streaming.md): u_gas.init onto the gas
    // state, rank 0 reading it slab by slab; per cell the thermal energy, w
    // (read_init_conditions) and e (init_unsplit_state). Collective.
    auto ics_to_level = [&]() {
      const std::string path = ics_dir + "/u_gas.init";
      std::ifstream f;
      if (io) {
        f.open(path, std::ios::binary | std::ios::ate);
        if (!f) throw std::runtime_error("cannot open " + path);
        const std::int64_t want = 5 * cfg.Ncell * static_cast<std::int64_t>(sizeof(float));
        if (static_cast<std::int64_t>(f.tellg()) < want)
          throw std::runtime_error(path + " is shorter than 5 N^3 floats");
      }
      std::vector<float> ug;
      amrex_state_fill_gas(s, 0, 10, [&](const StreamSlab& sl, float* d) {
        const std::int64_t n = sl.ncell();
        ug.resize(static_cast<std::size_t>(5 * n));
        f.seekg(static_cast<std::streamoff>(5 * sl.c0() * static_cast<std::int64_t>(sizeof(float))));
        f.read(reinterpret_cast<char*>(ug.data()), static_cast<std::streamsize>(5 * n * sizeof(float)));
        if (!f) throw std::runtime_error("cannot read " + path);
#pragma omp parallel for schedule(static)
        for (std::int64_t c = 0; c < n; ++c) {
          float* u = &ug[static_cast<std::size_t>(5 * c)];
          float w[4];
          pm::ics_gas_cell(u, w, ics_Tinit);
          float* o = d + 10 * c;
          for (int m = 0; m < 5; ++m) o[kStreamGasU + m] = u[m];
          o[kStreamGasE] = unsplit_init_e(u, s.units);
          for (int m = 0; m < 4; ++m) o[kStreamGasW + m] = w[m];
        }
      });
    };

    if (out)
      *out << " radhydro: N = " << cfg.Ngas << ", Ndm = " << cfg.Ndm
           << (cfg.pm_refine > 1 ? ", pm_refine = " + std::to_string(cfg.pm_refine) +
                                       " (PM mesh " + std::to_string(pm.pc.Nm) + "^3" +
                                       (fft::fft_threaded() ? ", threaded FFTW)" : ", serial FFTW)")
                                 : std::string())
           << ", box = " << cfg.box
           << " Mpc/h, z_init = " << P.z_init << ", Nray_max = " << cfg.Nray_max
           << "\n legacy_rt = " << P.legacy_rt << ", legacy_physics = " << P.legacy_physics << ", c_reduced = " << P.c_reduced
           << (P.c_sched_c.empty() ? std::string() : " (overridden by c_reduced_schedule)")
           << (P.redshift_per_distance ? ", redshift_per_distance = 1" : "")
           << (P.redshift_per_distance && !(P.c_reduced < 1) ? " (no effect: c_reduced >= 1)" : "")
           << (P.rt_cells_per_substep != 1
                   ? ", rt_cells_per_substep = " + std::to_string(P.rt_cells_per_substep)
                   : std::string())
           << (P.photon_conserve >= 0 ? ", photon_conserve = " + std::to_string(P.photon_conserve) +
                                            (P.photon_conserve == int(!cfg.legacy_physics) ? "" : " (set explicitly)")
                                      : std::string())
           << (P.ray_deposit ? ", ray_deposit = exposure" : "")
           << (P.cool_suppress_overdensity != 200.0
                   ? ", cool_suppress_overdensity = " + std::to_string(P.cool_suppress_overdensity)
                   : std::string())
           << (plane_summary(cfg).empty() ? std::string() : "\n " + plane_summary(cfg))
           << (cfg.thermal_init == 1 || cfg.thermal_floor != 1.0
                   ? std::string(", thermal_init = ") + (cfg.thermal_init == 1 ? "recfast" : "floor") +
                         ", thermal_floor = " + std::to_string(cfg.thermal_floor) + " K"
                   : std::string())
           << "\n " << hydro_scheme_summary(cfg) << ", " << rt_method_summary(cfg)
           << ", " << rt_backend_summary(cfg) << ", " << hydro_backend_summary(cfg)
           << ", " << pm_backend_summary(cfg) << (on_level ? ", " + amrex_state_summary(cfg) : std::string())
           << ", ranks = " << backend_nranks()
           << ", deterministic_parallel = " << P.deterministic_parallel
           << ", threads = " << par::max_threads() << ", sources = "
           << (P.halo.mode == HaloSourceParams::Mode::Legacy ? "legacy" : "onthefly")
           << (P.halo.fof_dt_myr > 0 ? ", fof_dt_myr = " + std::to_string(P.halo.fof_dt_myr)
                                     : std::string())
           << ", zdata1: " << Ndata1 << " redshifts down to " << zdata1.back() << "\n"
           << (cosmology_summary(cfg).empty() ? std::string() : " " + cosmology_summary(cfg) + "\n");

    if (io) plane_write_spectrum(odir, s);  // plane_parallel = 1: the bins and fluxes
    Timers T;
    std::int64_t& Nhalo = dst.Nhalo;
    // amrex_balance = 1: the balance efficiency of every RT step (mean over
    // max rank cost) and whether the mapping was changed
    std::int64_t balance_moves = 0;
    auto balance_log = [&](bool moved, double e0, double e1) {
      if (moved) ++balance_moves;
      if (out)
        *out << " balance: efficiency " << static_cast<float>(e0) << (moved ? " -> " : " (proposed ")
             << static_cast<float>(e1) << (moved ? ", remapped" : ", kept)") << "\n";
    };
    auto summaries = [&]() {
      const auto ts = Clock::now();
      if (on_level) {
        sync_xv();
        if (io) pm::summarize_dm(pm, s.units, out);
      } else {
        pm::summarize_dm(pm, s.units, out);
      }
      if (on_level) {
        // the clamps on the level, the sums on rank 0 over the gas streamed
        // slab by slab in the serial order (docs/rank0_streaming.md; the
        // grid velocity is zero with the unsplit hydro); the maxima of the
        // timestep onto every rank
        amrex_state_clamp_gas(s);
        GasSummer gsr(cfg, s.units);
        std::vector<float> su, sw, se;
        amrex_state_stream_gas(s, 0, 10, [&](const StreamSlab& sl) {
          const std::int64_t n = sl.ncell();
          su.resize(static_cast<std::size_t>(5 * n));
          sw.resize(static_cast<std::size_t>(4 * n));
          se.resize(static_cast<std::size_t>(n));
          const float* a = sl.data;
#pragma omp parallel for schedule(static)
          for (std::int64_t c = 0; c < n; ++c) {
            for (int m = 0; m < 5; ++m) su[5 * c + m] = a[10 * c + kStreamGasU + m];
            se[c] = a[10 * c + kStreamGasE];
            for (int m = 0; m < 4; ++m) sw[4 * c + m] = a[10 * c + kStreamGasW + m];
          }
          gsr.add_planes(sl.nk, su.data(), sw.data(), se.data(), nullptr);
        });
        double v[3] = {0, 0, 0};
        if (io) {
          const GasSummary gsum = gsr.finish(s.clk, out);
          hydrolog << format_gas_record(gsum.record) << "\n";
          v[0] = s.clk.cmax_gas;
          v[1] = s.clk.vmax_gas;
          v[2] = s.clk.vrms_gas;
        }
        backend_bcast(v, 3);
        s.clk.cmax_gas = v[0];
        s.clk.vmax_gas = v[1];
        s.clk.vrms_gas = v[2];
      } else {
        const GasSummary gsum = summarize_gas(cfg, s.units, s.clk, s.gas, hs, out);
        hydrolog << format_gas_record(gsum.record) << "\n";
      }
      summarize_rt(s);
      T.summaries += std::chrono::duration<double>(Clock::now() - ts).count();
    };

    // The sums of write_stats over the RT grid (k,j,i order, real(8)), also
    // those of update_sched_x: from rank 0's host grid, or (amrex_state = 1)
    // streamed from the level slab by slab in the same order
    // (docs/rank0_streaming.md). Collective with amrex_state = 1; the sums
    // are rank 0's.
    struct GridSums {
      double xv = 0, xm = 0, mass = 0, tv = 0, tm = 0, xhe3 = 0;
      bool filled = false;  // the RT has written the grid (allocate zero-fills it)
      void add(double fHI, double rho, double T_, float fHeI, float fHeII) {
        xv = xv + (1 - fHI);
        xm = xm + rho * (1 - fHI);
        mass = mass + rho;
        tv = tv + T_;
        tm = tm + rho * T_;
        xhe3 = xhe3 + (1 - static_cast<double>(fHeI) - static_cast<double>(fHeII));
        filled = filled || fHI != 0;
      }
    };
    auto grid_sums = [&]() {
      GridSums g;
      if (on_level) {
        // rhog, fHI, fHeI, fHeII, temp (the first 5 components of a cell)
        amrex_state_stream_rt(s, 0, 5, [&](const StreamSlab& sl) {
          const std::int64_t n = sl.ncell();
          for (std::int64_t c = 0; c < n; ++c) {
            const float* a = sl.data + 5 * c;
            g.add(a[1], a[0], a[4], a[2], a[3]);
          }
        });
      } else if (io) {
        const RTGrid& G = s.rt;
        for (std::int64_t c = 0; c < G.ncell; ++c) g.add(G.fHI[c], G.rhog[c], G.temp[c], G.fHeI[c], G.fHeII[c]);
      }
      return g;
    };
    // write_stats' sums, kept for update_sched_x right after it
    GridSums stats_sums;
    bool stats_fresh = false;

    auto write_stats = [&](const std::string& d) {
      stats_sums = grid_sums();
      stats_fresh = true;
      if (!io) return;
      // global diagnostics over the RT grid, k,j,i order, real(8) sums
      const double xv = stats_sums.xv, xm = stats_sums.xm, mass = stats_sums.mass, tv = stats_sums.tv,
                   tm = stats_sums.tm, xhe3 = stats_sums.xhe3;
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
      if (!P.c_sched_c.empty()) d1("c_reduced", s.cfg.c_reduced);  // this step's (schedule only)
      if (s.plane.on) {  // plane_parallel = 1: the step's photon bookkeeping (plane_parallel.hpp)
        d1("plane_emitted", s.plane.emitted);
        d1("plane_absorbed", s.plane.absorbed());
        d1("plane_escaped", s.plane.escaped());
      }
      std::fclose(fp);
    };

    // c_reduced_schedule: x_HII,V of the box at the end of the last step
    // (the RT grid, summed as in write_stats: write_stats' sums when it has
    // just run, else a pass of its own, e.g. after reading a checkpoint) and
    // the speed it selects for the next step
    double sched_x = 0;
    auto update_sched_x = [&]() {
      const bool fresh = stats_fresh;
      stats_fresh = false;
      if (P.c_sched_c.empty()) return;
      const GridSums g = fresh ? stats_sums : grid_sums();
      double v = 0;
      if (io) v = g.filled ? g.xv / cfg.nrtcr : 0;  // before the first RT step: c0
      if (on_level) backend_bcast(&v, 1);
      sched_x = v;
    };
    bool sched_announced = false;
    auto apply_sched = [&]() {
      if (P.c_sched_c.empty()) return;
      double c = P.c_sched_c[0];
      if (P.c_sched_interp && !P.c_sched_x.empty()) {
        // linear between the plateau midpoints, constant outside them
        const std::size_t n = P.c_sched_c.size();
        std::vector<double> xm(n);
        for (std::size_t k = 0; k < n; ++k) {
          const double lo = k == 0 ? 0.0 : P.c_sched_x[k - 1];
          const double hi = k + 1 == n ? 1.0 : P.c_sched_x[k];
          xm[k] = 0.5 * (lo + hi);
        }
        if (sched_x >= xm[n - 1]) {
          c = P.c_sched_c[n - 1];
        } else {
          for (std::size_t k = 0; k + 1 < n; ++k)
            if (sched_x >= xm[k] && sched_x < xm[k + 1]) {
              c = P.c_sched_c[k] +
                  (P.c_sched_c[k + 1] - P.c_sched_c[k]) * (sched_x - xm[k]) / (xm[k + 1] - xm[k]);
              break;
            }
        }
      } else {
        for (std::size_t i = 0; i < P.c_sched_x.size(); ++i)
          if (sched_x >= P.c_sched_x[i]) c = P.c_sched_c[i + 1];
      }
      if (c != s.cfg.c_reduced || !sched_announced) {
        if (io) {
          std::printf(" c_reduced_schedule: c_reduced = %g from step %d (z = %.5f, x_HII,V = %.6f)\n", c,
                      s.clk.istep + 1, s.clk.redshift, sched_x);
          std::fflush(stdout);
        }
        s.cfg.c_reduced = c;
        sched_announced = true;
      }
    };

    // The fields of fields_every / fields_at_data1 (fields: the list; empty:
    // all but the optional ones)
    auto field_listed = [&](const char* name) {
      return P.fields.empty() || ("," + P.fields + ",").find(std::string(",") + name + ",") != std::string::npos;
    };
    auto field_opted = [&](const char* name) {  // only when listed explicitly
      return ("," + P.fields + ",").find(std::string(",") + name + ",") != std::string::npos;
    };
    // amrex_state = 1: every file streamed from the levels slab by slab
    // (docs/rank0_streaming.md), the same files in the same order as the
    // host path below. Collective.
    auto write_fields_streamed = [&](const std::string& d) {
      auto SR = [&](const char* name, int c0, int nc) {  // RT cell components
        if (!field_listed(name)) return;
        F32Out f(d + "/" + name + ".bin", io);
        amrex_state_stream_rt(s, c0, nc, [&](const StreamSlab& sl) { f.write(sl.data, sl.ncell() * nc); });
        f.close();
      };
      auto SG = [&](const char* name, int c0, int nc) {  // gas state components
        if (!field_listed(name)) return;
        F32Out f(d + "/" + name + ".bin", io);
        amrex_state_stream_gas(s, c0, nc, [&](const StreamSlab& sl) { f.write(sl.data, sl.ncell() * nc); });
        f.close();
      };
      using CD = rrk::CellData;
      constexpr int fl = sizeof(float);
      SR("fHI", offsetof(CD, fHI) / fl, 1);
      SR("fHeI", offsetof(CD, fHeI) / fl, 1);
      SR("fHeII", offsetof(CD, fHeII) / fl, 1);
      SR("temp", offsetof(CD, temp) / fl, 1);
      SR("rhog", offsetof(CD, rhog) / fl, 1);
      SR("ngam", offsetof(CD, ngam) / fl, kNfreq);
      SR("nu", offsetof(CD, nu) / fl, kNfreq);
      SR("fgam", offsetof(CD, fgam) / fl, kNfreq);
      SR("tau", offsetof(CD, tau) / fl, kNfreq);
      SR("gamma", offsetof(CD, gamma) / fl, 3);
      if (field_opted("gammaHI")) SR("gammaHI", offsetof(CD, gamma) / fl, 1);
      if (field_opted("gammaHI_avg")) {
        F32Out f(d + "/gammaHI_avg.bin", io);
        amrex_gamma_avg_stream(s, [&](const StreamSlab& sl) { f.write(sl.data, sl.ncell()); });
        f.close();
      }
      if (field_opted("Luv")) SR("Luv", offsetof(CD, Luv) / fl, 1);
      SR("lambda", offsetof(CD, lambda) / fl, 3);
      SR("dx", offsetof(CD, dx) / fl, 3);
      SG("u_gas", kStreamGasU, 5);
      SG("w_gas", kStreamGasW, 4);
      if (field_listed("v_gas")) {  // the grid velocity: zero with the unsplit hydro
        F32Out f(d + "/v_gas.bin", io);
        f.zeros(static_cast<std::size_t>(3 * cfg.Ncell));
        f.close();
      }
      SG("e_gas", kStreamGasE, 1);
      if (field_listed("xv_dm")) {
        sync_xv();
        if (io) write_f32(d + "/xv_dm.bin", pm.xv.data(), pm.xv.size());
      }
    };

    // write_data_files for amrex_state = 1: the grids streamed from the
    // levels slab by slab (docs/rank0_streaming.md), the particles from
    // rank 0's copy, the Ray files from the host rays; the same files in
    // the same order. Collective.
    auto write_data_files_streamed = [&](const std::string& tag, bool data1, const std::string& which) {
      auto SG = [&](const std::string& path, int c0, int nc) {
        F32Out f(path, io);
        amrex_state_stream_gas(s, c0, nc, [&](const StreamSlab& sl) { f.write(sl.data, sl.ncell() * nc); });
        f.close();
      };
      if (which == "gas") {
        if (io) std::filesystem::create_directories(odir + "/Gas");
        SG(odir + "/Gas/u_gas." + tag, kStreamGasU, 5);
        return;
      }
      if (which != "nodm") {
        sync_xv();
        if (io) write_f32(odir + "/DM/xv_dm." + tag, pm.xv.data(), pm.xv.size());
      }
      if (which == "dm") return;
      if (io) std::filesystem::create_directories(odir + "/Gas");
      SG(odir + "/Gas/u_gas." + tag, kStreamGasU, 5);
      {  // the grid velocity: zero with the unsplit hydro
        F32Out f(odir + "/Gas/v_gas." + tag, io);
        f.zeros(static_cast<std::size_t>(3 * cfg.Ncell));
        f.close();
      }
      SG(odir + "/Gas/w_gas." + tag, kStreamGasW, 4);
      SG(odir + "/Gas/e_gas." + tag, kStreamGasE, 1);  // hydro_scheme = unsplit
      if (!data1) return;
      if (io) std::filesystem::create_directories(odir + "/RT");
      {  // the rt_grid record: the first 36 floats of a cell (rrk::CellData)
        F32Out f(odir + "/RT/rt_grid." + tag, io);
        amrex_state_stream_rt(s, 0, 36, [&](const StreamSlab& sl) { f.write(sl.data, sl.ncell() * 36); });
        f.close();
      }
      if (which == "norays" || which == "nodm") return;
      if (io) write_ray_files(odir, tag, s);
    };

    auto write_fields = [&](const std::string& d) {
      if (on_level) {
        write_fields_streamed(d);
        return;
      }
      // gammaHI_avg on rt_backend = amrex: the level's sums onto rank 0 (collective)
      std::vector<float> ga;
      if (gavg_fields && use_amrex(cfg)) amrex_gamma_avg_field(s, ga);
      if (!io) return;
      const RTGrid& G = s.rt;
      auto W = [&](const char* name, const std::vector<float>& v) {
        if (!field_listed(name)) return;
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
      // gammaHI: only the HI rate, gamma(1) (written only when listed in fields)
      if (("," + P.fields + ",").find(",gammaHI,") != std::string::npos) {
        std::vector<float> g1(static_cast<std::size_t>(G.ncell));
        for (std::int64_t c = 0; c < G.ncell; ++c) g1[c] = G.gamma[c * 3];
        W("gammaHI", g1);
      }
      // gammaHI_avg: the substep-averaged HI rate of the last hydro step with
      // RT (written only when listed in fields)
      if (("," + P.fields + ",").find(",gammaHI_avg,") != std::string::npos) {
        if (!use_amrex(cfg)) gamma_avg_field(s, ga);
        W("gammaHI_avg", ga);
      }
      // Luv: the source luminosity per RT cell (written only when listed)
      if (("," + P.fields + ",").find(",Luv,") != std::string::npos) W("Luv", G.Luv);
      W("lambda", G.lambda);
      W("dx", G.dx);
      W("u_gas", s.gas.u);
      W("w_gas", s.gas.w);
      W("v_gas", hs.v);
      if (!s.gas.e.empty()) W("e_gas", s.gas.e);
      W("xv_dm", pm.xv);
    };

    // fields_levels = 1: the plane at the finest resolution (collective)
    auto write_levels_slice = [&](const std::string& d) {
      if (!P.fields_levels) return;
      const int axis = P.fields_levels_axis;
      const int idx = P.fields_levels_index >= 0 ? P.fields_levels_index : cfg.Nrt / 2;
      std::vector<float> out;
      int nf = 0, finest = 0, plane = 0;
      if (use_amrex(cfg)) {
        amrex_rt_level_slice(s, axis, idx, out, nf, finest, plane);
      } else if (io) {  // the serial RT: the base plane of rank 0's grid
        const RTGrid& G = s.rt;
        const int N = cfg.Nrt;
        if (idx < 0 || idx >= N) throw std::runtime_error("fields_levels_index must be within 0..N-1");
        nf = N;
        plane = idx;
        out.assign(static_cast<std::size_t>(kLevelSliceComps) * N * N, 0.0f);
        const int ia = axis == 0 ? 1 : 0, ib = axis == 2 ? 1 : 2;
        for (int b = 0; b < N; ++b)
          for (int a = 0; a < N; ++a) {
            int ijk[3];
            ijk[axis] = idx;
            ijk[ia] = a;
            ijk[ib] = b;
            const std::int64_t c = ijk[0] + static_cast<std::int64_t>(N) * (ijk[1] + static_cast<std::int64_t>(N) * ijk[2]);
            float* o = out.data() + (static_cast<std::size_t>(a) + static_cast<std::size_t>(N) * b) * kLevelSliceComps;
            o[0] = G.rhog[c];
            o[1] = G.fHI[c];
            o[2] = G.fHeI[c];
            o[3] = G.fHeII[c];
            o[4] = G.temp[c];
          }
      }
      if (!io) return;
      write_f32(d + "/levels_slice.bin", out.data(), out.size());
      std::ofstream t(d + "/levels_slice.txt");
      const char* names = "xyz";
      const int ia = axis == 0 ? 1 : 0, ib = axis == 2 ? 1 : 2;
      t << "# fields_levels: one plane of the RT cells at the finest resolution\n"
        << "# levels_slice.bin: float32 (ncomp, n, n) in Fortran order, the first in-plane axis fastest\n"
        << "ncomp = " << kLevelSliceComps << "\ncomponents = rhog fHI fHeI fHeII temp level\n"
        << "n = " << nf << "\nfinest = " << finest << "\ncells_per_base_cell = " << (1 << finest)
        << "\naxis = " << names[axis] << "\nin_plane = " << names[ia] << " " << names[ib]
        << "\nbase_index = " << idx << "\nplane = " << plane
        << "\n# plane: the fine-cell index along the axis; its centre is at (plane + 0.5) / cells_per_base_cell base cells\n";
    };

    auto step_dir = [&]() {
      char name[32];
      std::snprintf(name, sizeof name, "step_%04d", s.clk.istep);
      const std::string d = odir + "/" + name;
      if (io) std::filesystem::create_directories(d);
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
    // movie_slab > 0: the slabs of movie_slab cells normal to each axis in
    // movie_slab_axes, one after the other in one file. Each block is the
    // sub-array with that axis's index restricted to s0 .. s0+t-1 (0-based,
    // s0 = Nrt/2 - t/2), in Fortran order (ncomp, nx, ny, nz) with
    // n_axis = t and the other two = Nrt; ncomp = 1 except gamma (3) and
    // ngam (Nfreq).
    const int mN = cfg.Nrt;
    const int mt = std::min(P.movie_slab, mN);
    const int ms0 = mN / 2 - mt / 2;
    std::vector<float> slabbuf;
    auto write_slabs = [&](const std::string& path, const std::vector<float>& v) {
      const std::int64_t nc = static_cast<std::int64_t>(mN) * mN * mN;
      const std::int64_t ncomp = nc > 0 ? static_cast<std::int64_t>(v.size()) / nc : 0;
      slabbuf.clear();
      for (char a : P.movie_slab_axes) {
        int lo[3] = {0, 0, 0}, hi[3] = {mN, mN, mN};
        const int d = a - 'x';
        lo[d] = ms0;
        hi[d] = ms0 + mt;
        for (int k = lo[2]; k < hi[2]; ++k)
          for (int j = lo[1]; j < hi[1]; ++j) {
            const std::int64_t c0 =
                lo[0] + static_cast<std::int64_t>(mN) * (j + static_cast<std::int64_t>(mN) * k);
            const float* p = v.data() + c0 * ncomp;
            slabbuf.insert(slabbuf.end(), p, p + (hi[0] - lo[0]) * ncomp);
          }
      }
      write_f32(path, slabbuf.data(), slabbuf.size());
    };
    // gammaHI: the HI component of rt_grid%gamma alone (one float per cell)
    std::vector<float> gammahi_buf, gion, gavg_buf;
    auto gamma_hi = [&](const RTGrid& G) -> const std::vector<float>& {
      gammahi_buf.resize(static_cast<std::size_t>(G.ncell));
      for (std::int64_t c = 0; c < G.ncell; ++c) gammahi_buf[c] = G.gamma[c * 3];
      return gammahi_buf;
    };
    // A frame's line of movie/frames.txt (rank 0): the whole-box x_HII sums
    // and the Gamma_HI statistics (gamma(1) of the last RT substep, s^-1):
    // volume mean, median and mean over the ionised cells (x_HII > 0.5),
    // volume fraction with Gamma_HI > 1e-14 (serial sums: thread-independent)
    struct FrameSums {
      double xv = 0, xm = 0, mass = 0, gsum = 0, gisum = 0;
      std::int64_t nhi = 0;
      void add(float fHI_, float rho_) {
        xv = xv + (1 - static_cast<double>(fHI_));
        xm = xm + static_cast<double>(rho_) * (1 - static_cast<double>(fHI_));
        mass = mass + static_cast<double>(rho_);
      }
    };
    auto frame_line = [&](const Simulation& ss, const FrameSums& fs) {
      const std::int64_t nion = static_cast<std::int64_t>(gion.size());
      double gmed = 0;
      if (nion > 0) {
        auto mid = gion.begin() + nion / 2;
        std::nth_element(gion.begin(), mid, gion.end());
        gmed = *mid;
        if (nion % 2 == 0) gmed = 0.5 * (gmed + *std::max_element(gion.begin(), mid));
      }
      char line[384];
      std::snprintf(line, sizeof line, "%6d %6d %6d %.8f %.8e %.8e %.6e %.6e %.6e %.8e\n",
                    dst.movie_frame, ss.clk.istep, ss.clk.istep_rt, ss.clk.redshift,
                    fs.xv / cfg.nrtcr, fs.xm / fs.mass, fs.gsum / cfg.nrtcr, gmed,
                    nion > 0 ? fs.gisum / static_cast<double>(nion) : 0.0,
                    static_cast<double>(fs.nhi) / cfg.nrtcr);
      movielog << line << std::flush;
    };
    // the Gamma_HI part of a cell's statistics (gion: the ionised cells' values)
    auto frame_gamma = [&](FrameSums& fs, float fHI_, float g) {
      fs.gsum = fs.gsum + static_cast<double>(g);
      if (g > 1e-14) fs.nhi += 1;
      if (1 - static_cast<double>(fHI_) > 0.5) {
        fs.gisum = fs.gisum + static_cast<double>(g);
        gion.push_back(g);
      }
    };
    // movie_slab > 0: the cells of a streamed slab (ncomp floats per cell)
    // that fall in each axis block, appended in the order of write_slabs
    auto slab_add = [&](std::vector<std::vector<float>>& blocks, const float* v, int k0, int nk, int ncomp) {
      blocks.resize(P.movie_slab_axes.size());
      for (std::size_t ai = 0; ai < P.movie_slab_axes.size(); ++ai) {
        int lo[3] = {0, 0, 0}, hi[3] = {mN, mN, mN};
        const int d = P.movie_slab_axes[ai] - 'x';
        lo[d] = ms0;
        hi[d] = ms0 + mt;
        std::vector<float>& out = blocks[ai];
        for (int k = std::max(lo[2], k0); k < std::min(hi[2], k0 + nk); ++k)
          for (int j = lo[1]; j < hi[1]; ++j) {
            const std::int64_t c0 =
                lo[0] + static_cast<std::int64_t>(mN) * (j + static_cast<std::int64_t>(mN) * (k - k0));
            const float* p = v + c0 * ncomp;
            out.insert(out.end(), p, p + static_cast<std::int64_t>(hi[0] - lo[0]) * ncomp);
          }
      }
    };

    // amrex_state = 1: the frame streamed from the level slab by slab
    // (docs/rank0_streaming.md): one pass over the components the frame
    // needs (rhog, fHI and gamma(1) for the statistics, and the movie
    // fields'), and one for gammaHI_avg. The same files and line as the
    // host path. Collective.
    auto write_frame_streamed = [&](const Simulation& ss) {
      using CD = rrk::CellData;
      constexpr int fl = sizeof(float);
      const int cRho = offsetof(CD, rhog) / fl, cHI = offsetof(CD, fHI) / fl, cG = offsetof(CD, gamma) / fl;
      std::vector<int> comps = {cRho, cHI, cG};
      auto pos_of = [&](int c) {
        for (std::size_t m = 0; m < comps.size(); ++m)
          if (comps[m] == c) return static_cast<int>(m);
        comps.push_back(c);
        return static_cast<int>(comps.size()) - 1;
      };
      struct MF {
        std::string name;
        std::vector<int> pos;  // the field's components in the streamed slab
        bool avg = false;      // gammaHI_avg (a stream of its own)
        std::vector<std::vector<float>> blocks;
        std::unique_ptr<F32Out> cube;
      };
      std::vector<MF> mf(mfields.size());
      for (std::size_t q = 0; q < mfields.size(); ++q) {
        const std::string& fld = mfields[q];
        MF& m = mf[q];
        m.name = fld;
        int c0 = -1, nc = 1;
        if (fld == "fHI") c0 = cHI;
        else if (fld == "fHeI") c0 = offsetof(CD, fHeI) / fl;
        else if (fld == "fHeII") c0 = offsetof(CD, fHeII) / fl;
        else if (fld == "temp") c0 = offsetof(CD, temp) / fl;
        else if (fld == "rhog") c0 = cRho;
        else if (fld == "ngam") c0 = offsetof(CD, ngam) / fl, nc = kNfreq;
        else if (fld == "gamma") c0 = cG, nc = 3;
        else if (fld == "gammaHI") c0 = cG;
        else if (fld == "gammaHI_avg") m.avg = true;
        else throw std::runtime_error("movie_fields: unknown field " + fld);
        for (int c = 0; c < nc && c0 >= 0; ++c) m.pos.push_back(pos_of(c0 + c));
        char name[64];
        if (P.movie_slab <= 0) {
          std::snprintf(name, sizeof name, "/movie/frame_%05d_%s.bin", dst.movie_frame, fld.c_str());
          m.cube = std::make_unique<F32Out>(odir + name, io);
        }
      }
      auto feed = [&](MF& m, const float* v, const StreamSlab& sl, int ncomp) {
        if (P.movie_slab > 0)
          slab_add(m.blocks, v, sl.k0, sl.nk, ncomp);
        else
          m.cube->write(v, static_cast<std::size_t>(sl.ncell() * ncomp));
      };
      FrameSums fs;
      gion.clear();
      std::vector<float> tmp;
      const int nc = static_cast<int>(comps.size());
      amrex_state_stream_rt(s, comps, [&](const StreamSlab& sl) {
        const std::int64_t n = sl.ncell();
        const float* a = sl.data;
        for (std::int64_t c = 0; c < n; ++c) {
          const float* p = a + c * nc;
          fs.add(p[1], p[0]);
          frame_gamma(fs, p[1], p[2]);
        }
        for (MF& m : mf) {
          if (m.avg) continue;
          const int fc = static_cast<int>(m.pos.size());
          tmp.resize(static_cast<std::size_t>(n * fc));
          for (std::int64_t c = 0; c < n; ++c)
            for (int q = 0; q < fc; ++q) tmp[static_cast<std::size_t>(c * fc + q)] = a[c * nc + m.pos[q]];
          feed(m, tmp.data(), sl, fc);
        }
      });
      for (MF& m : mf)
        if (m.avg) amrex_gamma_avg_stream(s, [&](const StreamSlab& sl) { feed(m, sl.data, sl, 1); });
      if (io) {
        for (MF& m : mf) {
          if (P.movie_slab > 0) {
            char name[64];
            std::snprintf(name, sizeof name, "/movie/frame_%05d_%s.slab", dst.movie_frame, m.name.c_str());
            F32Out f(odir + name, true);
            for (const std::vector<float>& b : m.blocks) f.write(b.data(), b.size());
            f.close();
          } else {
            m.cube->close();
          }
        }
        frame_line(ss, fs);
      }
      dst.movie_frame += 1;
      dst.movie_next -= P.movie_dz;
    };

    auto write_frame = [&](const Simulation& ss) {
      if (on_level) {
        write_frame_streamed(ss);
        return;
      }
      // gammaHI_avg on rt_backend = amrex: the level's running sums onto
      // rank 0 (collective)
      if (gavg_movie && use_amrex(cfg)) amrex_gamma_avg_field(s, gavg_buf);
      if (!io) {
        dst.movie_frame += 1;
        dst.movie_next -= P.movie_dz;
        return;
      }
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
        else if (fld == "gammaHI") v = &gamma_hi(G);
        else if (fld == "gammaHI_avg") {  // the running substep average of this hydro step
          if (!use_amrex(cfg)) gamma_avg_field(ss, gavg_buf);
          v = &gavg_buf;
        } else throw std::runtime_error("movie_fields: unknown field " + fld);
        if (P.movie_slab > 0) {
          std::snprintf(name, sizeof name, "/movie/frame_%05d_%s.slab", dst.movie_frame, fld.c_str());
          if (io) write_slabs(odir + name, *v);
        } else {
          std::snprintf(name, sizeof name, "/movie/frame_%05d_%s.bin", dst.movie_frame, fld.c_str());
          if (io) write_f32(odir + name, v->data(), v->size());
        }
      }
      FrameSums fs;
      gion.clear();
      for (std::int64_t c = 0; c < G.ncell; ++c) fs.add(G.fHI[c], G.rhog[c]);
      for (std::int64_t c = 0; c < G.ncell; ++c) frame_gamma(fs, G.fHI[c], G.gamma[c * 3]);
      frame_line(ss, fs);
      dst.movie_frame += 1;
      dst.movie_next -= P.movie_dz;
    };
    if (P.movie_dz > 0) {
      if (io) {
        std::filesystem::create_directories(odir + "/movie");
        if (P.movie_slab > 0) {
          if (P.movie_slab > mN)
            throw std::runtime_error("movie_slab = " + std::to_string(P.movie_slab) + " > N = " +
                                     std::to_string(mN));
          // the slab geometry, for the tools (tools/movies/movie_frames.py)
          std::ostringstream sd;
          sd << "# movie slab frames (movie_slab, docs/movies.md): frame_NNNNN_<field>.slab holds,\n"
                "# for each axis in 'axes' in turn, the cells start .. start+thick-1 (0-based) along\n"
                "# that axis and all N along the others: float32, Fortran order (ncomp, nx, ny, nz)\n"
             << "N = " << mN << "\nthick = " << mt << "\naxes = " << P.movie_slab_axes
             << "\nstart = " << ms0 << "\n";
          const std::string sp = odir + "/movie/slab.txt";
          if (!restart_from.empty() && std::filesystem::exists(sp)) {
            std::ifstream in(sp);
            std::stringstream prev;
            prev << in.rdbuf();
            if (prev.str() != sd.str())
              throw std::runtime_error(sp + " differs from movie_slab/movie_slab_axes of this run: "
                                       "keep them the same at every restart");
          } else {
            std::ofstream sout(sp);
            sout << sd.str();
          }
        }
        movielog.open(odir + "/movie/frames.txt", logmode);
      }
      if (restart_from.empty()) {
        if (P.movie_slab > 0)
          movielog << "# frame istep istep_rt redshift xHII_V xHII_M GammaHI_V GammaHI_ion_med "
                      "GammaHI_ion_mean fV_GammaHI_gt_1e-14 (frame_NNNNN_<field>.slab: "
                      "float32 slabs of rt_grid fields, geometry in slab.txt; xHII and GammaHI "
                      "[s^-1, last RT substep; ion: x_HII > 0.5] over the whole box)\n";
        else
          movielog << "# frame istep istep_rt redshift xHII_V xHII_M GammaHI_V GammaHI_ion_med "
                      "GammaHI_ion_mean fV_GammaHI_gt_1e-14 (frame_NNNNN_<field>.bin: "
                      "float32 (N,N,N) rt_grid fields; GammaHI [s^-1, last RT substep; ion: "
                      "x_HII > 0.5])\n";
        dst.movie_next = P.movie_z0 > 0 ? P.movie_z0 : zdata1.front();
      }
      s.substep_hook = [&](Simulation& ss) {
        if (P.movie_every_substep) {  // a frame after every RT substep
          if (use_amrex(cfg) && !on_level) amrex_rt_sync_fields(ss);
          write_frame(ss);
          return;
        }
        if (!(ss.clk.redshift <= dst.movie_next)) return;
        // rt_backend = amrex: the fields live on the level inside a hydro step
        // (amrex_state = 1: write_frame gathers them onto rank 0 itself)
        if (use_amrex(cfg) && !on_level) amrex_rt_sync_fields(ss);
        while (ss.clk.redshift <= dst.movie_next) write_frame(ss);
      };
    }

    // Restart: the complete state of the end of a step (checkpoint.hpp)
    if (!restart_from.empty()) {
      read_checkpoint(restart_from, s, hs, pm, dst);
      if (out) *out << " radhydro: restarted from " << restart_from << " at step " << s.clk.istep << "\n";
      calc_sim_var(cfg, s.clk, s.units);
      init_unsplit_state(cfg, s.units, s.gas, hs);  // no-op unless a legacy checkpoint
      if (on_level) {  // the grids from the file and rank 0's particles onto the levels
        read_checkpoint_grids(restart_from, s);
        amrex_pm_load_particles(pm);
        if (cfg.sf_overdensity > 0) amrex_pm_load_stars(pm);  // the sink's stars (version 8)
      }
      // Stage F: the refined hydro levels of the checkpoint (a checkpoint
      // without them starts with level 0 only and regrids at the first step)
      if (use_amrex_hydro(cfg) && cfg.amrex_max_level > 0 && std::filesystem::is_directory(restart_from + ".amr"))
        amrex_hydro_read_levels(cfg, cfg.Ngas, restart_from + ".amr");
      update_sched_x();  // the checkpoint's RT grid: the last step's x_HII,V
    } else {
      // Summarize the initial state
      calc_sim_var(cfg, s.clk, s.units);
      init_unsplit_state(cfg, s.units, s.gas, hs);  // hydro_scheme = unsplit: gas.e
      if (on_level) {  // the initial conditions onto the levels
        ics_to_level();
        amrex_pm_load_particles(pm);
      }
      summaries();
      const std::string d = step_dir();
      write_stats(d);
      update_sched_x();
      if (P.fields_every > 0) {
        write_fields(d);
        if (!use_amrex(cfg) || s.amrex) write_levels_slice(d);  // (the RT level is made at the first RT step)
      }
    }

    std::ofstream timing;
    if (io) timing.open(odir + "/timing.txt", logmode);
    if (restart_from.empty())
      timing << "# istep redshift nstep_rt Nsource nray_active t_hydro t_pm t_halos t_rt "
                "t_summaries t_io t_total [s, this step]\n";

    // ------------------------------------------------------------------
    // Time stepping
    // ------------------------------------------------------------------
    std::vector<double> pk_dm, pk_gas;
    int steps_run = 0;
    for (;;) {
      if (s.clk.istep >= P.max_steps)
        throw std::runtime_error("max_steps reached before the last zdata1 redshift");
      const Timers T0 = T;
      const auto ts0 = Clock::now();

      calc_sim_var(cfg, s.clk, s.units);
      apply_sched();
      // land_outputs = 1: the step ends on the next zdata1 output (config.hpp)
      s.cfg.land_a = (cfg.land_outputs && idata1 <= Ndata1) ? adata1[idata1 - 1] : 0.0;
      const TimestepInfo ti =
          calc_timestep(cfg, s.units, TimestepMode::Original, 0.0, s.clk, out);
      // Set IO flags (calc_timestep)
      if (ti.a2 >= adata1[idata1 - 1]) data1_out = true;
      if (idata2 <= Ndata2)
        if (ti.a2 >= adata2[idata2 - 1]) data2_out = true;

      // radhydro: forward hydro step with transposes in between (legacy),
      // or one unsplit step over dt_sim/2 (hydro_scheme = unsplit)
      auto th = Clock::now();
      nstar_step = 0;
      hydro_half_step(cfg, s.units, s.clk, HydroHalf::First, s.gas, hs);
      sink_to_stars();
      T.hydro += std::chrono::duration<double>(Clock::now() - th).count();

      s.clk.year = s.clk.year + s.clk.dt_yr / 2;
      s.clk.scalefactor = s.clk.scalefactor + s.clk.da_sim / 2;
      s.clk.redshift = 1 / s.clk.scalefactor - 1;
      calc_sim_var(cfg, s.clk, s.units);

      // Gravity for gas and dark matter (gas ordered x,y,z, same as DM)
      th = Clock::now();
      if (cfg.hydro_scheme == HydroScheme::Unsplit)
        pm::particle_mesh_unsplit(pm, s.gas.u.data(), s.gas.e.data(),
                                  cfg.dual_energy ? cfg.dual_eta1 : 0.0, s.clk, s.units);
      else
        pm::particle_mesh(pm, s.gas.u.data(), hs.v.data(), s.clk, s.units);
      T.pm += std::chrono::duration<double>(Clock::now() - th).count();
      xv_fresh = false;

      // Reverse hydro step with transposes in between (legacy), or the
      // second unsplit step over dt_sim/2
      th = Clock::now();
      hydro_half_step(cfg, s.units, s.clk, HydroHalf::Second, s.gas, hs);
      sink_to_stars();
      T.hydro += std::chrono::duration<double>(Clock::now() - th).count();

      s.clk.year = s.clk.year + s.clk.dt_yr / 2;
      s.clk.scalefactor = s.clk.scalefactor + s.clk.da_sim / 2;
      s.clk.redshift = 1 / s.clk.scalefactor - 1;
      calc_sim_var(cfg, s.clk, s.units);

      // thermal_init = recfast (docs/thermal_history.md): the CMB coupling
      // of the residual electrons over the step, at its middle
      if (cfg.thermal_init == 1) {
        th = Clock::now();
        Units um;
        SimClock cm = s.clk;
        cm.scalefactor = s.clk.scalefactor - s.clk.da_sim / 2;
        cm.redshift = 1 / cm.scalefactor - 1;
        calc_sim_var(cfg, cm, um);
        compton_residual_step(um, cm.scalefactor, s.clk.dt_sec, s.rates.coef_cool(14, 1), s.gas);
        T.hydro += std::chrono::duration<double>(Clock::now() - th).count();
      }

      // Ray tracing radiative transfer (after the first zdata1 output);
      // plane_parallel = 1 (plane_parallel.hpp): from plane_z_re instead,
      // the faces of the RT domains as the only sources (no halo finder)
      if (cfg.plane_parallel) {
        // before the turn-on the RT runs only for the X-ray background
        // (xray_heating), with no face sources
        s.plane.active = plane_rt_active(cfg, s.clk.scalefactor - s.clk.da_sim);
        if (s.plane.active || s.xray.on) {
          th = Clock::now();
          if (on_level && cfg.amrex_balance) {
            double e0 = 0, e1 = 0;
            const bool moved = amrex_balance_step(s, pm, e0, e1);
            balance_log(moved, e0, e1);
          }
          radiative_transfer(s);
          T.rt += std::chrono::duration<double>(Clock::now() - th).count();
        }
      } else if (idata1 > 1) {
        th = Clock::now();
        SourceCatalog cat;
        if (on_level) {
          // the halo finder on rank 0's particle mirror, the catalogue onto
          // every rank, the sources applied everywhere (halo_sources.hpp)
          sync_xv();
          if (io) cat = source_catalog(s, pm.xv.data(), pm.pc.Npart, P.halo, out);
          bcast_catalog(cat);
          apply_sources(s, cat, P.halo);
        } else {
          cat = update_sources(s, pm.xv.data(), pm.pc.Npart, P.halo, out);
        }
        Nhalo = static_cast<std::int64_t>(cat.halo.size());
        T.halos += std::chrono::duration<double>(Clock::now() - th).count();
        if (cat.Nsource > 0) {
          th = Clock::now();
          if (on_level && cfg.amrex_balance) {
            // amrex_balance = 1: the shared decomposition on this step's
            // cost (docs/stage_e_state.md section 5)
            double e0 = 0, e1 = 0;
            const bool moved = amrex_balance_step(s, pm, e0, e1);
            balance_log(moved, e0, e1);
          }
          if (P.tensor_dump && data1_out) tensor_diag_begin(s);  // diagnostic
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
        sync_mirror();
        if (on_level)
          write_data_files_streamed(pm::ztag(s.clk.redshift), false, P.data_files);
        else if (io)
          write_data_files(odir, pm::ztag(s.clk.redshift), s, pm.xv, hs.v, false, P.data_files);
        if (cfg.sf_overdensity > 0 && io) {
          write_stars(odir + "/DM/stars." + pm::ztag(s.clk.redshift), pm.stars);
          write_sink_summary(odir + "/sink.txt", pm);
        }
        idata2 = idata2 + 1;
        data2_out = false;
      }
      if (data1_out) {
        const auto tp = Clock::now();
        sync_mirror();
        // (amrex_state = 1: the host P(k) on rank 0's PM mesh; its FFT and
        // the DM linking are the host's, redundant until Stage E's P(k) on
        // the level)
        if (on_level) {
          // the gas density streamed onto rank 0 (its only grid array here;
          // the DM part still uses rank 0's host PM mesh,
          // docs/rank0_streaming.md step 11)
          sync_xv();
          std::vector<float> rho(io ? static_cast<std::size_t>(cfg.Ncell) : 0);
          amrex_state_stream_gas(s, kStreamGasU, 1, [&](const StreamSlab& sl) {
            std::copy(sl.data, sl.data + sl.ncell(), rho.begin() + sl.c0());
          });
          if (io) pm::calc_power_spectrum(pm, rho.data(), cfg.box, pk_dm, pk_gas, 1);
        } else {
          pm::calc_power_spectrum(pm, s.gas.u.data(), cfg.box, pk_dm, pk_gas);
        }
        if (io)
          pm::write_power_spectrum(odir + "/Pk/pk." + pm::ztag(s.clk.redshift), pm, pk_dm,
                                   pk_gas);
        T.pk += std::chrono::duration<double>(Clock::now() - tp).count();
        // write_data_1 (F1 fix: all files)
        if (on_level)
          write_data_files_streamed(pm::ztag(s.clk.redshift), true, P.data_files);
        else if (io)
          write_data_files(odir, pm::ztag(s.clk.redshift), s, pm.xv, hs.v, true, P.data_files);
        if (cfg.sf_overdensity > 0 && io) {
          write_stars(odir + "/DM/stars." + pm::ztag(s.clk.redshift), pm.stars);
          write_sink_summary(odir + "/sink.txt", pm);
        }
        idata1 = idata1 + 1;
        data1_out = false;
        if (P.fields_at_data1) fields = true;
        if (s.tdiag && s.tdiag->active) {
          if (io) std::filesystem::create_directories(d);
          tensor_diag_write(s, d, io);
        }
      }
      const bool last = s.clk.redshift <= zdata1[Ndata1 - 1];
      write_stats(d);
      update_sched_x();
      if (fields || last) {
        sync_mirror();
        write_fields(d);
        write_levels_slice(d);
      }
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
        if (quiet && io) {
          std::printf(" step %4d  z = %9.5f  nstep_rt = %5d  Nsource = %5lld  %.3f s\n",
                      s.clk.istep, s.clk.redshift, s.clk.nstep_rt,
                      static_cast<long long>(s.Nsource), tstep);
          std::fflush(stdout);
        }
        // the backends' phase timers so far (cumulative; per step by difference)
        if (P.step_timers && io) {
          std::printf(" step %4d timers:", s.clk.istep);
          if (use_amrex_pm(cfg)) std::printf(" %s;", amrex_pm_timers(pm).c_str());
          if (use_amrex_hydro(cfg)) std::printf(" %s;", amrex_hydro_timers().c_str());
          std::printf(" halos %.3f, P(k) %.3f\n", T.halos, T.pk);
          std::fflush(stdout);
        }
      }

      // Frames not reached inside RT substeps (steps without RT)
      if (P.movie_dz > 0)
        while (s.clk.redshift <= dst.movie_next) write_frame(s);

      // Checkpoint (complete state at the end of this step)
      if (P.restart_every > 0 && (s.clk.istep % P.restart_every == 0 || last)) {
        if (io) std::filesystem::create_directories(odir + "/restart");
        backend_barrier();
        sync_mirror();
        char name[64];
        if (P.restart_keep_latest && P.restart_keep == 0)
          std::snprintf(name, sizeof name, "/restart/checkpoint.bin");
        else
          std::snprintf(name, sizeof name, "/restart/checkpoint_%04d.bin", s.clk.istep);
        write_checkpoint(odir + name, s, hs, pm, dst);
        // Stage F: the refined hydro levels next to it (docs/stage_f.md)
        if (use_amrex_hydro(cfg) && cfg.amrex_max_level > 0) amrex_hydro_write_levels(odir + name + ".amr");
        if (P.restart_keep > 0 && io) {
          // keep the newest restart_keep numbered checkpoints (names sort by step)
          std::vector<std::string> ck;
          for (const auto& e : std::filesystem::directory_iterator(odir + "/restart")) {
            const std::string f = e.path().filename().string();
            if (f.size() == 19 && f.compare(0, 11, "checkpoint_") == 0 &&
                f.compare(15, 4, ".bin") == 0)
              ck.push_back(e.path().string());
          }
          std::sort(ck.begin(), ck.end());
          for (std::size_t i = 0; i + P.restart_keep < ck.size(); ++i) {
            std::filesystem::remove(ck[i]);
            std::filesystem::remove_all(ck[i] + ".amr");  // its refined levels, if any
          }
        }
      }

      // Stop
      if (last) break;
      if (P.run_steps > 0 && ++steps_run >= P.run_steps) break;
    }

    backend_barrier();
    if (cfg.sf_overdensity > 0) {  // the sink's final count and mass (sink.txt)
      amrex_pm_gather_stars(pm);
      if (io) write_sink_summary(odir + "/sink.txt", pm);
    }
    const double total = std::chrono::duration<double>(Clock::now() - t0).count();
    // Peak RSS of this rank, and its range over the ranks of an MPI run
    const double rss = peak_rss_mb();
    double rss_lo = rss, rss_hi = rss;
    backend_reduce_minmax(rss_lo, rss_hi);
    if (io) {
      std::printf(
          " radhydro: %d steps, total %.2f s: hydro %.2f, pm %.2f (+ P(k) %.2f), halos %.2f, "
          "rt %.2f, summaries %.2f, io %.2f; peak RSS %.1f MB\n",
          s.clk.istep, total, T.hydro, T.pm, T.pk, T.halos, T.rt, T.summaries, T.io, rss);
      if (backend_nranks() > 1)
        std::printf(" peak RSS over %d ranks: min %.1f MB, max %.1f MB\n", backend_nranks(), rss_lo,
                    rss_hi);
      std::printf(" %s\n ray pool guard: split limited in %lld substeps, merging forced in %lld\n",
                  format_rt_timers(s.timers).c_str(), static_cast<long long>(s.rays.npool_limited),
                  static_cast<long long>(s.rays.npool_forced));
      if (use_amrex(cfg)) std::printf(" %s\n", amrex_rt_timers(s).c_str());
      if (use_amrex_pm(cfg)) std::printf(" %s\n", amrex_pm_timers(pm).c_str());
      if (use_amrex_hydro(cfg)) std::printf(" %s\n", amrex_hydro_timers().c_str());
      if (cfg.amrex_balance) std::printf(" balance: remapped in %lld of %lld RT steps\n",
                                         static_cast<long long>(balance_moves), static_cast<long long>(s.clk.istep));
    }
    if (cfg.hydro_scheme == HydroScheme::Unsplit && io) {
      const UnsplitCounters& U = hs.ucounters;
      std::printf(" unsplit hydro: %lld calls, %lld substeps (max %lld per call), dt_cfl min %.4g; "
                  "cell-substeps with p from e %lld (E reset %lld), synced %lld; floors T %lld, "
                  "rho %lld\n",
                  static_cast<long long>(U.calls), static_cast<long long>(U.substeps),
                  static_cast<long long>(U.max_substeps), U.dt_cfl_min,
                  static_cast<long long>(U.cells_from_e), static_cast<long long>(U.cells_reset),
                  static_cast<long long>(U.cells_synced), static_cast<long long>(U.floor_T),
                  static_cast<long long>(U.floor_rho));
      timing << "# unsplit calls " << U.calls << " substeps " << U.substeps << " max_substeps "
             << U.max_substeps << " cells_from_e " << U.cells_from_e << " cells_reset "
             << U.cells_reset << " cells_synced " << U.cells_synced << " floor_T " << U.floor_T
             << " floor_rho " << U.floor_rho << "\n";
    }
    timing << "# " << format_rt_timers(s.timers) << "\n";
    if (use_amrex_hydro(cfg)) timing << "# " << amrex_hydro_timers() << "\n";
    timing << "# total " << total << " hydro " << T.hydro << " pm " << T.pm << " pk " << T.pk
           << " halos " << T.halos << " rt " << T.rt << " summaries " << T.summaries << " io "
           << T.io << " peak_rss_MB " << rss;
    if (backend_nranks() > 1)
      timing << " ranks " << backend_nranks() << " peak_rss_min_MB " << rss_lo << " peak_rss_max_MB "
             << rss_hi;
    timing << "\n";
    if (P.halo.fof_dt_myr > 0)
      timing << "# fof_dt_myr " << P.halo.fof_dt_myr << " fof_calls " << s.fof.nfof
             << " reused " << s.fof.nreuse << "\n";
  } catch (const std::exception& e) {
    std::fprintf(stderr, "radhydro: error: %s (peak RSS %.1f MB)\n", e.what(), peak_rss_mb());
    backend_abort(1);  // more than one rank: end the job (the others would wait)
    return 1;
  }
  return 0;
}
