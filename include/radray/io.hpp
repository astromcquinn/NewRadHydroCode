// Problem input and validation output (docs/validation_spec.md).
#pragma once

#include <map>
#include <string>
#include <vector>

#include "radray/rt.hpp"
#include "radray/simulation.hpp"

namespace radray {

struct ProblemParams {
  int N = 32;
  double z_init = 10;
  double da_sim = 1e-4;
  int nsteps = 10;
  int nray_max_factor = 64;
  int legacy_rt = 1;  // 1: original RT bookkeeping; 0: fixes B1-B4, B8, B9
  int legacy_physics = 1;  // 1: original atomic physics; 0: fixes P1-P6
  // Ray snapshots (ray_dump.hpp); off unless ray_dump_dir is set and
  // ray_dump_every > 0 or ray_dump_at lists global RT substeps.
  std::string ray_dump_dir;
  long long ray_dump_every = 0;
  std::string ray_dump_at;
  // Lossy ray-precision experiment (ray_quant.hpp); research only, off.
  int ray_quant_bits = 0;
  int ray_quant_dither = 1;
  int ray_quant_struct = 0;
  int ray_quant_offbits = 16;
  double c_reduced = 1.0;  // reduced speed of light (calc_timestep cscale)
  int deterministic_parallel = 0;  // 1: parallel ray tracing (docs/parallel.md)
  int threads = 0;                 // OpenMP threads (0: OpenMP default)
  int ray_rotate = 0;              // B16 fix: Config::ray_rotate (0 = original)
  long long ray_rotate_seed = 0;   // Config::ray_rotate_seed
  int redshift_per_distance = 0;   // Config::redshift_per_distance (0 = original)
  int rt_cells_per_substep = 1;    // Config::rt_cells_per_substep (1 = original)
  int rt_substep_carry = 0;        // Config::rt_substep_carry (0 = original; legacy_rt = 0)
  int photon_conserve = -1;        // Config::photon_conserve (B17): -1 = !legacy_physics, else 0 | 1
  int ray_deposit = 0;             // Config::ray_deposit (0 = share, the original; 1 = exposure)
  // rt_method = rays | m1 | otvet and the moment-solver parameters
  // (set_rt_method_param, rt_moments.hpp; docs/rt_moments.md)
  Config rtm;
};

// params.txt: `key = value` lines; '#' starts a comment.
ProblemParams read_params(const std::string& path);

// Sets one params.txt entry (also used for the driver's --set key=value).
// Unknown keys are ignored, as in the Fortran harness.
void set_param(ProblemParams& p, const std::string& key, std::string val);

// gas.bin: float32 (9,N,N,N) -> u_gas(1:5), w_gas(1:4)
void read_gas(const std::string& path, GasGrid& gas);

// sources.txt: `i j k Luv sfr` per line
std::vector<Source> read_sources(const std::string& path);

// Writes <dir>/<field>.bin for all spec fields and <dir>/stats.txt.
void write_step_output(const std::string& dir, const Simulation& s);

void write_raw_f32(const std::string& path, const std::vector<float>& v);

}  // namespace radray
