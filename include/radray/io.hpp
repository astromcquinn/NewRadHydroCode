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
