// Complete binary checkpoint of a radhydro run, for a bit-identical restart
// (docs/production_run.md). Unlike the original's write_data_1 files (which
// do not hold the ray free lists, the stale thresholds, the RT statistics or
// the clock), this holds every piece of state the time loop reads, so a run
// restarted from step n produces exactly the files of an uninterrupted run
// from step n+1 on, in every mode (legacy and fixed, deterministic_parallel
// or not, any thread count).
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "radray/hydro.hpp"
#include "radray/pm.hpp"
#include "radray/simulation.hpp"

namespace radray {

// Driver state of apps/radhydro that is not part of Simulation/HydroState/PM
struct DriverState {
  int idata1 = 1, idata2 = 1;
  bool data1_out = false, data2_out = false;
  std::int64_t Nhalo = 0;
  std::int32_t movie_frame = 0;
  double movie_next = 0;
};

void write_checkpoint(const std::string& path, const Simulation& s, const HydroState& hs,
                      const pm::PMState& pm, const DriverState& d);

// Reads a checkpoint into an initialised Simulation/HydroState/PMState of the
// same sizes (N, Ndm, Nray_max are checked). Returns its version.
int read_checkpoint(const std::string& path, Simulation& s, HydroState& hs, pm::PMState& pm,
                     DriverState& d);

}  // namespace radray
