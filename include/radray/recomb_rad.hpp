// Case A recombination and recombination radiation transported by the rays
// (docs/recomb_radiation.md). New options, not in the original: with the
// defaults (recomb_case = B, recomb_rad = 0) nothing here runs and the code
// is the original.
//
// recomb_rad = 1: photo_cell counts each cell's ground-state recombinations
// alpha_A - alpha_B (HII -> group 1, HeII -> group 3, HeIII -> group 5, the
// mapping of the original's disabled `if (.false.)` block) into a per-cell
// buffer. At the next RT substep (recomb_prepare, in rt_step_1) the buffers
// are summed over blocks of recomb_rad_block^3 cells. A block either
// launches its photons as 12 nside^2 rays from one of its cells
// (ray_source_recomb, in ray_tracing_1), or, when it is optically thick or
// the ray pool has no room, each cell absorbs its own photons on the spot:
// they are added to the cell's photon field before photo_cell, as the
// optically thick branch of ray_source does for stars.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "radray/config.hpp"

namespace radray {

struct Simulation;

// Parses recomb_case, recomb_rad, recomb_rad_nside, recomb_rad_block and
// recomb_rad_pool; false for any other key.
bool set_recomb_param(Config& c, const std::string& key, const std::string& val);
void copy_recomb_param(const Config& from, Config& to);
// Throws for an inconsistent combination (recomb_rad without case A, ...).
void check_recomb_param(const Config& c);

struct RecombState {
  // Ground-state recombination photons of the last substep per RT cell and
  // frequency group (kNfreq per cell, only groups 1, 3, 5 used), in the
  // units of rt_grid%ngam (per <nH> per RT cell volume), and the
  // photon-weighted energy ngam*nu [eV]. Carried from photo_cell to the next
  // rt_step_1; part of the checkpoint (version 4).
  std::vector<float> ngam, egam;
  // Photons to launch from this cell in the current substep (filled by
  // recomb_prepare, emptied by ray_source_recomb).
  std::vector<float> emit_ngam, emit_egam;
  std::vector<std::uint8_t> emit;  // 1: the cell launches rays this substep
  // Budget since the start, in photons per H atom of the box: produced in
  // photo_cell, launched as rays, absorbed on the spot (by optical depth or
  // for lack of pool room).
  double produced = 0, launched = 0, local_thick = 0, local_pool = 0;
  std::int64_t nemit_last = 0, nblock_last = 0;  // emitting / bright blocks, last substep
  bool active() const { return !ngam.empty(); }
};

void recomb_init(Simulation& s);
// rt_step_1, after calc_rt_1 and the ray histogram: chooses the emitting
// blocks within the pool budget (live rays L, stellar source cells nsrc,
// rays created by the expected stellar splits nsplit) and applies the
// on-the-spot deposits. Returns the number of rays that ray_tracing_1 will
// create for recombination radiation.
std::int64_t recomb_prepare(Simulation& s, std::int64_t L, std::int64_t nsrc, std::int64_t nsplit);
// ray_tracing_1: launches the recombination rays of cell (irt,jrt,krt).
void ray_source_recomb(Simulation& s, int irt, int jrt, int krt);
// photo_cell: adds the ground-state recombinations nrab[3] (per cm^3 of gas,
// this substep) of one gas cell of temperature T [K] to RT cell c.
void recomb_count(Simulation& s, std::int64_t c, const double nrab[3], double T,
                  double vol_gas);

}  // namespace radray
