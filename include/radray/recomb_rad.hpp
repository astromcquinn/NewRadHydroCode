// Case A recombination and recombination radiation transported by the rays
// (docs/recomb_radiation.md). New options, not in the original: with the
// defaults (recomb_case = B, recomb_rad = 0) nothing here runs and the code
// is the original.
//
// recomb_rad = 1: every RT cell is, for one substep, in one of two modes,
// chosen per block of recomb_rad_block^3 cells at the start of the substep
// (recomb_prepare, in rt_step_1):
//  * launch: photo_cell uses case A and counts the cell's ground-state
//    recombinations alpha_A - alpha_B (HII -> group 1, HeII -> group 3,
//    HeIII -> group 5, the mapping of the original's disabled
//    `if (.false.)` block) into a buffer; at the next substep the block
//    launches them as 12 nside^2 rays from one of its cells, rotating every
//    substep (ray_source_recomb, in ray_tracing_1);
//  * on the spot: photo_cell uses case B, i.e. the ground-state
//    recombination photons are absorbed where they are made (counted only).
// A block is put on the spot when its launching cell is optically thick,
// when its photons per ray would be below ngam_min, or when the ray pool has
// no room for it (the brightest blocks, by their production in the last
// substep, get the room first).
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
  // Ground-state recombination photons of the cells in launch mode, made in
  // the last substep, per RT cell and frequency group (kNfreq per cell, only
  // groups 1, 3, 5 used), in the units of rt_grid%ngam (per <nH> per RT cell
  // volume), and the photon-weighted energy ngam*nu [eV]. Carried from
  // photo_cell to the next rt_step_1.
  std::vector<float> ngam, egam;
  // The same count for the cells on the spot (all groups summed), for the
  // budget and the next choice of modes.
  std::vector<float> ots;
  // Mode of each cell in the current substep: 1 launch (case A), 0 on the
  // spot (case B).
  std::vector<std::uint8_t> mode;
  // Photons to launch from this cell in the current substep (filled by
  // recomb_prepare, emptied by ray_source_recomb).
  std::vector<float> emit_ngam, emit_egam;
  std::vector<std::uint8_t> emit;  // 1: the cell launches rays this substep
  // Budget since the start, in photons per H atom of the box: launched as
  // rays, dropped (a launch-mode block whose photons per ray came out below
  // ngam_min), and absorbed on the spot.
  double launched = 0, dropped = 0, onspot = 0;
  std::int64_t nemit_last = 0, nlaunch_last = 0;  // emitting blocks / blocks in launch mode
  bool active() const { return !ngam.empty(); }
};

void recomb_init(Simulation& s);
// rt_step_1, after calc_rt_1 and the ray histogram: moves the buffers of the
// last substep's launch-mode blocks to their emitting cells, and chooses the
// modes of this substep within the pool budget (live rays L, stellar source
// cells nsrc, rays of the expected stellar splits nsplit). Returns the number
// of rays that ray_tracing_1 will create for recombination radiation.
std::int64_t recomb_prepare(Simulation& s, std::int64_t L, std::int64_t nsrc, std::int64_t nsplit);
// ray_tracing_1: launches the recombination rays of cell (irt,jrt,krt).
void ray_source_recomb(Simulation& s, int irt, int jrt, int krt);
// photo_cell: whether RT cell c uses case B this substep (on the spot).
inline bool recomb_onspot(const RecombState& Rc, std::int64_t c) { return Rc.mode[c] == 0; }
// photo_cell: adds the ground-state recombinations nrab[3] (per cm^3 of gas,
// this substep) of one gas cell of temperature T [K] to RT cell c (to the
// launch buffer or to the on-the-spot count, by the cell's mode).
void recomb_count(Simulation& s, std::int64_t c, const double nrab[3], double T,
                  double vol_gas);

}  // namespace radray
