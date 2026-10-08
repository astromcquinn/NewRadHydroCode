// Optional ray snapshots for the ray-compression study (docs/ray_compression.md).
//
// Off by default. When enabled, the dump
//  * keeps a provenance side array (RayState::prov, one RayProv per storage
//    slot) that ray_source / ray_cast_split / ray_cast_merge fill in, and
//  * after selected RT substeps (after rt_step_3, i.e. the state that would
//    be stored between substeps) writes every live ray (cell chains, in
//    cell order then chain order) with cell context to
//    <dir>/rays_s<substep>.npy plus a JSON sidecar <dir>/rays_s<substep>.json.
//
// Neither touches anything the RT reads, so results are bit-identical with
// the dump on or off (checked by the ray_dump_* CTests).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace radray {

struct Simulation;

// One record of the .npy file (numpy structured dtype in ray_dump.cpp).
struct RayDumpRecord {
  std::int32_t cell;      // 0-based linear RT cell of the chain holding the ray
  std::int32_t chainpos;  // position in that cell's chain (0 = head)
  // Raw ray fields (bit copies of radray::Ray, without the links)
  std::int32_t nside, ipix;
  float r, dr;
  float x[3], dx[3];
  float nu[5], ngam[5];
  // Provenance (RayProv)
  std::int32_t src;     // emitting source cell (0-based), -1 if mixed
  std::int32_t birth;   // global RT substep at emission (of the oldest ancestor)
  std::uint16_t nmerge; // number of merges in the ancestry (saturating)
  std::uint8_t nclamp;  // splits that clamped r (moved the effective origin)
  std::uint8_t flags;   // RayProv flags
  // Cell context (the ray's chain cell, after rt_step_3)
  float c_tau[5], c_fgam[5], c_ngam[5];
  float c_fHI, c_fHeI, c_fHeII, c_temp, c_rhog;
};
static_assert(sizeof(RayDumpRecord) == 4 * 45, "RayDumpRecord must be unpadded");

struct RayDumpConfig {
  std::string dir;                // empty: dumping off
  std::int64_t every = 0;         // dump every n-th global RT substep (0: never)
  std::vector<std::int64_t> at;   // and at these global RT substeps (1-based)
  bool enabled() const { return !dir.empty() && (every > 0 || !at.empty()); }
};

// Parses a list "10 20,30" of substep numbers.
std::vector<std::int64_t> parse_substep_list(const std::string& s);

// Turns on provenance tracking; call after init_simulation and before the
// first RT step. Stores cfg in s.raydump.
void ray_dump_enable(Simulation& s, const RayDumpConfig& cfg);

// Called by radiative_transfer after each rt_step_3 when s.raydump is set.
void ray_dump_after_substep(Simulation& s);

// Writes one snapshot of the current live rays (used by the above).
void ray_dump_write(const Simulation& s, const std::string& stem);

struct RayDump {
  RayDumpConfig cfg;
  std::int64_t substep = 0;  // global RT substep counter (1-based after the first)
  std::int64_t ndumped = 0;
};

}  // namespace radray
