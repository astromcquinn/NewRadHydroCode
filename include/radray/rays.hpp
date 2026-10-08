// Ray storage: ray(Nray_max), ray_index(Nray_max), ray_stat/ray_save
// (radhydro_global.f90: ray_type, ray_stat_type) and the free-list buffer
// logic of ray_tracing_tools.f90.
//
// The Fortran allocator is reproduced exactly (ray_index permutation,
// per-thread free buffer ifree, Nray/Nrbuf bookkeeping). Results only depend
// on chain ORDER, not on storage indices, except through one quirk that we
// must keep: rt_step_1/calc_ray_1 histograms every storage slot with
// nside > 0, which includes freed-but-not-yet-zeroed rays sitting in the
// free buffer. Reproducing the allocator makes that quirk (and therefore
// the ngam_ray thresholds) identical, and gives identical storage indices,
// which is handy for debugging against the Fortran.
//
// All ray storage indices `ir` are Fortran 1-based; slot 0 is a dummy
// (ir == 0 means "no ray" in the chains).
//
// Memory: the Fortran allocates Nray_max = 64*N^3 rays up front. We reserve
// that capacity (virtual memory only) and grow the touched part lazily; it
// can be shown that only slots <= max(Nray) ever get handed out, because
// ray_index(1:hwm) is always a permutation of 1..hwm.
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "radray/constants.hpp"

namespace radray {

// radhydro_global.f90: ray_type
struct Ray {
  std::int32_t nside = 0, ipix = 0;
  std::int64_t next = 0, prev = 0;
  float r = 0, dr = 0;
  float x[3] = {0, 0, 0};
  float dx[3] = {0, 0, 0};
  float nu[kNfreq] = {0, 0, 0, 0, 0};
  float ngam[kNfreq] = {0, 0, 0, 0, 0};
};

// radhydro_global.f90: ray_stat_type
struct RayStat {
  std::int64_t ncst1 = 0, ncst2 = 0;
  std::int64_t nmrg1 = 0, nmrg2 = 0;
  std::int64_t nspl1 = 0, nspl2 = 0;
  std::int64_t ndel1 = 0;
  std::int64_t nfree = 0;
  // ifree(1:nfree_max); element 0 unused (1-based like Fortran)
  std::vector<std::int64_t> ifree;
};

struct RayState {
  std::int64_t Nray_max = 0;
  std::int64_t Nfree_max = 0;

  std::int64_t Nray = 0;   // number of ray_index entries handed out
  std::int64_t Nrbuf = 0;  // buffer size set in prep_rt

  std::vector<Ray> ray;                 // ray(ir), ir >= 1
  std::vector<std::int64_t> ray_index_;  // ray_index(i), i >= 1

  RayStat stat;  // ray_stat(1)   (Ncpu = 1)
  RayStat save;  // ray_save

  std::array<double, 64> ngam_ray{};  // ngam_ray(1:64) -> [0..63]

  void init(std::int64_t nray_max, std::int64_t nfree_max);

  // Highest storage slot materialised so far.
  std::int64_t hwm() const { return static_cast<std::int64_t>(ray.size()) - 1; }

  // Grow storage so that slots/positions 1..i exist (identity ray_index,
  // zeroed rays). Never reallocates (capacity reserved in init), so
  // references to rays stay valid.
  void ensure(std::int64_t i);

  // ray_index(i) (Fortran 1-based). Untouched positions hold i.
  std::int64_t ray_index(std::int64_t i) {
    ensure(i);
    return ray_index_[i];
  }
  void set_ray_index(std::int64_t i, std::int64_t v) {
    ensure(i);
    ray_index_[i] = v;
  }

  Ray& at(std::int64_t ir) { return ray[ir]; }
  const Ray& at(std::int64_t ir) const { return ray[ir]; }
};

}  // namespace radray
