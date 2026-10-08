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
// Memory: the Fortran allocates Nray_max = 64*N^3 rays up front. We map
// that capacity (calloc: virtual memory only, zero pages) and touch it
// lazily; only slots <= max(Nray) ever get handed out, because
// ray_index(1:hwm) is always a permutation of 1..hwm.
#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <vector>

#include "radray/constants.hpp"

namespace radray {

// radhydro_global.f90: ray_type. The links are int32 here (int8 in the
// original): Config::make refuses Nray_max >= 2^31. This saves 8 B per ray
// (88 B instead of 96 B) and changes no result; see docs/parallel.md.
struct Ray {
  std::int32_t nside = 0, ipix = 0;
  std::int32_t next = 0, prev = 0;
  float r = 0, dr = 0;
  float x[3] = {0, 0, 0};
  float dx[3] = {0, 0, 0};
  float nu[kNfreq] = {0, 0, 0, 0, 0};
  float ngam[kNfreq] = {0, 0, 0, 0, 0};
};
static_assert(sizeof(Ray) == 88, "Ray layout");

// radhydro_global.f90: ray_stat_type
struct RayStat {
  std::int64_t ncst1 = 0, ncst2 = 0;
  std::int64_t nmrg1 = 0, nmrg2 = 0;
  std::int64_t nspl1 = 0, nspl2 = 0;
  std::int64_t ndel1 = 0;
  std::int64_t nfree = 0;
  // ifree(1:nfree_max); element 0 unused (1-based like Fortran). int32
  // storage of the int8 Fortran indices (values < 2^31).
  std::vector<std::int32_t> ifree;
};

// Per-thread ray bookkeeping of the deterministic parallel ray tracing
// (deterministic_parallel = 1): a local free list refilled from / spilled to
// the shared pool, and the ray_stat counters of the thread.
struct ThreadRays {
  std::vector<std::int32_t> free;
  std::int64_t ncst1 = 0, ncst2 = 0, nmrg1 = 0, nmrg2 = 0, nspl1 = 0, nspl2 = 0;
  std::int64_t ndel1 = 0, ncreated = 0;
  char pad[64];  // keep the counters of different threads on different lines
};

// Zero-initialised storage for Nray_max+1 rays. calloc of a large block
// maps zero pages lazily, so only the touched part costs memory (as the
// lazily grown std::vector did before).
struct RayBuffer {
  Ray* p = nullptr;
  std::int64_t cap = 0;
  void reset(std::int64_t n);
  ~RayBuffer() { std::free(p); }
  RayBuffer() = default;
  RayBuffer(const RayBuffer&) = delete;
  RayBuffer& operator=(const RayBuffer&) = delete;
  Ray& operator[](std::int64_t i) { return p[i]; }
  const Ray& operator[](std::int64_t i) const { return p[i]; }
};

struct RayState {
  std::int64_t Nray_max = 0;
  std::int64_t Nfree_max = 0;

  // Fortran-emulating pool (deterministic_parallel = 0): number of ray_index
  // entries handed out. With deterministic_parallel = 1: the number of live
  // rays (updated after each sweep).
  std::int64_t Nray = 0;
  std::int64_t Nrbuf = 0;  // buffer size set in prep_rt

  RayBuffer ray;                         // ray(ir), ir >= 1; slot 0 dummy
  std::vector<std::int32_t> ray_index_;  // ray_index(i), i >= 1

  RayStat stat;  // ray_stat(1)   (Ncpu = 1)
  RayStat save;  // ray_save

  std::array<double, 64> ngam_ray{};  // ngam_ray(1:64) -> [0..63]

  // F3 pool guard (legacy_rt = 0): substeps in which the split threshold was
  // raised / splitting stopped and merging forced (cumulative)
  std::int64_t npool_limited = 0, npool_forced = 0;

  // deterministic_parallel = 1: shared pool of free slots (beyond the
  // per-thread lists) and the next never-used slot.
  std::mutex pool_mtx;
  std::vector<std::int32_t> pool_free;
  std::int64_t next_fresh = 1;
  std::vector<ThreadRays> tl;

  void init(std::int64_t nray_max, std::int64_t nfree_max);

  // Highest storage slot materialised so far (slots above it are zero).
  std::int64_t hwm() const { return hwm_; }

  // Make positions/slots 1..i exist (identity ray_index, zeroed rays).
  void ensure(std::int64_t i);

  // ray_index(i) (Fortran 1-based). Untouched positions hold i.
  std::int64_t ray_index(std::int64_t i) {
    ensure(i);
    return ray_index_[i];
  }
  void set_ray_index(std::int64_t i, std::int64_t v) {
    ensure(i);
    ray_index_[i] = static_cast<std::int32_t>(v);
  }

  Ray& at(std::int64_t ir) { return ray[ir]; }
  const Ray& at(std::int64_t ir) const { return ray[ir]; }

  // deterministic_parallel: slot management (thread-safe)
  void det_take(ThreadRays& t, int nr, std::int64_t* indx);
  void det_give(ThreadRays& t, std::int64_t ir);
  void det_set_hwm(std::int64_t h) { hwm_ = h; }

 private:
  std::int64_t hwm_ = 0;
};

}  // namespace radray
