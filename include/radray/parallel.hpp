// OpenMP helpers. The code builds and runs without OpenMP (RADRAY_OPENMP=OFF
// or a compiler without it); every parallel loop is then serial.
//
// Determinism rules used throughout (docs/parallel.md):
//   * a parallel loop may only write data owned by its iteration (a pencil, a
//     plane, a cell, a domain of independent cells), or per-thread buffers;
//   * floating-point reductions are never done with OpenMP reduction
//     clauses: the addends are recorded per iteration and summed afterwards
//     serially, in the serial loop order, so the result is the same bits for
//     any thread count (and the same as the original serial code);
//   * integer sums and max/min reductions are order independent and may be
//     combined from per-thread partial results.
#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

namespace radray {
namespace par {

inline int max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

// Threads of the current team (1 outside a parallel region)
inline int num_threads() {
#ifdef _OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}

inline int thread_id() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

// Sets the number of threads of later parallel regions (n <= 0: leave the
// OpenMP default, i.e. OMP_NUM_THREADS or the number of cores).
inline void set_threads(int n) {
#ifdef _OPENMP
  if (n > 0) omp_set_num_threads(n);
#else
  (void)n;
#endif
}

// Whether a loop over `cells` grid cells is worth a parallel region (small
// grids: the fork/join cost of thousands of RT substeps dominates). Only
// affects speed, never results.
inline bool worth(long long cells) { return cells >= 32LL * 32 * 32; }

inline bool enabled() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

}  // namespace par
}  // namespace radray
