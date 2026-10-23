// Exact sums across ranks for the AMReX backends (internal to src/amrex):
// 128-bit fixed point, whose sums are associative, so a statistic does not
// depend on the decomposition or the thread count (rt_level.cpp, the PM
// level's kick statistics; docs/stage_d.md, docs/stage_e_moments.md).
#pragma once
#include <cmath>
#include <cstddef>
#include <vector>

#include <AMReX_INT.H>
#include <AMReX_ParallelDescriptor.H>

namespace radray {
namespace xs {

using i128 = __int128;

// x * 2^bits rounded to the nearest integer (resolution 2^-bits)
inline i128 to_fx128(double x, int bits) { return static_cast<i128>(std::floor(std::ldexp(x, bits) + 0.5)); }
inline double from_fx128(i128 v, int bits) { return std::ldexp(static_cast<double>(v), -bits); }

// MPI sum of 128-bit integers as four 32-bit limbs in int64 (each limb sum
// fits for < 2^31 ranks); two's complement makes the limb sum of signed
// values the signed sum modulo 2^128.
inline void allreduce_i128(i128* v, int n) {
  namespace ax = amrex;
  if (ax::ParallelDescriptor::NProcs() == 1) return;
  std::vector<ax::Long> limbs(static_cast<std::size_t>(4 * n));
  for (int i = 0; i < n; ++i) {
    const unsigned __int128 u = static_cast<unsigned __int128>(v[i]);
    for (int k = 0; k < 4; ++k) limbs[4 * i + k] = static_cast<ax::Long>((u >> (32 * k)) & 0xffffffffULL);
  }
  ax::ParallelDescriptor::ReduceLongSum(limbs.data(), 4 * n);
  for (int i = 0; i < n; ++i) {
    unsigned __int128 r = 0;
    for (int k = 0; k < 4; ++k)
      r += static_cast<unsigned __int128>(static_cast<unsigned long long>(limbs[4 * i + k])) << (32 * k);
    v[i] = static_cast<i128>(r);
  }
}

}  // namespace xs
}  // namespace radray
