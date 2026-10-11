// Device-usable copies of the array constants of radray/constants.hpp
// (moved from spike/common; this copy is the maintained one).
// nvcc cannot read namespace-scope constexpr std::arrays from device code;
// function-local constexpr arrays initialised from them can be.
#pragma once

#include "radray/constants.hpp"
#include "radray/kernels/rr_hd.hpp"

namespace rrk {

RR_HD double freq_init(int n) {
  using radray::cp::freq_init;
  constexpr double v[5] = {freq_init[0], freq_init[1], freq_init[2], freq_init[3], freq_init[4]};
  return v[n];
}
RR_HD double pop2_init(int n) {
  using radray::cp::pop2_init;
  constexpr double v[5] = {pop2_init[0], pop2_init[1], pop2_init[2], pop2_init[3], pop2_init[4]};
  return v[n];
}

}  // namespace rrk
