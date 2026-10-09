// Host/device annotation and the Fortran-semantics helpers the kernels in
// radray/kernels/ need, in a form that compiles for CUDA/HIP as well as the
// host. They reproduce include/radray/fortran.hpp exactly on the host (same
// operations), so kernels built from them stay bit-identical to the port.
// Header-only and framework-independent on purpose: the same file is
// compiled by the AMReX port (it was spike/common/rr_hd.hpp on the
// framework-spike branch; this copy is the one maintained from now on).
#pragma once

#include <cmath>
#include <cstdint>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define RR_HD __host__ __device__ inline
#else
#define RR_HD inline
#endif

namespace rrf {

RR_HD std::int32_t int4(double x) {
  if (!(x > -2147483649.0 && x < 2147483648.0)) return INT32_MIN;
  return static_cast<std::int32_t>(x);
}
RR_HD std::int32_t floor4(double x) { return int4(std::floor(x)); }
RR_HD std::int32_t ceiling4(double x) { return int4(std::ceil(x)); }
RR_HD double mod(double a, double p) { return std::fmod(a, p); }
RR_HD std::int32_t mod(std::int32_t a, std::int32_t p) { return a % p; }
template <class T>
RR_HD T max(T a, T b) { return (b > a) ? b : a; }
template <class T>
RR_HD T min(T a, T b) { return (b < a) ? b : a; }
template <class T>
RR_HD T pow2(T x) { return x * x; }
template <class T>
RR_HD T pow4(T x) { T t = x * x; return t * t; }
template <int N>
RR_HD float sum4(const float* a) {
  float s = 0.0f;
  for (int i = 0; i < N; ++i) s += a[i];
  return s;
}
template <int N>
RR_HD double sum8(const double* a) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) s += a[i];
  return s;
}

}  // namespace rrf
