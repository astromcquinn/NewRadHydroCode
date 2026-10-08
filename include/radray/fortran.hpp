// Helpers reproducing Fortran (gfortran, x86-64) semantics exactly.
//
// Precision rules used throughout the port (see docs/port_notes.md):
//   real(4) -> float, real(8) -> double.
//   Mixed float/double arithmetic promotes to double (same as C++).
//   Single-precision Fortran literals (e.g. 0.99, 1E-6) are written as
//   float literals (0.99f, 1e-6f) and promoted where Fortran promotes them.
//   Integer powers x**2, x**3, x**4 are expanded like gfortran does
//   (x*x, (x*x)*x, (x*x)*(x*x)).
//   sum() over a real(4) array accumulates in float, left to right from 0.
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

namespace radray {
namespace f {

// Fortran INT(x) for a real -> default integer(4): truncation toward zero.
// gfortran emits cvttsd2si/cvttss2si, which returns INT32_MIN (the
// "integer indefinite" value) for NaN and out-of-range inputs. We reproduce
// that instead of invoking C++ undefined behaviour.
inline std::int32_t int4(double x) {
  if (!(x > -2147483649.0 && x < 2147483648.0)) return INT32_MIN;
  return static_cast<std::int32_t>(x);
}

// Fortran INT(x, kind=8).
inline std::int64_t int8(double x) {
  if (!(x >= -9223372036854775808.0 && x < 9223372036854775808.0))
    return INT64_MIN;
  return static_cast<std::int64_t>(x);
}

// Fortran NINT(x) -> integer(4): round half away from zero (gfortran lround).
inline std::int32_t nint(double x) { return int4(std::round(x)); }

// Fortran FLOOR / CEILING -> default integer(4).
inline std::int32_t floor4(double x) { return int4(std::floor(x)); }
inline std::int32_t ceiling4(double x) { return int4(std::ceil(x)); }

// Fortran MOD for integers: result has the sign of the first argument,
// which is exactly C++'s % (truncating division).
inline std::int32_t mod(std::int32_t a, std::int32_t p) { return a % p; }
inline std::int64_t mod(std::int64_t a, std::int64_t p) { return a % p; }

// Fortran MOD for reals: a - int(a/p)*p computed exactly; gfortran uses
// fmod, which is exact.
inline double mod(double a, double p) { return std::fmod(a, p); }

// Fortran MAX/MIN (two arguments, same kind after promotion).
// NaN behaviour is processor dependent in gfortran >= 9; inputs here are
// never NaN in normal operation.
template <class T>
inline T max(T a, T b) { return (b > a) ? b : a; }
template <class T>
inline T min(T a, T b) { return (b < a) ? b : a; }

// Integer powers as expanded by gfortran (powi addition chains).
template <class T>
inline T pow2(T x) { return x * x; }
template <class T>
inline T pow3(T x) { return (x * x) * x; }
template <class T>
inline T pow4(T x) { T t = x * x; return t * t; }

// sum() of a real(4) array: float accumulator, sequential.
template <int N>
inline float sum4(const float* a) {
  float s = 0.0f;
  for (int i = 0; i < N; ++i) s += a[i];
  return s;
}
inline float sum4(const float* a, int n) {
  float s = 0.0f;
  for (int i = 0; i < n; ++i) s += a[i];
  return s;
}

// sum() of a real(8) array: double accumulator, sequential.
template <int N>
inline double sum8(const double* a) {
  double s = 0.0;
  for (int i = 0; i < N; ++i) s += a[i];
  return s;
}
inline double sum8(const double* a, int n) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += a[i];
  return s;
}

// huge(0.) and huge(0) of default kinds.
constexpr float huge_r4 = std::numeric_limits<float>::max();
constexpr std::int32_t huge_i4 = std::numeric_limits<std::int32_t>::max();
// epsilon(real(0.))
constexpr float epsilon_r4 = std::numeric_limits<float>::epsilon();

}  // namespace f
}  // namespace radray
