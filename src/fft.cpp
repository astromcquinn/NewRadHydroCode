// fft_tools.f90 with FFTW3f in place of MKL DFTI. Every call below has an
// exact counterpart in reference/harness_pm/src/fft_tools.f90.
#include "radray/fft.hpp"

#include <fftw3.h>

#include <cstdint>
#include <stdexcept>

namespace radray {
namespace fft {

namespace {
// FFTW_ESTIMATE never touches the arrays and is deterministic;
// FFTW_UNALIGNED makes the plan independent of the array address, so the
// Fortran harness and this port (different allocators) get the same plan.
constexpr unsigned kFlags = FFTW_ESTIMATE | FFTW_UNALIGNED;
}  // namespace

void fft_3d(char c, int n1, int n2, int n3, float* a) {
  fftwf_complex* ac = reinterpret_cast<fftwf_complex*>(a);
  switch (c) {
    case 'f': {
      // Row-major dimensions (n3,n2,n1) = Fortran a(n1,n2,n3)
      fftwf_plan p = fftwf_plan_dft_r2c_3d(n3, n2, n1, a, ac, kFlags);
      if (!p) throw std::runtime_error("fftwf_plan_dft_r2c_3d failed");
      fftwf_execute_dft_r2c(p, a, ac);
      fftwf_destroy_plan(p);
      break;
    }
    case 'b': {
      // bscale = 1E0/(int(n1,8)*int(n2,8)*int(n3,8)): real(4)/integer(8),
      // the integer is converted to real(4) first.
      const std::int64_t ntot = static_cast<std::int64_t>(n1) * n2 * n3;
      const float bscale = 1.0f / static_cast<float>(ntot);
      fftwf_plan p = fftwf_plan_dft_c2r_3d(n3, n2, n1, ac, a, kFlags);
      if (!p) throw std::runtime_error("fftwf_plan_dft_c2r_3d failed");
      fftwf_execute_dft_c2r(p, ac, a);
      fftwf_destroy_plan(p);
      // DFTI_BACKWARD_SCALE on the n1 real outputs of every row
      const std::int64_t ld = n1 + 2;
      for (std::int64_t k = 0; k < n3; ++k)
        for (std::int64_t j = 0; j < n2; ++j) {
          float* row = a + ld * (j + static_cast<std::int64_t>(n2) * k);
          for (int i = 0; i < n1; ++i) row[i] = row[i] * bscale;
        }
      break;
    }
    default:
      break;  // the Fortran select case ignores other values
  }
}

void fft_1d(char c, int n, float* a) {
  fftwf_complex* ac = reinterpret_cast<fftwf_complex*>(a);
  switch (c) {
    case 'f': {
      fftwf_plan p = fftwf_plan_dft_r2c_1d(n, a, ac, kFlags);
      if (!p) throw std::runtime_error("fftwf_plan_dft_r2c_1d failed");
      fftwf_execute_dft_r2c(p, a, ac);
      fftwf_destroy_plan(p);
      break;
    }
    case 'b': {
      const float bscale = 1.0f / static_cast<float>(n);  // 1E0/n
      fftwf_plan p = fftwf_plan_dft_c2r_1d(n, ac, a, kFlags);
      if (!p) throw std::runtime_error("fftwf_plan_dft_c2r_1d failed");
      fftwf_execute_dft_c2r(p, ac, a);
      fftwf_destroy_plan(p);
      for (int i = 0; i < n; ++i) a[i] = a[i] * bscale;
      break;
    }
    default:
      break;
  }
}

}  // namespace fft
}  // namespace radray
