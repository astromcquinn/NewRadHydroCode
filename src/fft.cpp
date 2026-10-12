// MKL DFTI semantics on FFTW3f (see fft.hpp). Every call below has an exact
// counterpart in reference/harness_common/mkl_dfti_shim.f90.
#include "radray/fft.hpp"

#include <fftw3.h>

#include <cstdint>
#include <cstring>
#include <new>
#include <stdexcept>

namespace radray {
namespace fft {

namespace {
// FFTW_ESTIMATE never touches the arrays and is deterministic;
// FFTW_UNALIGNED makes the plan independent of the array address, so the
// Fortran harness and this port (different allocators) get the same plan.
constexpr unsigned kFlags = FFTW_ESTIMATE | FFTW_UNALIGNED;

// Plans made while this is alive use n threads (threaded FFTW only); the
// planner is reset to one thread afterwards, so every other plan (the
// default PM, the hydro 1-D transforms) stays the serial one. FFTW's
// threads are initialised only on the first n > 1.
struct PlanThreads {
  bool on = false;
  explicit PlanThreads(int n) {
#ifdef RADRAY_FFTW_OMP
    if (n > 1) {
      static const bool init = fftwf_init_threads() != 0;
      if (!init) throw std::runtime_error("fftwf_init_threads failed");
      fftwf_plan_with_nthreads(n);
      on = true;
    }
#else
    (void)n;
#endif
  }
  ~PlanThreads() {
#ifdef RADRAY_FFTW_OMP
    if (on) fftwf_plan_with_nthreads(1);
#endif
  }
};
}  // namespace

bool fft_threaded() {
#ifdef RADRAY_FFTW_OMP
  return true;
#else
  return false;
#endif
}

void fft_3d(char c, int n1, int n2, int n3, float* a, int nthreads) {
  fftwf_complex* ac = reinterpret_cast<fftwf_complex*>(a);
  const PlanThreads pt(nthreads);
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
#pragma omp parallel for schedule(static) if (pt.on)
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

Fft1d::Fft1d(int n, float bscale) : n_(n), bscale_(bscale) {
  buf_ = fftwf_alloc_real(static_cast<std::size_t>(n + 2));
  if (!buf_) throw std::bad_alloc();
  auto* c = reinterpret_cast<fftwf_complex*>(buf_);
  fftwf_plan pf = fftwf_plan_dft_r2c_1d(n, buf_, c, FFTW_ESTIMATE);
  fftwf_plan pb = fftwf_plan_dft_c2r_1d(n, c, buf_, FFTW_ESTIMATE);
  if (!pf || !pb) throw std::runtime_error("FFTW planning failed");
  pf_ = pf;
  pb_ = pb;
}

Fft1d::~Fft1d() {
  fftwf_destroy_plan(static_cast<fftwf_plan>(pf_));
  fftwf_destroy_plan(static_cast<fftwf_plan>(pb_));
  fftwf_free(buf_);
}

void Fft1d::forward(float* x) {
  std::memcpy(buf_, x, sizeof(float) * n_);
  fftwf_execute_dft_r2c(static_cast<fftwf_plan>(pf_), buf_,
                        reinterpret_cast<fftwf_complex*>(buf_));
  std::memcpy(x, buf_, sizeof(float) * (n_ + 2));
}

void Fft1d::backward(float* x) {
  std::memcpy(buf_, x, sizeof(float) * (n_ + 2));
  fftwf_execute_dft_c2r(static_cast<fftwf_plan>(pb_), reinterpret_cast<fftwf_complex*>(buf_),
                        buf_);
  for (int i = 0; i < n_; ++i) x[i] = buf_[i] * bscale_;
}

void fft_1d(char c, int n, float* a) {
  switch (c) {
    case 'f': {
      Fft1d d(n, 1.0f);
      d.forward(a);
      break;
    }
    case 'b': {
      Fft1d d(n, 1.0f / static_cast<float>(n));  // bscale = 1E0/n
      d.backward(a);
      break;
    }
    default:
      break;
  }
}

}  // namespace fft
}  // namespace radray
