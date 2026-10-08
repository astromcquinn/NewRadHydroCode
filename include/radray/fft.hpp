// The one FFT layer of the port: what the original gets from Intel MKL DFTI,
// in hydro_tools.f90 (1-D transforms in smooth) and in fft_tools.f90
// (fft_3d for the PM, fft_1d unused).
//
// MKL is not available, so this port and the Fortran reference harnesses
// both call FFTW3 single precision with identical plans, which makes them
// agree bit for bit. On the Fortran side every harness uses one shim,
// reference/harness_common/mkl_dfti_shim.f90 (module MKL_DFTI), under the
// unmodified hydro_tools.f90 and fft_tools.f90; each routine below has its
// exact counterpart there. The MKL semantics the code relies on:
//
//   * 1-D real descriptor of length n (Fft1d; fft_1d): in place on a(n+2),
//     MKL's default CCS packed format = FFTW's r2c layout (n/2+1 interleaved
//     complex numbers). Forward unnormalised; backward = c2r times the
//     real(4) DFTI_BACKWARD_SCALE. Plans: FFTW_ESTIMATE, made once per
//     descriptor on an FFTW-allocated (aligned) n+2 buffer; data are copied
//     in and out of that buffer.
//   * 3-D real descriptor (fft_3d): array a(n1+2,n2,n3) (Fortran order, i
//     fastest). Forward: real input in a(1:n1,:,:), conjugate-even output in
//     the same memory as complex c(n1/2+1,n2,n3) (DFTI_COMPLEX_COMPLEX with
//     complex strides (1, n1/2+1, (n1/2+1)*n2)). Unnormalised, exp(-i k x).
//     Backward: the inverse, scaled by bscale = 1E0/(n1*n2*n3) (real(4)),
//     applied to the n1 real outputs of each row. Plans: FFTW_ESTIMATE |
//     FFTW_UNALIGNED, in place on the array itself, one per call.
//
// Results differ from an MKL build at round-off level.
#pragma once

namespace radray {
namespace fft {

// MKL DFTI 1-D single-precision real descriptor (DftiCreateDescriptor +
// DFTI_BACKWARD_SCALE + DftiCommitDescriptor ... DftiFreeDescriptor).
class Fft1d {
 public:
  Fft1d(int n, float bscale);
  ~Fft1d();
  Fft1d(const Fft1d&) = delete;
  Fft1d& operator=(const Fft1d&) = delete;

  // DftiComputeForward: x(1:n) -> x(1:n+2)
  void forward(float* x);
  // DftiComputeBackward: x(1:n+2) -> x(1:n) = c2r(x) * bscale
  void backward(float* x);

 private:
  int n_;
  float bscale_;
  float* buf_ = nullptr;
  void* pf_ = nullptr;  // fftwf_plan
  void* pb_ = nullptr;
};

// fft_tools.f90: fft_3d(c,n1,n2,n3,a). c = 'f' (forward) or 'b' (backward).
void fft_3d(char c, int n1, int n2, int n3, float* a);

// fft_tools.f90: fft_1d(c,n,a). Unused by the original (the hydro calls are
// commented out); one Fft1d descriptor per call, as the original creates
// one MKL descriptor per call.
void fft_1d(char c, int n, float* a);

}  // namespace fft
}  // namespace radray
