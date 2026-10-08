// fft_tools.f90: real <-> complex FFTs on padded in-place arrays.
//
// The original uses Intel MKL DFTI. MKL is not available, so this port and
// the Fortran reference harness (reference/harness_pm/src/fft_tools.f90)
// both call FFTW3 single precision with identical plans, which makes the two
// agree bit for bit. The MKL semantics that the code relies on are kept:
//
//   * fft_3d: array a(n1+2,n2,n3) (Fortran order, i fastest). Forward:
//     real input in a(1:n1,:,:), conjugate-even output in the same memory as
//     complex c(n1/2+1,n2,n3) (DFTI_COMPLEX_COMPLEX storage with complex
//     strides (1, n1/2+1, (n1/2+1)*n2)). Unnormalised, exp(-i k x).
//     Backward: the inverse, scaled by bscale = 1E0/(n1*n2*n3) (real(4)),
//     applied to the n1 real outputs of each row.
//   * fft_1d: a(n+2), the same layout in 1D (MKL CCS format), backward scale
//     1E0/n.
//
// See reference/harness_pm/PATCHES.md ("FFT") for the full mapping. Results
// differ from an MKL build at round-off level.
#pragma once

namespace radray {
namespace fft {

// fft_tools.f90: fft_3d(c,n1,n2,n3,a). c = 'f' (forward) or 'b' (backward).
void fft_3d(char c, int n1, int n2, int n3, float* a);

// fft_tools.f90: fft_1d(c,n,a). Unused by the original (the hydro calls are
// commented out) but ported for completeness.
void fft_1d(char c, int n, float* a);

}  // namespace fft
}  // namespace radray
