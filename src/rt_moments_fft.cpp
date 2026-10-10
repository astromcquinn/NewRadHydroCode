// FFT convolution of the OTVET tensor (rt_moments.hpp), FFTW3 single
// precision when fftw3.h is available (CMakeLists.txt then links FFTW3f);
// otherwise a stub, and the OTVET tensor is always summed directly.
#include <stdexcept>

#include "radray/rt_moments.hpp"

#if __has_include(<fftw3.h>)
#include <fftw3.h>

namespace radray {

bool otvet_fft_available() { return true; }

void otvet_fft_convolve(int N, const std::vector<double>& rho,
                        std::array<std::vector<double>, 6>& P) {
  const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
  const std::int64_t nk = static_cast<std::int64_t>(N) * N * (N / 2 + 1);
  float* r = static_cast<float*>(fftwf_malloc(sizeof(float) * nc));
  fftwf_complex* rk = static_cast<fftwf_complex*>(fftwf_malloc(sizeof(fftwf_complex) * nk));
  fftwf_complex* kk = static_cast<fftwf_complex*>(fftwf_malloc(sizeof(fftwf_complex) * nk));
  if (!r || !rk || !kk) throw std::runtime_error("otvet_fft_convolve: out of memory");
  // Row-major (k, j, i): the RT grid's linear index i + N (j + N k)
  fftwf_plan pf = fftwf_plan_dft_r2c_3d(N, N, N, r, rk, FFTW_ESTIMATE);
  fftwf_plan pb = fftwf_plan_dft_c2r_3d(N, N, N, kk, r, FFTW_ESTIMATE);

  for (std::int64_t c = 0; c < nc; ++c) r[c] = static_cast<float>(rho[c]);
  fftwf_execute_dft_r2c(pf, r, rk);

  auto mi = [N](int d) { return 2 * d > N ? d - N : d; };
  for (int m = 0; m < 6; ++m) {
    for (int k = 0; k < N; ++k)
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
          double K[6];
          otvet_kernel(N, mi(i), mi(j), mi(k), K);
          r[i + static_cast<std::int64_t>(N) * (j + static_cast<std::int64_t>(N) * k)] =
              static_cast<float>(K[m]);
        }
    fftwf_execute_dft_r2c(pf, r, kk);
    for (std::int64_t q = 0; q < nk; ++q) {
      const float ar = rk[q][0], ai = rk[q][1], br = kk[q][0], bi = kk[q][1];
      kk[q][0] = ar * br - ai * bi;
      kk[q][1] = ar * bi + ai * br;
    }
    fftwf_execute_dft_c2r(pb, kk, r);
    P[m].resize(nc);
    for (std::int64_t c = 0; c < nc; ++c) P[m][c] = static_cast<double>(r[c]) / nc;
  }
  fftwf_destroy_plan(pf);
  fftwf_destroy_plan(pb);
  fftwf_free(r);
  fftwf_free(rk);
  fftwf_free(kk);
}

}  // namespace radray

#else

namespace radray {

bool otvet_fft_available() { return false; }

void otvet_fft_convolve(int, const std::vector<double>&, std::array<std::vector<double>, 6>&) {
  throw std::runtime_error("otvet_fft_convolve: built without FFTW3f");
}

}  // namespace radray

#endif
