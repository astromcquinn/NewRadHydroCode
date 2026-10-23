// The guard of the moment level's threaded OTVET FFT plans
// (docs/stage_e_moments.md 4.1): OtvetConvolver with plans made for n
// threads gives bitwise the convolution of the serial plans
// (convolve_checked), for the thin kernel on a sparse source field, and
// equals otvet_fft_convolve (the serial code) component by component.
//
// usage: test_otvet_fft_threads N [threads ...]
#include <fftw3.h>

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "check.hpp"
#include "radray/otvet_closure.hpp"
#include "radray/rt_moments.hpp"

int main(int argc, char** argv) {
  using namespace radray;
  const int N = argc > 1 ? std::atoi(argv[1]) : 64;
  std::vector<int> threads;
  for (int a = 2; a < argc; ++a) threads.push_back(std::atoi(argv[a]));
  if (threads.empty()) threads = {2, 4};
  const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
  // a sparse source field (1 % of the cells), as rho = L / Lmax
  std::vector<double> rho(nc, 0.0);
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> U(0, 1);
  for (std::int64_t c = 0; c < nc; ++c)
    if (U(rng) < 0.01) rho[c] = U(rng);
  // the serial code
  std::array<std::vector<double>, 6> P;
  otvet_fft_convolve(N, rho, P);
  if (fftwf_init_threads() == 0) {
    std::printf("fftwf_init_threads failed\n");
    return 1;
  }
  for (int nt : threads) {
    fftwf_plan_with_nthreads(nt);
    OtvetConvolver cv(N);
    fftwf_plan_with_nthreads(1);
    for (std::int64_t c = 0; c < nc; ++c) cv.r()[c] = static_cast<float>(rho[c]);
    cv.forward();
    for (int m = 0; m < 6; ++m) {
      const void* ks = cv.kernel(m, otvet_thin_fill, nullptr, nullptr, 0.0, cv.out());
      if (m == 0) {
        CHECK(cv.convolve_checked(ks, cv.out()));
      } else {
        cv.convolve(ks, cv.out());
      }
      std::int64_t ndiff = 0;
      for (std::int64_t c = 0; c < nc; ++c)
        if (static_cast<double>(cv.out()[c]) / nc != P[m][c]) ++ndiff;
      std::printf("N = %d, %d threads, component %d: %lld cells differ from the serial code\n", N, nt, m,
                  static_cast<long long>(ndiff));
      CHECK(ndiff == 0);
    }
  }
  TEST_MAIN_END();
}
