// Threaded FFTW (fft_3d with nthreads > 1, used by pm_refine > 1,
// docs/pm_refine.md) must not depend on the thread count: the forward and
// backward transforms of a random field are byte-identical for every count.
//
// usage: test_fft_threads [N] [threads ...]   (default 64, 1 2 3 4 8)
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#include "check.hpp"
#include "radray/fft.hpp"

int main(int argc, char** argv) {
  using namespace radray;
  const int N = argc > 1 ? std::atoi(argv[1]) : 64;
  std::vector<int> threads;
  for (int a = 2; a < argc; ++a) threads.push_back(std::atoi(argv[a]));
  if (threads.empty()) threads = {1, 2, 3, 4, 8};
  std::printf("fft_3d %d^3, threaded FFTW: %s\n", N, fft::fft_threaded() ? "yes" : "no");
  const std::size_t n = static_cast<std::size_t>(N + 2) * N * N;
  std::vector<float> a0(n, 0.0f);
  {
    std::mt19937 rng(1);
    std::uniform_real_distribution<float> U(-1, 1);
    for (std::size_t k = 0; k < static_cast<std::size_t>(N) * N; ++k)
      for (int i = 0; i < N; ++i) a0[k * (N + 2) + i] = U(rng);
  }
  std::vector<float> f1, b1, a(n);
  for (int t : threads) {
    std::memcpy(a.data(), a0.data(), n * sizeof(float));
    fft::fft_3d('f', N, N, N, a.data(), t);
    const std::vector<float> f = a;
    fft::fft_3d('b', N, N, N, a.data(), t);
    if (f1.empty()) {
      f1 = f;
      b1 = a;
    }
    const bool sf = std::memcmp(f.data(), f1.data(), n * sizeof(float)) == 0;
    const bool sb = std::memcmp(a.data(), b1.data(), n * sizeof(float)) == 0;
    std::printf("  %3d threads: forward %s, backward %s\n", t, sf ? "identical" : "DIFFERENT",
                sb ? "identical" : "DIFFERENT");
    CHECK(sf && sb);
  }
  // nthreads = 1 after threaded plans is still the serial plan
  std::memcpy(a.data(), a0.data(), n * sizeof(float));
  fft::fft_3d('f', N, N, N, a.data());
  CHECK(std::memcmp(a.data(), f1.data(), n * sizeof(float)) == 0);
  TEST_MAIN_END();
}
