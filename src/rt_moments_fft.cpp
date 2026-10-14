// FFT convolution of the OTVET tensor (rt_moments.hpp), FFTW3 single
// precision when fftw3.h is available (CMakeLists.txt then links FFTW3f);
// otherwise a stub, and the OTVET tensor is always summed directly.
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>

#include "radray/otvet_closure.hpp"
#include "radray/parallel.hpp"
#include "radray/rt_moments.hpp"

#if __has_include(<fftw3.h>)
#include <fftw3.h>

namespace radray {

bool otvet_fft_available() { return true; }

// The original OTVET kernel (otvet_kernel), component m
static void thin_fill(int N, int m, const void*, float* K) {
  auto mi = [N](int d) { return 2 * d > N ? d - N : d; };
#pragma omp parallel for schedule(static) if (par::worth(static_cast<long long>(N) * N * N))
  for (int k = 0; k < N; ++k)
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i) {
        double Kd[6];
        otvet_kernel(N, mi(i), mi(j), mi(k), Kd);
        K[i + static_cast<std::int64_t>(N) * (j + static_cast<std::int64_t>(N) * k)] =
            static_cast<float>(Kd[m]);
      }
}

void otvet_fft_convolve(int N, const std::vector<double>& rho,
                        std::array<std::vector<double>, 6>& P) {
  otvet_fft_convolve_kernel(N, rho, thin_fill, nullptr, P, "thin", 0.0);
}

namespace {

// Kernel spectra kept across calls (they depend only on N, the kernel and
// its parameter: the thin kernel and the cones' near field are the same
// every hydro step). A few entries, least recently used first out.
struct KernelSpectra {
  std::string tag;
  int N = 0;
  double param = 0;
  std::array<fftwf_complex*, 6> k{};
  ~KernelSpectra() {
    for (auto* a : k) fftwf_free(a);
  }
};
std::vector<std::unique_ptr<KernelSpectra>>& kernel_cache() {
  static std::vector<std::unique_ptr<KernelSpectra>> c;
  return c;
}
constexpr std::size_t kCacheEntries = 3;

template <class T>
T* fft_alloc(std::int64_t n) {
  T* p = static_cast<T*>(fftwf_malloc(sizeof(T) * static_cast<std::size_t>(n)));
  if (!p) throw std::runtime_error("otvet_fft_convolve: out of memory");
  return p;
}

}  // namespace

// rho -> its spectrum (one serial FFT); the six kernel spectra (cached, or
// filled and transformed); then per component, in parallel: product and
// inverse FFT. Each transform is a serial FFTW execution of the same plans
// on arrays of the same alignment (fftwf_malloc), so the result does not
// depend on the thread count, and equals the one-component-at-a-time loop
// it replaces bit for bit.
void otvet_fft_convolve_kernel(int N, const std::vector<double>& rho, OtvetKernelFill fill,
                               const void* ctx, std::array<std::vector<double>, 6>& P,
                               const char* cache_tag, double cache_param) {
  const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
  const std::int64_t nk = static_cast<std::int64_t>(N) * N * (N / 2 + 1);
  const bool big = par::worth(nc);
  float* r = fft_alloc<float>(nc);
  fftwf_complex* rk = fft_alloc<fftwf_complex>(nk);
  std::array<float*, 6> out{};
  std::array<fftwf_complex*, 6> scr{};
  for (int m = 0; m < 6; ++m) {
    out[m] = fft_alloc<float>(nc);
    scr[m] = fft_alloc<fftwf_complex>(nk);
  }
  // Row-major (k, j, i): the RT grid's linear index i + N (j + N k)
  fftwf_plan pf = fftwf_plan_dft_r2c_3d(N, N, N, r, rk, FFTW_ESTIMATE);
  fftwf_plan pb = fftwf_plan_dft_c2r_3d(N, N, N, scr[0], out[0], FFTW_ESTIMATE);

#pragma omp parallel for schedule(static) if (big)
  for (std::int64_t c = 0; c < nc; ++c) r[c] = static_cast<float>(rho[c]);
  fftwf_execute_dft_r2c(pf, r, rk);

  // kernel spectra
  auto& cache = kernel_cache();
  KernelSpectra* ks = nullptr;
  if (cache_tag) {
    for (std::size_t e = 0; e < cache.size(); ++e)
      if (cache[e]->tag == cache_tag && cache[e]->N == N && cache[e]->param == cache_param) {
        std::rotate(cache.begin() + e, cache.begin() + e + 1, cache.end());  // most recent last
        ks = cache.back().get();
        break;
      }
  }
  std::unique_ptr<KernelSpectra> fresh;
  if (!ks) {
    fresh = std::make_unique<KernelSpectra>();
    fresh->tag = cache_tag ? cache_tag : "";
    fresh->N = N;
    fresh->param = cache_param;
    for (int m = 0; m < 6; ++m) {
      fresh->k[m] = fft_alloc<fftwf_complex>(nk);
      fill(N, m, ctx, out[m]);  // parallel over cells inside
    }
#pragma omp parallel for schedule(static, 1) if (big)
    for (int m = 0; m < 6; ++m) fftwf_execute_dft_r2c(pf, out[m], fresh->k[m]);
    ks = fresh.get();
    if (cache_tag) {
      if (cache.size() >= kCacheEntries) cache.erase(cache.begin());
      cache.push_back(std::move(fresh));
    }
  }

  // product and inverse transform, one component per thread
#pragma omp parallel for schedule(static, 1) if (big)
  for (int m = 0; m < 6; ++m) {
    const fftwf_complex* kk = ks->k[m];
    fftwf_complex* sc = scr[m];
    for (std::int64_t q = 0; q < nk; ++q) {
      const float ar = rk[q][0], ai = rk[q][1], br = kk[q][0], bi = kk[q][1];
      sc[q][0] = ar * br - ai * bi;
      sc[q][1] = ar * bi + ai * br;
    }
    fftwf_execute_dft_c2r(pb, sc, out[m]);
  }
  for (int m = 0; m < 6; ++m) {
    P[m].resize(nc);
    const float* o = out[m];
    double* pm = P[m].data();
#pragma omp parallel for schedule(static) if (big)
    for (std::int64_t c = 0; c < nc; ++c) pm[c] = static_cast<double>(o[c]) / nc;
  }
  fftwf_destroy_plan(pf);
  fftwf_destroy_plan(pb);
  fftwf_free(r);
  fftwf_free(rk);
  for (int m = 0; m < 6; ++m) {
    fftwf_free(out[m]);
    fftwf_free(scr[m]);
  }
}

}  // namespace radray

#else

namespace radray {

bool otvet_fft_available() { return false; }

void otvet_fft_convolve(int, const std::vector<double>&, std::array<std::vector<double>, 6>&) {
  throw std::runtime_error("otvet_fft_convolve: built without FFTW3f");
}

void otvet_fft_convolve_kernel(int, const std::vector<double>&, OtvetKernelFill, const void*,
                               std::array<std::vector<double>, 6>&, const char*, double) {
  throw std::runtime_error("otvet_fft_convolve: built without FFTW3f");
}

}  // namespace radray

#endif
