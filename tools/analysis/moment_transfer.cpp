// Constant-opacity transfer function of the moment RT (docs/rt_moments.md
// section 7): the steady-state photon density N of the M1 and OTVET
// transport for sources L in a medium of uniform mean free path lambda, to
// compare T(k) = (P_N / P_L)^(1/2) with the exact arctan(k lambda)/(k lambda)
// (Gnedin 2025, arXiv:2504.10571, eq. 2) and with the linear theory of the
// closures.
//
// The transport is the solver's own: mk::flux_cell and mk::glf_cell of
// kernels/moment_cell.hpp in the loop of transport_group (rt_moments.cpp),
// and the OTVET tensor is otvet_fft_convolve with the normalisation of
// compute_otvet_tensor. Only the physics is replaced: every RT substep
// (moment_nsub subcycles, one cell of light travel) ends with N, F *=
// exp(-1/lambda) instead of photo_cell's fgam, and the sources inject after
// the transport (moment_inject = 1) or before it (0).
//
//   moment_transfer N lambda nsrc m1|otvet|iso nsteps out_prefix [seed] [inject] [nsub] [order]
//
// order = 2: moment_order = 2 (moment_transport2 of rt_moments.cpp; nsub >= 6).
//
// lambda in cells; nsrc = 0: every cell a source, L = 1 + 0.3 g (g a unit
// Gaussian, L >= 0.05), the linear regime; nsrc > 0: that many random cells
// with L = 1. Writes out_prefix.L and out_prefix.N (float32, N^3, x
// fastest): N as radiative_physics would see it. iso = P = N I/3 (the
// linearised M1).
//
// Build (after the library): g++ -O2 -fopenmp -ffp-contract=off -Iinclude
//   tools/analysis/moment_transfer.cpp build/libradray.a -lfftw3f -o moment_transfer
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "radray/kernels/moment_cell.hpp"
#include "radray/rt_moments.hpp"

using namespace radray;

int main(int argc, char** argv) {
  if (argc < 7) {
    std::fprintf(stderr, "usage: %s N lambda nsrc m1|otvet|iso nsteps out [seed] [inject] [nsub] [order]\n", argv[0]);
    return 1;
  }
  const int N = std::atoi(argv[1]);
  const double lambda = std::atof(argv[2]);
  const long nsrc = std::atol(argv[3]);
  const std::string method = argv[4];
  const int nsteps = std::atoi(argv[5]);
  const std::string out = argv[6];
  const unsigned seed = argc > 7 ? std::atoi(argv[7]) : 1;
  const int inject_after = argc > 8 ? std::atoi(argv[8]) : 1;
  const int nsub = argc > 9 ? std::atoi(argv[9]) : 4;
  const int order = argc > 10 ? std::atoi(argv[10]) : 1;
  if (order == 2 && nsub < 6) {
    std::fprintf(stderr, "order 2 needs nsub >= 6\n");
    return 2;
  }
  std::vector<float> work2;
  const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
  const mk::Closure cl =
      method == "m1" ? mk::Closure::M1 : method == "otvet" ? mk::Closure::Otvet : mk::Closure::Iso;

  // sources
  std::vector<float> L(nc, 0.0f);
  std::mt19937_64 rng(seed);
  if (nsrc == 0) {
    std::normal_distribution<double> g(0.0, 1.0);
    for (std::int64_t c = 0; c < nc; ++c) L[c] = static_cast<float>(std::max(0.05, 1.0 + 0.3 * g(rng)));
  } else {
    std::uniform_int_distribution<std::int64_t> u(0, nc - 1);
    for (long s = 0; s < nsrc; ++s) L[u(rng)] += 1.0f;
  }

  // OTVET tensor: h = P / tr P, P = L convolved with the kernel
  std::vector<float> h;
  if (cl == mk::Closure::Otvet) {
    std::vector<double> rho(L.begin(), L.end());
    std::array<std::vector<double>, 6> P;
    otvet_fft_convolve(N, rho, P);
    h.assign(6 * nc, 0.0f);
    for (std::int64_t c = 0; c < nc; ++c) {
      const double tr = (P[0][c] + P[1][c]) + P[2][c];
      for (int m = 0; m < 6; ++m)
        h[6 * c + m] = tr > 0 ? static_cast<float>(P[m][c] / tr) : (m < 3 ? 1.0f / 3 : 0.0f);
    }
  }

  std::vector<float> n(nc, 0.0f), fx(nc, 0.0f), fy(nc, 0.0f), fz(nc, 0.0f), e(nc, 0.0f);
  std::vector<float> nn(nc), nfx(nc), nfy(nc), nfz(nc), ne(nc), P(9 * nc);
  const double C = 1.0 / nsub;
  const float surv = static_cast<float>(std::exp(-1.0 / lambda));
  const std::int64_t sj = N, sk = static_cast<std::int64_t>(N) * N;

  auto inject = [&]() {
#pragma omp parallel for schedule(static)
    for (std::int64_t c = 0; c < nc; ++c) {
      n[c] += L[c];
      e[c] += L[c];
    }
  };

  for (int step = 0; step < nsteps; ++step) {
    if (!inject_after) inject();
    for (int it = 0; it < nsub; ++it) {
      if (order == 2) {
        float* U[5] = {n.data(), fx.data(), fy.data(), fz.data(), e.data()};
        moment_transport2(N, U, static_cast<int>(cl), h.empty() ? nullptr : h.data(), nullptr, 0.0, C, work2);
        continue;
      }
#pragma omp parallel for schedule(static)
      for (std::int64_t c = 0; c < nc; ++c)
        mk::flux_cell(cl, n[c], fx[c], fy[c], fz[c], e[c], h.empty() ? nullptr : &h[6 * c], 0.0, &P[9 * c]);
#pragma omp parallel for schedule(static)
      for (int k = 0; k < N; ++k) {
        const std::int64_t kl = (k == 0 ? N - 1 : k - 1) * sk, kr = (k == N - 1 ? 0 : k + 1) * sk;
        for (int j = 0; j < N; ++j) {
          const std::int64_t jl = (j == 0 ? N - 1 : j - 1) * sj, jr = (j == N - 1 ? 0 : j + 1) * sj;
          const std::int64_t row = j * sj + k * sk;
          for (int i = 0; i < N; ++i) {
            const std::int64_t c = row + i;
            const std::int64_t xl = row + (i == 0 ? N - 1 : i - 1), xr = row + (i == N - 1 ? 0 : i + 1);
            const std::int64_t yl = c - j * sj + jl, yr = c - j * sj + jr;
            const std::int64_t zl = c - k * sk + kl, zr = c - k * sk + kr;
            const std::int64_t st[7] = {c, xl, xr, yl, yr, zl, zr};
            const float* v[5] = {n.data(), fx.data(), fy.data(), fz.data(), e.data()};
            float o[5];
            mk::glf_cell([&](int q, int p) { return v[q][st[p]]; },
                         [&](int p, int m) { return P[9 * st[p] + m]; }, C, o);
            nn[c] = o[0];
            nfx[c] = o[1];
            nfy[c] = o[2];
            nfz[c] = o[3];
            ne[c] = o[4];
          }
        }
      }
      n.swap(nn);
      fx.swap(nfx);
      fy.swap(nfy);
      fz.swap(nfz);
      e.swap(ne);
    }
    if (inject_after) inject();
    if (step == nsteps - 1) break;  // n is now what radiative_physics sees
#pragma omp parallel for schedule(static)
    for (std::int64_t c = 0; c < nc; ++c) {
      n[c] *= surv;
      e[c] *= surv;
      fx[c] *= surv;
      fy[c] *= surv;
      fz[c] *= surv;
    }
    if (step % 50 == 0 || step == nsteps - 2) {
      double t = 0;
      for (std::int64_t c = 0; c < nc; ++c) t += n[c];
      std::fprintf(stderr, "step %d photons %.6g\n", step, t);
    }
  }

  auto dump = [&](const std::string& f, const std::vector<float>& v) {
    FILE* fp = std::fopen(f.c_str(), "wb");
    std::fwrite(v.data(), sizeof(float), v.size(), fp);
    std::fclose(fp);
  };
  dump(out + ".L", L);
  dump(out + ".N", n);
  return 0;
}
