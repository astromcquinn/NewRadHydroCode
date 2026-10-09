// Serial/OpenMP reference of the spike hydro kernel (unsplit
// MUSCL-Hancock + HLLC, spike/common/hydro_mh.hpp) on a periodic n^3 grid.
// The framework prototypes must reproduce its output bit for bit.
//
// usage: hydro_serial [--n 64] [--steps 20] [--ic blast|sod] [--cfl 0.4]
//                     [--out U.bin]
// Output: raw float64, (5, n, n, n), component fastest, i fastest.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "hydro_mh.hpp"

int main(int argc, char** argv) {
  int n = 64, steps = 20, ic = 0;
  double cfl = 0.4;
  std::string out;
  for (int a = 1; a < argc; ++a) {
    const std::string s = argv[a];
    if (s == "--n") n = std::atoi(argv[++a]);
    else if (s == "--steps") steps = std::atoi(argv[++a]);
    else if (s == "--ic") ic = std::string(argv[++a]) == "sod" ? 1 : 0;
    else if (s == "--cfl") cfl = std::atof(argv[++a]);
    else if (s == "--out") out = argv[++a];
    else { std::fprintf(stderr, "unknown argument %s\n", s.c_str()); return 2; }
  }
  constexpr int G = rrh::NG, NV = rrh::NV;
  const int m = n + 2 * G;
  const long long nc = static_cast<long long>(m) * m * m;
  std::vector<double> U(nc * NV), Q(nc * rrh::NQ);
  std::vector<double> F[3];
  for (auto& f : F) f.assign(nc * NV, 0.0);
  auto at = [m](int i, int j, int k) {  // i,j,k in [-G, n+G)
    return ((i + G) + static_cast<long long>(m) * ((j + G) + static_cast<long long>(m) * (k + G)));
  };
  for (int k = 0; k < n; ++k)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) rrh::initial_state(ic, i, j, k, n, &U[at(i, j, k) * NV]);
  auto fill_ghosts = [&]() {
#pragma omp parallel for schedule(static)
    for (int k = -G; k < n + G; ++k)
      for (int j = -G; j < n + G; ++j)
        for (int i = -G; i < n + G; ++i) {
          if (i >= 0 && i < n && j >= 0 && j < n && k >= 0 && k < n) continue;
          const int ii = (i + n) % n, jj = (j + n) % n, kk = (k + n) % n;
          std::memcpy(&U[at(i, j, k) * NV], &U[at(ii, jj, kk) * NV], NV * sizeof(double));
        }
  };
  auto Ua = [&](int i, int j, int k, int c) { return U[at(i, j, k) * NV + c]; };
  const double dx = 1.0 / n;
  double t = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (int s = 0; s < steps; ++s) {
    fill_ghosts();
    double smax = 0;
#pragma omp parallel for schedule(static) reduction(max : smax)
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) smax = rrf::max(smax, rrh::max_speed(&U[at(i, j, k) * NV]));
    const double dt = cfl * dx / smax;
    const double dtdx = dt / dx;
#pragma omp parallel for schedule(static)
    for (int k = -1; k < n + 1; ++k)
      for (int j = -1; j < n + 1; ++j)
        for (int i = -1; i < n + 1; ++i) rrh::predict_cell(Ua, i, j, k, dtdx, &Q[at(i, j, k) * rrh::NQ]);
    for (int d = 0; d < 3; ++d) {
      const int ei = d == 0, ej = d == 1, ek = d == 2;
#pragma omp parallel for schedule(static)
      for (int k = 0; k < n + ek; ++k)
        for (int j = 0; j < n + ej; ++j)
          for (int i = 0; i < n + ei; ++i)
            rrh::face_flux(&Q[at(i - ei, j - ej, k - ek) * rrh::NQ], &Q[at(i, j, k) * rrh::NQ], d,
                           &F[d][at(i, j, k) * NV]);
    }
#pragma omp parallel for schedule(static)
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
          rrh::update_cell(&U[at(i, j, k) * NV], &F[0][at(i, j, k) * NV], &F[0][at(i + 1, j, k) * NV],
                           &F[1][at(i, j, k) * NV], &F[1][at(i, j + 1, k) * NV], &F[2][at(i, j, k) * NV],
                           &F[2][at(i, j, k + 1) * NV], dtdx);
    t += dt;
  }
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  int nth = 1;
#ifdef _OPENMP
  nth = omp_get_max_threads();
#endif
  double mass = 0, etot = 0;
  for (int k = 0; k < n; ++k)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        mass += U[at(i, j, k) * NV];
        etot += U[at(i, j, k) * NV + 4];
      }
  std::printf(" n %d steps %d t %.6e mass %.17g energy %.17g\n", n, steps, t, mass * dx * dx * dx,
              etot * dx * dx * dx);
  std::printf(" threads %d wall %.3f s  ns/cell-update %.2f (core-ns %.2f)\n", nth, wall,
              1e9 * wall / (double(n) * n * n * steps), 1e9 * wall * nth / (double(n) * n * n * steps));
  if (!out.empty()) {
    std::vector<double> o(static_cast<long long>(n) * n * n * NV);
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
          for (int c = 0; c < NV; ++c)
            o[((i + static_cast<long long>(n) * (j + static_cast<long long>(n) * k)) * NV) + c] =
                U[at(i, j, k) * NV + c];
    std::FILE* fp = std::fopen(out.c_str(), "wb");
    std::fwrite(o.data(), sizeof(double), o.size(), fp);
    std::fclose(fp);
  }
  return 0;
}
