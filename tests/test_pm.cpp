// Unit tests of the PM port: FFT layout/scaling (MKL DFTI semantics on top
// of FFTW3f), domain decomposition, linked lists, CIC mass conservation and
// the potential of a single Fourier mode.
#include <cmath>
#include <complex>
#include <cstdint>
#include <random>
#include <vector>

#include "check.hpp"
#include "radray/constants.hpp"
#include "radray/fft.hpp"
#include "radray/hydro.hpp"
#include "radray/pm.hpp"

using namespace radray;

// Naive DFT of a(n1,n2,n3) (Fortran order) at frequency (m1,m2,m3), sign -1
static std::complex<double> dft(const std::vector<float>& a, int n1, int n2, int n3, int m1,
                                int m2, int m3) {
  std::complex<double> s = 0;
  const int ld = n1 + 2;
  for (int k = 0; k < n3; ++k)
    for (int j = 0; j < n2; ++j)
      for (int i = 0; i < n1; ++i) {
        const double ph = -2 * pc::pi * (double(m1) * i / n1 + double(m2) * j / n2 + double(m3) * k / n3);
        s += double(a[i + ld * (j + n2 * k)]) * std::complex<double>(std::cos(ph), std::sin(ph));
      }
  return s;
}

static void test_fft_3d() {
  const int n1 = 8, n2 = 6, n3 = 4, ld = n1 + 2;
  std::mt19937 rng(7);
  std::uniform_real_distribution<float> U(-1, 1);
  std::vector<float> a(ld * n2 * n3, 0.0f);
  for (int k = 0; k < n3; ++k)
    for (int j = 0; j < n2; ++j)
      for (int i = 0; i < n1; ++i) a[i + ld * (j + n2 * k)] = U(rng);
  const std::vector<float> a0 = a;

  fft::fft_3d('f', n1, n2, n3, a.data());
  // conjugate-even storage: complex c(n1/2+1,n2,n3), strides 1, n1/2+1, ...
  double maxerr = 0;
  for (int m3 = 0; m3 < n3; ++m3)
    for (int m2 = 0; m2 < n2; ++m2)
      for (int m1 = 0; m1 <= n1 / 2; ++m1) {
        const std::complex<double> ref = dft(a0, n1, n2, n3, m1, m2, m3);
        const std::size_t c = m1 + (n1 / 2 + 1) * (m2 + std::size_t(n2) * m3);
        const std::complex<double> got(a[2 * c], a[2 * c + 1]);
        maxerr = std::fmax(maxerr, std::abs(got - ref));
      }
  CHECK(maxerr < 1e-4);

  // backward with DFTI_BACKWARD_SCALE = 1/(n1 n2 n3): round trip
  fft::fft_3d('b', n1, n2, n3, a.data());
  double rt = 0;
  for (int k = 0; k < n3; ++k)
    for (int j = 0; j < n2; ++j)
      for (int i = 0; i < n1; ++i) {
        const int p = i + ld * (j + n2 * k);
        rt = std::fmax(rt, std::fabs(a[p] - a0[p]));
      }
  CHECK(rt < 1e-6);
}

static void test_fft_1d() {
  const int n = 12;
  std::vector<float> a(n + 2, 0.0f);
  for (int i = 0; i < n; ++i) a[i] = std::sin(0.3f * i) + 0.1f * i;
  const std::vector<float> a0 = a;
  fft::fft_1d('f', n, a.data());
  // CCS: a(2m+1), a(2m+2) = Re, Im of X_m, m = 0..n/2
  for (int m = 0; m <= n / 2; ++m) {
    std::complex<double> s = 0;
    for (int i = 0; i < n; ++i) s += double(a0[i]) * std::polar(1.0, -2 * pc::pi * m * i / n);
    CHECK(std::abs(std::complex<double>(a[2 * m], a[2 * m + 1]) - s) < 1e-4);
  }
  fft::fft_1d('b', n, a.data());
  for (int i = 0; i < n; ++i) CHECK(std::fabs(a[i] - a0[i]) < 1e-6);
}

static void test_domains() {
  Domains d;
  init_domain_gas(Config::make(32), d);
  CHECK(d.Ndomain == 8);
  CHECK(d.Nx_dom == 16);
  // domain 6 (1-based): i = 2, j = 1, k = 2
  CHECK(d.box[5][0] == 17 && d.box[5][1] == 32);
  CHECK(d.box[5][2] == 1 && d.box[5][3] == 16);
  CHECK(d.box[5][4] == 17 && d.box[5][5] == 32);
}

static void test_link_and_cic() {
  pm::PMState s;
  {
    Domains d;
    init_domain_gas(Config::make(8), d);
    pm::init_pm(s, pm::PMConfig::make(8, 4), d);
  }
  std::mt19937 rng(3);
  std::uniform_real_distribution<float> U(0, 4);
  for (std::int64_t ip = 1; ip <= s.pc.Npart; ++ip)
    for (int c = 0; c < 3; ++c) s.xv_p(ip)[c] = U(rng);
  s.xv_p(5)[0] = 4.0f;  // a coordinate rounded up to exactly Ndm: maps to cell 1
  pm::link_dm(s);
  // every particle exactly once, in its gas cell, chains in increasing ip
  std::vector<int> seen(s.pc.Npart + 1, 0);
  for (int k = 1; k <= 8; ++k)
    for (int j = 1; j <= 8; ++j)
      for (int i = 1; i <= 8; ++i) {
        std::int64_t prev = 0;
        for (std::int64_t ip = s.hoc[s.gcell(i, j, k)]; ip > 0; ip = s.ll[ip]) {
          ++seen[ip];
          CHECK(ip > prev);
          prev = ip;
          const float* p = s.xv_p(ip);
          CHECK(1 + (int(p[0] * 2.0) % 8) == i);
          CHECK(1 + (int(p[1] * 2.0) % 8) == j);
          CHECK(1 + (int(p[2] * 2.0) % 8) == k);
        }
      }
  for (std::int64_t ip = 1; ip <= s.pc.Npart; ++ip) CHECK(seen[ip] == 1);

  // CIC conserves mass: sum = Npart*(1-fb)*8 on top of the gas
  std::vector<float> u(5 * 512, 0.0f);
  for (int c = 0; c < 512; ++c) u[5 * c] = 0.25f;
  pm::calc_density_field(s, u.data());
  double sum = 0;
  for (int k = 1; k <= 8; ++k)
    for (int j = 1; j <= 8; ++j)
      for (int i = 1; i <= 8; ++i) sum += s.density[s.dcell(i, j, k)];
  const double expect = 512 * 0.25 + s.pc.Npart * (1 - cp::omegab / cp::omegam) * 8;
  CHECK_CLOSE(sum, expect, 1e-6);
}

static void test_potential_single_mode() {
  // gas only, rho = 1 + eps cos(2 pi x/N): phi = Aphi*rho_k/K^2 with the
  // finite-difference K^2 = (2 sin(k/2))^2, Aphi = -2a/3
  const int N = 16;
  pm::PMState s;
  {
    Domains d;
    init_domain_gas(Config::make(N), d);
    pm::init_pm(s, pm::PMConfig::make(N, 8), d);
  }
  std::vector<float> u(5 * N * N * N, 0.0f);
  const double eps = 0.1, kk = 2 * pc::pi / N, a = 0.5;
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j)
      for (int i = 1; i <= N; ++i)
        u[5 * s.gcell(i, j, k)] = float(1 + eps * std::cos(kk * (i - 1)));
  pm::calc_density_field(s, u.data());  // no particles linked: gas only
  pm::calc_gravity_field(s, a);
  const double K2 = std::pow(2 * std::sin(kk / 2), 2);
  for (int i = 1; i <= N; ++i) {
    const double expect = -(2 * a / 3) * eps * std::cos(kk * (i - 1)) / K2;
    CHECK(std::fabs(s.density[s.dcell(i, 3, 5)] - expect) < 1e-5);
  }
}

int main() {
  test_fft_3d();
  test_fft_1d();
  test_domains();
  test_link_and_cic();
  test_potential_single_mode();
  TEST_MAIN_END();
}
