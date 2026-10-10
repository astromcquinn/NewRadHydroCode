// Unit tests of the moment RT (rt_moments.hpp, docs/rt_moments.md):
// transport conservation and realizability, the M1 free-streaming beam, the
// OTVET tensor (direct sum vs FFT, single source), thread independence and
// the redshift.
#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "check.hpp"
#include "radray/constants.hpp"
#include "radray/parallel.hpp"
#include "radray/rt.hpp"
#include "radray/rt_moments.hpp"

using namespace radray;

namespace {

std::unique_ptr<Simulation> make_sim(int N, RTMethod m) {
  auto s = std::make_unique<Simulation>();
  init_simulation(*s, N, 4, 10.0);
  s->cfg.rt_method = m;
  s->clk.dtsec_rt = 1e12;
  moment_init(*s);
  return s;
}

double total(const std::vector<float>& v) {
  double t = 0;
  for (float x : v) t += x;
  return t;
}

bool realizable(const MomentState& M) {
  for (std::size_t i = 0; i < M.n.size(); ++i) {
    if (!(M.n[i] >= 0) || !(M.e[i] >= 0)) return false;
    const double f = std::sqrt(double(M.fx[i]) * M.fx[i] + double(M.fy[i]) * M.fy[i] +
                               double(M.fz[i]) * M.fz[i]);
    if (f > 1.000001 * M.n[i]) return false;
  }
  return true;
}

// Point "source" put directly into N (no sources in the grid): transport
// only, no physics. N, E conserved; realizable; the photons stay within the
// light cone (|x| <= substeps + numerical diffusion)
void test_vacuum(RTMethod m) {
  const int N = 16;
  auto sp = make_sim(N, m);
  Simulation& s = *sp;
  MomentState& M = *s.moments;
  const std::int64_t c0 = s.rt.cell(8, 8, 8);
  for (int g = 0; g < kNfreq; ++g) {
    M.n[g * M.ncell + c0] = 1.0f;
    M.e[g * M.ncell + c0] = static_cast<float>(cp::freq_init[g]);
  }
  moment_prep(s);
  const double n0 = total(M.n), e0 = total(M.e);
  for (int it = 0; it < 6; ++it) {
    moment_step_1(s);
    CHECK_CLOSE(total(M.n), n0, 2e-6);
    CHECK_CLOSE(total(M.e), e0, 2e-6);
    CHECK(realizable(M));
  }
  // the deposit for radiative_physics is N
  CHECK(s.rt.ngam[c0 * kNfreq] == M.n[c0]);
  // mean radius of the photons: free streaming at <= one cell per substep
  double r1 = 0;
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j)
      for (int i = 1; i <= N; ++i) {
        const double r = std::sqrt(double((i - 8) * (i - 8) + (j - 8) * (j - 8) + (k - 8) * (k - 8)));
        r1 += r * M.n[s.rt.cell(i, j, k)];
      }
  r1 /= total(std::vector<float>(M.n.begin(), M.n.begin() + M.ncell));
  std::printf("  %s point pulse: mean radius after 6 substeps %.3f cells\n",
              rt_method_name(m), r1);
  // a shell at c t = 6 broadened by the numerical diffusion (sigma ~ 2
  // cells after 24 subcycles of C = 1/4) has <|r|> ~ c t + sigma^2/r
  CHECK(r1 <= 7.0);
  CHECK(r1 > 5.0);
}

// M1: a plane-parallel beam (F = N along +x) moves one cell per substep
// without spreading sideways (the free-streaming limit is exact in M1)
void test_m1_beam() {
  const int N = 32;
  auto sp = make_sim(N, RTMethod::M1);
  Simulation& s = *sp;
  MomentState& M = *s.moments;
  for (int k = 1; k <= N; ++k)
    for (int j = 1; j <= N; ++j) {
      const std::int64_t c = s.rt.cell(3, j, k);
      M.n[c] = 1.0f;
      M.fx[c] = 1.0f;
      M.e[c] = static_cast<float>(cp::freq_init[0]);
    }
  moment_prep(s);
  const int nsub = 5;
  for (int it = 0; it < nsub; ++it) moment_step_1(s);
  double sn = 0, sx = 0;
  bool flat = true;
  for (int i = 1; i <= N; ++i) {
    const float v = M.n[s.rt.cell(i, 1, 1)];
    for (int k = 1; k <= N; ++k)
      for (int j = 1; j <= N; ++j) {
        const std::int64_t c = s.rt.cell(i, j, k);
        flat = flat && M.n[c] == v && M.fy[c] == 0 && M.fz[c] == 0;
        sn += M.n[c];
        sx += M.n[c] * i;
      }
  }
  CHECK(flat);
  CHECK_CLOSE(sn, N * N, 1e-6);
  CHECK_CLOSE(sx / sn, 3.0 + nsub, 1e-5);
  // nothing upstream: pure upwind for a free-streaming beam
  CHECK(M.n[s.rt.cell(2, 1, 1)] == 0);
}

void test_otvet_tensor() {
  const int N = 16;
  auto sp = make_sim(N, RTMethod::OTVET);
  Simulation& s = *sp;
  RTGrid& G = s.rt;
  // one source: h = d d^T / |d|^2 (rank 1), isotropic in the source cell
  G.sfr[G.cell(5, 6, 7)] = 1;
  G.Luv[G.cell(5, 6, 7)] = 3;
  std::vector<float> h;
  CHECK(compute_otvet_tensor(G, s.units, 1, h) == 1);
  const float* a = &h[6 * G.cell(9, 6, 7)];
  CHECK_CLOSE(a[0], 1.0, 1e-7);
  CHECK(a[1] == 0 && a[2] == 0 && a[3] == 0);
  const float* b = &h[6 * G.cell(6, 7, 8)];  // d = (1,1,1)
  for (int m = 0; m < 6; ++m) CHECK_CLOSE(b[m], 1.0 / 3, 1e-6);
  const float* c = &h[6 * G.cell(5, 6, 7)];
  CHECK_CLOSE(c[0], 1.0 / 3, 1e-7);
  CHECK(c[3] == 0);
  // the periodic image at exactly N/2: no off-diagonal term
  const float* d = &h[6 * G.cell(5 + N / 2, 7, 7)];
  CHECK(d[3] == 0);
  // random sources: trace 1; direct sum and FFT agree
  std::mt19937 rng(3);
  std::uniform_int_distribution<int> U(1, N);
  std::uniform_real_distribution<double> L(0.1, 10);
  for (int n = 0; n < 20; ++n) {
    const std::int64_t q = G.cell(U(rng), U(rng), U(rng));
    G.sfr[q] = 1;
    G.Luv[q] = static_cast<float>(L(rng));
  }
  std::vector<float> h1, h2;
  compute_otvet_tensor(G, s.units, 1, h1);
  double emax = 0;
  for (std::int64_t q = 0; q < G.ncell; ++q)
    CHECK_CLOSE(double(h1[6 * q]) + h1[6 * q + 1] + h1[6 * q + 2], 1.0, 1e-6);
  if (otvet_fft_available()) {
    CHECK(compute_otvet_tensor(G, s.units, 2, h2) == 2);
    for (std::size_t i = 0; i < h1.size(); ++i) emax = std::fmax(emax, std::fabs(h1[i] - h2[i]));
    std::printf("  OTVET tensor: max |direct - FFT| = %.2e\n", emax);
    CHECK(emax < 1e-4);
  }
}

// The transport is a Jacobi stencil: identical for any thread count
void test_threads(RTMethod m) {
  std::vector<float> ref;
  for (int nt : {1, 3}) {
    par::set_threads(nt);
    auto sp = make_sim(32, m);  // >= 32^3: the parallel loops are on
    Simulation& s = *sp;
    MomentState& M = *s.moments;
    std::mt19937 rng(11);
    std::uniform_real_distribution<float> U(0, 1);
    for (std::size_t i = 0; i < M.n.size(); ++i) {
      M.n[i] = U(rng) < 0.05f ? U(rng) : 0.0f;
      M.fx[i] = 0.5f * M.n[i] * (U(rng) - 0.5f);
      M.e[i] = 20 * M.n[i];
    }
    s.rt.sfr[s.rt.cell(3, 4, 5)] = 1;
    s.rt.Luv[s.rt.cell(3, 4, 5)] = 1;
    moment_prep(s);
    for (int it = 0; it < 3; ++it) moment_step_1(s);
    std::vector<float> all = M.n;
    all.insert(all.end(), M.fx.begin(), M.fx.end());
    all.insert(all.end(), M.fz.begin(), M.fz.end());
    all.insert(all.end(), M.e.begin(), M.e.end());
    if (ref.empty())
      ref = all;
    else
      CHECK(all == ref);
  }
  par::set_threads(1);
}

// Redshift: the photons move down in frequency (those of the lowest group
// out of the band, counted as lost), F follows the photons, realizability
// is kept
void test_redshift() {
  auto sp = make_sim(8, RTMethod::M1);
  Simulation& s = *sp;
  MomentState& M = *s.moments;
  for (int g = 0; g < kNfreq; ++g) {
    M.n[g * M.ncell] = 1.0f;
    M.fx[g * M.ncell] = 0.5f;
    M.e[g * M.ncell] = static_cast<float>(cp::freq_init[g]);
  }
  s.clk.scalefactor = 0.1;
  s.clk.da_sim = 0.01;
  const double n0 = total(M.n);
  moment_redshift(s);
  const double n1 = total(M.n);
  CHECK(n1 < n0 && n1 > n0 - 1);
  CHECK_CLOSE(M.lost * s.cfg.nrtcr, n0 - n1, 1e-6);
  double fx = 0;
  for (float v : M.fx) fx += v;
  CHECK_CLOSE(fx, 0.5 * n1, 1e-6);
  CHECK(realizable(M));
  // group 1 keeps its frequency range
  const double nu1 = M.e[0] / M.n[0];
  CHECK(nu1 >= pc::HI_eV && nu1 <= std::sqrt(pc::HI_eV * pc::HeI_eV));
}

}  // namespace

int main() {
  test_vacuum(RTMethod::M1);
  test_vacuum(RTMethod::OTVET);
  test_m1_beam();
  test_otvet_tensor();
  test_threads(RTMethod::M1);
  test_threads(RTMethod::OTVET);
  test_redshift();
  TEST_MAIN_END();
}
