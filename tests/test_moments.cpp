// Unit tests of the moment RT (rt_moments.hpp, docs/rt_moments.md):
// transport conservation and realizability, the M1 free-streaming beam, the
// OTVET tensor (direct sum vs FFT, single source), thread independence and
// the redshift; the alternative OTVET tensors (otvet_closure.hpp).
#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "check.hpp"
#include "radray/constants.hpp"
#include "radray/otvet_closure.hpp"
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

// otvet_closure.hpp: the pieces, the split's direct light, the screened
// kernel's budget, the cone's shadows, thread independence
OtvetClosureIn closure_input(int N) {
  OtvetClosureIn in;
  in.N = N;
  const SigmaTable sig = photoion_sigma(cp::freq_init);
  for (int g = 0; g < kNfreq; ++g) {
    in.q[g] = cp::pop2_init[g] / cp::pop2_init[0];
    in.sigma18[g] = sig[0][g] * 1e18;
  }
  return in;
}

void test_otvet_closures() {
  // HEALPix ring centres: 12 nside^2 unit vectors with sum n n^T -> npix/3 I
  // (pixel centres: <z^2> = 0.296 at nside 1)
  for (int ns : {1, 2, 4}) {
    const auto d = otvet_healpix_dirs(ns);
    CHECK(static_cast<int>(d.size()) == 12 * ns * ns);
    double S[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (const auto& n : d) {
      CHECK_CLOSE(n[0] * n[0] + n[1] * n[1] + n[2] * n[2], 1.0, 1e-12);
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) S[a][b] += n[a] * n[b];
    }
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        CHECK(std::fabs(S[a][b] / d.size() - (a == b ? 1.0 / 3 : 0.0)) < 0.04 / (ns * ns));
  }
  // smallest eigenvalue
  {
    const double a[6] = {0.5, 0.3, 0.2, 0, 0, 0};
    CHECK_CLOSE(otvet_min_eigenvalue(a), 0.2, 1e-12);
    const double b[6] = {0.5, 0.5, 0, 0.5, 0, 0};  // eigenvalues 1, 0, 0
    CHECK(std::fabs(otvet_min_eigenvalue(b)) < 1e-7);  // acos at a double root
    const double c[6] = {0.4, 0.4, 0.2, 0.6, 0, 0};  // 1.0, -0.2, 0.2
    CHECK_CLOSE(otvet_min_eigenvalue(c), -0.2, 1e-12);
  }
  // thin + split: the original h, and N_dir = sum L/(4 pi r^2)
  {
    const int N = 16;
    auto sp = make_sim(N, RTMethod::OTVET);
    Simulation& s = *sp;
    RTGrid& G = s.rt;
    G.sfr[G.cell(5, 6, 7)] = 1;
    G.Luv[G.cell(5, 6, 7)] = 3;
    G.sfr[G.cell(12, 2, 9)] = 1;
    G.Luv[G.cell(12, 2, 9)] = 1.5f;
    std::vector<float> h0;
    compute_otvet_tensor(G, s.units, 1, h0);
    OtvetClosureIn in = closure_input(N);
    in.src = {{4, 5, 6, 3.0}, {11, 1, 8, 1.5}};
    OtvetClosureParams p;
    p.split = 1;
    OtvetClosureOut out;
    otvet_closure_tensor(p, 1, in, out);
    CHECK(out.h == h0);
    const std::int64_t c = G.cell(9, 6, 7);  // 4 cells from the first source
    const double d2 = 3 * 3 + 4 * 4 + 2 * 2;  // from the second
    CHECK_CLOSE(out.ndir[c], (3.0 / 16 + 1.5 / d2) / (4 * M_PI), 1e-6);
  }
  if (!otvet_fft_available()) return;
  // screen: sum of the screened light = 4 pi lambda L (+ the source cell)
  {
    const int N = 32;
    OtvetClosureIn in = closure_input(N);
    in.src = {{3, 7, 11, 2.0}};
    OtvetClosureParams p;
    p.closure = 1;
    p.split = 1;
    p.screen_lambda = 3;
    OtvetClosureOut out;
    otvet_closure_tensor(p, 0, in, out);
    double t = 0;
    for (float v : out.ndir) t += v;
    // the lattice sum over the 27 images, the source cell and the continuum
    double lat = 4;
    for (int k = 0; k < N; ++k)
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i)
          for (int sc = -1; sc <= 1; ++sc)
            for (int sb = -1; sb <= 1; ++sb)
              for (int sa = -1; sa <= 1; ++sa) {
                auto m = [N](int d) { return 2 * d > N ? d - N : d; };
                const double x = m(i) + sa * N, y = m(j) + sb * N, z = m(k) + sc * N;
                const double r2 = x * x + y * y + z * z;
                if (r2 > 0) lat += std::exp(-std::sqrt(r2) / 3) / r2;
              }
    lat += 4 * M_PI * 3 * std::exp(-3.0 * N * std::cbrt(3 / (4 * M_PI)) / 3);
    const double expect = lat * 2.0 / (4 * M_PI);
    std::printf("  screen lambda 3: sum N_dir / lattice sum = %.6f (continuum 4 pi lambda: %.3f)\n",
                t / expect, (4 * M_PI * 3 + 4) / lat);
    CHECK(std::fabs(t / expect - 1) < 1e-4);
    // budget: lambda = sum N / (q_w sum L)
    std::vector<double> nw(static_cast<std::size_t>(N) * N * N, 1e-3);
    in.nw = nw.data();
    double qw = 0;
    for (int g = 0; g < kNfreq; ++g) qw += in.q[g] * in.sigma18[g];
    p.screen_lambda = 0;
    otvet_closure_tensor(p, 0, in, out);
    CHECK_CLOSE(out.lambda, 1e-3 * nw.size() / (qw * 2.0), 1e-9);
  }
  // cone: two sources 11 cells from a receiver along -x and -y; an HI
  // block in front of the first one shadows it
  {
    const int N = 32;
    OtvetClosureIn in = closure_input(N);
    in.src = {{4, 15, 15, 1.0}, {15, 4, 15, 1.0}};
    std::vector<float> col(static_cast<std::size_t>(N) * N * N, 0.0f);
    in.col = col.data();
    OtvetClosureParams p;
    p.closure = 2;
    p.cone_nside = 2;
    p.cone_coarse = 2;
    p.cone_rc = 4;
    p.cone_budget = 0;
    p.cone_blend = 0;
    OtvetClosureOut open, shadow;
    otvet_closure_tensor(p, 0, in, open);
    for (int k = 12; k <= 18; ++k)
      for (int j = 12; j <= 18; ++j)
        for (int i = 7; i <= 12; ++i) col[i + N * (j + static_cast<std::int64_t>(N) * k)] = 1e3f;
    otvet_closure_tensor(p, 0, in, shadow);
    const std::int64_t c = 15 + N * (15 + static_cast<std::int64_t>(N) * 15);
    const float* a = &open.h[6 * c];
    const float* b = &shadow.h[6 * c];
    std::printf("  cone: h_xx, h_yy open %.3f %.3f, shadowed %.3f %.3f\n", a[0], a[1], b[0], b[1]);
    CHECK(std::fabs(a[0] - a[1]) < 0.05);
    CHECK(b[1] > b[0] + 0.1);
    for (std::int64_t q = 0; q < static_cast<std::int64_t>(col.size()); ++q) {
      const float* h = &shadow.h[6 * q];
      for (int m = 0; m < 6; ++m) CHECK(std::isfinite(h[m]));
      CHECK_CLOSE(double(h[0]) + h[1] + h[2], 1.0, 1e-5);
      const double hd[6] = {h[0], h[1], h[2], h[3], h[4], h[5]};
      CHECK(otvet_min_eigenvalue(hd) > -1e-5);
    }
    // the budget brings the direct light to the photons held
    std::vector<double> nw(col.size(), 0.03);
    in.nw = nw.data();
    p.cone_budget = 1;
    p.cone_blend = 2;
    otvet_closure_tensor(p, 0, in, shadow);
    std::printf("  cone budget: lambda_x = %.2f, N_dir/N = %.4f\n", shadow.lam_x, shadow.ratio);
    CHECK(shadow.lam_x > 1 && std::fabs(shadow.ratio - 1) < 1e-6);
  }
  // screen and cone: identical for any thread count
  for (int closure : {1, 2}) {
    const int N = 32;
    std::vector<float> ref;
    for (int nt : {1, 3}) {
      par::set_threads(nt);
      OtvetClosureIn in = closure_input(N);
      std::mt19937 rng(5);
      std::uniform_int_distribution<int> U(0, N - 1);
      std::uniform_real_distribution<double> L(0.1, 10);
      for (int n = 0; n < 12; ++n) in.src.push_back({U(rng), U(rng), U(rng), L(rng)});
      std::vector<float> col(static_cast<std::size_t>(N) * N * N);
      std::vector<double> nw(col.size());
      std::uniform_real_distribution<double> V(0, 1);
      for (std::size_t q = 0; q < col.size(); ++q) {
        col[q] = V(rng) < 0.1 ? 3.0f : 0.0f;
        nw[q] = 1e-3 * V(rng);
      }
      in.col = col.data();
      in.nw = nw.data();
      OtvetClosureParams p;
      p.closure = closure;
      p.cone_nside = 2;
      p.cone_coarse = 2;
      p.cone_rc = 4;
      OtvetClosureOut out;
      otvet_closure_tensor(p, 0, in, out);
      if (ref.empty())
        ref = out.h;
      else
        CHECK(out.h == ref);
    }
    par::set_threads(1);
  }
}

// otvet_split = 1 in the solver: thread independence, realizability
void test_split_threads() {
  std::vector<float> ref;
  for (int nt : {1, 3}) {
    par::set_threads(nt);
    auto sp = make_sim(32, RTMethod::OTVET);
    Simulation& s = *sp;
    s.cfg.otvet_split = 1;
    MomentState& M = *s.moments;
    std::mt19937 rng(11);
    std::uniform_real_distribution<float> U(0, 1);
    for (std::size_t i = 0; i < M.n.size(); ++i) {
      M.n[i] = U(rng) < 0.05f ? U(rng) : 0.0f;
      M.e[i] = 20 * M.n[i];
    }
    s.rt.sfr[s.rt.cell(3, 4, 5)] = 1;
    s.rt.Luv[s.rt.cell(3, 4, 5)] = 1;
    moment_prep(s);
    CHECK(!M.ndir.empty());
    for (int it = 0; it < 3; ++it) moment_step_1(s);
    CHECK(realizable(M));
    std::vector<float> all = M.n;
    all.insert(all.end(), M.fx.begin(), M.fx.end());
    all.insert(all.end(), M.e.begin(), M.e.end());
    if (ref.empty())
      ref = all;
    else
      CHECK(all == ref);
  }
  par::set_threads(1);
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
  test_otvet_closures();
  test_split_threads();
  TEST_MAIN_END();
}
