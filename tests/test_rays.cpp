// Unit tests: ray_segment, nside/radius relations, chain operations,
// ray_cast_split / ray_cast_merge conservation.
#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "check.hpp"
#include "radray/constants.hpp"
#include "radray/fortran.hpp"
#include "radray/healpix.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/rt.hpp"

using namespace radray;

static void test_segment() {
  const double seg_min = static_cast<double>(1e-4f);
  // Axis-aligned from cell centre: 0.5 to the face, then 0.5 in the next cell
  {
    const double x[3] = {2.5, 2.5, 2.5}, dx[3] = {1, 0, 0};
    const auto s = ray_segment(x, dx, seg_min);
    CHECK(s[0] == 0.5);
    CHECK(s[1] == 0.5);
    CHECK(s[2] == 0.0);
  }
  // Ray on a face moving in -x: full unit step inside cell [1,2)
  {
    const double x[3] = {2.0, 2.5, 2.5}, dx[3] = {-1, 0, 0};
    const auto s = ray_segment(x, dx, seg_min);
    CHECK(s[0] == 1.0);
  }
  // Random directions: segments non-negative, total <= 1, each segment ends
  // on a cell face (or at total length 1)
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> U(0, 8);
  std::normal_distribution<double> G;
  int n_trunc = 0;
  for (int n = 0; n < 200000; ++n) {
    double x[3] = {U(rng), U(rng), U(rng)}, dx[3] = {G(rng), G(rng), G(rng)};
    const double nr = std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    for (double& d : dx) d /= nr;
    const auto s = ray_segment(x, dx, seg_min);
    const double tot = s[0] + s[1] + s[2];
    CHECK(s[0] >= 0 && s[1] >= 0 && s[2] >= 0);
    CHECK(tot <= 1 + 1e-12);
    if (tot < 1 - 1e-9) ++n_trunc;  // 4th crossing dropped (original behaviour)
  }
  CHECK(n_trunc > 0);  // the 3-segment limit does truncate occasionally
}

static void test_nside_radius() {
  using namespace rayp;
  CHECK(rmin_from_nside(nside_min) == 0);
  CHECK(rmax_from_nside(4) == 4 * std::sqrt(3 / pc::pi));
  CHECK(rmin_from_nside(8) == 4 * std::sqrt(3 / pc::pi));
  CHECK(nray_from_nside(2) == nray_init);
  CHECK(nray_from_nside(16) == 12 * 256);
  CHECK(lside_from_radius(0.5) == lside_min);
  CHECK(lside_from_radius(1e9) == lside_max);
  // Consistency: rays per solid angle ~ ray_per_cell * 4 pi r^2
  for (int l = lside_min; l < lside_max; ++l) {
    const int ns = 1 << l;
    const double rmax = rmax_from_nside(ns);
    // just below rmax the radius maps to nside, just above to 2*nside
    CHECK(nside_from_radius(rmax * (1 - 1e-6)) == ns);
    CHECK(nside_from_radius(rmax * (1 + 1e-6)) == 2 * ns);
    // between rmin and rmax of nside (nside > nside_min) -> nside
    if (ns > nside_min) CHECK(nside_from_radius(0.5 * (rmin_from_nside(ns) + rmax)) == ns);
  }
  // monotone
  int prev = 0;
  for (double r = 0; r < 5000; r += 0.37) {
    const int l = lside_from_radius(r);
    CHECK(l >= prev);
    prev = l;
  }
}

static std::unique_ptr<Simulation> make_sim(int N) {
  auto s = std::make_unique<Simulation>();
  init_simulation(*s, N, 64, 10.0);
  s->clk.dtsec_rt = 1e12;
  // neutral uniform gas fields are irrelevant here (tau = 0 everywhere)
  s->rays.Nrbuf = s->rays.Nfree_max / 4;
  ray_buffer_allocate(*s);
  return s;
}

static double chain_count(Simulation& s, int i, int j, int k) {
  std::int64_t ir = s.rt.htoc_ray(3, s.rt.cell(i, j, k));
  int n = 0;
  std::int64_t prev = 0;
  while (ir > 0) {
    CHECK(s.rays.at(ir).prev == prev);
    prev = ir;
    ir = s.rays.at(ir).next;
    ++n;
  }
  CHECK(s.rt.htoc_ray(4, s.rt.cell(i, j, k)) == prev);
  return n;
}

static void test_chain_ops() {
  auto sp = make_sim(8);
  Simulation& s = *sp;
  std::int64_t ix[5];
  ray_create(s, 5, ix);
  for (int n = 0; n < 5; ++n) s.rays.at(ix[n]).nside = 4;
  ray_insert_tail(s, ix[0], 1, 1, 1);
  ray_insert_tail(s, ix[1], 1, 1, 1);
  ray_insert_head(s, ix[2], 1, 1, 1);
  ray_insert_tail(s, ix[3], 1, 1, 1);
  // chain: 2,0,1,3
  const std::int64_t c = s.rt.cell(1, 1, 1);
  std::vector<std::int64_t> order;
  for (std::int64_t ir = s.rt.htoc_ray(3, c); ir > 0; ir = s.rays.at(ir).next) order.push_back(ir);
  CHECK(order.size() == 4);
  CHECK(order[0] == ix[2] && order[1] == ix[0] && order[2] == ix[1] && order[3] == ix[3]);
  CHECK(chain_count(s, 1, 1, 1) == 4);
  const std::int64_t nfree0 = s.rays.stat.nfree;
  ray_remove(s, ix[2], 1, 1, 1, true);   // head
  ray_remove(s, ix[3], 1, 1, 1, false);  // tail, not freed
  CHECK(chain_count(s, 1, 1, 1) == 2);
  CHECK(s.rays.stat.nfree == nfree0 + 1);
  CHECK(s.rays.stat.ifree[s.rays.stat.nfree] == ix[2]);  // LIFO
  CHECK(s.rays.stat.ndel1 == 1);
  // freed ray keeps its data (zombie) until reused / deallocated
  CHECK(s.rays.at(ix[2]).nside == 4);
  std::int64_t again;
  ray_create(s, 1, &again);
  CHECK(again == ix[2]);
}

static double sum_grid_ngam(const Simulation& s) {
  double t = 0;
  for (float v : s.rt.ngam) t += v;
  return t;
}

static void test_split_conserves() {
  auto sp = make_sim(16);
  Simulation& s = *sp;
  std::mt19937_64 rng(3);
  std::normal_distribution<double> G;
  for (int trial = 0; trial < 200; ++trial) {
    std::int64_t ir;
    ray_create(s, 1, &ir);
    Ray& r = s.rays.at(ir);
    const int nside = 4 << (trial % 4);
    double v[3] = {G(rng), G(rng), G(rng)};
    const int ipix = healpix::vec2pix_nest(nside, v);
    const double nr = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    r.nside = nside;
    r.ipix = ipix;
    r.r = static_cast<float>(rmax_from_nside(nside) + 0.3);
    r.dr = -0.3f;
    for (int d = 0; d < 3; ++d) {
      r.x[d] = static_cast<float>(8.5 + 0.2 * d);
      r.dx[d] = static_cast<float>(v[d] / nr);
    }
    double ng0[kNfreq];
    for (int n = 0; n < kNfreq; ++n) {
      r.ngam[n] = 0.001f * (n + 1) * (1 + trial);
      r.nu[n] = static_cast<float>(cp::freq_init[n]);
      ng0[n] = r.ngam[n];
    }
    ray_insert_tail(s, ir, 9, 9, 9);
    std::fill(s.rt.ngam.begin(), s.rt.ngam.end(), 0.0f);
    const std::int64_t nfree_before = s.rays.stat.nfree;
    ray_cast_split(s, ir, 9, 9, 9);
    CHECK(s.rays.stat.nfree == nfree_before - 3);
    // collect the 4 daughters: parent slot + 3 most recently handed out
    std::vector<std::int64_t> d = {ir};
    for (int m = 1; m <= 3; ++m) d.push_back(s.rays.stat.ifree[nfree_before - 3 + m]);
    double tot[kNfreq] = {0, 0, 0, 0, 0};
    for (auto k : d) {
      const Ray& q = s.rays.at(k);
      CHECK(q.nside == std::min(2 * nside, rayp::nside_max));
      const double len = std::sqrt(double(q.dx[0]) * q.dx[0] + double(q.dx[1]) * q.dx[1] +
                                   double(q.dx[2]) * q.dx[2]);
      CHECK(std::fabs(len - 1) < 1e-6);
      for (int n = 0; n < kNfreq; ++n) tot[n] += q.ngam[n];
    }
    for (int n = 0; n < kNfreq; ++n) CHECK_CLOSE(tot[n], ng0[n], 1e-6);
    // tau = 0 -> each daughter deposits photons in proportion to path length
    CHECK(sum_grid_ngam(s) > 0);
    // daughters are linked somewhere: remove them for the next trial
    for (auto k : d) {
      const Ray& q = s.rays.at(k);
      const int i = 1 + f::mod(f::int4(q.x[0]), 16), j = 1 + f::mod(f::int4(q.x[1]), 16),
                kk = 1 + f::mod(f::int4(q.x[2]), 16);
      ray_remove(s, k, i, j, kk, true);
    }
  }
  // all chains empty again
  for (std::int64_t c = 0; c < s.rt.ncell; ++c) {
    CHECK(s.rt.htoc_ray(3, c) == 0);
    CHECK(s.rt.htoc_ray(4, c) == 0);
  }
}

static void test_merge_conserves() {
  auto sp = make_sim(16);
  Simulation& s = *sp;
  std::mt19937_64 rng(11);
  std::normal_distribution<double> G;
  std::uniform_real_distribution<double> U(0, 1);
  for (int trial = 0; trial < 200; ++trial) {
    const int nr = 2 + trial % 7;
    std::vector<std::int64_t> irc(nr);
    double ng[kNfreq] = {0, 0, 0, 0, 0}, en[kNfreq] = {0, 0, 0, 0, 0};
    // rays in one cell, roughly parallel (same cone)
    double base[3] = {G(rng), G(rng), G(rng)};
    for (int m = 0; m < nr; ++m) {
      ray_create(s, 1, &irc[m]);
      Ray& r = s.rays.at(irc[m]);
      double v[3];
      for (int d = 0; d < 3; ++d) v[d] = base[d] + 0.05 * G(rng);
      const double l = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      r.nside = 16;
      r.r = static_cast<float>(10 + U(rng));
      r.dr = static_cast<float>(rmax_from_nside(16) - r.r);
      for (int d = 0; d < 3; ++d) {
        r.x[d] = static_cast<float>(4 + U(rng));
        r.dx[d] = static_cast<float>(v[d] / l);
      }
      double vv[3] = {r.dx[0], r.dx[1], r.dx[2]};
      r.ipix = healpix::vec2pix_nest(16, vv);
      for (int n = 0; n < kNfreq; ++n) {
        r.ngam[n] = static_cast<float>(1e-3 * (1 + U(rng)));
        r.nu[n] = static_cast<float>(cp::freq_init[n] * (1 + 0.1 * U(rng)));
        ng[n] += r.ngam[n];
        en[n] += static_cast<double>(r.ngam[n] * r.nu[n]);
      }
      ray_insert_tail(s, irc[m], 5, 5, 5);
    }
    const std::int64_t ndel0 = s.rays.stat.ndel1;
    ray_cast_merge(s, nr, irc.data(), 5, 5, 5);
    CHECK(s.rays.stat.ndel1 == ndel0 + nr - 1);
    const Ray& q = s.rays.at(irc[0]);
    for (int n = 0; n < kNfreq; ++n) {
      CHECK_CLOSE(q.ngam[n], ng[n], 2e-7);
      CHECK_CLOSE(static_cast<double>(q.ngam[n]) * q.nu[n], en[n], 5e-7);
    }
    const int i = 1 + f::mod(f::int4(q.x[0]), 16), j = 1 + f::mod(f::int4(q.x[1]), 16),
              k = 1 + f::mod(f::int4(q.x[2]), 16);
    CHECK(chain_count(s, i, j, k) >= 1);
    ray_remove(s, irc[0], i, j, k, true);
    CHECK(chain_count(s, 5, 5, 5) == 0);
  }
}

int main() {
  test_segment();
  test_nside_radius();
  test_chain_ops();
  test_split_conserves();
  test_merge_conserves();
  TEST_MAIN_END();
}
