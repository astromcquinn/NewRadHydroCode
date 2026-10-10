// Unit test of the AMReX backend's ray encodings (src/amrex/ray_encoding.hpp)
// on a fake particle tile: the float layout round-trips exactly, the
// R-struct within its quantisation steps (docs/stage_d.md section 3), with
// the photon total exact and never-merged directions on the pixel centre.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include "check.hpp"
#include "radray/kernels/healpix_hd.hpp"
#include "../src/amrex/ray_encoding.hpp"

using namespace radray;
using rrk::RayP;

template <int NR, int NI>
struct FakeTile {
  std::vector<float> r[NR];
  std::vector<int> i[NI];
  struct Data {
    float* m_rdata[NR];
    int* m_idata[NI];
  };
  explicit FakeTile(int n) {
    for (auto& v : r) v.assign(n, 0.0f);
    for (auto& v : i) v.assign(n, 0);
  }
  Data data() {
    Data d;
    for (int k = 0; k < NR; ++k) d.m_rdata[k] = r[k].data();
    for (int k = 0; k < NI; ++k) d.m_idata[k] = i[k].data();
    return d;
  }
};

static RayP random_ray(std::mt19937_64& rng, int N, bool merged) {
  std::uniform_real_distribution<double> U(0, 1);
  std::normal_distribution<double> G;
  RayP r;
  const int lside = 2 + static_cast<int>(U(rng) * 8);  // nside 4..512
  r.nside = 1 << lside;
  double v[3] = {G(rng), G(rng), G(rng)};
  r.ipix = rrhp::vec2pix_nest(r.nside, v);
  double c[3];
  rrhp::pix2vec_nest(r.nside, r.ipix, c);
  if (merged) {
    // a deviation inside the pixel: up to ~0.6 pixel radius
    const double ps = std::sqrt(4 * pc::pi / (12.0 * r.nside * r.nside));
    double w[3] = {c[0] + 0.6 * ps * G(rng) / 2, c[1] + 0.6 * ps * G(rng) / 2, c[2] + 0.6 * ps * G(rng) / 2};
    const double n = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    for (int d = 0; d < 3; ++d) r.dx[d] = static_cast<float>(w[d] / n);
    r.ipix = rrhp::vec2pix_nest(r.nside, w);
    rrhp::pix2vec_nest(r.nside, r.ipix, c);
  } else {
    for (int d = 0; d < 3; ++d) r.dx[d] = static_cast<float>(c[d]);
  }
  for (int d = 0; d < 3; ++d) r.x[d] = static_cast<float>(U(rng) * N);
  r.r = static_cast<float>(U(rng) * 100);
  r.dr = static_cast<float>(rrk::rmax_from_nside(r.nside) - r.r);
  const double edge[6] = {pc::HI_eV, std::sqrt(pc::HI_eV * pc::HeI_eV), pc::HeI_eV,
                          std::sqrt(pc::HeI_eV * pc::HeII_eV), pc::HeII_eV, 100.0};
  for (int g = 0; g < 5; ++g) {
    r.nu[g] = static_cast<float>(edge[g] + U(rng) * (edge[g + 1] - edge[g]));
    r.ngam[g] = static_cast<float>(std::pow(10.0, -4 + 4 * U(rng)) * (U(rng) < 0.1 ? 0 : 1));
  }
  if (rrf::sum4<5>(r.ngam) <= 0) r.ngam[0] = 1e-3f;
  return r;
}

int main() {
  const int N = 32, n = 2000;
  std::mt19937_64 rng(7);
  // float layout: exact
  {
    using E = amrex_enc::FloatEnc;
    FakeTile<E::kNR, E::kNI> t(n);
    auto p = t.data();
    const auto ctx = E::make_ctx(N, 1);
    std::vector<RayP> in;
    for (int k = 0; k < n; ++k) {
      in.push_back(random_ray(rng, N, k % 2));
      E::store(p, k, in[k], ctx);
    }
    for (int k = 0; k < n; ++k) {
      const RayP o = E::load(p, k, ctx);
      CHECK(std::memcmp(&o, &in[k], sizeof o) == 0);
    }
  }
  // R-struct: within the quantisation steps
  {
    using E = amrex_enc::RStructEnc;
    FakeTile<E::kNR, E::kNI> t(n);
    auto p = t.data();
    const auto ctx = E::make_ctx(N, 3);
    std::vector<RayP> in;
    double worst_pos = 0, worst_r = 0, worst_dir = 0, worst_nu = 0, worst_shape = 0;
    int on_centre_exact = 0, merged = 0;
    for (int k = 0; k < n; ++k) {
      in.push_back(random_ray(rng, N, k % 2));
      E::store(p, k, in[k], ctx);
      // the backend sets the particle position to the home-cell centre
      for (int d = 0; d < 3; ++d) p.m_rdata[d][k] = static_cast<float>(rrk::cell_of(in[k].x[d], N)) + 0.5f;
    }
    for (int k = 0; k < n; ++k) {
      const RayP a = in[k], o = E::load(p, k, ctx);
      CHECK(o.nside == a.nside);
      CHECK(o.ipix == a.ipix);
      for (int d = 0; d < 3; ++d) {
        CHECK(rrk::cell_of(o.x[d], N) == rrk::cell_of(a.x[d], N));
        worst_pos = std::max(worst_pos, std::fabs(static_cast<double>(o.x[d]) - a.x[d]));
      }
      worst_r = std::max(worst_r, std::fabs(static_cast<double>(o.r) - a.r));
      CHECK(std::fabs(static_cast<double>(o.dr) - (rrk::rmax_from_nside(o.nside) - o.r)) < 1e-4);
      double dd = 0;
      for (int d = 0; d < 3; ++d) dd += (static_cast<double>(o.dx[d]) - a.dx[d]) * (static_cast<double>(o.dx[d]) - a.dx[d]);
      dd = std::sqrt(dd);
      const double ps = E::pixel_radius(a.nside);
      if (k % 2) {
        merged += 1;
        worst_dir = std::max(worst_dir, dd / ps);
      } else {
        double c[3];
        rrhp::pix2vec_nest(a.nside, a.ipix, c);
        bool exact = true;
        for (int d = 0; d < 3; ++d) exact = exact && o.dx[d] == static_cast<float>(c[d]);
        on_centre_exact += exact;
      }
      // photons: the total is exact, each group within 1/16 octave (dithered)
      CHECK(rrf::sum4<5>(o.ngam) == rrf::sum4<5>(a.ngam) || std::fabs(rrf::sum4<5>(o.ngam) / rrf::sum4<5>(a.ngam) - 1) < 1e-5);
      const double tot = rrf::sum4<5>(a.ngam);
      for (int g = 0; g < 5; ++g) {
        const double fa = a.ngam[g] / tot;
        if (fa > std::exp2(-15.5)) {
          const double fo = o.ngam[g] / rrf::sum4<5>(o.ngam);
          worst_shape = std::max(worst_shape, std::fabs(std::log2(fo / fa)));
        } else {
          CHECK(o.ngam[g] <= a.ngam[g] * 2 + 1e-30f || fa <= std::exp2(-15.5));
        }
        const double edge[6] = {pc::HI_eV, std::sqrt(pc::HI_eV * pc::HeI_eV), pc::HeI_eV,
                                std::sqrt(pc::HeI_eV * pc::HeII_eV), pc::HeII_eV, 100.0};
        worst_nu = std::max(worst_nu, std::fabs(static_cast<double>(o.nu[g]) - a.nu[g]) / (edge[g + 1] - edge[g]));
      }
    }
    std::printf("rstruct: worst position %.2e cell, r %.3f cell, direction %.4f pixel radii (merged), "
                "nu %.4f band, shape %.3f octave; on-centre exact %d/%d\n",
                worst_pos, worst_r, worst_dir, worst_nu, worst_shape, on_centre_exact, n - merged);
    CHECK(worst_pos <= 1.0 / 65536 + 1e-5);
    CHECK(worst_r <= 1.0 / 8 + 1e-6);
    CHECK(worst_dir <= 1.5 / 128 * 1.5 + 0.02);  // dithered 8-bit offsets, two axes
    CHECK(worst_nu <= 1.0 / 255 + 1e-6);
    // dithered rounding errs by up to one 1/16-octave step and the
    // renormalisation by the groups' sum adds up to another
    CHECK(worst_shape <= 2.0 / 16 + 0.01);
    CHECK(on_centre_exact == n - merged);
  }
  TEST_MAIN_END();
}
