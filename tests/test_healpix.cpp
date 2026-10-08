// Self-consistency tests of pix2vec_nest / vec2pix_nest.
// (Bit-level comparison with healpy is in healpix_vs_healpy.py.)
#include <cmath>
#include <cstdint>
#include <random>

#include "check.hpp"
#include "radray/healpix.hpp"

using namespace radray;

int main() {
  // Pixel centres map back to themselves, unit vectors
  for (int l = 0; l <= 6; ++l) {
    const int nside = 1 << l;
    const int npix = 12 * nside * nside;
    for (int p = 0; p < npix; ++p) {
      double v[3];
      healpix::pix2vec_nest(nside, p, v);
      CHECK(std::fabs(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] - 1) < 1e-14);
      CHECK(healpix::vec2pix_nest(nside, v) == p);
    }
  }
  // Large nside: sample pixels
  std::mt19937_64 rng(1);
  for (int l = 7; l <= 13; ++l) {
    const int nside = 1 << l;
    std::uniform_int_distribution<std::int32_t> dist(0, 12 * nside * nside - 1);
    for (int n = 0; n < 20000; ++n) {
      const int p = dist(rng);
      double v[3];
      healpix::pix2vec_nest(nside, p, v);
      CHECK(healpix::vec2pix_nest(nside, v) == p);
    }
  }
  // Nested hierarchy: parent of child = p/4
  std::normal_distribution<double> g;
  for (int n = 0; n < 100000; ++n) {
    double v[3] = {g(rng), g(rng), g(rng)};
    const int p8 = healpix::vec2pix_nest(8, v);
    const int p16 = healpix::vec2pix_nest(16, v);
    CHECK(p16 / 4 == p8);
  }
  // Poles and axis vectors
  {
    const double npole[3] = {0, 0, 1}, spole[3] = {0, 0, -1};
    for (int l = 0; l <= 13; ++l) {
      const int nside = 1 << l;
      const int np = healpix::vec2pix_nest(nside, npole);
      const int sp = healpix::vec2pix_nest(nside, spole);
      CHECK(np >= 0 && np < 12 * nside * nside);
      CHECK(sp >= 0 && sp < 12 * nside * nside);
    }
  }
  TEST_MAIN_END();
}
