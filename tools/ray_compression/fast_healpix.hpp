// Fast HEALPix NESTED routines for the ray codec (docs/ray_compression.md).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>

// Not the Fortran-faithful radray::healpix used by the RT: pix2vec agrees
// with it to ~1e-15 and vec2pix except at exact pixel-boundary ties.
// Encoder and decoder use the same functions, so the codec is exact anyway.
namespace fasthp {
inline std::uint32_t compact_bits(std::uint32_t v) {
  v &= 0x55555555u;
  v = (v | (v >> 1)) & 0x33333333u;
  v = (v | (v >> 2)) & 0x0f0f0f0fu;
  v = (v | (v >> 4)) & 0x00ff00ffu;
  v = (v | (v >> 8)) & 0x0000ffffu;
  return v;
}
inline std::uint32_t spread_bits(std::uint32_t v) {
  v &= 0x0000ffffu;
  v = (v | (v << 8)) & 0x00ff00ffu;
  v = (v | (v << 4)) & 0x0f0f0f0fu;
  v = (v | (v << 2)) & 0x33333333u;
  v = (v | (v << 1)) & 0x55555555u;
  return v;
}
// sin/cos to ~1e-15 (Cody-Waite reduction to |r| <= pi/4, Taylor to r^15)
inline void fsincos(double x, double& s, double& c) {
  const double k = std::nearbyint(x * (2 / M_PI));
  const double r = (x - k * 1.57079632679489655800e+00) - k * 6.12323399573676603587e-17;
  const double r2 = r * r;
  const double sr = r * (1 + r2 * (-1.0 / 6 + r2 * (1.0 / 120 + r2 * (-1.0 / 5040 + r2 * (1.0 / 362880 +
                    r2 * (-1.0 / 39916800 + r2 * (1.0 / 6227020800.0 + r2 * (-1.0 / 1307674368000.0))))))));
  const double cr = 1 + r2 * (-0.5 + r2 * (1.0 / 24 + r2 * (-1.0 / 720 + r2 * (1.0 / 40320 +
                    r2 * (-1.0 / 3628800 + r2 * (1.0 / 479001600.0 + r2 * (-1.0 / 87178291200.0)))))));
  switch (static_cast<long>(k) & 3) {
    case 0: s = sr; c = cr; break;
    case 1: s = cr; c = -sr; break;
    case 2: s = -sr; c = -cr; break;
    default: s = -cr; c = sr; break;
  }
}
// atan2 in [0, 2 pi) to ~1e-8 (only used to predict pixels)
inline double fatan2_pos(double y, double x) {
  const double ax = std::fabs(x), ay = std::fabs(y);
  const double mx = std::max(ax, ay), mn = std::min(ax, ay);
  double a = mx > 0 ? mn / mx : 0;
  double base = 0;
  if (a > 0.41421356237309503) {  // tan(pi/8)
    a = (a - 1) / (a + 1);
    base = M_PI / 4;
  }
  const double a2 = a * a;
  double t = base + a * (1 + a2 * (-1.0 / 3 + a2 * (1.0 / 5 + a2 * (-1.0 / 7 + a2 * (1.0 / 9 +
             a2 * (-1.0 / 11 + a2 * (1.0 / 13 + a2 * (-1.0 / 15 + a2 * (1.0 / 17)))))))));
  if (ay > ax) t = M_PI / 2 - t;
  if (x < 0) t = M_PI - t;
  if (y < 0) t = 2 * M_PI - t;
  if (t >= 2 * M_PI) t -= 2 * M_PI;
  return t;
}
constexpr int kJrll[12] = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
constexpr int kJpll[12] = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};

inline void pix2vec(int lside, std::uint32_t ipix, double v[3]) {
  const long nside = 1L << lside;
  const long npface = nside * nside;
  const int face = static_cast<int>(ipix >> (2 * lside));
  const std::uint32_t ip = ipix & static_cast<std::uint32_t>(npface - 1);
  const long ix = compact_bits(ip), iy = compact_bits(ip >> 1);
  const long jr = kJrll[face] * nside - ix - iy - 1;
  const double fact2 = 1.0 / (3.0 * npface);
  long nr, kshift;
  double z, sth;
  if (jr < nside) {
    nr = jr;
    const double tmp = nr * nr * fact2;
    z = 1 - tmp;
    sth = std::sqrt(tmp * (2 - tmp));
    kshift = 0;
  } else if (jr > 3 * nside) {
    nr = 4 * nside - jr;
    const double tmp = nr * nr * fact2;
    z = tmp - 1;
    sth = std::sqrt(tmp * (2 - tmp));
    kshift = 0;
  } else {
    nr = nside;
    z = (2 * nside - jr) * (2.0 / (3.0 * nside));
    sth = std::sqrt((1 - z) * (1 + z));
    kshift = (jr - nside) & 1;
  }
  long jp = (kJpll[face] * nr + ix - iy + 1 + kshift) / 2;
  if (jp > 4 * nside) jp -= 4 * nside;
  if (jp < 1) jp += 4 * nside;
  const double phi = (jp - (kshift + 1) * 0.5) * (M_PI / 2 / nr);
  double sp, cp;
  fsincos(phi, sp, cp);
  v[0] = sth * cp;
  v[1] = sth * sp;
  v[2] = z;
}

inline std::uint32_t vec2pix(int lside, const double v[3]) {
  const long nside = 1L << lside;
  const double rr = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  const double z = v[2] / rr, za = std::fabs(z);
  const double phi = fatan2_pos(v[1], v[0]);
  double tt = phi * (2 / M_PI);
  if (tt >= 4) tt -= 4;
  int face;
  long ix, iy;
  if (za <= 2.0 / 3.0) {
    const double t1 = nside * (0.5 + tt), t2 = nside * z * 0.75;
    const long jp = static_cast<long>(t1 - t2), jm = static_cast<long>(t1 + t2);
    const long ifp = jp >> lside, ifm = jm >> lside;
    face = static_cast<int>(ifp == ifm ? (ifp | 4) : (ifp < ifm ? ifp : ifm + 8));
    ix = jm & (nside - 1);
    iy = nside - (jp & (nside - 1)) - 1;
  } else {
    const int ntt = std::min(3, static_cast<int>(tt));
    const double tp = tt - ntt;
    const double tmp = nside * std::sqrt(3 * (1 - za));
    long jp = static_cast<long>(tp * tmp), jm = static_cast<long>((1 - tp) * tmp);
    jp = std::min(jp, nside - 1);
    jm = std::min(jm, nside - 1);
    if (z >= 0) {
      face = ntt;
      ix = nside - jm - 1;
      iy = nside - jp - 1;
    } else {
      face = ntt + 8;
      ix = jp;
      iy = jm;
    }
  }
  return (static_cast<std::uint32_t>(face) << (2 * lside)) + spread_bits(static_cast<std::uint32_t>(ix)) +
         (spread_bits(static_cast<std::uint32_t>(iy)) << 1);
}
}  // namespace fasthp
