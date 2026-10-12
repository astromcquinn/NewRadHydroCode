// HEALPix NESTED pix2vec/vec2pix for the kernels (host and device): the same
// arithmetic as src/healpix.cpp (HEALPix 3.11 Fortran), but table-free
// (bit (de)interleaving computed inline) and without exceptions, so it can
// run on a GPU. Integer results are identical; the vector is identical on
// the host (same libm calls in the same order).
#pragma once

#include "radray/kernels/rr_hd.hpp"

namespace rrhp {

// pix2x/pix2y of a 10-bit number: de-interleave even/odd bits
RR_HD void deinterleave10(std::int32_t kpix, std::int32_t& ix, std::int32_t& iy) {
  ix = 0;
  iy = 0;
  for (int b = 0; b < 5; ++b) {
    ix |= ((kpix >> (2 * b)) & 1) << b;
    iy |= ((kpix >> (2 * b + 1)) & 1) << b;
  }
}
// x2pix1 of a 7-bit number: spread bits to even positions
RR_HD std::int32_t spread7(std::int32_t i) {
  std::int32_t k = 0;
  for (int b = 0; b < 7; ++b) k |= ((i >> b) & 1) << (2 * b);
  return k;
}

RR_HD void pix2vec_nest(std::int32_t nside, std::int32_t ipix, double vec[3]) {
  constexpr double HALFPI = 1.570796326794896619231321691639751442099;
  const std::int32_t jrll[12] = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
  const std::int32_t jpll[12] = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};
  const std::int32_t npface = nside * nside;
  const std::int32_t nl4 = 4 * nside;
  const std::int32_t face_num = ipix / npface;
  const std::int32_t ipf = ipix % npface;
  const double fn = static_cast<double>(nside);
  const double fact1 = 1.0 / (3.0 * fn * fn);
  const double fact2 = 2.0 / (3.0 * fn);
  const std::int32_t ip_low = ipf & 1023;
  const std::int32_t ip_trunc = ipf / 1024;
  const std::int32_t ip_med = ip_trunc & 1023;
  const std::int32_t ip_hi = ip_trunc / 1024;
  std::int32_t xl, yl, xm, ym, xh, yh;
  deinterleave10(ip_low, xl, yl);
  deinterleave10(ip_med, xm, ym);
  deinterleave10(ip_hi, xh, yh);
  const std::int32_t ix = 1024 * xh + 32 * xm + xl;
  const std::int32_t iy = 1024 * yh + 32 * ym + yl;
  const std::int32_t jrt = ix + iy;
  const std::int32_t jpt = ix - iy;
  const std::int32_t jr = jrll[face_num] * nside - jrt - 1;
  std::int32_t nr = 0, kshift = 0;
  double z = 0, sth = 0;
  if (jr < nside) {
    nr = jr;
    z = 1.0 - nr * fact1 * nr;
    sth = nr * std::sqrt(fact1 * (1.0 + z));
    kshift = 0;
  } else if (jr <= 3 * nside) {
    nr = nside;
    z = (2 * nside - jr) * fact2;
    sth = std::sqrt((1.0 - z) * (1.0 + z));
    kshift = (jr - nside) & 1;
  } else {
    nr = nl4 - jr;
    z = -1.0 + nr * fact1 * nr;
    sth = nr * std::sqrt(fact1 * (1.0 - z));
    kshift = 0;
  }
  std::int32_t jp = (jpll[face_num] * nr + jpt + 1 + kshift) / 2;
  if (jp > nl4) jp = jp - nl4;
  if (jp < 1) jp = jp + nl4;
  const double phi = (jp - (kshift + 1) * 0.5) * (HALFPI / nr);
  const double cos_phi = std::cos(phi);
  const double sin_phi = std::sin(phi);
  vec[0] = sth * cos_phi;
  vec[1] = sth * sin_phi;
  vec[2] = z;
}

RR_HD std::int32_t vec2pix_nest(std::int32_t nside, const double vec[3]) {
  constexpr double HALFPI = 1.570796326794896619231321691639751442099;
  constexpr double TWOPI = 6.283185307179586476925286766559005768394;
  constexpr double TWOTHIRD = 0.6666666666666666666666666666666666666666;
  const double dnorm = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  const double z = vec[2] / dnorm;
  double phi = 0.0;
  if (vec[0] != 0.0 || vec[1] != 0.0) phi = std::atan2(vec[1], vec[0]);
  const double za = std::fabs(z);
  if (phi < 0.0) phi = phi + TWOPI;
  const double tt = phi / HALFPI;
  std::int32_t face_num, ix, iy;
  if (za <= TWOTHIRD) {
    const std::int32_t jp = rrf::int4(nside * (0.5 + tt - z * 0.75));
    const std::int32_t jm = rrf::int4(nside * (0.5 + tt + z * 0.75));
    const std::int32_t ifp = jp / nside;
    const std::int32_t ifm = jm / nside;
    if (ifp == ifm) {
      face_num = (ifp & 3) + 4;
    } else if (ifp < ifm) {
      face_num = ifp & 3;
    } else {
      face_num = (ifm & 3) + 8;
    }
    ix = jm & (nside - 1);
    iy = nside - (jp & (nside - 1)) - 1;
  } else {
    std::int32_t ntt = rrf::int4(tt);
    if (ntt >= 4) ntt = 3;
    const double tp = tt - ntt;
    double tmp = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1]) / dnorm;
    tmp = tmp * std::sqrt(3.0 / (1.0 + za));
    std::int32_t jp = rrf::int4(nside * tp * tmp);
    std::int32_t jm = rrf::int4(nside * (1.0 - tp) * tmp);
    jp = (nside - 1 < jp) ? nside - 1 : jp;
    jm = (nside - 1 < jm) ? nside - 1 : jm;
    if (z >= 0) {
      face_num = ntt;
      ix = nside - jm - 1;
      iy = nside - jp - 1;
    } else {
      face_num = ntt + 8;
      ix = jp;
      iy = jm;
    }
  }
  const std::int32_t ix_low = ix & 127;
  const std::int32_t iy_low = iy & 127;
  const std::int32_t ipf = spread7(ix_low) + 2 * spread7(iy_low) +
                           (spread7(ix / 128) + 2 * spread7(iy / 128)) * 16384;
  return ipf + face_num * nside * nside;
}

}  // namespace rrhp
