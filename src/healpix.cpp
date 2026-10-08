// HEALPix 3.11 (Fortran) NESTED pixel routines, see healpix.hpp.
#include "radray/healpix.hpp"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "radray/fortran.hpp"

namespace radray {
namespace healpix {

namespace {

// healpix_types.F90
constexpr double HALFPI = 1.570796326794896619231321691639751442099;
constexpr double TWOPI = 6.283185307179586476925286766559005768394;
constexpr double TWOTHIRD = 0.6666666666666666666666666666666666666666;

// pix_tools.F90: mk_pix2xy -- pix2x/pix2y for 10-bit pixel numbers
// (de-interleave even/odd bits). mk_xy2pix1 -- x2pix1/y2pix1 for 7-bit
// coordinates (interleave). Built once.
struct Tables {
  std::int32_t pix2x[1024];
  std::int32_t pix2y[1024];
  std::int32_t x2pix1[128];
  std::int32_t y2pix1[128];
  Tables() {
    for (int kpix = 0; kpix < 1024; ++kpix) {
      int jpix = kpix, ix = 0, iy = 0, ip = 1;
      while (jpix != 0) {
        int id = jpix % 2;
        jpix /= 2;
        ix += id * ip;
        id = jpix % 2;
        jpix /= 2;
        iy += id * ip;
        ip *= 2;
      }
      pix2x[kpix] = ix;
      pix2y[kpix] = iy;
    }
    for (int i = 0; i < 128; ++i) {
      int j = i, k = 0, ip = 1;
      while (j != 0) {
        int id = j % 2;
        j /= 2;
        k += ip * id;
        ip *= 4;
      }
      x2pix1[i] = k;
      y2pix1[i] = 2 * k;
    }
  }
};

const Tables& tables() {
  static const Tables t;
  return t;
}

// Fortran jrll/jpll (1-based in Fortran, indexed by face_num+1).
constexpr std::int32_t jrll[12] = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
constexpr std::int32_t jpll[12] = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};

}  // namespace

void pix2vec_nest(std::int32_t nside, std::int32_t ipix, double vec[3]) {
  if (nside < 1 || nside > ns_max4)
    throw std::runtime_error("pix2vec_nest> nside out of range");
  const std::int32_t npix = 12 * nside * nside;
  if (ipix < 0 || ipix > npix - 1)
    throw std::runtime_error("pix2vec_nest> ipix out of range");
  const Tables& T = tables();

  const std::int32_t npface = nside * nside;
  const std::int32_t nl4 = 4 * nside;

  const std::int32_t face_num = ipix / npface;
  // MODULO(ipix,npface) with ipix >= 0
  const std::int32_t ipf = ipix % npface;

  const double fn = static_cast<double>(nside);
  const double fact1 = 1.0 / (3.0 * fn * fn);
  const double fact2 = 2.0 / (3.0 * fn);

  const std::int32_t ip_low = ipf & 1023;
  const std::int32_t ip_trunc = ipf / 1024;
  const std::int32_t ip_med = ip_trunc & 1023;
  const std::int32_t ip_hi = ip_trunc / 1024;

  const std::int32_t ix =
      1024 * T.pix2x[ip_hi] + 32 * T.pix2x[ip_med] + T.pix2x[ip_low];
  const std::int32_t iy =
      1024 * T.pix2y[ip_hi] + 32 * T.pix2y[ip_med] + T.pix2y[ip_low];

  const std::int32_t jrt = ix + iy;
  const std::int32_t jpt = ix - iy;

  const std::int32_t jr = jrll[face_num] * nside - jrt - 1;

  std::int32_t nr = 0, kshift = 0;
  double z = 0, sth = 0;
  if (jr < nside) {  // north pole region
    nr = jr;
    z = 1.0 - nr * fact1 * nr;
    sth = nr * std::sqrt(fact1 * (1.0 + z));
    kshift = 0;
  } else if (jr <= 3 * nside) {  // equatorial region
    nr = nside;
    z = (2 * nside - jr) * fact2;
    sth = std::sqrt((1.0 - z) * (1.0 + z));
    kshift = (jr - nside) & 1;
  } else {  // south pole region
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

std::int32_t vec2pix_nest(std::int32_t nside, const double vec[3]) {
  if (nside < 1 || nside > ns_max4)
    throw std::runtime_error("vec2pix_nest> nside out of range");
  const Tables& T = tables();

  const double dnorm =
      std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  const double z = vec[2] / dnorm;
  double phi = 0.0;
  if (vec[0] != 0.0 || vec[1] != 0.0) phi = std::atan2(vec[1], vec[0]);

  const double za = std::fabs(z);
  if (phi < 0.0) phi = phi + TWOPI;
  const double tt = phi / HALFPI;

  std::int32_t face_num, ix, iy;
  if (za <= TWOTHIRD) {  // equatorial region
    const std::int32_t jp = f::int4(nside * (0.5 + tt - z * 0.75));
    const std::int32_t jm = f::int4(nside * (0.5 + tt + z * 0.75));

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
  } else {  // polar region
    std::int32_t ntt = f::int4(tt);
    if (ntt >= 4) ntt = 3;
    const double tp = tt - ntt;
    double tmp = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1]) / dnorm;
    tmp = tmp * std::sqrt(3.0 / (1.0 + za));

    std::int32_t jp = f::int4(nside * tp * tmp);
    std::int32_t jm = f::int4(nside * (1.0 - tp) * tmp);
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
  const std::int32_t ipf = T.x2pix1[ix_low] + T.y2pix1[iy_low] +
                           (T.x2pix1[ix / 128] + T.y2pix1[iy / 128]) * 16384;
  return ipf + face_num * nside * nside;
}

}  // namespace healpix
}  // namespace radray
