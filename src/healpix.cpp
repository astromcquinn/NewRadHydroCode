// HEALPix 3.11 (Fortran) NESTED pixel routines, see healpix.hpp. The
// arithmetic lives in include/radray/kernels/healpix_hd.hpp (rrhp::, table-
// free and device-capable, the same libm calls in the same order); this
// file adds the argument checks of the library (exceptions).
#include "radray/healpix.hpp"

#include <stdexcept>

#include "radray/kernels/healpix_hd.hpp"

namespace radray {
namespace healpix {

void pix2vec_nest(std::int32_t nside, std::int32_t ipix, double vec[3]) {
  if (nside < 1 || nside > ns_max4)
    throw std::runtime_error("pix2vec_nest> nside out of range");
  const std::int32_t npix = 12 * nside * nside;
  if (ipix < 0 || ipix > npix - 1)
    throw std::runtime_error("pix2vec_nest> ipix out of range");
  rrhp::pix2vec_nest(nside, ipix, vec);
}

std::int32_t vec2pix_nest(std::int32_t nside, const double vec[3]) {
  if (nside < 1 || nside > ns_max4)
    throw std::runtime_error("vec2pix_nest> nside out of range");
  return rrhp::vec2pix_nest(nside, vec);
}

}  // namespace healpix
}  // namespace radray
