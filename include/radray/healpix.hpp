// HEALPix NESTED pixelisation: pix2vec_nest / vec2pix_nest.
//
// Faithful C++ transcription of the Fortran HEALPix 3.11 routines
// (src/f90/mod/pixel_routines.F90, 32-bit integer path, nside <= 8192),
// which is what the reference code links against (Makefile: Healpix_3.11).
// Floating-point expressions are evaluated in exactly the same order as
// the Fortran so that boundary cases round identically.
#pragma once

#include <cstdint>

namespace radray {
namespace healpix {

constexpr std::int32_t ns_max4 = 8192;

// pixel_routines.F90: pix2vec_nest (no vertex output).
// ipix in [0, 12*nside^2-1]; vec receives the unit vector of the pixel centre.
void pix2vec_nest(std::int32_t nside, std::int32_t ipix, double vec[3]);

// pixel_routines.F90: vec2pix_nest. vec need not be normalised.
std::int32_t vec2pix_nest(std::int32_t nside, const double vec[3]);

}  // namespace healpix
}  // namespace radray
