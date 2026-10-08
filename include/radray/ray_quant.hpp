// Lossy ray-precision experiment for the compression study
// (docs/ray_compression.md, step 6). RESEARCH ONLY, off by default.
//
// When enabled, after every RT substep (after rt_step_3: the point where a
// compressed code would store the rays) each live ray is replaced by a
// quantised copy:
//
//  * bits > 0: every photon group ngam(i) is rounded to `bits` mantissa bits
//    (relative step 2^-bits), stochastically (dither = true: unbiased,
//    E[q] = ngam) or to nearest (dither = false). The dither is a hash of
//    (hydro step, RT substep, storage slot, group), so runs are
//    deterministic.
//  * structural = true: the "R-struct" round trip of
//    tools/ray_compression/ray_codec.cpp: position as a 16-bit offset in its
//    cell, r to 1/8 cell (dr = rmax(nside) - r), direction as the HEALPix
//    centre plus an offset quantised to 1/128 pixel, nu to 8 bits in its
//    band, and the photon groups as the exact float total times log2
//    fractions quantised to 1/128 octave (dithered if dither = true).
//
// This deliberately changes results. It is a measurement instrument, not a
// correction of the original, so it has no Fortran counterpart; with the
// defaults (bits = 0, structural = false) nothing is touched.
#pragma once

namespace radray {

struct Simulation;

struct RayQuantConfig {
  int bits = 0;             // photon mantissa bits (0: off)
  bool dither = true;       // stochastic (unbiased) rounding
  bool structural = false;  // R-struct round trip of geometry, nu and shape
  bool enabled() const { return bits > 0 || structural; }
};

// Quantises every live ray (cell chains) in place.
void ray_quantise(Simulation& s, const RayQuantConfig& q);

}  // namespace radray
