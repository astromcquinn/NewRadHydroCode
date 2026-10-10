// Ray kernels of the adaptive ray tracing (ray_tracing_tools.f90), written
// once for every host: the serial port (src/ray_tracing.cpp, per-cell ray
// chains) and the sorted-ray-array RT of the AMReX backend
// (src/amrex/rt_level.cpp). Header-only, host/device (rr_hd.hpp), no
// library state, no exceptions.
//
// Provenance: spike/common/ray_kernels.hpp of the framework-spike branch,
// which was a line-by-line copy of the serial port's arithmetic. This copy
// is the maintained one (docs/stage_d.md, milestone 1); the serial port
// calls these functions, so cast, advance, split, merge, emission and
// redshifting exist once. What stays with the host is only bookkeeping and
// ORDER: the serial port deposits into float32 cells in chain order and
// merges in chain order (bit-identical to the Fortran); the sorted-array
// path deposits into int64 fixed-point accumulators and merges in the
// canonical order (cone, payload hash, payload bytes), which makes its
// result independent of thread count, rank count and decomposition.
//
// Precision: ray storage is float (real(4)), arithmetic is double
// (real(8)), exactly as the Fortran. Cells are 0-based global indices
// (i,j,k), i fastest, periodic box of N cells; a ray's home cell is
// cell_of(x). The serial port's 1-based (irt,jrt,krt) is (i,j,k) + 1.
#pragma once

#include <cstdint>
#include <cstring>

#include "radray/constants.hpp"
#include "radray/kernels/consts_hd.hpp"
#include "radray/kernels/healpix_hd.hpp"
#include "radray/kernels/rr_hd.hpp"

namespace rrk {

constexpr int NF = radray::kNfreq;  // 5 frequency groups
constexpr int kNdep = 13;           // ngam(5), ngam*nu(5), dx*sum(ngam)(3)

// Ray payload: radhydro_global.f90 ray_type without the chain links.
// 2 ints + 18 floats = 80 bytes, no padding. The serial port's radray::Ray
// derives from it (adding the links), so a Ray is a RayP.
struct RayP {
  std::int32_t nside = 0, ipix = 0;
  float r = 0, dr = 0;
  float x[3] = {0, 0, 0};
  float dx[3] = {0, 0, 0};
  float nu[NF] = {0, 0, 0, 0, 0};
  float ngam[NF] = {0, 0, 0, 0, 0};
};
static_assert(sizeof(RayP) == 80, "RayP must be 80 bytes without padding");

// Fixed-point scales of the deposit accumulators of the sorted-array path
// (powers of two, so the products v*S are exact). Resolution 5.7e-14 for
// ngam and dx*ngam (photons per H atom), 1.5e-11 for ngam*nu (nu in eV,
// <= 100). RRK_FX_BITS (default 44) can be lowered to test the sensitivity
// of the results to round-off (spike/serial rt_serial_fx36).
#ifndef RRK_FX_BITS
#define RRK_FX_BITS 44
#endif
constexpr double kScaleNgam = static_cast<double>(1ULL << RRK_FX_BITS);
constexpr double kScaleNu = static_cast<double>(1ULL << (RRK_FX_BITS - 8));
constexpr double kScaleTot = 1099511627776.0;  // 2^40 (diagnostic sums)
RR_HD double dep_scale(int m) { return (m >= 5 && m < 10) ? kScaleNu : kScaleNgam; }
RR_HD long long to_fx(double v, double S) {
  return static_cast<long long>(std::floor(v * S + 0.5));
}
RR_HD double from_fx(long long a, double S) { return static_cast<double>(a) / S; }

// Run-time parameters of one RT substep (from Config/Units/SimClock).
struct RTP {
  int N = 0;
  double nrtr = 0, hrtr = 0;
  double seg_min = 0, ngam_min = 0, fgam_min = 0, sfr_min = 0;
  double vol_rt = 0, nH_avg = 0, dtsec_rt = 0;
  int nseg = 4;        // 3 with legacy_rt = 1 (B8), 4 otherwise
  int nside_merge = 2;  // HEALPix nside of the merge cones (12 nside^2 cones)
  double ngam_ray16 = 0, ngam_ray32 = 0;  // ngam_ray(16), ngam_ray(32)
};
constexpr int kNconeMax = 192;  // 12 * 4^2: nside_merge <= 4

// ---------------------------------------------------------------------------
// Geometry helpers (ray_tracing_tools.f90)
// ---------------------------------------------------------------------------

RR_HD int cell_of(double x, int N) { return rrf::mod(rrf::int4(x), N); }

// ray_tracing_tools.f90: lside_from_radius. log(4.) is real(4).
RR_HD int lside_from_radius(double r) {
  using namespace radray;
  if (r < 1) return rayp::lside_min;
  const int lside = rrf::ceiling4(
      std::log(pc::pi * rayp::ray_per_cell / 3.0 * rrf::pow2(r)) / folded::log4_r4);
  return rrf::min(rrf::max(lside, rayp::lside_min), rayp::lside_max);
}
RR_HD int nside_from_radius(double r) { return 1 << lside_from_radius(r); }
RR_HD int nray_from_nside(int nside) {
  using namespace radray;
  if (nside < rayp::nside_min) return rayp::nray_init;
  return 12 * nside * nside;
}
RR_HD double rmin_from_nside(int nside) {
  using namespace radray;
  if (nside == rayp::nside_min) return 0;
  return (nside / 2) * std::sqrt(3 / (pc::pi * rayp::ray_per_cell));
}
RR_HD double rmax_from_nside(int nside) {
  using namespace radray;
  return nside * std::sqrt(3 / (pc::pi * rayp::ray_per_cell));
}

// ray_tracing_tools.f90: ray_segment -- lengths of the cell crossings of a
// unit-length step from x0 along dx. nseg = 3 (the original; can truncate
// the step, B8) or 4 (legacy_rt = 0: always covers it). Elements beyond
// nseg are 0.
RR_HD void ray_segment(const double x0[3], const double dx[3], double seg_min, int nseg,
                       double seg[4]) {
  double l = 0;
  double x[3] = {x0[0], x0[1], x0[2]};
  const double eps = static_cast<double>(1e-6f);  // 1E-6 is real(4)
  for (int n = 0; n < 4; ++n) seg[n] = 0;
  for (int n = 0; n < nseg; ++n) {
    double dl = 1 - l;
    if (dl > seg_min) {
      for (int i = 0; i < 3; ++i) {
        if (std::fabs(dx[i]) > eps) {
          double ds;
          if (dx[i] < 0) {
            ds = rrf::floor4(x[i]) - x[i];
            if (ds > -eps) ds = -1;
          } else {
            ds = rrf::ceiling4(x[i]) - x[i];
            if (ds < +eps) ds = +1;
          }
          dl = rrf::min(dl, ds / dx[i]);
        }
      }
    }
    l = l + dl;
    for (int i = 0; i < 3; ++i) x[i] = x[i] + dl * dx[i];
    seg[n] = dl;
  }
}

// sum(ray%ngam): real(4) sum, then widened
RR_HD double sum_ngam(const RayP& r) { return static_cast<double>(rrf::sum4<NF>(r.ngam)); }

// ---------------------------------------------------------------------------
// ray_cast_1: deposit all photons of a ray along its <= nseg segments.
// tau(i,j,k) -> const float* (5 values); dep(i,j,k, const double v[13])
// receives the 13 addends (dngam(5), dngam*nu(5), dx*sum(dngam)(3)) of one
// segment. The host decides how to accumulate them (float32 cells in order
// in the serial port, int64 fixed point in the sorted-array path).
// ---------------------------------------------------------------------------
template <class TauF, class DepF>
RR_HD void cast_1(const RayP& ry, const RTP& P, const TauF& tau, const DepF& dep) {
  double x[3], dx[3], nu[NF], ngam[NF];
  for (int d = 0; d < 3; ++d) {
    x[d] = ry.x[d];
    dx[d] = ry.dx[d];
  }
  for (int n = 0; n < NF; ++n) {
    nu[n] = ry.nu[n];
    ngam[n] = ry.ngam[n];
  }
  double dl[4];
  ray_segment(x, dx, P.seg_min, P.nseg, dl);
  double l = 0;
  for (int n = 0; n < P.nseg; ++n) {
    l = l + dl[n];
    if (dl[n] >= P.seg_min) {
      double x1[3];
      for (int d = 0; d < 3; ++d)
        x1[d] = rrf::mod(x[d] + (l - dl[n] / 2) * dx[d] + P.nrtr, P.nrtr);
      const int i = cell_of(x1[0], P.N), j = cell_of(x1[1], P.N), k = cell_of(x1[2], P.N);
      // fdl = dl(n)/sum(dl(n:))
      double sdl = 0;
      for (int m = n; m < P.nseg; ++m) sdl += dl[m];
      const double fdl = dl[n] / sdl;
      const float* t = tau(i, j, k);
      double dng[NF];
      for (int q = 0; q < NF; ++q) {
        const double fabs_ = 1 - std::exp(static_cast<double>(-t[q]) * dl[n]);
        if (fabs_ < fdl) {
          dng[q] = ngam[q] * fdl;
          ngam[q] = ngam[q] - dng[q];
        } else {
          dng[q] = ngam[q] * fabs_;
          ngam[q] = ngam[q] - dng[q];
        }
      }
      double v[kNdep];
      for (int q = 0; q < NF; ++q) v[q] = dng[q];
      for (int q = 0; q < NF; ++q) v[5 + q] = dng[q] * nu[q];
      const double sdng = rrf::sum8<NF>(dng);
      for (int d = 0; d < 3; ++d) v[10 + d] = dx[d] * sdng;
      dep(i, j, k, v);
    }
  }
}

// ---------------------------------------------------------------------------
// ray_cast_3: attenuate by fgam, advance, delete faint rays.
// tf(i,j,k, const float*& tau, const float*& fgam). Returns false if the
// ray is deleted (its fields are then unchanged); the host relinks it to
// cell_of(ry.x) otherwise.
// ---------------------------------------------------------------------------
template <class TFF>
RR_HD bool cast_3(RayP& ry, const RTP& P, const TFF& tf) {
  double r = ry.r, dr = ry.dr, x[3], dx[3], ngam[NF];
  for (int d = 0; d < 3; ++d) {
    x[d] = ry.x[d];
    dx[d] = ry.dx[d];
  }
  for (int n = 0; n < NF; ++n) ngam[n] = ry.ngam[n];
  double dl[4];
  ray_segment(x, dx, P.seg_min, P.nseg, dl);
  // Save original
  const double r0 = r, dr0 = dr;
  const double x0[3] = {x[0], x[1], x[2]};
  double ng1[NF], ng2[NF];
  for (int n = 0; n < NF; ++n) {
    ng1[n] = ngam[n];
    ng2[n] = ngam[n];
  }
  const double ngam_cut = rrf::max(P.ngam_min, rrf::sum8<NF>(ngam) * P.fgam_min);
  double l = 0;
  for (int n = 0; n < NF; ++n) ngam[n] = 0;
  for (int n = 0; n < P.nseg; ++n) {
    l = l + dl[n];
    if (dl[n] >= P.seg_min && rrf::sum8<NF>(ng2) > ngam_cut) {
      double x1[3];
      for (int d = 0; d < 3; ++d)
        x1[d] = rrf::mod(x0[d] + (l - dl[n] / 2) * dx[d] + P.nrtr, P.nrtr);
      const int i = cell_of(x1[0], P.N), j = cell_of(x1[1], P.N), k = cell_of(x1[2], P.N);
      double sdl = 0;
      for (int m = n; m < P.nseg; ++m) sdl += dl[m];
      const double fdl = dl[n] / sdl;
      const float* tau;
      const float* fgam;
      tf(i, j, k, tau, fgam);
      double dng[NF];
      for (int q = 0; q < NF; ++q) {
        const double fabs_ = 1 - std::exp(static_cast<double>(-tau[q]) * dl[n]);
        if (fabs_ < fdl) {
          dng[q] = ng1[q] * fdl;
          ng1[q] = ng1[q] - dng[q];
        } else {
          dng[q] = ng1[q] * fabs_;
          ng1[q] = ng1[q] - dng[q];
        }
      }
      // Consume photons along ray segment: (1-fgam) is real(4) arithmetic
      for (int q = 0; q < NF; ++q) ng2[q] = ng2[q] - dng[q] * static_cast<double>(1.0f - fgam[q]);
      if (rrf::sum8<NF>(ng2) > ngam_cut) {
        // Advance ray to end of segment (and slightly beyond); 1E-6 real(4)
        for (int d = 0; d < 3; ++d)
          x[d] = rrf::mod(x0[d] + (l + static_cast<double>(1e-6f)) * dx[d] + P.nrtr, P.nrtr);
        r = r0 + l;
        dr = dr0 - l;
        for (int q = 0; q < NF; ++q) ngam[q] = ngam[q] + dng[q] * static_cast<double>(fgam[q]);
      }
    }
  }
  if (rrf::sum8<NF>(ngam) < ngam_cut) return false;
  ry.r = static_cast<float>(r);
  ry.dr = static_cast<float>(dr);
  for (int d = 0; d < 3; ++d) ry.x[d] = static_cast<float>(x[d]);
  for (int q = 0; q < NF; ++q) ry.ngam[q] = static_cast<float>(ngam[q]);
  return true;
}

// ---------------------------------------------------------------------------
// ray_cast_split: 4 daughters (child[0] corresponds to the parent's slot).
// clamped (optional) reports whether the daughters' radius was clamped to
// [rmin, rmin+1] of the new level, i.e. the effective origin moved
// (ray-dump provenance of the serial port).
// ---------------------------------------------------------------------------
RR_HD void split(const RayP& p, const RTP& P, RayP child[4], bool* clamped = nullptr) {
  using namespace radray;
  const double nrtr = P.nrtr;
  int nside = p.nside, ipix;
  double r = p.r, dr, x[3], dx[3], nu[NF], ngam[NF];
  for (int d = 0; d < 3; ++d) {
    x[d] = p.x[d];
    dx[d] = p.dx[d];
  }
  for (int n = 0; n < NF; ++n) {
    nu[n] = p.nu[n];
    ngam[n] = p.ngam[n];
  }
  // Current indices: deviation of the ray direction from its pixel centre
  const int ipixc = p.ipix;
  double dc[3];
  rrhp::pix2vec_nest(nside, ipixc, dc);
  for (int d = 0; d < 3; ++d) dc[d] = dx[d] - dc[d];
  // Daughter rays
  nside = rrf::min(2 * nside, rayp::nside_max);
  const double rmin = rmin_from_nside(nside);
  const double rmax = rmax_from_nside(nside);
  const double r_parent = r;
  r = rrf::min(rrf::max(r, rmin), rmin + 1);
  dr = rmax - r;
  if (clamped) *clamped = (r != r_parent);
  // Effective origin, position corresponding to r = 0
  double x0[3];
  for (int d = 0; d < 3; ++d) x0[d] = rrf::mod(x[d] - rrf::mod(r * dx[d], nrtr) + nrtr, nrtr);
  // Divide photon count
  for (int n = 0; n < NF; ++n) ngam[n] = ngam[n] / 4;
  for (int n = 1; n <= 4; ++n) {
    ipix = 4 * ipixc + (n - 1);
    // New direction: default direction of new pixel + parent deviation
    rrhp::pix2vec_nest(nside, ipix, dx);
    for (int d = 0; d < 3; ++d) dx[d] = dx[d] + dc[d];
    const double sq[3] = {dx[0] * dx[0], dx[1] * dx[1], dx[2] * dx[2]};
    const double nrm = std::sqrt(rrf::sum8<3>(sq));
    for (int d = 0; d < 3; ++d) dx[d] = dx[d] / nrm;
    ipix = rrhp::vec2pix_nest(nside, dx);
    // New position
    for (int d = 0; d < 3; ++d) x[d] = rrf::mod(x0[d] + rrf::mod(r * dx[d], nrtr) + nrtr, nrtr);
    RayP& c = child[n - 1];
    c.nside = nside;
    c.ipix = ipix;
    c.r = static_cast<float>(r);
    c.dr = static_cast<float>(dr);
    for (int d = 0; d < 3; ++d) {
      c.x[d] = static_cast<float>(x[d]);
      c.dx[d] = static_cast<float>(dx[d]);
    }
    for (int k = 0; k < NF; ++k) {
      c.nu[k] = static_cast<float>(nu[k]);
      c.ngam[k] = static_cast<float>(ngam[k]);
    }
  }
}

// ---------------------------------------------------------------------------
// ray_cast_merge as a reduction: add() each ray of the cone in the host's
// order (chain order in the serial port, canonical order in the sorted-array
// path), then finish(). cc = (/irt,jrt,krt/) - 0.5 in real(4): the centre
// of the home cell.
// ---------------------------------------------------------------------------
struct MergeAcc {
  double r = 0, dr = 0, ntot = 0;
  double x[3] = {0, 0, 0}, dx[3] = {0, 0, 0};
  double nu[NF] = {0, 0, 0, 0, 0}, ngam[NF] = {0, 0, 0, 0, 0};
};

RR_HD void cell_centre(int i, int j, int k, float cc[3]) {
  // (/irt,jrt,krt/) - 0.5 in real(4), irt = i + 1
  cc[0] = static_cast<float>(i + 1) - 0.5f;
  cc[1] = static_cast<float>(j + 1) - 0.5f;
  cc[2] = static_cast<float>(k + 1) - 0.5f;
}

RR_HD void merge_add(MergeAcc& a, const RayP& ry, const float cc[3], const RTP& P) {
  double xd[3];
  // xd = ray(ir)%x - ((/irt,jrt,krt/)-0.5): real(4) subtraction
  for (int d = 0; d < 3; ++d) xd[d] = static_cast<double>(ry.x[d] - cc[d]);
  for (int d = 0; d < 3; ++d) {
    if (std::fabs(xd[d]) > P.hrtr) {
      if (xd[d] > 0)
        xd[d] = xd[d] - P.nrtr;
      else
        xd[d] = xd[d] + P.nrtr;
    }
  }
  const double w = sum_ngam(ry);
  a.ntot = a.ntot + w;
  for (int d = 0; d < 3; ++d) a.x[d] = a.x[d] + w * xd[d];
  // ray(ir)%r**2 is a real(4) product
  a.r = a.r + w * static_cast<double>(ry.r * ry.r);
  a.dr = a.dr + w * static_cast<double>(ry.dr);
  for (int d = 0; d < 3; ++d) a.dx[d] = a.dx[d] + w * static_cast<double>(ry.dx[d]);
  // ray(ir)%ngam*ray(ir)%nu is a real(4) product
  for (int n = 0; n < NF; ++n) a.nu[n] = a.nu[n] + static_cast<double>(ry.ngam[n] * ry.nu[n]);
  for (int n = 0; n < NF; ++n) a.ngam[n] = a.ngam[n] + static_cast<double>(ry.ngam[n]);
}

RR_HD RayP merge_finish(const MergeAcc& a, const float cc[3], const RTP& P) {
  using namespace radray;
  double r = a.r / a.ntot;
  const double dr = a.dr / a.ntot;
  double xd[3], dx[3], x[3], nu[NF], ngam[NF];
  for (int d = 0; d < 3; ++d) xd[d] = a.x[d] / a.ntot;
  for (int d = 0; d < 3; ++d) dx[d] = a.dx[d] / a.ntot;
  r = std::sqrt(r);
  {
    const double nrm = std::sqrt((0.0 + dx[0] * dx[0]) + dx[1] * dx[1] + dx[2] * dx[2]);
    for (int d = 0; d < 3; ++d) dx[d] = dx[d] / nrm;
  }
  for (int d = 0; d < 3; ++d) x[d] = rrf::mod(static_cast<double>(cc[d]) + xd[d] + P.nrtr, P.nrtr);
  // Conserve photon number; conserve energy via effective frequency
  for (int n = 0; n < NF; ++n) {
    if (a.ngam[n] < P.ngam_min) {
      nu[n] = rrk::freq_init(n);
      ngam[n] = 0;
    } else {
      nu[n] = a.nu[n] / a.ngam[n];
      ngam[n] = a.ngam[n];
    }
  }
  // nside from mean radius, then distance from splitting, final radius
  const int nside = nside_from_radius(r);
  const double rmin = rmin_from_nside(nside);
  const double rmax = rmax_from_nside(nside);
  r = rrf::max(rmax - dr, rmin);
  const double drn = rmax - r;
  RayP o;
  o.nside = nside;
  o.ipix = rrhp::vec2pix_nest(nside, dx);
  o.r = static_cast<float>(r);
  o.dr = static_cast<float>(drn);
  for (int d = 0; d < 3; ++d) {
    o.x[d] = static_cast<float>(x[d]);
    o.dx[d] = static_cast<float>(dx[d]);
  }
  for (int n = 0; n < NF; ++n) {
    o.nu[n] = static_cast<float>(nu[n]);
    o.ngam[n] = static_cast<float>(ngam[n]);
  }
  return o;
}

// ---------------------------------------------------------------------------
// ray_source
// ---------------------------------------------------------------------------

// exp(-minval(tau)) < 1E-2, all in real(4): optically thick source cell
RR_HD bool source_is_thick(const float* tau) {
  float tmin = tau[0];
  for (int n = 1; n < NF; ++n)
    if (tau[n] < tmin) tmin = tau[n];
  return std::exp(-tmin) < 1e-2f;
}

// The rays of one emission: 12 nside^2 rays from the centre of cell
// (i,j,k) at r = 0, with the same spectrum (nu, ngam per ray).
struct SourceSpec {
  int nside = 0;
  float x[3] = {0, 0, 0};
  float dr = 0;
  float nu[NF] = {0, 0, 0, 0, 0}, ngam[NF] = {0, 0, 0, 0, 0};
};

// Photons per group of a stellar source cell (ray_source):
// fResc*pop2_init*(Luv*dtsec_rt)/vol_rt/nH_avg, per <nH> and RT cell volume.
RR_HD void source_photons(float Luv, const RTP& P, double ngam[NF]) {
  using namespace radray;
  const double Ldt = static_cast<double>(Luv) * P.dtsec_rt;
  for (int n = 0; n < NF; ++n) {
    ngam[n] = cp::fResc * rrk::pop2_init(n) * Ldt;
    ngam[n] = ngam[n] / P.vol_rt / P.nH_avg;
  }
}

// Direct deposit of a thick source cell (ngam, ngam*nu; no direction)
RR_HD void source_thick_deposit(float Luv, const RTP& P, double v[kNdep]) {
  double ngam[NF];
  source_photons(Luv, P, ngam);
  for (int n = 0; n < NF; ++n) {
    v[n] = ngam[n];
    v[5 + n] = ngam[n] * rrk::freq_init(n);
  }
  for (int d = 0; d < 3; ++d) v[10 + d] = 0;
}

// Emission of nside_init^2*12 rays of a stellar source cell (ray_source);
// nu = freq_init, ngam = the cell's photons / nray_init per ray.
RR_HD SourceSpec source_spec(int i, int j, int k, float Luv, const RTP& P) {
  using namespace radray;
  SourceSpec s;
  s.nside = rayp::nside_init;
  cell_centre(i, j, k, s.x);
  s.dr = static_cast<float>(rmax_from_nside(rayp::nside_init));
  double ngam[NF];
  source_photons(Luv, P, ngam);
  for (int n = 0; n < NF; ++n) {
    s.nu[n] = static_cast<float>(rrk::freq_init(n));
    s.ngam[n] = static_cast<float>(ngam[n] / rayp::nray_init);
  }
  return s;
}

// Emission of 12 nside^2 rays with a given spectrum (ray_source_recomb of
// the recombination radiation, recomb_rad.hpp): nu[n], ngam[n] per ray.
RR_HD SourceSpec source_spec_spectrum(int i, int j, int k, int nside, const double nu[NF],
                                      const double ngam[NF]) {
  SourceSpec s;
  s.nside = nside;
  cell_centre(i, j, k, s.x);
  s.dr = static_cast<float>(rmax_from_nside(nside));
  for (int n = 0; n < NF; ++n) {
    s.nu[n] = static_cast<float>(nu[n]);
    s.ngam[n] = static_cast<float>(ngam[n]);
  }
  return s;
}

// Ray ipix (0 .. 12 nside^2 - 1) of an emission
RR_HD RayP source_ray(const SourceSpec& s, int ipix) {
  RayP ry;
  ry.nside = s.nside;
  ry.ipix = ipix;
  ry.r = 0;
  ry.dr = s.dr;
  double dx[3];
  rrhp::pix2vec_nest(s.nside, ipix, dx);
  for (int d = 0; d < 3; ++d) {
    ry.x[d] = s.x[d];
    ry.dx[d] = static_cast<float>(dx[d]);
  }
  for (int n = 0; n < NF; ++n) {
    ry.nu[n] = s.nu[n];
    ry.ngam[n] = s.ngam[n];
  }
  return ry;
}
// Source ray ipix (0..191) of the stellar source in cell (i,j,k)
RR_HD RayP source_ray(int ipix, int i, int j, int k, float Luv, const RTP& P) {
  return source_ray(source_spec(i, j, k, Luv, P), ipix);
}

// ---------------------------------------------------------------------------
// Classification (ray_tracing_1; last match wins)
// ---------------------------------------------------------------------------
enum RayCase : int { kMerge = 0, kCast = 1, kSplit = 2, kRemove = 3 };
RR_HD int classify(const RayP& ry, const RTP& P) {
  const double ntot = sum_ngam(ry);
  int c = kMerge;
  if (ntot > P.ngam_ray16) c = kCast;
  if (ry.r == 0) c = kCast;
  if (ry.dr < 0 && ntot > P.ngam_ray32) c = kSplit;
  if (ntot < P.ngam_min) c = kRemove;
  if (ry.nside == radray::rayp::nside_max) c = kRemove;
  return c;
}
// Merge cone of a ray: HEALPix pixel of its direction at nside_merge
// (2 in the original: 48 cones; recomb_prepare may lower it to 1)
RR_HD int merge_cone(const RayP& ry, int nside_merge) {
  const double dxv[3] = {static_cast<double>(ry.dx[0]), static_cast<double>(ry.dx[1]),
                         static_cast<double>(ry.dx[2])};
  return rrhp::vec2pix_nest(nside_merge, dxv);
}

// Canonical order of the sorted-array path: FNV-1a hash of the payload
// bytes, then the bytes themselves.
RR_HD std::uint32_t payload_hash(const RayP& r) {
  const unsigned char* b = reinterpret_cast<const unsigned char*>(&r);
  std::uint32_t h = 2166136261u;
  for (int n = 0; n < static_cast<int>(sizeof(RayP)); ++n) {
    h ^= b[n];
    h *= 16777619u;
  }
  return h;
}
RR_HD int payload_cmp(const RayP& a, const RayP& b) {
  const unsigned char* pa = reinterpret_cast<const unsigned char*>(&a);
  const unsigned char* pb = reinterpret_cast<const unsigned char*>(&b);
  for (int n = 0; n < static_cast<int>(sizeof(RayP)); ++n)
    if (pa[n] != pb[n]) return pa[n] < pb[n] ? -1 : 1;
  return 0;
}

// Photon-count histogram of rt_step_1/calc_ray_1 (bins 1..N_pdf over
// log10(sum ngam) in [-6, 3) with lgn_del = log10(1.01))
constexpr int kNpdfMax = 2100;
constexpr double kLgnMin = -6.0, kLgnMax = 3.0;
RR_HD int pdf_nbins() {
  return static_cast<int>(std::round((kLgnMax - kLgnMin) / radray::folded::lgn_del));
}
RR_HD int pdf_bin(double ntot, int N_pdf) {
  int n = 1 + rrf::int4((std::log10(ntot) - kLgnMin) / radray::folded::lgn_del);
  return rrf::min(rrf::max(n, 1), N_pdf);
}

// ---------------------------------------------------------------------------
// Sorted-array path: per-cell substep-1 work on the rays of one cell
// (cell-sorted input). Two passes: count (exact number of output rays),
// then write. The thick-source deposit is the host's (source_thick_deposit).
// ---------------------------------------------------------------------------
struct CellStats {
  long long ncst = 0, nspl = 0, nmrg = 0, ndel = 0;
};

// Count pass. `in` holds the n rays whose home is cell (i,j,k). `src`
// (optional) is a stellar emission that is not thick (nullptr: none).
RR_HD int step1_count(const RayP* in, int n, const SourceSpec* src, const RTP& P, CellStats& st) {
  int nout = 0;
  if (src) {
    const RayP s0 = source_ray(*src, 0);  // all rays of an emission have equal ngam
    const int nr = nray_from_nside(src->nside);
    if (classify(s0, P) == kRemove) {
      st.ndel += nr;
    } else {
      nout += nr;
      st.ncst += nr;
    }
  }
  const int ncone = 12 * P.nside_merge * P.nside_merge;
  int cone_n[kNconeMax];
  for (int c = 0; c < ncone; ++c) cone_n[c] = 0;
  for (int m = 0; m < n; ++m) {
    switch (classify(in[m], P)) {
      case kCast: nout += 1; st.ncst += 1; break;
      case kSplit: nout += 4; st.ncst += 4; st.nspl += 4; break;
      case kRemove: st.ndel += 1; break;
      default: cone_n[merge_cone(in[m], P.nside_merge)] += 1; break;
    }
  }
  for (int c = 0; c < ncone; ++c) {
    if (cone_n[c] >= 1) {
      nout += 1;
      st.ncst += 1;
      st.nmrg += cone_n[c] - 1;
      st.ndel += cone_n[c] - 1;
    }
  }
  return nout;
}

// Heap sort of the merge-class indices idx[0..m) by (cone, hash, bytes)
struct MergeKey {
  int cone;
  std::uint32_t hash;
};
RR_HD bool key_less(const RayP* in, const MergeKey* key, int a, int b) {
  if (key[a].cone != key[b].cone) return key[a].cone < key[b].cone;
  if (key[a].hash != key[b].hash) return key[a].hash < key[b].hash;
  return payload_cmp(in[a], in[b]) < 0;
}
RR_HD void sift_down(int* idx, int start, int end, const RayP* in, const MergeKey* key) {
  int root = start;
  while (2 * root + 1 < end) {
    int child = 2 * root + 1;
    if (child + 1 < end && key_less(in, key, idx[child], idx[child + 1])) child += 1;
    if (key_less(in, key, idx[root], idx[child])) {
      const int t = idx[root];
      idx[root] = idx[child];
      idx[child] = t;
      root = child;
    } else {
      return;
    }
  }
}
RR_HD void heap_sort(int* idx, int m, const RayP* in, const MergeKey* key) {
  for (int s = m / 2 - 1; s >= 0; --s) sift_down(idx, s, m, in, key);
  for (int e = m - 1; e > 0; --e) {
    const int t = idx[0];
    idx[0] = idx[e];
    idx[e] = t;
    sift_down(idx, 0, e, in, key);
  }
}

// Write pass. idx/key: scratch of length n. out: room for the count of
// step1_count. Returns the number of rays written.
RR_HD int step1_write(const RayP* in, int n, const SourceSpec* src, int i, int j, int k,
                      const RTP& P, int* idx, MergeKey* key, RayP* out) {
  int w = 0;
  if (src) {
    const RayP s0 = source_ray(*src, 0);
    if (classify(s0, P) != kRemove) {
      const int nr = nray_from_nside(src->nside);
      for (int ip = 0; ip < nr; ++ip) out[w++] = source_ray(*src, ip);
    }
  }
  int nm = 0;
  for (int m = 0; m < n; ++m) {
    const int c = classify(in[m], P);
    if (c == kCast) {
      out[w++] = in[m];
    } else if (c == kSplit) {
      split(in[m], P, out + w);
      w += 4;
    } else if (c == kMerge) {
      key[m].cone = merge_cone(in[m], P.nside_merge);
      key[m].hash = payload_hash(in[m]);
      idx[nm++] = m;
    }
  }
  if (nm > 0) {
    heap_sort(idx, nm, in, key);
    float cc[3];
    cell_centre(i, j, k, cc);
    int a = 0;
    while (a < nm) {
      int b = a + 1;
      while (b < nm && key[idx[b]].cone == key[idx[a]].cone) ++b;
      if (b - a == 1) {
        out[w++] = in[idx[a]];
      } else {
        MergeAcc acc;
        for (int q = a; q < b; ++q) merge_add(acc, in[idx[q]], cc, P);
        out[w++] = merge_finish(acc, cc, P);
      }
      a = b;
    }
  }
  return w;
}

// ---------------------------------------------------------------------------
// ray_redshifting (non-FREQ1): the tables and the per-spectrum update,
// shared by the rays (serial port and sorted-array path) and the moment RT.
// ---------------------------------------------------------------------------
constexpr int kZbin = 121;  // 1 + nint((bmax-bmin)/bdel)
struct RedshiftTable {
  double bmin = -6, bmax = 6;
  double bdel = static_cast<double>(0.1f);
  int Nbin = kZbin;
  double freq[NF + 1] = {};
  // indx(2,Nbin,Nfreq-1) stored as tab[j][2*(i-1)+(row-1)]: (value, b) pairs
  double tab[NF - 1][2 * kZbin] = {};
  double rmin = 0, rmax = 0;
};

// general_tools.f90 interpolate(x,f,a,b,'lin') for f(2,n) stored as pairs
// f[2*i+0] = f(1,i+1), f[2*i+1] = f(2,i+1); a,b in {1,2}.
RR_HD double interp_lin(double x, const double* fa, int n, int a = 1, int b = 2) {
  auto F = [&](int row, int i) { return fa[2 * (i - 1) + (row - 1)]; };
  const double dxmin = static_cast<double>(1.1920928955078125e-07f);  // epsilon(real(0.))
  int i1 = 1, i2 = n;
  if (F(a, 1) < F(a, n)) {
    while (i2 - i1 > 1) {
      const int i = (i1 + i2) / 2;
      if (x > F(a, i)) i1 = i; else i2 = i;
    }
  } else {
    while (i2 - i1 > 1) {
      const int i = (i1 + i2) / 2;
      if (x < F(a, i)) i1 = i; else i2 = i;
    }
  }
  const double x1 = F(a, i1), y1 = F(b, i1), x2 = F(a, i2), y2 = F(b, i2);
  const double ddx = x2 - x1;
  if (std::fabs(ddx) < dxmin) return (y1 + y2) / 2;
  return y1 + (y2 - y1) * ((x - x1) / (x2 - x1));
}

// ray_tracing_tools.f90: the tables of ray_redshifting (Nfreq = 5 case).
// Host only (built once per hydro step).
inline RedshiftTable build_redshift_table() {
  using namespace radray;
  RedshiftTable T;
  double* freq = T.freq;
  freq[0] = pc::HI_eV;
  freq[1] = std::sqrt(pc::HI_eV * pc::HeI_eV);
  freq[2] = pc::HeI_eV;
  freq[3] = std::sqrt(pc::HeI_eV * pc::HeII_eV);
  freq[4] = pc::HeII_eV;
  freq[5] = static_cast<double>(1e2f);
  double rmax = 0., rmin = static_cast<double>(3.40282346638528859812e+38F);  // huge(0.)
  for (int j = 1; j <= NF - 1; ++j) {
    const double fj = freq[j - 1], fj1 = freq[j], fj2 = freq[j + 1];
    for (int i = 1; i <= T.Nbin; ++i) {
      const double b = (i - 1) * T.bdel + T.bmin;
      const double c = b + 1;
      double v;
      if (std::fabs(c) < T.bdel / 1e2)
        v = std::log(fj2 / fj1) / std::log(fj1 / fj);
      else
        v = (std::pow(fj2, c) - std::pow(fj1, c)) / (std::pow(fj1, c) - std::pow(fj, c));
      T.tab[j - 1][2 * (i - 1) + 1] = b;
      T.tab[j - 1][2 * (i - 1) + 0] = v;
      rmax = rrf::max(rmax, v);
      rmin = rrf::min(rmin, v);
    }
  }
  T.rmin = rmin;
  T.rmax = rmax;
  return T;
}

// redshift_ray (#else FREQ1) of one 5-group spectrum without the rays'
// final minimum: nu, ngam (per group) are updated in place. keep[g] and
// down[g] (optional) are the fractions of the photons of group g that
// stay in g and move to g-1 (both 0 for a group at or below ngam_min).
RR_HD void redshift_spectrum(const RedshiftTable& T, double fscale, double ngam_min, double* nu,
                             double* ngam, double* keep = nullptr, double* down = nullptr) {
  using namespace radray;
  const double* freq = T.freq;
  double egam[NF], ngam1[NF], egam1[NF];
  for (int i = 0; i < NF; ++i) {
    egam[i] = nu[i] * ngam[i];
    ngam1[i] = 0;
    egam1[i] = 0;
    if (keep) keep[i] = 0;
    if (down) down[i] = 0;
  }
  for (int i = 1; i <= NF; ++i) {
    const int im = i - 1, ip = i + 1;
    if (ngam[i - 1] > ngam_min) {
      const double fq = freq[i - 1] / fscale;
      double rN, rE, b;
      if (i < NF) {
        rN = ngam[ip - 1] / ngam[i - 1];
        rN = rrf::min(rrf::max(rN, T.rmin), T.rmax);
        b = interp_lin(rN, T.tab[i - 1], T.Nbin);
        b = rrf::min(rrf::max(b, T.bmin), T.bmax);
      } else {
        rN = ngam[i - 1] / rrf::max(ngam[im - 1], ngam_min);
        rN = rrf::min(rrf::max(rN, T.rmin), T.rmax);
        b = interp_lin(rN, T.tab[im - 1], T.Nbin);
        b = rrf::min(rrf::max(b, T.bmin), T.bmax);
      }
      const double c = b + 1;
      const double d = b + 2;
      const double fi = freq[i - 1], fip = freq[ip - 1];
      if (std::fabs(c) < T.bdel / 1e2)
        rN = std::log(fq / fi) / std::log(fip / fi);
      else
        rN = (std::pow(fq, c) - std::pow(fi, c)) / (std::pow(fip, c) - std::pow(fi, c));
      if (std::fabs(d) < T.bdel / 1e2)
        rE = std::log(fq / fi) / std::log(fip / fi);
      else
        rE = (std::pow(fq, d) - std::pow(fi, d)) / (std::pow(fip, d) - std::pow(fi, d));
      rN = rrf::min(rrf::max(rN, 0.0), 1.0);
      rE = rrf::min(rrf::max(rE, 0.0), 1.0);
      ngam1[i - 1] = (1 - rN) * ngam[i - 1];
      egam1[i - 1] = (1 - rE) * egam[i - 1];
      if (keep) keep[i - 1] = 1 - rN;
      if (i > 1) {
        ngam1[im - 1] = ngam1[im - 1] + rN * ngam[i - 1];
        egam1[im - 1] = egam1[im - 1] + rE * egam[i - 1];
        if (down) down[i - 1] = rN;
      }
    }
  }
  for (int i = 1; i <= NF; ++i) {
    if (ngam1[i - 1] < ngam_min) {
      nu[i - 1] = rrk::freq_init(i - 1);
      ngam[i - 1] = 0;
    } else {
      nu[i - 1] = egam1[i - 1] / ngam1[i - 1] * fscale;
      nu[i - 1] = rrf::min(rrf::max(nu[i - 1], freq[i - 1]), freq[i]);
      ngam[i - 1] = ngam1[i - 1];
    }
  }
}

// ray_redshifting: one ray (float storage, double arithmetic), with the
// minimum for stability (real(4) sum)
RR_HD void redshift_ray(RayP& ry, const RedshiftTable& T, double fscale, double ngam_min) {
  double nu[NF], ngam[NF];
  for (int i = 0; i < NF; ++i) {
    nu[i] = ry.nu[i];
    ngam[i] = ry.ngam[i];
  }
  redshift_spectrum(T, fscale, ngam_min, nu, ngam);
  for (int i = 0; i < NF; ++i) {
    ry.nu[i] = static_cast<float>(nu[i]);
    ry.ngam[i] = static_cast<float>(ngam[i]);
  }
  if (static_cast<double>(rrf::sum4<NF>(ry.ngam)) < ngam_min)
    ry.ngam[0] = static_cast<float>(ngam_min);
}

}  // namespace rrk
