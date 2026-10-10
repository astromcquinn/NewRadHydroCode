// Ray storage layouts of the AMReX backend (ray_encoding = float | rstruct;
// docs/stage_d.md milestone 3, docs/ray_compression.md section 4 and 8).
// Each encoding maps an rrk::RayP to and from the SoA components of a
// particle tile (the 3 ParticleReal positions are the backend's: the centre
// of the home cell, which Redistribute uses). Header-only, host/device.
//
//   FloatEnc    the original's float fields: 18 float + 2 int components
//               (80 B of ray data; the bitwise-validation layout).
//   RStructEnc  the R-struct: 7 int32 components (28 B of ray data):
//               position as 3 x 16-bit offsets in the home cell (1/65536
//               cell), a 30-bit HEALPix id (4 nside^2 + ipix) with the
//               on-centre flag, the direction as an 8+8-bit tangent-plane
//               offset from the pixel centre (1/128 of the pixel radius;
//               absent for never-merged rays, whose direction is the
//               centre), r at 1/8 cell (dr = rmax(nside) - r), the photon
//               total as a float, the spectral shape as 5 x 8-bit
//               log-quantised fractions (1/16 octave over 16 octaves) and nu
//               as 5 x 8 bits linear in the group's band. Every lossy
//               field is rounded with deterministic dithering keyed on
//               (substep, the ray's float payload, field), so the
//               quantisation is unbiased and the result does not depend
//               on the thread or rank layout.
#pragma once

#include <cstdint>
#include <cstring>

#include "radray/constants.hpp"
#include "radray/kernels/ray_kernels.hpp"

// Usable without AMReX (tests/test_ray_encoding.cpp): the tile data type is
// anything with m_rdata[comp][n] and m_idata[comp][n]
#ifndef AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_HOST_DEVICE
#endif
#ifndef AMREX_FORCE_INLINE
#define AMREX_FORCE_INLINE inline
#endif

namespace radray {
namespace amrex_enc {

struct EncCtx {
  int N = 0;
  std::uint64_t substep = 0;
};

// ---------------------------------------------------------------------------
// float: the original's fields
// ---------------------------------------------------------------------------
struct FloatEnc {
  // real comps: 0-2 position (home-cell centre), 3 r, 4 dr, 5-7 x, 8-10 dx,
  // 11-15 nu, 16-20 ngam; int comps: 0 nside, 1 ipix.
  static constexpr int kNR = 3 + 18, kNI = 2;
  static constexpr int kBytes = 80;
  static constexpr int kR = 3, kDR = 4, kX = 5, kDX = 8, kNU = 11, kNGAM = 16;
  using Ctx = EncCtx;
  static Ctx make_ctx(int N, std::uint64_t substep) { return Ctx{N, substep}; }

  template <class PTD>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static rrk::RayP load(PTD const& p, int n, const Ctx&) {
    rrk::RayP r;
    r.nside = p.m_idata[0][n];
    r.ipix = p.m_idata[1][n];
    r.r = p.m_rdata[kR][n];
    r.dr = p.m_rdata[kDR][n];
    for (int d = 0; d < 3; ++d) {
      r.x[d] = p.m_rdata[kX + d][n];
      r.dx[d] = p.m_rdata[kDX + d][n];
    }
    for (int q = 0; q < rrk::NF; ++q) {
      r.nu[q] = p.m_rdata[kNU + q][n];
      r.ngam[q] = p.m_rdata[kNGAM + q][n];
    }
    return r;
  }
  template <class PTD>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void store(PTD const& p, int n, const rrk::RayP& r, const Ctx&) {
    p.m_idata[0][n] = r.nside;
    p.m_idata[1][n] = r.ipix;
    p.m_rdata[kR][n] = r.r;
    p.m_rdata[kDR][n] = r.dr;
    for (int d = 0; d < 3; ++d) {
      p.m_rdata[kX + d][n] = r.x[d];
      p.m_rdata[kDX + d][n] = r.dx[d];
    }
    for (int q = 0; q < rrk::NF; ++q) {
      p.m_rdata[kNU + q][n] = r.nu[q];
      p.m_rdata[kNGAM + q][n] = r.ngam[q];
    }
  }
};

// ---------------------------------------------------------------------------
// rstruct: the compact R-struct
// ---------------------------------------------------------------------------

// Deterministic dithering: u in [0,1) from a hash of (key, field)
RR_HD std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
RR_HD double dither_u(std::uint64_t key, int field) {
  return static_cast<double>(splitmix64(key ^ (0x5bd1e995ULL * static_cast<std::uint64_t>(field + 1))) >> 11) *
         (1.0 / 9007199254740992.0);
}
// floor(v + u): an unbiased integer estimate of v
RR_HD long long round_dither(double v, double u) { return static_cast<long long>(std::floor(v + u)); }
RR_HD int clampi(long long v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : static_cast<int>(v)); }

// Tangent basis (e_theta, e_phi) at the unit vector c
RR_HD void tangent_basis(const double c[3], double et[3], double ep[3]) {
  const double rxy = std::sqrt(c[0] * c[0] + c[1] * c[1]);
  if (rxy > 1e-12) {
    ep[0] = -c[1] / rxy;
    ep[1] = c[0] / rxy;
    ep[2] = 0;
  } else {
    ep[0] = 0;
    ep[1] = 1;
    ep[2] = 0;
  }
  // e_theta = e_phi x c
  et[0] = ep[1] * c[2] - ep[2] * c[1];
  et[1] = ep[2] * c[0] - ep[0] * c[2];
  et[2] = ep[0] * c[1] - ep[1] * c[0];
}

struct RStructEnc {
  static constexpr int kNR = 3, kNI = 7;
  static constexpr int kBytes = 28;
  static constexpr double kOctaves = 16.0;     // range of the log shape
  static constexpr double kShapeStep = 16.0;   // codes per octave (8 bits: 255 = zero)
  static constexpr double kPosM = 65536.0;     // position quantum 1/65536 cell
  static constexpr double kRQ = 8.0;           // r quantum 1/8 cell
  static constexpr double kDirQ = 128.0;       // direction quantum: pixel radius / 128
  static constexpr double kOnCentre = 2e-6;    // |dx - centre| below: never merged
  struct Ctx {
    int N = 0;
    std::uint64_t substep = 0;
    double edge[rrk::NF + 1];  // the frequency bands (redshift table edges)
  };
  static Ctx make_ctx(int N, std::uint64_t substep) {
    Ctx c;
    c.N = N;
    c.substep = substep;
    c.edge[0] = pc::HI_eV;
    c.edge[1] = std::sqrt(pc::HI_eV * pc::HeI_eV);
    c.edge[2] = pc::HeI_eV;
    c.edge[3] = std::sqrt(pc::HeI_eV * pc::HeII_eV);
    c.edge[4] = pc::HeII_eV;
    c.edge[5] = static_cast<double>(1e2f);
    return c;
  }
  RR_HD static double pixel_radius(int nside) {
    return std::sqrt(4 * pc::pi / (12.0 * nside * nside));
  }

  template <class PTD>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void store(PTD const& p, int n, const rrk::RayP& r, const Ctx& c) {
    const std::uint64_t key = splitmix64((c.substep << 32) ^ static_cast<std::uint64_t>(rrk::payload_hash(r)));
    std::uint32_t w[7] = {0, 0, 0, 0, 0, 0, 0};
    // position: offset in the home cell
    std::uint32_t q[3];
    for (int d = 0; d < 3; ++d) {
      const int cell = rrk::cell_of(r.x[d], c.N);
      const double frac = static_cast<double>(r.x[d]) - cell;  // [0,1)
      q[d] = static_cast<std::uint32_t>(clampi(round_dither(frac * kPosM, dither_u(key, d)), 0, 65535));
    }
    // r at 1/8 cell
    const std::uint32_t qr = static_cast<std::uint32_t>(clampi(round_dither(static_cast<double>(r.r) * kRQ, dither_u(key, 3)), 0, 65535));
    w[0] = q[0] | (q[1] << 16);
    w[1] = q[2] | (qr << 16);
    // HEALPix id and the direction offset
    double cen[3];
    rrhp::pix2vec_nest(r.nside, r.ipix, cen);
    bool on_centre = true;
    for (int d = 0; d < 3; ++d)
      if (std::fabs(static_cast<double>(r.dx[d]) - static_cast<double>(static_cast<float>(cen[d]))) > kOnCentre)
        on_centre = false;
    const std::uint32_t id = static_cast<std::uint32_t>(4 * r.nside) * static_cast<std::uint32_t>(r.nside) +
                             static_cast<std::uint32_t>(r.ipix);
    w[2] = id | (on_centre ? (1u << 30) : 0u);
    std::uint32_t qt = 128, qp = 128;
    if (!on_centre) {
      double et[3], ep[3];
      tangent_basis(cen, et, ep);
      double t = 0, f = 0;
      for (int d = 0; d < 3; ++d) {
        t += (static_cast<double>(r.dx[d]) - cen[d]) * et[d];
        f += (static_cast<double>(r.dx[d]) - cen[d]) * ep[d];
      }
      const double sc = kDirQ / pixel_radius(r.nside);
      qt = static_cast<std::uint32_t>(clampi(round_dither(t * sc, dither_u(key, 4)) + 128, 1, 255));
      qp = static_cast<std::uint32_t>(clampi(round_dither(f * sc, dither_u(key, 5)) + 128, 1, 255));
    }
    // photons: float total, log-quantised shape
    const float tot = rrf::sum4<rrk::NF>(r.ngam);
    std::uint32_t qs[rrk::NF];
    for (int g = 0; g < rrk::NF; ++g) {
      const double fr = tot > 0 ? static_cast<double>(r.ngam[g]) / static_cast<double>(tot) : 0.0;
      if (fr > 0) {
        const double v = -std::log2(fr) * kShapeStep;  // >= 0
        const long long s = round_dither(v, dither_u(key, 6 + g));
        qs[g] = static_cast<std::uint32_t>(s < 0 ? 0 : (s >= 255 ? 255 : s));
      } else {
        qs[g] = 255;
      }
    }
    // nu: linear in the band
    std::uint32_t qn[rrk::NF];
    for (int g = 0; g < rrk::NF; ++g) {
      const double lo = c.edge[g], hi = c.edge[g + 1];
      const double pos = (static_cast<double>(r.nu[g]) - lo) / (hi - lo);
      qn[g] = static_cast<std::uint32_t>(clampi(round_dither(pos * 255.0, dither_u(key, 11 + g)), 0, 255));
    }
    w[3] = qt | (qp << 8) | (qs[0] << 16) | (qs[1] << 24);
    w[4] = qs[2] | (qs[3] << 8) | (qs[4] << 16) | (qn[0] << 24);
    w[5] = qn[1] | (qn[2] << 8) | (qn[3] << 16) | (qn[4] << 24);
    std::memcpy(&w[6], &tot, sizeof tot);
    for (int m = 0; m < 7; ++m) {
      std::int32_t v;
      std::memcpy(&v, &w[m], sizeof v);
      p.m_idata[m][n] = v;
    }
  }

  template <class PTD>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static rrk::RayP load(PTD const& p, int n, const Ctx& c) {
    std::uint32_t w[7];
    for (int m = 0; m < 7; ++m) {
      const std::int32_t v = p.m_idata[m][n];
      std::memcpy(&w[m], &v, sizeof v);
    }
    rrk::RayP r;
    // home cell from the particle position (inside the domain after
    // Redistribute), the offset from the code
    const std::uint32_t q[3] = {w[0] & 0xffffu, w[0] >> 16, w[1] & 0xffffu};
    for (int d = 0; d < 3; ++d) {
      const int cell = rrk::cell_of(static_cast<double>(p.m_rdata[d][n]), c.N);
      r.x[d] = static_cast<float>(cell + q[d] / kPosM);
    }
    // nside, ipix from the id: 4 nside^2 + ipix with ipix < 12 nside^2
    const std::uint32_t id = w[2] & 0x3fffffffu;
    const bool on_centre = (w[2] >> 30) & 1u;
    int nside = rayp::nside_min;
    while (static_cast<std::uint32_t>(16) * nside * nside <= id) nside *= 2;
    r.nside = nside;
    r.ipix = static_cast<std::int32_t>(id - static_cast<std::uint32_t>(4 * nside) * nside);
    // r, dr
    const double rr = (w[1] >> 16) / kRQ;
    r.r = static_cast<float>(rr);
    r.dr = static_cast<float>(rrk::rmax_from_nside(nside) - rr);
    // direction
    double cen[3];
    rrhp::pix2vec_nest(nside, r.ipix, cen);
    if (on_centre) {
      for (int d = 0; d < 3; ++d) r.dx[d] = static_cast<float>(cen[d]);
    } else {
      const int qt = static_cast<int>(w[3] & 0xffu) - 128, qp = static_cast<int>((w[3] >> 8) & 0xffu) - 128;
      double et[3], ep[3];
      tangent_basis(cen, et, ep);
      const double sc = pixel_radius(nside) / kDirQ;
      double v[3], nrm = 0;
      for (int d = 0; d < 3; ++d) {
        v[d] = cen[d] + qt * sc * et[d] + qp * sc * ep[d];
        nrm += v[d] * v[d];
      }
      nrm = std::sqrt(nrm);
      for (int d = 0; d < 3; ++d) r.dx[d] = static_cast<float>(v[d] / nrm);
    }
    // photons
    float tot;
    std::memcpy(&tot, &w[6], sizeof tot);
    const std::uint32_t qs[rrk::NF] = {(w[3] >> 16) & 0xffu, w[3] >> 24, w[4] & 0xffu, (w[4] >> 8) & 0xffu,
                                       (w[4] >> 16) & 0xffu};
    double f[rrk::NF], fs = 0;
    for (int g = 0; g < rrk::NF; ++g) {
      f[g] = qs[g] < 255 ? std::exp2(-static_cast<double>(qs[g]) / kShapeStep) : 0.0;
      fs += f[g];
    }
    for (int g = 0; g < rrk::NF; ++g)
      r.ngam[g] = fs > 0 ? static_cast<float>(static_cast<double>(tot) * f[g] / fs) : 0.0f;
    // nu
    const std::uint32_t qn[rrk::NF] = {w[4] >> 24, w[5] & 0xffu, (w[5] >> 8) & 0xffu, (w[5] >> 16) & 0xffu, w[5] >> 24};
    for (int g = 0; g < rrk::NF; ++g) {
      const double lo = c.edge[g], hi = c.edge[g + 1];
      r.nu[g] = static_cast<float>(lo + qn[g] / 255.0 * (hi - lo));
    }
    return r;
  }
};

}  // namespace amrex_enc
}  // namespace radray
