// Lossy ray-precision experiment (include/radray/ray_quant.hpp).
#include "radray/ray_quant.hpp"

#include <cmath>
#include <cstdint>

#include "radray/constants.hpp"
#include "radray/healpix.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/simulation.hpp"

namespace radray {

namespace {

inline std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ull;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
  return x ^ (x >> 31);
}

// Uniform in [0,1) from a deterministic key.
inline double uniform(std::uint64_t key) {
  return static_cast<double>(splitmix64(key) >> 11) * (1.0 / 9007199254740992.0);
}

// v > 0 rounded to `bits` mantissa bits; stochastic if u >= 0.
inline double round_mantissa(double v, int bits, double u) {
  if (!(v > 0)) return v;
  int e;
  std::frexp(v, &e);  // v = m 2^e, m in [0.5,1)
  const double ulp = std::ldexp(1.0, e - bits);
  const double lo = std::floor(v / ulp) * ulp;
  if (u < 0) return (v - lo >= 0.5 * ulp) ? lo + ulp : lo;
  return (u < (v - lo) / ulp) ? lo + ulp : lo;
}

// x rounded to a multiple of step; stochastic if u >= 0.
inline double round_step(double x, double step, double u) {
  const double lo = std::floor(x / step) * step;
  if (u < 0) return (x - lo >= 0.5 * step) ? lo + step : lo;
  return (u < (x - lo) / step) ? lo + step : lo;
}

// ray_redshifting's band edges (Nfreq = 5)
const double kFreqEdges[6] = {pc::HI_eV, std::sqrt(pc::HI_eV * pc::HeI_eV), pc::HeI_eV,
                              std::sqrt(pc::HeI_eV * pc::HeII_eV), pc::HeII_eV,
                              static_cast<double>(1e2f)};

void tangent_basis(const double c[3], double eth[3], double eph[3]) {
  const double rho = std::sqrt(c[0] * c[0] + c[1] * c[1]);
  if (rho > 0) {
    const double cp = c[0] / rho, sp = c[1] / rho;
    eth[0] = c[2] * cp; eth[1] = c[2] * sp; eth[2] = -rho;
    eph[0] = -sp; eph[1] = cp; eph[2] = 0;
  } else {
    eth[0] = c[2]; eth[1] = 0; eth[2] = 0;
    eph[0] = 0; eph[1] = 1; eph[2] = 0;
  }
}

void structural_round_trip(Ray& ry, int N, bool dither, std::uint64_t key) {
  // position: 16-bit offset inside the cell (the cell does not change)
  for (int d = 0; d < 3; ++d) {
    const double x = ry.x[d];
    const double c = std::floor(x);
    double q = std::floor((x - c) * 65536);
    q = std::fmin(std::fmax(q, 0.0), 65535.0);
    double xn = c + (q + 0.5) / 65536;
    if (xn >= N) xn -= N;
    ry.x[d] = static_cast<float>(xn);
  }
  // r to 1/8 cell, dr from nside
  const double r = std::round(static_cast<double>(ry.r) * 8) / 8;
  const double rmax = rmax_from_nside(ry.nside);
  ry.r = static_cast<float>(r);
  ry.dr = static_cast<float>(rmax - r);
  // direction: centre + offset quantised to 1/128 pixel
  double c[3], eth[3], eph[3];
  healpix::pix2vec_nest(ry.nside, ry.ipix, c);
  tangent_basis(c, eth, eph);
  const double ps = std::sqrt(4 * 3.14159265358979323846 / (12.0 * ry.nside * ry.nside));
  double t = 0, p = 0;
  for (int d = 0; d < 3; ++d) {
    t += (ry.dx[d] - c[d]) * eth[d];
    p += (ry.dx[d] - c[d]) * eph[d];
  }
  const double qt = std::fmax(-127.0, std::fmin(127.0, std::round(t / ps * 128))) * ps / 128;
  const double qp = std::fmax(-127.0, std::fmin(127.0, std::round(p / ps * 128))) * ps / 128;
  double v[3], nrm = 0;
  for (int d = 0; d < 3; ++d) {
    v[d] = c[d] + qt * eth[d] + qp * eph[d];
    nrm += v[d] * v[d];
  }
  nrm = std::sqrt(nrm);
  for (int d = 0; d < 3; ++d) ry.dx[d] = static_cast<float>(v[d] / nrm);
  // nu: 8 bits inside the band
  for (int n = 0; n < kNfreq; ++n) {
    const double lo = kFreqEdges[n], hi = kFreqEdges[n + 1];
    const double pos = std::fmin(std::fmax((ry.nu[n] - lo) / (hi - lo), 0.0), 1.0);
    ry.nu[n] = static_cast<float>(lo + std::round(pos * 255) / 255 * (hi - lo));
  }
  // photons: exact float total, log2 fractions to 1/128 octave
  float tot = 0;
  double dtot = 0;
  for (int n = 0; n < kNfreq; ++n) {
    tot += ry.ngam[n];
    dtot += ry.ngam[n];
  }
  if (!(dtot > 0)) return;
  double f[kNfreq], fs = 0;
  for (int n = 0; n < kNfreq; ++n) {
    const double fr = ry.ngam[n] / dtot;
    if (fr > 0 && fr >= std::exp2(-4095.0 / 128)) {
      const double u = dither ? uniform(key ^ (0x51ull + n)) : -1;
      f[n] = std::exp2(round_step(std::log2(fr), 1.0 / 128, u));
    } else {
      f[n] = 0;
    }
    fs += f[n];
  }
  for (int n = 0; n < kNfreq; ++n) ry.ngam[n] = static_cast<float>(tot * f[n] / fs);
}

}  // namespace

void ray_quantise(Simulation& s, const RayQuantConfig& qc) {
  if (!qc.enabled()) return;
  RTGrid& G = s.rt;
  RayState& R = s.rays;
  const std::uint64_t step_key =
      (static_cast<std::uint64_t>(s.clk.istep) << 40) ^ (static_cast<std::uint64_t>(s.clk.istep_rt) << 20);
  for (std::int64_t c = 0; c < G.ncell; ++c) {
    std::int64_t ir = G.htoc_ray(3, c);
    const std::int64_t tail = G.htoc_ray(4, c);
    while (ir > 0) {
      Ray& ry = R.at(ir);
      const std::uint64_t key = splitmix64(step_key ^ (static_cast<std::uint64_t>(ir) * 0x100000001b3ull));
      if (qc.structural) structural_round_trip(ry, s.cfg.Nrt, qc.dither, key);
      if (qc.bits > 0) {
        for (int n = 0; n < kNfreq; ++n) {
          const double u = qc.dither ? uniform(key ^ (0x9e37ull + n)) : -1;
          ry.ngam[n] = static_cast<float>(round_mantissa(ry.ngam[n], qc.bits, u));
        }
      }
      if (ir == tail) break;
      ir = ry.next;
    }
  }
}

}  // namespace radray
