// Ray-compression prototype (docs/ray_compression.md, steps 4-5).
//
// usage: ray_codec TRAIN_STEM TEST_STEM [--mbits m] [--block n] [--reps k]
//                  [--no-origin] [--json]
//
// Reads two ray snapshots written by the RT's optional ray dump
// (include/radray/ray_dump.hpp): models are fitted on TRAIN and used to
// code TEST (pass the same stem twice for the "oracle" model).
//
// 1. Structural quantisation "Q" (step 4): position as a 16-bit offset in
//    the chain cell, direction as (nside, ipix) + 8-bit tangent-plane offset
//    for merged rays only, r to 1/8 cell (dr dropped), photons as a float32
//    total (optionally m mantissa bits) plus 5 log2 shape fractions to
//    1/128 octave, nu to 8 bits inside its band.  Reported: fixed-width bits
//    per ray and the reconstruction error.
// 2. Entropy coding (step 5) of Q with rANS (32-bit state, byte renorm) and
//    static conditional histograms. Rays whose origin is exact (never
//    merged, split never clamped r, single source) are coded through their
//    origin: source id, r given the source distance, nside given r, and ipix
//    given the pixel predicted from the cell centre; position and direction
//    are then reconstructed as x_src + r*centre(ipix).  Other rays code the
//    Q fields with conditional models.  Rays are coded in blocks of whole
//    cells (>= --block rays) so blocks decode independently.
//
// Output: bytes per ray by field group, encode/decode ns per ray (decode
// includes reconstruction of the float ray fields), and reconstruction
// errors. Lossless with respect to Q except for origin-coded rays, whose
// reconstructed position differs from the stored float by < 1e-4 cell.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "fast_healpix.hpp"
#include "radray/ray_dump.hpp"

using namespace radray;

namespace {

// ---------------------------------------------------------------------------
// snapshot input
// ---------------------------------------------------------------------------

struct Snapshot {
  std::vector<RayDumpRecord> rec;
  int N = 0;
  long long substep = 0;
  std::vector<long long> src_cells;
};

std::string slurp(const std::string& p) {
  std::ifstream in(p, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + p);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

long long json_int(const std::string& js, const std::string& key) {
  const auto p = js.find("\"" + key + "\"");
  if (p == std::string::npos) throw std::runtime_error("json: no key " + key);
  return std::atoll(js.c_str() + js.find(':', p) + 1);
}

Snapshot load(const std::string& stem) {
  Snapshot s;
  const std::string npy = slurp(stem + ".npy");
  if (npy.compare(0, 6, "\x93NUMPY") != 0) throw std::runtime_error("not npy");
  const std::size_t hl = static_cast<unsigned char>(npy[8]) | (static_cast<unsigned char>(npy[9]) << 8);
  const std::size_t off = 10 + hl;
  const std::size_t n = (npy.size() - off) / sizeof(RayDumpRecord);
  s.rec.resize(n);
  std::memcpy(s.rec.data(), npy.data() + off, n * sizeof(RayDumpRecord));
  const std::string js = slurp(stem + ".json");
  s.N = static_cast<int>(json_int(js, "N"));
  s.substep = json_int(js, "substep");
  const auto p = js.find("\"source_cells\"");
  const auto a = js.find('[', p), b = js.find(']', p);
  std::istringstream in(js.substr(a + 1, b - a - 1));
  std::string tok;
  while (std::getline(in, tok, ','))
    if (tok.find_first_not_of(" \n") != std::string::npos) s.src_cells.push_back(std::atoll(tok.c_str()));
  return s;
}

// ---------------------------------------------------------------------------
// helpers shared by encoder and decoder
// ---------------------------------------------------------------------------

constexpr double kSqrt3OverPi = 0.97720502380583984317;  // sqrt(3/pi), ray_per_cell = 1
inline double rmax_of(int nside) { return nside * kSqrt3OverPi; }
inline double rmin_of(int nside) { return nside == 4 ? 0.0 : (nside / 2) * kSqrt3OverPi; }
inline int lside_of_r(double r) {
  if (r < 1) return 2;
  const int l = static_cast<int>(std::ceil(std::log(M_PI / 3.0 * r * r) / std::log(4.0)));
  return std::min(std::max(l, 2), 13);
}
inline int bitlen(std::uint32_t v) { return v ? 32 - __builtin_clz(v) : 0; }
inline std::uint32_t zigzag(int v) { return v >= 0 ? 2u * v : 2u * static_cast<std::uint32_t>(-v) - 1; }
inline int unzigzag(std::uint32_t u) { return (u & 1) ? -static_cast<int>((u + 1) / 2) : static_cast<int>(u / 2); }

// nested ipix <-> (face, ix, iy)
using fasthp::compact_bits;
using fasthp::spread_bits;
inline void nest2xyf(int lside, std::uint32_t ipix, int& face, int& ix, int& iy) {
  const std::uint32_t npf = 1u << (2 * lside);
  face = static_cast<int>(ipix / npf);
  const std::uint32_t ip = ipix & (npf - 1);
  ix = static_cast<int>(compact_bits(ip));
  iy = static_cast<int>(compact_bits(ip >> 1));
}
inline std::uint32_t xyf2nest(int lside, int face, int ix, int iy) {
  return (static_cast<std::uint32_t>(face) << (2 * lside)) + spread_bits(ix) + (spread_bits(iy) << 1);
}

// Pixel symbol: 0..8 = (dix,diy) in {-1,0,1}^2 on the predicted pixel's
// face, 9 = escape (raw ipix).
inline int pix_symbol(int lside, std::uint32_t pred, std::uint32_t ipix) {
  int f1, x1, y1, f2, x2, y2;
  nest2xyf(lside, pred, f1, x1, y1);
  nest2xyf(lside, ipix, f2, x2, y2);
  if (f1 != f2) return 9;
  const int dx = x2 - x1, dy = y2 - y1;
  if (dx < -1 || dx > 1 || dy < -1 || dy > 1) return 9;
  return (dx + 1) * 3 + (dy + 1);
}
inline std::uint32_t pix_from_symbol(int lside, std::uint32_t pred, int sym) {
  int f, x, y;
  nest2xyf(lside, pred, f, x, y);
  return xyf2nest(lside, f, x + sym / 3 - 1, y + sym % 3 - 1);
}

inline void tangent_basis(const double c[3], double eth[3], double eph[3]) {
  // e_theta = (z cos(phi), z sin(phi), -sin(theta)), e_phi = (-sin(phi), cos(phi), 0)
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
inline double pixel_size(int nside) { return std::sqrt(4 * M_PI / (12.0 * nside * nside)); }

const double kFreqEdges[6] = {13.6, std::sqrt(13.6 * 24.59), 24.59, std::sqrt(24.59 * 54.42), 54.42, 100.0};

// ---------------------------------------------------------------------------
// Quantised ray (the structural encoding, step 4)
// ---------------------------------------------------------------------------

enum Cls { kOrigin = 0, kGenUnmerged = 1, kGenMerged = 2 };
constexpr int kShapeZero = -32768;

struct QRay {
  int cls = 0;
  int src = -1;  // index in the snapshot's source list (or -1)
  std::uint16_t off[3] = {0, 0, 0};
  int rq = 0;  // round(8 r)
  int lside = 2;
  std::uint32_t ipix = 0;
  int dth = 0, dph = 0;  // tangent-plane offset, pixel size / 128
  std::uint32_t ntot_bits = 0;
  int s[5] = {0, 0, 0, 0, 0};  // round(128 log2 f_i), kShapeZero for 0
  std::uint8_t nq[5] = {0, 0, 0, 0, 0};
};

struct Recon {  // decoded float fields (what the RT works with)
  int nside, ipix;
  float r, dr, x[3], dx[3], nu[5], ngam[5];
};

struct Geometry {
  int N;
  std::vector<std::array<double, 3>> src_x;  // source cell centres (0-based coords)
  std::map<long long, int> src_index;
};

Geometry make_geometry(const Snapshot& s) {
  Geometry g;
  g.N = s.N;
  for (std::size_t i = 0; i < s.src_cells.size(); ++i) {
    const long long c = s.src_cells[i];
    // (irt,jrt,krt) - 0.5 in real(4), 1-based -> ijk + 0.5
    g.src_x.push_back({static_cast<double>(static_cast<float>(c % s.N + 1) - 0.5f),
                       static_cast<double>(static_cast<float>((c / s.N) % s.N + 1) - 0.5f),
                       static_cast<double>(static_cast<float>(c / (1LL * s.N * s.N) + 1) - 0.5f)});
    g.src_index[c] = static_cast<int>(i);
  }
  return g;
}

inline void cell_ijk(long long c, int N, int ijk[3]) {
  ijk[0] = static_cast<int>(c % N);
  ijk[1] = static_cast<int>((c / N) % N);
  ijk[2] = static_cast<int>(c / (1LL * N * N));
}

// Vector from the source to point p (min image). Only used for rays with
// r < N/2 - 1, which cannot have wrapped around the periodic box, so the
// min image is the true displacement; farther rays are coded without a
// geometric prediction.
inline bool can_predict(double r, int N) { return r < 0.5 * N - 1; }
inline void image_vector(const double p[3], const std::array<double, 3>& xs, int N, double v[3]) {
  for (int d = 0; d < 3; ++d) {
    v[d] = p[d] - xs[d];
    v[d] -= N * std::floor(v[d] / N + 0.5);
  }
}

inline std::uint32_t predict_pix(int lside, const double v[3]) {
  double w[3] = {v[0], v[1], v[2]};
  if (w[0] == 0 && w[1] == 0 && w[2] == 0) w[2] = 1;
  return fasthp::vec2pix(lside, w);
}

struct QuantOpts {
  int mbits = 23;       // mantissa bits kept in the float32 photon total
  bool origin = true;   // use the origin class
};

QRay quantise(const RayDumpRecord& a, const Geometry& g, const QuantOpts& o) {
  QRay q;
  const int N = g.N;
  int ijk[3];
  cell_ijk(a.cell, N, ijk);
  for (int d = 0; d < 3; ++d) {
    const double off = static_cast<double>(a.x[d]) - ijk[d];
    q.off[d] = static_cast<std::uint16_t>(std::min(std::max(std::floor(off * 65536), 0.0), 65535.0));
  }
  q.lside = bitlen(static_cast<std::uint32_t>(a.nside)) - 1;
  q.ipix = static_cast<std::uint32_t>(a.ipix);
  q.rq = static_cast<int>(std::lround(static_cast<double>(a.r) * 8));
  const bool merged = (a.flags & 1) != 0;
  const bool exact = !merged && !(a.flags & 4) && a.src >= 0 && (q.rq % 8 == 0);
  auto it = g.src_index.find(a.src);
  q.src = (a.src >= 0 && it != g.src_index.end()) ? it->second : -1;
  if (o.origin && exact && q.src >= 0)
    q.cls = kOrigin;
  else
    q.cls = merged ? kGenMerged : kGenUnmerged;
  if (q.cls == kGenMerged) {
    double c[3], eth[3], eph[3];
    fasthp::pix2vec(bitlen(static_cast<std::uint32_t>(a.nside)) - 1, static_cast<std::uint32_t>(a.ipix), c);
    tangent_basis(c, eth, eph);
    const double ps = pixel_size(a.nside);
    double t = 0, p = 0;
    for (int d = 0; d < 3; ++d) {
      t += (a.dx[d] - c[d]) * eth[d];
      p += (a.dx[d] - c[d]) * eph[d];
    }
    q.dth = static_cast<int>(std::max(-127L, std::min(127L, std::lround(t / ps * 128))));
    q.dph = static_cast<int>(std::max(-127L, std::min(127L, std::lround(p / ps * 128))));
  }
  // photons: float32 total (optionally rounded to mbits), log2 shape
  float tot = 0;
  for (int n = 0; n < 5; ++n) tot += a.ngam[n];
  std::uint32_t tb;
  std::memcpy(&tb, &tot, 4);
  if (o.mbits < 23) {
    const std::uint32_t drop = 23 - o.mbits;
    tb = ((tb + (1u << (drop - 1))) >> drop) << drop;  // round to nearest
  }
  q.ntot_bits = tb;
  double dtot = 0;
  for (int n = 0; n < 5; ++n) dtot += a.ngam[n];
  for (int n = 0; n < 5; ++n) {
    const double f = dtot > 0 ? a.ngam[n] / dtot : 0;
    if (f <= 0) {
      q.s[n] = kShapeZero;
    } else {
      const long v = std::lround(std::log2(f) * 128);
      q.s[n] = (v < -4095) ? kShapeZero : static_cast<int>(v);
    }
  }
  for (int n = 0; n < 5; ++n) {
    const double pos = (a.nu[n] - kFreqEdges[n]) / (kFreqEdges[n + 1] - kFreqEdges[n]);
    q.nq[n] = static_cast<std::uint8_t>(std::lround(std::min(std::max(pos, 0.0), 1.0) * 255));
  }
  if (q.cls == kOrigin) {
    for (int d = 0; d < 3; ++d) q.off[d] = 0;  // not stored
  }
  return q;
}

// Fixed-width size of Q in bits (step 4 accounting)
int structural_bits(const QRay& q) {
  // off 3x16, uniq (nside+ipix) 30, merged flag 1, dev 2x8 if merged, r 16,
  // float total 32 (or 9+mbits), shape 5x13, nu 5x8
  return 48 + 30 + 1 + (q.cls == kGenMerged ? 16 : 0) + 16 + 32 + 65 + 40;
}

void reconstruct(const QRay& q, long long cell, const Geometry& g, Recon& out) {
  const int N = g.N;
  int ijk[3];
  cell_ijk(cell, N, ijk);
  const int nside = 1 << q.lside;
  out.nside = nside;
  out.ipix = static_cast<int>(q.ipix);
  const double r = q.rq / 8.0;
  out.r = static_cast<float>(r);
  out.dr = static_cast<float>(rmax_of(nside) - r);
  double c[3];
  fasthp::pix2vec(q.lside, q.ipix, c);
  if (q.cls == kGenMerged && (q.dth || q.dph)) {
    double eth[3], eph[3];
    tangent_basis(c, eth, eph);
    const double ps = pixel_size(nside) / 128;
    double v[3], nrm = 0;
    for (int d = 0; d < 3; ++d) {
      v[d] = c[d] + ps * (q.dth * eth[d] + q.dph * eph[d]);
      nrm += v[d] * v[d];
    }
    nrm = std::sqrt(nrm);
    for (int d = 0; d < 3; ++d) out.dx[d] = static_cast<float>(v[d] / nrm);
  } else {
    for (int d = 0; d < 3; ++d) out.dx[d] = static_cast<float>(c[d]);
  }
  if (q.cls == kOrigin) {
    const auto& xs = g.src_x[q.src];
    for (int d = 0; d < 3; ++d) {
      double x = xs[d] + r * c[d];
      x -= N * std::floor(x / N);
      // keep it inside the chain cell
      const double lo = ijk[d], hi = ijk[d] + 1;
      if (x < lo || x >= hi) {
        // wrap-aware clamp
        double xx = x;
        if (xx - lo > N / 2) xx -= N;
        if (lo - xx > N / 2) xx += N;
        xx = std::min(std::max(xx, lo + 1e-6), hi - 1e-6);
        x = xx;
      }
      out.x[d] = static_cast<float>(x);
    }
  } else {
    for (int d = 0; d < 3; ++d) out.x[d] = static_cast<float>(ijk[d] + (q.off[d] + 0.5) / 65536.0);
  }
  float tot;
  std::memcpy(&tot, &q.ntot_bits, 4);
  double f[5], fs = 0;
  static const std::array<double, 128> kExp2Frac = [] {
    std::array<double, 128> t{};
    for (int i = 0; i < 128; ++i) t[i] = std::exp2(-i / 128.0);
    return t;
  }();
  static const std::array<double, 33> kExp2Int = [] {
    std::array<double, 33> t{};
    for (int i = 0; i < 33; ++i) t[i] = std::exp2(-i);
    return t;
  }();
  for (int n = 0; n < 5; ++n) {
    // 2^(s/128) for -4095 <= s <= 0 from two tables
    f[n] = q.s[n] == kShapeZero ? 0 : kExp2Frac[(-q.s[n]) & 127] * kExp2Int[(-q.s[n]) >> 7];
    fs += f[n];
  }
  for (int n = 0; n < 5; ++n) out.ngam[n] = static_cast<float>(fs > 0 ? tot * f[n] / fs : 0);
  for (int n = 0; n < 5; ++n)
    out.nu[n] = static_cast<float>(kFreqEdges[n] + q.nq[n] / 255.0 * (kFreqEdges[n + 1] - kFreqEdges[n]));
}

// ---------------------------------------------------------------------------
// rANS (ryg_rans style, 32-bit state, byte-wise renormalisation)
// ---------------------------------------------------------------------------

constexpr std::uint32_t kRansL = 1u << 23;
constexpr int kScale = 12;  // probability precision of every model

struct Model {
  int nsym = 0;
  int scale = 12;
  std::vector<std::uint32_t> freq, cum;
  // decode table: slot -> symbol, and per-symbol (freq, cum) packed together
  std::vector<std::uint8_t> lut;
  std::vector<std::uint32_t> fc;  // freq | cum << 16
  void fit(const std::vector<double>& counts, int scale_bits) {
    nsym = static_cast<int>(counts.size());
    scale = scale_bits;
    if (nsym > 256 || scale != kScale) throw std::runtime_error("model: alphabet or scale too large");
    const std::uint32_t M = 1u << scale;
    freq.assign(nsym, 1);
    double tot = 0;
    for (double c : counts) tot += c;
    std::uint32_t used = nsym;
    if (tot > 0 && M > used) {
      const double avail = M - used;
      for (int s = 0; s < nsym; ++s) {
        const std::uint32_t extra = static_cast<std::uint32_t>(counts[s] / tot * avail);
        freq[s] += extra;
        used += extra;
      }
    }
    // hand out the rounding remainder to the most frequent symbol
    int best = static_cast<int>(std::max_element(counts.begin(), counts.end()) - counts.begin());
    if (used > M) throw std::runtime_error("model: alphabet larger than scale");
    freq[best] += M - used;
    cum.assign(nsym + 1, 0);
    for (int s = 0; s < nsym; ++s) cum[s + 1] = cum[s] + freq[s];
    lut.assign(M, 0);
    fc.assign(nsym, 0);
    for (int s = 0; s < nsym; ++s) {
      fc[s] = freq[s] | (cum[s] << 16);
      for (std::uint32_t k = cum[s]; k < cum[s + 1]; ++k) lut[k] = static_cast<std::uint8_t>(s);
    }
  }
  double bits(int s) const { return scale - std::log2(static_cast<double>(freq[s])); }
};

// A context-indexed family of models with counting for training.
struct Family {
  int nsym, nctx, scale;
  std::vector<double> counts;
  std::vector<Model> m;
  Family(int nsym_, int nctx_, int scale_) : nsym(nsym_), nctx(nctx_), scale(scale_),
      counts(static_cast<std::size_t>(nsym_) * nctx_, 0.0) {}
  void count(int ctx, int s) { counts[static_cast<std::size_t>(ctx) * nsym + s] += 1; }
  void fit() {
    m.resize(nctx);
    // pooled marginal as a prior, so sparse contexts are not degenerate
    std::vector<double> pool(nsym, 0.0);
    for (int c = 0; c < nctx; ++c)
      for (int s = 0; s < nsym; ++s) pool[s] += counts[static_cast<std::size_t>(c) * nsym + s];
    double ptot = 0;
    for (double v : pool) ptot += v;
    for (int c = 0; c < nctx; ++c) {
      std::vector<double> cc(counts.begin() + static_cast<std::ptrdiff_t>(c) * nsym,
                             counts.begin() + static_cast<std::ptrdiff_t>(c + 1) * nsym);
      if (ptot > 0)
        for (int s = 0; s < nsym; ++s) cc[s] += 2.0 * pool[s] / ptot;
      m[c].fit(cc, scale);
    }
  }
};

// Two interleaved rANS states sharing one byte stream: symbol k uses state
// k & 1, which halves the serial dependency chain in the decoder.

struct RansEncoder {
  // symbols are pushed forward and encoded in reverse at flush()
  struct Item { std::uint16_t start, freq; };
  std::vector<Item> items;
  void put(const Model& m, int s) {
    items.push_back({static_cast<std::uint16_t>(m.cum[s]), static_cast<std::uint16_t>(m.freq[s])});
  }
  void flush(std::vector<std::uint8_t>& dst) {
    std::uint32_t x[2] = {kRansL, kRansL};
    std::vector<std::uint8_t> rev;
    rev.reserve(items.size());
    for (std::size_t i = items.size(); i-- > 0;) {
      const Item& it = items[i];
      std::uint32_t& xs = x[i & 1];
      const std::uint32_t xmax = ((kRansL >> kScale) << 8) * it.freq;
      while (xs >= xmax) {
        rev.push_back(static_cast<std::uint8_t>(xs & 0xff));
        xs >>= 8;
      }
      xs = ((xs / it.freq) << kScale) + (xs % it.freq) + it.start;
    }
    for (int st = 1; st >= 0; --st)
      for (int k = 0; k < 4; ++k) {
        rev.push_back(static_cast<std::uint8_t>(x[st] & 0xff));
        x[st] >>= 8;
      }
    dst.insert(dst.end(), rev.rbegin(), rev.rend());
    items.clear();
  }
};

struct RansDecoder {
  const std::uint8_t* p;
  std::uint32_t x[2];
  unsigned k = 0;
  explicit RansDecoder(const std::uint8_t* src) : p(src) {
    for (int st = 0; st < 2; ++st) {
      x[st] = (static_cast<std::uint32_t>(p[0]) << 24) | (static_cast<std::uint32_t>(p[1]) << 16) |
              (static_cast<std::uint32_t>(p[2]) << 8) | p[3];
      p += 4;
    }
  }
  inline int get(const Model& m) {
    std::uint32_t& xs = x[k++ & 1];
    const std::uint32_t slot = xs & ((1u << kScale) - 1);
    const int s = m.lut[slot];
    const std::uint32_t fc = m.fc[s];
    xs = (fc & 0xffff) * (xs >> kScale) + slot - (fc >> 16);
    while (xs < kRansL) xs = (xs << 8) | *p++;
    return s;
  }
};

struct BitWriter {
  std::vector<std::uint8_t> buf;
  std::uint64_t acc = 0;
  int n = 0;
  void put(std::uint32_t v, int nb) {
    if (nb == 0) return;
    acc |= static_cast<std::uint64_t>(v & ((nb == 32) ? 0xffffffffu : ((1u << nb) - 1))) << n;
    n += nb;
    while (n >= 8) {
      buf.push_back(static_cast<std::uint8_t>(acc & 0xff));
      acc >>= 8;
      n -= 8;
    }
  }
  void flush() {
    while (n > 0) {
      buf.push_back(static_cast<std::uint8_t>(acc & 0xff));
      acc >>= 8;
      n -= 8;
    }
    n = 0;
    acc = 0;
  }
};

struct BitReader {
  const std::uint8_t* p;
  std::uint64_t acc = 0;
  int n = 0;
  explicit BitReader(const std::uint8_t* src) : p(src) {}
  inline std::uint32_t get(int nb) {
    if (nb == 0) return 0;
    while (n < nb) {
      acc |= static_cast<std::uint64_t>(*p++) << n;
      n += 8;
    }
    const std::uint32_t v = static_cast<std::uint32_t>(acc & ((nb == 32) ? 0xffffffffull : ((1ull << nb) - 1)));
    acc >>= nb;
    n -= nb;
    return v;
  }
};

// ---------------------------------------------------------------------------
// The models and the symmetric field traversal
// ---------------------------------------------------------------------------

enum Group { gCount, gClass, gSrc, gR, gLside, gPix, gOff, gDev, gNtot, gShape, gNu, kNGroups };
const char* kGroupName[kNGroups] = {"count", "class", "src", "r", "lside", "ipix", "offset", "dev", "ntot", "shape", "nu"};

struct Models {
  // All models use 12-bit probabilities and alphabets <= 256, so every
  // decode table is <= 4 KB + 1 KB and the whole set stays in L2.
  Family count{64, 16, 12};       // rays in cell (63 = escape), ctx: prev count bucket
  Family cls{3, 4, 12};           // ctx: class of previous ray in cell (3 = first)
  Family srcsame{2, 2, 12};       // same source as previous in block; ctx O / G
  Family srchint{2, 1, 12};       // G: has a source hint
  Family rdelta{64, 1, 12};       // O: r - round(|v|) + 32, 63 = escape
  Family lsideO{4, 1, 12};        // O: lside(r) - lside, 3 = escape
  Family pixO{10, 12, 12};        // O: pixel symbol, ctx lside-2
  Family lsideG{12, 3, 12};       // G: lside-2, ctx class
  Family rG{20, 24, 12};          // G: bitlen(zigzag(rq - 8 rmin)), ctx (lside-2, class)
  Family pixG{10, 12, 12};        // G with source: pixel symbol from position
  Family coneG{48, 1, 12};        // G without source: nside=2 cone
  Family dth{256, 1, 12};
  Family dph{256, 4, 12};         // ctx |dth| bucket
  Family nexp{64, 12, 12};        // float exponent - base(lside) + 32, 63 = escape
  Family s0{256, 1, 12};          // first shape value / 8 (0 = zero, 255 = escape)
  Family sd{256, 8, 12};          // shape step / 4, ctx (group-1, prev zero)
  Family nudict{256, 1, 12};      // nu tuple dictionary index, 255 = escape
  Family nuraw{256, 5, 12};
  std::map<std::uint64_t, int> nu_index;  // packed nq tuple -> dictionary index
  std::vector<std::uint64_t> nu_tuple;
  int exp_base[12] = {127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127};
  std::vector<Family*> all() {
    return {&count, &cls, &srcsame, &srchint, &rdelta, &lsideO, &pixO, &lsideG, &rG, &pixG,
            &coneG, &dth, &dph, &nexp, &s0, &sd, &nudict, &nuraw};
  }
  void fit() {
    for (Family* f : all()) f->fit();
  }
};

inline std::uint64_t pack_nu(const std::uint8_t nq[5]) {
  std::uint64_t v = 0;
  for (int n = 0; n < 5; ++n) v |= static_cast<std::uint64_t>(nq[n]) << (8 * n);
  return v;
}

// Per-block coding state (reset at block start)
struct BlockState {
  int last_src = -1;
};

// IO policies: Count (training), Enc (encode + bit accounting), Dec.
struct CountIO {
  static constexpr bool decoding = false;
  int sym(Family& f, int ctx, int s, int /*grp*/) { f.count(ctx, s); return s; }
  std::uint32_t raw(std::uint32_t v, int /*nb*/, int /*grp*/) { return v; }
};
struct EncIO {
  static constexpr bool decoding = false;
  RansEncoder* rans;
  BitWriter* bits;
  double* acct;  // bits per group
  int sym(Family& f, int ctx, int s, int grp) {
    const Model& m = f.m[ctx];
    rans->put(m, s);
    acct[grp] += m.bits(s);
    return s;
  }
  std::uint32_t raw(std::uint32_t v, int nb, int grp) {
    bits->put(v, nb);
    acct[grp] += nb;
    return v;
  }
};
struct DecIO {
  static constexpr bool decoding = true;
  RansDecoder* rans;
  BitReader* bits;
  int sym(Family& f, int ctx, int /*s*/, int /*grp*/) { return rans->get(f.m[ctx]); }
  std::uint32_t raw(std::uint32_t /*v*/, int nb, int /*grp*/) { return bits->get(nb); }
};

// ipix without a geometric prediction: nside=2 cone, then raw bits
template <class IO>
void code_cone(IO& io, Models& M, QRay& q) {
  const int sh = 2 * (q.lside - 1);
  const int cone = io.sym(M.coneG, 0, static_cast<int>(q.ipix >> sh), gPix);
  const std::uint32_t low = io.raw(q.ipix & ((1u << sh) - 1), sh, gPix);
  q.ipix = (static_cast<std::uint32_t>(cone) << sh) | low;
}

// Codes one ray (values in q on encode; filled in q on decode). `ctx_rec`
// supplies the cell context (tau), which the decoder also has.
template <class IO>
void code_ray(IO& io, Models& M, QRay& q, const RayDumpRecord& ctx_rec, long long cell,
              int prev_cls, BlockState& bs, const Geometry& g, int nsrc_bits) {
  const int N = g.N;
  q.cls = io.sym(M.cls, prev_cls, q.cls, gClass);
  int ijk[3];
  cell_ijk(cell, N, ijk);
  if (q.cls == kOrigin) {
    // source
    {
      const int same = io.sym(M.srcsame, 0, q.src == bs.last_src ? 1 : 0, gSrc);
      if (same)
        q.src = bs.last_src;
      else
        q.src = static_cast<int>(io.raw(static_cast<std::uint32_t>(q.src), nsrc_bits, gSrc));
      bs.last_src = q.src;
    }
    const auto& xs = g.src_x[q.src];
    const double cc[3] = {ijk[0] + 0.5, ijk[1] + 0.5, ijk[2] + 0.5};
    // r given the min-image distance from the source to the cell centre
    double b[3], d0 = 0;
    for (int d = 0; d < 3; ++d) {
      b[d] = cc[d] - xs[d];
      b[d] -= N * std::floor(b[d] / N + 0.5);
      d0 += b[d] * b[d];
    }
    d0 = std::sqrt(d0);
    {
      const int rint = q.rq / 8;
      const int delta = rint - static_cast<int>(std::lround(d0));
      int s = (delta >= -32 && delta < 31) ? delta + 32 : 63;
      s = io.sym(M.rdelta, 0, s, gR);
      if (s == 63)
        q.rq = 8 * static_cast<int>(io.raw(static_cast<std::uint32_t>(rint), 16, gR));
      else
        q.rq = 8 * (s - 32 + static_cast<int>(std::lround(d0)));
    }
    const double r = q.rq / 8.0;
    {
      const int lr = lside_of_r(r);
      int s = lr - q.lside;
      s = (s >= 0 && s < 3) ? s : 3;
      s = io.sym(M.lsideO, 0, s, gLside);
      if (s == 3)
        q.lside = static_cast<int>(io.raw(static_cast<std::uint32_t>(q.lside), 4, gLside));
      else
        q.lside = lr - s;
    }
    if (can_predict(r, N)) {
      double v[3];
      image_vector(cc, xs, N, v);
      const std::uint32_t pred = predict_pix(q.lside, v);
      int s = IO::decoding ? 0 : pix_symbol(q.lside, pred, q.ipix);
      s = io.sym(M.pixO, q.lside - 2, s, gPix);
      if (s == 9)
        q.ipix = io.raw(q.ipix, 2 * q.lside + 4, gPix);
      else
        q.ipix = pix_from_symbol(q.lside, pred, s);
    } else {
      code_cone(io, M, q);
    }
    q.dth = q.dph = 0;
  } else {
    // offsets
    // (nearly uniform: stored raw)
    for (int d = 0; d < 3; ++d) q.off[d] = static_cast<std::uint16_t>(io.raw(q.off[d], 16, gOff));
    // nside, r
    q.lside = 2 + io.sym(M.lsideG, q.cls, q.lside - 2, gLside);
    {
      const int base = static_cast<int>(std::lround(8 * rmin_of(1 << q.lside)));
      const std::uint32_t z = zigzag(q.rq - base);
      const int bl = io.sym(M.rG, (q.lside - 2) * 2 + (q.cls == kGenMerged), bitlen(z), gR);
      std::uint32_t zz = 0;
      if (bl > 0) zz = (1u << (bl - 1)) | io.raw(z & ((1u << (bl - 1)) - 1), bl - 1, gR);
      q.rq = base + unzigzag(zz);
    }
    // source hint and ipix
    {
      const int has = io.sym(M.srchint, 0, q.src >= 0 ? 1 : 0, gSrc);
      if (has) {
        const int same = io.sym(M.srcsame, 1, q.src == bs.last_src ? 1 : 0, gSrc);
        if (same)
          q.src = bs.last_src;
        else
          q.src = static_cast<int>(io.raw(static_cast<std::uint32_t>(q.src), nsrc_bits, gSrc));
        bs.last_src = q.src;
        if (can_predict(q.rq / 8.0, N)) {
          const double p[3] = {ijk[0] + (q.off[0] + 0.5) / 65536.0, ijk[1] + (q.off[1] + 0.5) / 65536.0,
                               ijk[2] + (q.off[2] + 0.5) / 65536.0};
          double v[3];
          image_vector(p, g.src_x[q.src], N, v);
          const std::uint32_t pred = predict_pix(q.lside, v);
          int s = IO::decoding ? 0 : pix_symbol(q.lside, pred, q.ipix);
          s = io.sym(M.pixG, q.lside - 2, s, gPix);
          if (s == 9)
            q.ipix = io.raw(q.ipix, 2 * q.lside + 4, gPix);
          else
            q.ipix = pix_from_symbol(q.lside, pred, s);
        } else {
          code_cone(io, M, q);
        }
      } else {
        q.src = -1;
        code_cone(io, M, q);
      }
    }
    if (q.cls == kGenMerged) {
      q.dth = io.sym(M.dth, 0, q.dth + 128, gDev) - 128;
      const int b = std::min(std::abs(q.dth) / 16, 3);
      q.dph = io.sym(M.dph, b, q.dph + 128, gDev) - 128;
    } else {
      q.dth = q.dph = 0;
    }
  }
  // photons: float total (exponent modelled around a per-lside base)
  {
    const int base = M.exp_base[q.lside - 2];
    const int e = static_cast<int>((q.ntot_bits >> 23) & 0xff);
    int s = (e - base >= -32 && e - base < 31) ? e - base + 32 : 63;
    s = io.sym(M.nexp, q.lside - 2, s, gNtot);
    const int ee = (s == 63) ? static_cast<int>(io.raw(static_cast<std::uint32_t>(e), 8, gNtot)) : s - 32 + base;
    const std::uint32_t man = io.raw(q.ntot_bits & 0x7fffff, 23, gNtot);
    q.ntot_bits = (static_cast<std::uint32_t>(ee) << 23) | man;
  }
  // shape: s0, then steps between adjacent groups
  {
    // absolute value v = -s in [0, 4095] or zero: 0 = zero, 1..254 = v/8 + 1
    // (+3 raw bits), 255 = escape (12 raw bits)
    auto code_abs = [&](Family& f, int ctx, int v) {
      int s = (v == kShapeZero) ? 0 : ((-v) / 8 < 254 ? (-v) / 8 + 1 : 255);
      s = io.sym(f, ctx, s, gShape);
      if (s == 0) return kShapeZero;
      if (s == 255) return -static_cast<int>(io.raw(static_cast<std::uint32_t>(-v), 12, gShape));
      return -static_cast<int>(8 * (s - 1) + io.raw(static_cast<std::uint32_t>((-v) & 7), 3, gShape));
    };
    q.s[0] = code_abs(M.s0, 0, q.s[0]);
    for (int n = 1; n < 5; ++n) {
      const bool pz = q.s[n - 1] == kShapeZero;
      const int ctx = (n - 1) * 2 + (pz ? 1 : 0);
      if (pz) {
        q.s[n] = code_abs(M.sd, ctx, q.s[n]);
      } else {
        // step d = s_n - s_{n-1}: 0 = zero, 1..253 = floor(d/4) + 127 (+2 raw
        // bits), 255 = escape (12 raw bits of -s_n)
        const int v = q.s[n];
        int s = 0;
        int dlt = 0;
        if (v != kShapeZero) {
          dlt = v - q.s[n - 1];
          const int c = (dlt >= 0) ? dlt / 4 : -((-dlt + 3) / 4);  // floor(dlt/4)
          s = (c >= -126 && c <= 126) ? c + 127 : 255;
        }
        s = io.sym(M.sd, ctx, s, gShape);
        if (s == 0)
          q.s[n] = kShapeZero;
        else if (s == 255)
          q.s[n] = -static_cast<int>(io.raw(static_cast<std::uint32_t>(-v), 12, gShape));
        else
          q.s[n] = q.s[n - 1] + 4 * (s - 127) + static_cast<int>(io.raw(static_cast<std::uint32_t>(dlt & 3), 2, gShape));
      }
    }
  }
  // nu
  {
    int idx = 255;
    if (!IO::decoding) {
      auto it = M.nu_index.find(pack_nu(q.nq));
      if (it != M.nu_index.end()) idx = it->second;
    }
    idx = io.sym(M.nudict, 0, idx, gNu);
    if (idx == 255) {
      for (int n = 0; n < 5; ++n) q.nq[n] = static_cast<std::uint8_t>(io.sym(M.nuraw, n, q.nq[n], gNu));
    } else if (IO::decoding) {
      const std::uint64_t t = M.nu_tuple[idx];
      for (int n = 0; n < 5; ++n) q.nq[n] = static_cast<std::uint8_t>((t >> (8 * n)) & 0xff);
    }
  }
}

// Block layout: per block, [u32 rans_bytes][u32 raw_bytes][rans][raw]
struct Encoded {
  std::vector<std::uint8_t> data;
  std::vector<std::size_t> block_offset;
  std::vector<std::pair<long long, long long>> block_cells;  // [c0, c1)
  std::vector<std::size_t> block_ray0;
};

// Cells: the snapshot lists rays grouped by cell in increasing cell order.
struct CellIndex {
  std::vector<long long> cells;     // cells with rays, increasing
  std::vector<std::size_t> start;   // ray index range per listed cell
};
CellIndex index_cells(const std::vector<RayDumpRecord>& r) {
  CellIndex ci;
  for (std::size_t i = 0; i < r.size(); ++i) {
    if (ci.cells.empty() || r[i].cell != ci.cells.back()) {
      ci.cells.push_back(r[i].cell);
      ci.start.push_back(i);
    }
  }
  ci.start.push_back(r.size());
  return ci;
}

template <class IO>
void code_count(IO& io, Models& M, int& prevc, int& n) {
  const int ctx = std::min(prevc, 15);
  int s = io.sym(M.count, ctx, n < 63 ? n : 63, gCount);
  if (s == 63) s = 63 + static_cast<int>(io.raw(static_cast<std::uint32_t>(n - 63), 24, gCount));
  n = s;
  prevc = n;
}

int nsrc_bits_of(const Geometry& g) { return std::max(1, bitlen(static_cast<std::uint32_t>(g.src_x.size()))); }

void train(Models& M, const Snapshot& s, const QuantOpts& o) {
  const Geometry g = make_geometry(s);
  CountIO io;
  const CellIndex ci = index_cells(s.rec);
  const long long ncell = 1LL * s.N * s.N * s.N;
  BlockState bs;
  int prevc = 0;
  std::size_t k = 0;
  std::map<std::uint64_t, double> nu_counts;
  std::vector<QRay> qs(s.rec.size());
  for (std::size_t i = 0; i < s.rec.size(); ++i) qs[i] = quantise(s.rec[i], g, o);
  for (const QRay& q : qs) nu_counts[pack_nu(q.nq)] += 1;
  {
    std::vector<std::array<int, 256>> h(12);
    for (auto& a : h) a.fill(0);
    for (const QRay& q : qs) ++h[q.lside - 2][(q.ntot_bits >> 23) & 0xff];
    for (int l = 0; l < 12; ++l) {
      const int m = static_cast<int>(std::max_element(h[l].begin(), h[l].end()) - h[l].begin());
      if (h[l][m] > 0) M.exp_base[l] = m;
    }
  }
  {
    std::vector<std::pair<double, std::uint64_t>> v;
    for (auto& kv : nu_counts) v.push_back({kv.second, kv.first});
    std::sort(v.rbegin(), v.rend());
    for (std::size_t i = 0; i < v.size() && i < 255; ++i) {
      M.nu_index[v[i].second] = static_cast<int>(i);
      M.nu_tuple.push_back(v[i].second);
    }
  }
  for (long long c = 0; c < ncell; ++c) {
    int n = 0;
    std::size_t i0 = 0;
    if (k < ci.cells.size() && ci.cells[k] == c) {
      i0 = ci.start[k];
      n = static_cast<int>(ci.start[k + 1] - i0);
      ++k;
    }
    code_count(io, M, prevc, n);
    int prev_cls = 3;
    for (int j = 0; j < n; ++j) {
      QRay q = qs[i0 + j];
      code_ray(io, M, q, s.rec[i0 + j], c, prev_cls, bs, g, nsrc_bits_of(g));
      prev_cls = q.cls;
    }
  }
  M.fit();
}

struct EncodeResult {
  Encoded enc;
  double acct[kNGroups] = {};
  double seconds = 0;
};

EncodeResult encode(Models& M, const Snapshot& s, const std::vector<QRay>& qs, std::size_t block_rays) {
  EncodeResult R;
  const Geometry g = make_geometry(s);
  const CellIndex ci = index_cells(s.rec);
  const long long ncell = 1LL * s.N * s.N * s.N;
  const auto t0 = std::chrono::steady_clock::now();
  RansEncoder re;
  BitWriter bw;
  EncIO io{&re, &bw, R.acct};
  BlockState bs;
  int prevc = 0;
  std::size_t k = 0, nblock = 0;
  long long c0 = 0;
  std::size_t ray0 = 0;
  auto close_block = [&](long long c1, std::size_t ray1) {
    std::vector<std::uint8_t> rb;
    re.flush(rb);
    bw.flush();
    R.enc.block_offset.push_back(R.enc.data.size());
    R.enc.block_cells.push_back({c0, c1});
    R.enc.block_ray0.push_back(ray0);
    const std::uint32_t h[2] = {static_cast<std::uint32_t>(rb.size()), static_cast<std::uint32_t>(bw.buf.size())};
    const std::uint8_t* hp = reinterpret_cast<const std::uint8_t*>(h);
    R.enc.data.insert(R.enc.data.end(), hp, hp + 8);
    R.enc.data.insert(R.enc.data.end(), rb.begin(), rb.end());
    R.enc.data.insert(R.enc.data.end(), bw.buf.begin(), bw.buf.end());
    bw.buf.clear();
    c0 = c1;
    ray0 = ray1;
    bs = BlockState{};
    prevc = 0;
  };
  for (long long c = 0; c < ncell; ++c) {
    int n = 0;
    std::size_t i0 = 0;
    if (k < ci.cells.size() && ci.cells[k] == c) {
      i0 = ci.start[k];
      n = static_cast<int>(ci.start[k + 1] - i0);
      ++k;
    }
    code_count(io, M, prevc, n);
    int prev_cls = 3;
    for (int j = 0; j < n; ++j) {
      QRay q = qs[i0 + j];
      code_ray(io, M, q, s.rec[i0 + j], c, prev_cls, bs, g, nsrc_bits_of(g));
      prev_cls = q.cls;
    }
    nblock += n;
    if (nblock >= block_rays) {
      close_block(c + 1, i0 + n);
      nblock = 0;
    }
  }
  close_block(ncell, s.rec.size());
  R.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return R;
}

// Decodes every block; returns seconds. Rays come out in snapshot order.
double decode(Models& M, const Snapshot& s, const Encoded& E, std::vector<QRay>& qout,
              std::vector<Recon>& rout) {
  const Geometry g = make_geometry(s);
  qout.assign(s.rec.size(), QRay{});
  rout.assign(s.rec.size(), Recon{});
  const auto t0 = std::chrono::steady_clock::now();
  for (std::size_t b = 0; b < E.block_offset.size(); ++b) {
    const std::uint8_t* p = E.data.data() + E.block_offset[b];
    std::uint32_t h[2];
    std::memcpy(h, p, 8);
    RansDecoder rd(p + 8);
    BitReader br(p + 8 + h[0]);
    DecIO io{&rd, &br};
    BlockState bs;
    int prevc = 0;
    std::size_t ir = E.block_ray0[b];
    for (long long c = E.block_cells[b].first; c < E.block_cells[b].second; ++c) {
      int n = 0;
      code_count(io, M, prevc, n);
      int prev_cls = 3;
      for (int j = 0; j < n; ++j) {
        QRay& q = qout[ir];
        // cell context (tau) comes from the grid; take it from the record
        code_ray(io, M, q, s.rec[ir], c, prev_cls, bs, g, nsrc_bits_of(g));
        prev_cls = q.cls;
        reconstruct(q, c, g, rout[ir]);
        ++ir;
      }
    }
  }
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

struct Err {
  double x_max = 0, x_rms = 0, ang_max = 0, ang_rms = 0, r_max = 0, dr_max = 0;
  double ng_rel_max = 0, ng_rel_rms = 0, ntot_rel_max = 0, nu_max = 0;
  double ng_tot_rel = 0;  // relative change in the total photon count
  long long n = 0, nside_mismatch = 0, ipix_mismatch = 0;
};

Err compare(const Snapshot& s, const std::vector<Recon>& rr) {
  Err e;
  const int N = s.N;
  double sum0 = 0, sum1 = 0, ngw = 0;
  for (std::size_t i = 0; i < s.rec.size(); ++i) {
    const RayDumpRecord& a = s.rec[i];
    const Recon& b = rr[i];
    double dxm = 0, dot = 0;
    for (int d = 0; d < 3; ++d) {
      double dd = std::fabs(static_cast<double>(a.x[d]) - b.x[d]);
      dd = std::min(dd, N - dd);
      dxm = std::max(dxm, dd);
      dot += static_cast<double>(a.dx[d]) * b.dx[d];
    }
    double na = 0, nb = 0;
    for (int d = 0; d < 3; ++d) {
      na += static_cast<double>(a.dx[d]) * a.dx[d];
      nb += static_cast<double>(b.dx[d]) * b.dx[d];
    }
    const double ang = std::acos(std::min(1.0, dot / std::sqrt(na * nb)));
    e.x_max = std::max(e.x_max, dxm);
    e.x_rms += dxm * dxm;
    e.ang_max = std::max(e.ang_max, ang);
    e.ang_rms += ang * ang;
    e.r_max = std::max(e.r_max, std::fabs(static_cast<double>(a.r) - b.r));
    e.dr_max = std::max(e.dr_max, std::fabs(static_cast<double>(a.dr) - b.dr));
    if (a.nside != b.nside) ++e.nside_mismatch;
    if (a.ipix != b.ipix) ++e.ipix_mismatch;
    double ta = 0, tb = 0;
    for (int n = 0; n < 5; ++n) {
      ta += a.ngam[n];
      tb += b.ngam[n];
    }
    for (int n = 0; n < 5; ++n) {
      // relative to the ray's total (groups far below the total do not matter)
      const double rel = std::fabs(static_cast<double>(a.ngam[n]) - b.ngam[n]) / std::max(ta, 1e-300);
      e.ng_rel_max = std::max(e.ng_rel_max, rel);
      e.ng_rel_rms += rel * rel * ta;
      e.nu_max = std::max(e.nu_max, std::fabs(static_cast<double>(a.nu[n]) - b.nu[n]));
    }
    ngw += ta;
    e.ntot_rel_max = std::max(e.ntot_rel_max, std::fabs(ta - tb) / std::max(ta, 1e-300));
    sum0 += ta;
    sum1 += tb;
    ++e.n;
  }
  e.x_rms = std::sqrt(e.x_rms / std::max<long long>(e.n, 1));
  e.ang_rms = std::sqrt(e.ang_rms / std::max<long long>(e.n, 1));
  e.ng_rel_rms = std::sqrt(e.ng_rel_rms / std::max(ngw, 1e-300));
  e.ng_tot_rel = (sum1 - sum0) / std::max(sum0, 1e-300);
  return e;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s TRAIN_STEM TEST_STEM [--mbits m] [--block n] [--reps k] [--no-origin]\n",
                 argv[0]);
    return 2;
  }
  QuantOpts o;
  std::size_t block = 4096;
  int reps = 3;
  for (int a = 3; a < argc; ++a) {
    const std::string s = argv[a];
    if (s == "--mbits" && a + 1 < argc)
      o.mbits = std::atoi(argv[++a]);
    else if (s == "--block" && a + 1 < argc)
      block = static_cast<std::size_t>(std::atoll(argv[++a]));
    else if (s == "--reps" && a + 1 < argc)
      reps = std::atoi(argv[++a]);
    else if (s == "--no-origin")
      o.origin = false;
    else {
      std::fprintf(stderr, "unknown argument %s\n", s.c_str());
      return 2;
    }
  }
  try {
    const Snapshot tr = load(argv[1]);
    const Snapshot te = load(argv[2]);
    Models M;
    train(M, tr, o);
    const Geometry g = make_geometry(te);
    const std::size_t n = te.rec.size();
    if (n == 0) throw std::runtime_error("empty snapshot");

    // Structural (step 4): quantise + fixed-width accounting + error
    std::vector<QRay> qs(n);
    QuantOpts os = o;
    os.origin = false;
    double sbits = 0;
    std::vector<Recon> rs(n);
    long long ncls[3] = {0, 0, 0};
    for (std::size_t i = 0; i < n; ++i) {
      const QRay q = quantise(te.rec[i], g, os);
      sbits += structural_bits(q);
      reconstruct(q, te.rec[i].cell, g, rs[i]);
    }
    const Err es = compare(te, rs);
    for (std::size_t i = 0; i < n; ++i) {
      qs[i] = quantise(te.rec[i], g, o);
      ++ncls[qs[i].cls];
    }

    // Entropy coding (step 5)
    EncodeResult best;
    double tenc = 1e300, tdec = 1e300;
    std::vector<QRay> qd;
    std::vector<Recon> rd;
    for (int k = 0; k < reps; ++k) {
      EncodeResult R = encode(M, te, qs, block);
      tenc = std::min(tenc, R.seconds);
      if (k == 0) best = R;
    }
    for (int k = 0; k < reps; ++k) tdec = std::min(tdec, decode(M, te, best.enc, qd, rd));
    // lossless w.r.t. Q?
    long long qmis = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const QRay& a = qs[i];
      const QRay& b = qd[i];
      bool same = a.cls == b.cls && a.rq == b.rq && a.lside == b.lside && a.ipix == b.ipix &&
                  a.dth == b.dth && a.dph == b.dph && a.ntot_bits == b.ntot_bits;
      for (int d = 0; d < 3 && same; ++d) same = a.off[d] == b.off[d];
      for (int m = 0; m < 5 && same; ++m) same = a.s[m] == b.s[m] && a.nq[m] == b.nq[m];
      if (!same) ++qmis;
    }
    const Err ee = compare(te, rd);
    // time the reconstruction alone (to separate entropy decode from pix2vec etc.)
    double trec = 1e300;
    for (int k = 0; k < reps; ++k) {
      const auto t0 = std::chrono::steady_clock::now();
      for (std::size_t i = 0; i < n; ++i) reconstruct(qd[i], te.rec[i].cell, g, rd[i]);
      trec = std::min(trec, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }

    std::printf("{\n  \"train\": \"%s\", \"test\": \"%s\", \"nrays\": %zu, \"N\": %d,\n", argv[1], argv[2], n, te.N);
    std::printf("  \"mbits\": %d, \"block\": %zu, \"origin\": %d,\n", o.mbits, block, o.origin ? 1 : 0);
    std::printf("  \"class_frac\": [%.4f, %.4f, %.4f],\n", ncls[0] / double(n), ncls[1] / double(n),
                ncls[2] / double(n));
    std::printf("  \"bytes_raw\": 80, \"bytes_structural\": %.3f,\n", sbits / 8 / n);
    std::printf("  \"bytes_entropy\": %.3f, \"bytes_entropy_model\": %.3f,\n", best.enc.data.size() / double(n),
                [&] { double t = 0; for (double v : best.acct) t += v; return t / 8 / n; }());
    std::printf("  \"bits_by_group\": {");
    for (int k = 0; k < kNGroups; ++k)
      std::printf("%s\"%s\": %.3f", k ? ", " : "", kGroupName[k], best.acct[k] / n);
    std::printf("},\n");
    std::printf("  \"nblocks\": %zu,\n", best.enc.block_offset.size());
    std::printf("  \"encode_ns_per_ray\": %.1f, \"decode_ns_per_ray\": %.1f, \"reconstruct_ns_per_ray\": %.1f,\n",
                tenc / n * 1e9, tdec / n * 1e9, trec / n * 1e9);
    std::printf("  \"q_mismatch\": %lld,\n", qmis);
    auto perr = [](const char* name, const Err& e, bool last) {
      std::printf("  \"%s\": {\"x_max\": %.3g, \"x_rms\": %.3g, \"ang_max\": %.3g, \"ang_rms\": %.3g, "
                  "\"r_max\": %.3g, \"dr_max\": %.3g, \"ng_rel_max\": %.3g, \"ng_rel_rms\": %.3g, "
                  "\"ntot_rel_max\": %.3g, \"nu_max\": %.3g, \"ng_tot_rel\": %.3g, \"nside_mismatch\": %lld, "
                  "\"ipix_mismatch\": %lld}%s\n",
                  name, e.x_max, e.x_rms, e.ang_max, e.ang_rms, e.r_max, e.dr_max, e.ng_rel_max,
                  e.ng_rel_rms, e.ntot_rel_max, e.nu_max, e.ng_tot_rel, e.nside_mismatch,
                  e.ipix_mismatch, last ? "" : ",");
    };
    perr("err_structural", es, false);
    perr("err_entropy", ee, true);
    std::printf("}\n");
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ray_codec: %s\n", e.what());
    return 1;
  }
  return 0;
}
