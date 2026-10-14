// Alternative OTVET Eddington tensors (otvet_closure.hpp,
// docs/otvet_closure.md). The numpy prototypes are
// tools/analysis/tensor_compare.py (split, screened kernel) and
// tools/analysis/cone_tensor.py (cones, budget cutoff).
#include "radray/otvet_closure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "radray/config.hpp"
#include "radray/parallel.hpp"
#include "radray/rt_moments.hpp"

namespace radray {

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

bool set_otvet_closure_param(Config& c, const std::string& key, const std::string& val) {
  auto pos_int = [&](int lo) {
    const int v = std::stoi(val);
    if (v < lo) throw std::runtime_error(key + " must be >= " + std::to_string(lo));
    return v;
  };
  auto pos_dbl = [&](double lo, bool strict) {
    const double v = std::stod(val);
    if (strict ? !(v > lo) : !(v >= lo))
      throw std::runtime_error(key + (strict ? " must be > " : " must be >= ") + std::to_string(lo));
    return v;
  };
  if (key == "otvet_closure") {
    if (val == "thin" || val == "0")
      c.otvet_closure = 0;
    else if (val == "screen" || val == "1")
      c.otvet_closure = 1;
    else if (val == "cone" || val == "2")
      c.otvet_closure = 2;
    else
      throw std::runtime_error("otvet_closure must be thin, screen or cone");
    return true;
  }
  if (key == "otvet_split") {
    c.otvet_split = std::stoi(val);
    if (c.otvet_split != 0 && c.otvet_split != 1) throw std::runtime_error("otvet_split must be 0 or 1");
    return true;
  }
  if (key == "otvet_screen_lambda") {
    c.otvet_screen_lambda = pos_dbl(0, false);
    return true;
  }
  if (key == "otvet_cone_nside") {
    c.otvet_cone_nside = pos_int(1);
    return true;
  }
  if (key == "otvet_cone_coarse") {
    c.otvet_cone_coarse = pos_int(1);
    return true;
  }
  if (key == "otvet_cone_rc") {
    c.otvet_cone_rc = pos_dbl(1, false);
    return true;
  }
  if (key == "otvet_cone_taucap") {
    c.otvet_cone_taucap = pos_dbl(0, true);
    return true;
  }
  if (key == "otvet_cone_smax") {
    c.otvet_cone_smax = pos_dbl(0, true);
    return true;
  }
  if (key == "otvet_cone_budget") {
    c.otvet_cone_budget = std::stoi(val);
    if (c.otvet_cone_budget != 0 && c.otvet_cone_budget != 1)
      throw std::runtime_error("otvet_cone_budget must be 0 or 1");
    return true;
  }
  if (key == "otvet_cone_blend") {
    c.otvet_cone_blend = pos_dbl(0, false);
    return true;
  }
  return false;
}

OtvetClosureParams otvet_closure_params(const Config& c) {
  OtvetClosureParams p;
  p.closure = c.otvet_closure;
  p.split = c.otvet_split;
  p.screen_lambda = c.otvet_screen_lambda;
  p.cone_nside = c.otvet_cone_nside;
  p.cone_coarse = c.otvet_cone_coarse;
  p.cone_rc = c.otvet_cone_rc;
  p.cone_taucap = c.otvet_cone_taucap;
  p.cone_smax = c.otvet_cone_smax;
  p.cone_budget = c.otvet_cone_budget;
  p.cone_blend = c.otvet_cone_blend;
  return p;
}

void check_otvet_closure_param(const Config& c) {
  const OtvetClosureParams p = otvet_closure_params(c);
  if (p.active() && c.rt_method != RTMethod::OTVET)
    throw std::runtime_error("otvet_closure / otvet_split require rt_method = otvet");
  if (p.closure == 2 && p.split)
    throw std::runtime_error(
        "otvet_split applies to otvet_closure = thin or screen; the cone uses otvet_cone_blend");
  if (p.closure == 2 && c.Nrt % p.cone_coarse != 0)
    throw std::runtime_error("otvet_cone_coarse must divide the RT grid size");
}

std::string otvet_closure_summary(const Config& c) {
  const OtvetClosureParams p = otvet_closure_params(c);
  if (!p.active()) return "";
  std::ostringstream o;
  o << "otvet_closure = " << (p.closure == 0 ? "thin" : p.closure == 1 ? "screen" : "cone");
  if (p.split) o << ", otvet_split = 1";
  if (p.closure == 1) {
    if (p.screen_lambda > 0)
      o << ", lambda = " << p.screen_lambda << " cells";
    else
      o << ", lambda from the photon budget";
  }
  if (p.closure == 2)
    o << " (nside " << p.cone_nside << ", coarse " << p.cone_coarse << ", r_c " << p.cone_rc
      << ", tau cap " << p.cone_taucap << ", " << p.cone_smax << " boxes, budget "
      << p.cone_budget << ", blend " << p.cone_blend << ")";
  return o.str();
}

// ---------------------------------------------------------------------------
// Pieces
// ---------------------------------------------------------------------------

std::vector<std::array<double, 3>> otvet_healpix_dirs(int nside) {
  // as tools/analysis/cone_tensor.py healpix_dirs (HEALPix ring centres)
  std::vector<std::array<double, 3>> out;
  const double ns = nside;
  for (int i = 1; i < 4 * nside; ++i) {
    double z, shift;
    int nr;
    if (i < nside) {
      z = 1 - i * static_cast<double>(i) / (3.0 * ns * ns);
      nr = 4 * i;
      shift = 0.5;
    } else if (i <= 3 * nside) {
      z = 4.0 / 3 - 2.0 * i / (3.0 * ns);
      nr = 4 * nside;
      shift = 0.5 * ((i - nside + 1) % 2);
    } else {
      const int ii = 4 * nside - i;
      z = -(1 - ii * static_cast<double>(ii) / (3.0 * ns * ns));
      nr = 4 * ii;
      shift = 0.5;
    }
    const double st = std::sqrt(std::max(0.0, 1 - z * z));
    for (int j = 0; j < nr; ++j) {
      const double phi = 2 * M_PI * (j + shift) / nr;
      out.push_back({st * std::cos(phi), st * std::sin(phi), z});
    }
  }
  return out;
}

double otvet_min_eigenvalue(const double h[6]) {
  const double a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5];
  const double p1 = (d * d + e * e) + f * f;
  if (p1 == 0) return std::min(a, std::min(b, c));
  const double q = ((a + b) + c) / 3;
  const double p2 = ((a - q) * (a - q) + (b - q) * (b - q)) + (c - q) * (c - q) + 2 * p1;
  const double p = std::sqrt(p2 / 6);
  const double A = (a - q) / p, B = (b - q) / p, C = (c - q) / p, D = d / p, E = e / p, F = f / p;
  const double det = A * (B * C - F * F) - D * (D * C - F * E) + E * (D * F - B * E);
  const double r = std::max(-1.0, std::min(1.0, det / 2));
  const double phi = std::acos(r) / 3;
  return q + 2 * p * std::cos(phi + 2 * M_PI / 3);
}

namespace {

inline double smoothstep(double r, double r0, double r1) {
  const double t = std::max(0.0, std::min(1.0, (r - r0) / (r1 - r0)));
  return t * t * (3 - 2 * t);
}

inline int mi(int d, int N) { return 2 * d > N ? d - N : d; }

inline std::int64_t idx(int i, int j, int k, int N) {
  return i + static_cast<std::int64_t>(N) * (j + static_cast<std::int64_t>(N) * k);
}

constexpr int kComp[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};

// Screened kernel over the 27 images of the minimum-image separation, plus
// the isotropic continuum of the images beyond (tensor_compare.py
// thin_kernel(images=True, lam)); the source cell is an isotropic 4/3 I.
struct ScreenCtx {
  double lam;
};
void screen_fill(int N, int m, const void* ctx, float* K) {
  const double lam = static_cast<const ScreenCtx*>(ctx)->lam;
  const double R = 3.0 * N * std::cbrt(3 / (4 * M_PI));
  const double bg = 4 * M_PI * lam * std::exp(-R / lam) / (static_cast<double>(N) * N * N);
  const int p = kComp[m][0], q = kComp[m][1];
#pragma omp parallel for schedule(static) if (par::worth(static_cast<long long>(N) * N * N))
  for (int k = 0; k < N; ++k)
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i) {
        const int d0[3] = {mi(i, N), mi(j, N), mi(k, N)};
        double v = 0;
        for (int sc = -1; sc <= 1; ++sc)
          for (int sb = -1; sb <= 1; ++sb)
            for (int sa = -1; sa <= 1; ++sa) {
              const double d[3] = {static_cast<double>(d0[0] + sa * N),
                                   static_cast<double>(d0[1] + sb * N),
                                   static_cast<double>(d0[2] + sc * N)};
              const double r2 = (d[0] * d[0] + d[1] * d[1]) + d[2] * d[2];
              if (r2 == 0) continue;
              v += std::exp(-std::sqrt(r2) / lam) / (r2 * r2) * d[p] * d[q];
            }
        if (p == q) {
          v += bg / 3;
          if (d0[0] == 0 && d0[1] == 0 && d0[2] == 0) v += 4.0 / 3;
        }
        K[idx(i, j, k, N)] = static_cast<float>(v);
      }
}

// Thin kernel (otvet_kernel) inside r < 2 r_c, weight 1 - smoothstep(r, r_c,
// 2 r_c) (the near field of the cones; cone_tensor.py Kn)
struct NearCtx {
  double rc;
};
void near_fill(int N, int m, const void* ctx, float* K) {
  const double rc = static_cast<const NearCtx*>(ctx)->rc;
#pragma omp parallel for schedule(static) if (par::worth(static_cast<long long>(N) * N * N))
  for (int k = 0; k < N; ++k)
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i) {
        const int dx = mi(i, N), dy = mi(j, N), dz = mi(k, N);
        const double r = std::sqrt((static_cast<double>(dx) * dx + static_cast<double>(dy) * dy) +
                                   static_cast<double>(dz) * dz);
        const double w = 1 - smoothstep(r, rc, 2 * rc);
        double Kd[6] = {0, 0, 0, 0, 0, 0};
        if (w > 0) otvet_kernel(N, dx, dy, dz, Kd);
        K[idx(i, j, k, N)] = static_cast<float>(w * Kd[m]);
      }
}

// Mip-map level: n^3 means of the level below
struct Level {
  int n = 0;
  std::vector<float> a;
};

std::vector<Level> pyramid(const std::vector<float>& fine, int N) {
  std::vector<Level> L;
  L.push_back({N, fine});
  while (L.back().n % 2 == 0 && L.back().n >= 2) {
    const Level& b = L.back();
    Level c;
    c.n = b.n / 2;
    c.a.assign(static_cast<std::size_t>(c.n) * c.n * c.n, 0.0f);
#pragma omp parallel for schedule(static) if (par::worth(static_cast<long long>(c.n) * c.n * c.n))
    for (int k = 0; k < c.n; ++k)
      for (int j = 0; j < c.n; ++j)
        for (int i = 0; i < c.n; ++i) {
          double s = 0;
          for (int dk = 0; dk < 2; ++dk)
            for (int dj = 0; dj < 2; ++dj)
              for (int di = 0; di < 2; ++di)
                s += b.a[idx(2 * i + di, 2 * j + dj, 2 * k + dk, b.n)];
          c.a[idx(i, j, k, c.n)] = static_cast<float>(s / 8);
        }
    L.push_back(std::move(c));
  }
  return L;
}

inline int wrapi(long long i, int n) {
  long long r = i % n;
  if (r < 0) r += n;
  return static_cast<int>(r);
}

// Trilinear periodic sample of level l (cell centres at (i + 0.5) 2^l in
// fine-cell units; scipy map_coordinates order 1, grid-wrap)
inline double sample(const Level& L, double f, const double pos[3]) {
  const int n = L.n;
  double t[3];
  int i0[3], i1[3];
  for (int d = 0; d < 3; ++d) {
    const double c = pos[d] / f - 0.5;
    const double fl = std::floor(c);
    t[d] = c - fl;
    i0[d] = wrapi(static_cast<long long>(fl), n);
    i1[d] = i0[d] + 1 == n ? 0 : i0[d] + 1;
  }
  const float* a = L.a.data();
  auto at = [&](int i, int j, int k) { return static_cast<double>(a[idx(i, j, k, n)]); };
  const double c00 = at(i0[0], i0[1], i0[2]) * (1 - t[0]) + at(i1[0], i0[1], i0[2]) * t[0];
  const double c10 = at(i0[0], i1[1], i0[2]) * (1 - t[0]) + at(i1[0], i1[1], i0[2]) * t[0];
  const double c01 = at(i0[0], i0[1], i1[2]) * (1 - t[0]) + at(i1[0], i0[1], i1[2]) * t[0];
  const double c11 = at(i0[0], i1[1], i1[2]) * (1 - t[0]) + at(i1[0], i1[1], i1[2]) * t[0];
  const double c0 = c00 * (1 - t[1]) + c10 * t[1];
  const double c1 = c01 * (1 - t[1]) + c11 * t[1];
  return c0 * (1 - t[2]) + c1 * t[2];
}

// Periodic separable Gaussian (sigma in cells, truncated at 4 sigma as
// scipy.ndimage.gaussian_filter) of one component stored with stride 6
void gauss_blur(std::vector<float>& a, int comp, int N, double sigma) {
  const int R = static_cast<int>(4 * sigma + 0.5);
  std::vector<double> w(2 * R + 1);
  double ws = 0;
  for (int k = -R; k <= R; ++k) ws += (w[k + R] = std::exp(-0.5 * k * k / (sigma * sigma)));
  for (double& x : w) x /= ws;
  for (int axis = 0; axis < 3; ++axis) {
#pragma omp parallel for schedule(static) if (par::worth(static_cast<long long>(N) * N * N))
    for (int b = 0; b < N; ++b) {
      std::vector<double> line(N), outl(N);
      for (int a2 = 0; a2 < N; ++a2) {
        auto cell = [&](int t) -> std::int64_t {
          if (axis == 0) return idx(t, a2, b, N);
          if (axis == 1) return idx(a2, t, b, N);
          return idx(a2, b, t, N);
        };
        for (int t = 0; t < N; ++t) line[t] = a[6 * cell(t) + comp];
        for (int t = 0; t < N; ++t) {
          double s = 0;
          for (int k = -R; k <= R; ++k) s += w[k + R] * line[wrapi(t + k, N)];
          outl[t] = s;
        }
        for (int t = 0; t < N; ++t) a[6 * cell(t) + comp] = static_cast<float>(outl[t]);
      }
    }
  }
}

// h = P / tr P (isotropic where tr P <= 0), as compute_otvet_tensor's store
inline void normalise(const double P[6], float* h) {
  const double tr = (P[0] + P[1]) + P[2];
  if (!(tr > 0)) {
    h[0] = h[1] = h[2] = static_cast<float>(1.0 / 3);
    h[3] = h[4] = h[5] = 0;
    return;
  }
  for (int m = 0; m < 6; ++m) h[m] = static_cast<float>(P[m] / tr);
}

// The same for the screened and cone tensors, whose P can fall to the noise
// of the single-precision FFT (exponential cutoffs): isotropic below trmin,
// otherwise h = P / tr P with its eigenvalues clamped to >= 0 (trace 1).
inline void clamp_eigen(double h[6]) {
  const double e = otvet_min_eigenvalue(h);
  if (e < 0) {
    for (int m = 0; m < 3; ++m) h[m] -= e;
    for (int m = 0; m < 6; ++m) h[m] /= 1 - 3 * e;
  }
}
inline void normalise_safe(const double P[6], double trmin, float* hout) {
  const double tr = (P[0] + P[1]) + P[2];
  if (!(tr > trmin)) {
    hout[0] = hout[1] = hout[2] = static_cast<float>(1.0 / 3);
    hout[3] = hout[4] = hout[5] = 0;
    return;
  }
  double h[6];
  for (int m = 0; m < 6; ++m) h[m] = P[m] / tr;
  clamp_eigen(h);
  for (int m = 0; m < 6; ++m) hout[m] = static_cast<float>(h[m]);
}
// FFT noise floor of tr P for sources rho (normalised units)
constexpr double kNoiseFloor = 1e-6;

// Sum in plane order (thread-count independent)
template <class F>
double plane_sum(int N, F&& plane) {
  std::vector<double> part(N, 0.0);
#pragma omp parallel for schedule(static) if (par::worth(static_cast<long long>(N) * N * N))
  for (int k = 0; k < N; ++k) part[k] = plane(k);
  double t = 0;
  for (int k = 0; k < N; ++k) t = t + part[k];
  return t;
}

double total(const double* a, int N) {
  const std::int64_t sk = static_cast<std::int64_t>(N) * N;
  return plane_sum(N, [&](int k) {
    double t = 0;
    for (std::int64_t c = k * sk; c < (k + 1) * sk; ++c) t += a[c];
    return t;
  });
}

// Gamma-weighted photons per substep per unit L
double qw_of(const OtvetClosureIn& in) {
  double q = 0;
  for (int g = 0; g < kNfreq; ++g) q += in.q[g] * in.sigma18[g];
  return q;
}

// Luminosity per cell, normalised by the largest source (single-precision
// FFT); returns the normalisation
double source_density(const OtvetClosureIn& in, std::vector<double>& rho) {
  const int N = in.N;
  double Lmax = 0;
  for (const OtvetSource& s : in.src) Lmax = std::max(Lmax, s.L);
  rho.assign(static_cast<std::size_t>(N) * N * N, 0.0);
  if (!(Lmax > 0)) return 0;
  for (const OtvetSource& s : in.src) rho[idx(s.i, s.j, s.k, N)] += s.L / Lmax;
  return Lmax;
}

// The thin tensor and N_dir = tr P/(4 pi) in L units (compute_otvet_tensor)
void thin_tensor(const OtvetClosureIn& in, int method, std::vector<float>& h,
                 std::vector<float>& ndir) {
  const int N = in.N;
  const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
  h.assign(6 * nc, 0.0f);
  ndir.assign(nc, 0.0f);
  if (in.src.empty()) {
    for (std::int64_t c = 0; c < nc; ++c) {
      h[6 * c] = h[6 * c + 1] = h[6 * c + 2] = static_cast<float>(1.0 / 3);
    }
    return;
  }
  if (method == 0)
    method = (!otvet_fft_available() ||
              static_cast<double>(in.src.size()) * static_cast<double>(nc) <= 2.7e8)
                 ? 1
                 : 2;
  if (method == 1) {
#pragma omp parallel for schedule(static) if (par::worth(nc))
    for (int k = 0; k < N; ++k)
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
          double P[6] = {0, 0, 0, 0, 0, 0}, K[6];
          for (const OtvetSource& s : in.src) {
            otvet_kernel(N, mi((i - s.i + N) % N, N), mi((j - s.j + N) % N, N),
                         mi((k - s.k + N) % N, N), K);
            for (int m = 0; m < 6; ++m) P[m] += s.L * K[m];
          }
          const std::int64_t c = idx(i, j, k, N);
          normalise(P, &h[6 * c]);
          ndir[c] = static_cast<float>(((P[0] + P[1]) + P[2]) / (4 * M_PI));
        }
    return;
  }
  std::vector<double> rho;
  const double Lmax = source_density(in, rho);
  std::array<std::vector<double>, 6> P;
  otvet_fft_convolve(N, rho, P);
#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (std::int64_t c = 0; c < nc; ++c) {
    const double Pc[6] = {P[0][c], P[1][c], P[2][c], P[3][c], P[4][c], P[5][c]};
    normalise(Pc, &h[6 * c]);
    ndir[c] = static_cast<float>(std::max(0.0, (Pc[0] + Pc[1]) + Pc[2]) * Lmax / (4 * M_PI));
  }
}

}  // namespace

// ---------------------------------------------------------------------------
// The closures
// ---------------------------------------------------------------------------

void otvet_closure_tensor(const OtvetClosureParams& p, int tensor_method, const OtvetClosureIn& in,
                          OtvetClosureOut& out) {
  const double t0 = wall_seconds();
  const int N = in.N;
  const std::int64_t nc = static_cast<std::int64_t>(N) * N * N;
  const double qw = qw_of(in);
  const double target = in.nw ? total(in.nw, N) : 0.0;
  double Lsum = 0;
  for (const OtvetSource& s : in.src) Lsum += s.L;
  out.lambda = out.lam_x = out.ratio = 0;
  out.ndir.clear();

  if (p.closure == 0) {
    // thin (+ split): the original tensor and its trace
    thin_tensor(in, tensor_method, out.h, out.ndir);
    if (!p.split) out.ndir.clear();
    if (target > 0 && !out.ndir.empty()) {
      std::vector<double> nd(out.ndir.begin(), out.ndir.end());
      out.ratio = qw * total(nd.data(), N) / target;
    }
    out.seconds = wall_seconds() - t0;
    return;
  }

  if (!otvet_fft_available()) throw std::runtime_error("otvet_closure = screen | cone need FFTW3f");
  out.h.assign(6 * nc, 0.0f);
  if (in.src.empty() || !(Lsum > 0)) {
    for (std::int64_t c = 0; c < nc; ++c) out.h[6 * c] = out.h[6 * c + 1] = out.h[6 * c + 2] =
        static_cast<float>(1.0 / 3);
    if (p.split) out.ndir.assign(nc, 0.0f);
    out.seconds = wall_seconds() - t0;
    return;
  }

  std::vector<double> rho;
  const double Lmax = source_density(in, rho);

  if (p.closure == 1) {
    // screen: lambda from the budget, sum q_w tr P/(4 pi) = q_w Lsum lambda = target
    double lam = p.screen_lambda;
    if (!(lam > 0)) lam = (target > 0 && qw > 0) ? target / (qw * Lsum) : 1.0;
    lam = std::max(lam, 1.0);
    out.lambda = lam;
    std::array<std::vector<double>, 6> P;
    const ScreenCtx ctx{lam};
    // cached only for a fixed lambda (the budget's changes every step)
    otvet_fft_convolve_kernel(N, rho, screen_fill, &ctx, P,
                              p.screen_lambda > 0 ? "screen" : nullptr, lam);
    if (p.split) out.ndir.assign(nc, 0.0f);
    const double trmin = kNoiseFloor * total(rho.data(), N);
#pragma omp parallel for schedule(static) if (par::worth(nc))
    for (std::int64_t c = 0; c < nc; ++c) {
      const double Pc[6] = {P[0][c], P[1][c], P[2][c], P[3][c], P[4][c], P[5][c]};
      normalise_safe(Pc, trmin, &out.h[6 * c]);
      if (p.split)
        out.ndir[c] = static_cast<float>(std::max(0.0, (Pc[0] + Pc[1]) + Pc[2]) * Lmax / (4 * M_PI));
    }
    if (target > 0) {
      const double tr = plane_sum(N, [&](int k) {
        double t = 0;
        for (std::int64_t c = k * static_cast<std::int64_t>(N) * N;
             c < (k + 1) * static_cast<std::int64_t>(N) * N; ++c)
          t += (P[0][c] + P[1][c]) + P[2][c];
        return t;
      });
      out.ratio = qw * tr * Lmax / (4 * M_PI) / target;
    }
    out.seconds = wall_seconds() - t0;
    return;
  }

  // ---- cone ---------------------------------------------------------------
  const int co = p.cone_coarse;
  const int ncr = N / co;
  const std::int64_t nrec = static_cast<std::int64_t>(ncr) * ncr * ncr;
  const double rc = p.cone_rc;
  const double smax = p.cone_smax * N;

  // h_s: the split thin tensor (blend), from the Gamma-weighted N
  std::vector<float> hs;
  if (p.cone_blend > 0) {
    std::vector<float> nd;
    thin_tensor(in, tensor_method, hs, nd);
#pragma omp parallel for schedule(static) if (par::worth(nc))
    for (std::int64_t c = 0; c < nc; ++c) {
      const double Nw = in.nw ? in.nw[c] : 0.0;
      if (!(Nw > 0)) continue;
      const double beta = std::min(1.0, qw * static_cast<double>(nd[c]) / Nw);
      float* h = &hs[6 * c];
      for (int m = 0; m < 6; ++m) {
        const double iso = m < 3 ? 1.0 / 3 : 0.0;
        h[m] = static_cast<float>(iso + beta * (static_cast<double>(h[m]) - iso));
      }
    }
  }

  // near field r < 2 r_c (FFT), in L/Lmax units
  std::array<std::vector<double>, 6> Pn;
  {
    const NearCtx ctx{rc};
    otvet_fft_convolve_kernel(N, rho, near_fill, &ctx, Pn, "near", rc);
  }
  // mip-maps of the emissivity (L/Lmax per cell) and of the capped HI column
  std::vector<Level> R, C;
  {
    std::vector<float> f(nc);
    for (std::int64_t c = 0; c < nc; ++c) f[c] = static_cast<float>(rho[c]);
    R = pyramid(f, N);
    const double cap = in.sigma18[0] > 0 ? p.cone_taucap / in.sigma18[0] : 0.0;
    for (std::int64_t c = 0; c < nc; ++c) {
      const double col = in.col ? static_cast<double>(in.col[c]) : 0.0;
      f[c] = static_cast<float>(std::min(col, cap));
    }
    C = pyramid(f, N);
  }
  const double trmin = kNoiseFloor * total(rho.data(), N);
  std::vector<double>().swap(rho);
  const int nlev = static_cast<int>(R.size());

  // transmission weights of the source spectrum (Gamma-weighted)
  std::array<double, kNfreq> W{};
  double Ws = 0;
  for (int g = 0; g < kNfreq; ++g) Ws += (W[g] = in.q[g] * in.sigma18[g]);
  for (int g = 0; g < kNfreq; ++g) W[g] = Ws > 0 ? W[g] / Ws : 0.0;
  auto trans = [&](double col) {
    double t = 0;
    for (int g = 0; g < kNfreq; ++g) t += W[g] * std::exp(-in.sigma18[g] * col);
    return t;
  };

  // half-octave distance bins from r_c (budget), else one bin
  std::vector<double> edges;
  if (p.cone_budget) {
    const int ne = static_cast<int>(std::ceil(2 * std::log2(smax / rc) + 2));
    for (int i = 0; i < ne; ++i) edges.push_back(rc * std::pow(2.0, i / 2.0));
  } else {
    edges = {rc, smax};
  }
  const int nb = static_cast<int>(edges.size()) - 1;

  const std::vector<std::array<double, 3>> dirs = otvet_healpix_dirs(p.cone_nside);
  const int npix = static_cast<int>(dirs.size());
  const double omega = 4 * M_PI / npix;
  const double tanth = std::tan(std::acos(1 - omega / (2 * M_PI)));
  constexpr double kStep = 0.5, kTmin = 1e-4;

  std::vector<float> Pb(static_cast<std::size_t>(nrec) * nb * 6, 0.0f);
#pragma omp parallel for schedule(dynamic, 64) if (nrec >= 4096)
  for (std::int64_t r = 0; r < nrec; ++r) {
    const int ix = static_cast<int>(r % ncr), iy = static_cast<int>((r / ncr) % ncr),
              iz = static_cast<int>(r / (static_cast<std::int64_t>(ncr) * ncr));
    const double X[3] = {(ix + 0.5) * co, (iy + 0.5) * co, (iz + 0.5) * co};
    std::vector<double> acc(static_cast<std::size_t>(nb) * 6, 0.0), Ib(nb);
    for (int pix = 0; pix < npix; ++pix) {
      const std::array<double, 3>& n = dirs[pix];
      const double nn[6] = {n[0] * n[0], n[1] * n[1], n[2] * n[2],
                            n[0] * n[1], n[0] * n[2], n[1] * n[2]};
      std::fill(Ib.begin(), Ib.end(), 0.0);
      double s = rc, colacc = 0;
      while (s < smax) {
        const double d = std::max(2 * s * tanth, 1.0);
        const double ds = kStep * d;
        const double sm = s + ds / 2;
        const int l = std::max(0, std::min(nlev - 1, static_cast<int>(std::floor(std::log2(d)))));
        const double f = std::ldexp(1.0, l);
        const double pos[3] = {X[0] + sm * n[0], X[1] + sm * n[1], X[2] + sm * n[2]};
        const double rv = sample(R[l], f, pos);
        const double cv = sample(C[l], f, pos);
        const double w = smoothstep(sm, rc, 2 * rc);
        const double T = trans(colacc + cv * ds / 2);
        int b = static_cast<int>(std::upper_bound(edges.begin(), edges.end(), sm) - edges.begin()) - 1;
        b = std::max(0, std::min(nb - 1, b));
        Ib[b] += w * rv * omega * ds * T;
        colacc += cv * ds;
        s += ds;
        if (T < kTmin) break;
      }
      for (int b = 0; b < nb; ++b)
        for (int m = 0; m < 6; ++m) acc[6 * b + m] += Ib[b] * nn[m];
    }
    float* o = &Pb[static_cast<std::size_t>(r) * nb * 6];
    for (int q = 0; q < nb * 6; ++q) o[q] = static_cast<float>(acc[q]);
  }

  // budget: lambda_x such that the direct light equals the photons held, both
  // summed over the ionised cells (x_HII > 0.5; all cells without xion):
  // before overlap the unattenuated near field puts light into the neutral
  // gas around small H II regions, which the solution never holds there
  std::vector<double> fac(nb, 1.0);
  auto ion = [&](std::int64_t c) { return !in.xion || in.xion[c] > 0.5f; };
  const double target_ion = in.nw ? plane_sum(N, [&](int k) {
    double t = 0;
    for (std::int64_t c = k * static_cast<std::int64_t>(N) * N;
         c < (k + 1) * static_cast<std::int64_t>(N) * N; ++c)
      if (ion(c)) t += in.nw[c];
    return t;
  }) : 0.0;
  if (target_ion > 0 && qw > 0) {
    const double unit = qw * Lmax / (4 * M_PI);  // photons per (L/Lmax) unit of tr P
    const double near_tr = plane_sum(N, [&](int k) {
      double t = 0;
      for (std::int64_t c = k * static_cast<std::int64_t>(N) * N;
           c < (k + 1) * static_cast<std::int64_t>(N) * N; ++c)
        if (ion(c)) t += (Pn[0][c] + Pn[1][c]) + Pn[2][c];
      return t;
    });
    // ionised fraction of each receiver's block of coarse^3 cells
    std::vector<double> wr(nrec, 1.0);
    if (in.xion) {
#pragma omp parallel for schedule(static) if (nrec >= 4096)
      for (std::int64_t r = 0; r < nrec; ++r) {
        const int ix = static_cast<int>(r % ncr), iy = static_cast<int>((r / ncr) % ncr),
                  iz = static_cast<int>(r / (static_cast<std::int64_t>(ncr) * ncr));
        int n = 0;
        for (int dz = 0; dz < co; ++dz)
          for (int dy = 0; dy < co; ++dy)
            for (int dx = 0; dx < co; ++dx) n += ion(idx(ix * co + dx, iy * co + dy, iz * co + dz, N));
        wr[r] = static_cast<double>(n) / (static_cast<double>(co) * co * co);
      }
    }
    std::vector<double> tb(nb, 0.0), sb(nb);
    for (int b = 0; b < nb; ++b) sb[b] = std::sqrt(edges[b] * edges[b + 1]);
    for (std::int64_t r = 0; r < nrec; ++r)
      for (int b = 0; b < nb; ++b) {
        const float* o = &Pb[(static_cast<std::size_t>(r) * nb + b) * 6];
        tb[b] += wr[r] * ((static_cast<double>(o[0]) + o[1]) + o[2]);
      }
    const double cvol = static_cast<double>(co) * co * co;
    auto light = [&](double lx) {
      double t = near_tr;
      for (int b = 0; b < nb; ++b) t += cvol * tb[b] * (lx > 0 ? std::exp(-sb[b] / lx) : 1.0);
      return t * unit;
    };
    const double target = target_ion;
    double lx = 0;
    if (p.cone_budget && light(0) > target) {
      if (near_tr * unit >= target) {
        lx = 0;  // the near field alone exceeds the budget: no cutoff
      } else {
        double lo = 0, hi = 12;  // log10 lambda_x (cells)
        for (int it = 0; it < 60; ++it) {
          const double mid = 0.5 * (lo + hi);
          if (light(std::pow(10.0, mid)) > target)
            hi = mid;
          else
            lo = mid;
        }
        lx = std::pow(10.0, 0.5 * (lo + hi));
      }
      if (lx > 0)
        for (int b = 0; b < nb; ++b) fac[b] = std::exp(-sb[b] / lx);
    }
    out.lam_x = lx;
    out.ratio = light(lx) / target;
  }

  // cone part per receiver, prolonged (trilinear, periodic), plus the near
  // field: h_cone; then the blend with h_s
  std::vector<double> Pr(static_cast<std::size_t>(nrec) * 6, 0.0);
#pragma omp parallel for schedule(static) if (nrec >= 4096)
  for (std::int64_t r = 0; r < nrec; ++r)
    for (int b = 0; b < nb; ++b) {
      const float* o = &Pb[(static_cast<std::size_t>(r) * nb + b) * 6];
      for (int m = 0; m < 6; ++m) Pr[6 * r + m] += fac[b] * static_cast<double>(o[m]);
    }
  std::vector<float>().swap(Pb);
#pragma omp parallel for schedule(static) if (par::worth(nc))
  for (int k = 0; k < N; ++k) {
    double t[3];
    int a0[3], a1[3];
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i) {
        const int fi[3] = {i, j, k};
        for (int d = 0; d < 3; ++d) {
          const double cc = (fi[d] + 0.5) / co - 0.5;
          const double fl = std::floor(cc);
          t[d] = cc - fl;
          a0[d] = wrapi(static_cast<long long>(fl), ncr);
          a1[d] = a0[d] + 1 == ncr ? 0 : a0[d] + 1;
        }
        const std::int64_t c = idx(i, j, k, N);
        double P[6];
        for (int m = 0; m < 6; ++m) P[m] = Pn[m][c];
        for (int dz = 0; dz < 2; ++dz)
          for (int dy = 0; dy < 2; ++dy)
            for (int dx = 0; dx < 2; ++dx) {
              const double wgt = (dx ? t[0] : 1 - t[0]) * (dy ? t[1] : 1 - t[1]) * (dz ? t[2] : 1 - t[2]);
              const std::int64_t rr = idx(dx ? a1[0] : a0[0], dy ? a1[1] : a0[1], dz ? a1[2] : a0[2], ncr);
              for (int m = 0; m < 6; ++m) P[m] += wgt * Pr[6 * rr + m];
            }
        normalise_safe(P, trmin, &out.h[6 * c]);
      }
  }
  for (auto& v : Pn) std::vector<double>().swap(v);

  if (p.cone_blend > 0) {
    // h = h_s + G * (h_cone - h_s), eigenvalues clamped to >= 0
    for (std::int64_t q = 0; q < 6 * nc; ++q) out.h[q] = out.h[q] - hs[q];
    for (int m = 0; m < 6; ++m) gauss_blur(out.h, m, N, p.cone_blend);
#pragma omp parallel for schedule(static) if (par::worth(nc))
    for (std::int64_t c = 0; c < nc; ++c) {
      double h[6];
      for (int m = 0; m < 6; ++m) h[m] = static_cast<double>(hs[6 * c + m]) + out.h[6 * c + m];
      const double tr = (h[0] + h[1]) + h[2];
      for (int m = 0; m < 6; ++m) h[m] /= tr;
      clamp_eigen(h);
      for (int m = 0; m < 6; ++m) out.h[6 * c + m] = static_cast<float>(h[m]);
    }
  }
  out.seconds = wall_seconds() - t0;
}

}  // namespace radray
