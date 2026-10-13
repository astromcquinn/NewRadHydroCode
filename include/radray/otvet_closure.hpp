// Alternative Eddington tensors for rt_method = otvet (docs/otvet_closure.md).
// Not in the original code; the default (otvet_closure = thin, otvet_split
// = 0) is the thin minimum-image tensor of compute_otvet_tensor, unchanged.
//
// The direct light of the sources, as a pressure tensor P (xx, yy, zz, xy,
// xz, yz) per RT cell, in units where a source of luminosity L contributes
// L d d^T / |d|^4 at separation d (cells), i.e. tr P = sum L / r^2:
//
//   thin    the minimum-image sum (compute_otvet_tensor);
//   screen  e^(-r/lambda) d d^T / r^4 over the 27 periodic images plus an
//           isotropic continuum for the images beyond, with lambda from the
//           photon budget (sum of the screened light = the photons the
//           solution holds) or fixed (otvet_screen_lambda);
//   cone    the thin near field r < r_c (FFT, handed over smoothly to the
//           cones between r_c and 2 r_c) plus 12 nside^2 HEALPix cones from
//           receivers on a grid otvet_cone_coarse times coarser, marched
//           through mip-maps of the emissivity and of the HI column (per-cell
//           group-1 tau capped at otvet_cone_taucap) out to otvet_cone_smax
//           boxes, with the Gamma-weighted transmission of the source
//           spectrum, and an extra e^(-s/lambda_x) with lambda_x from the
//           photon budget (otvet_cone_budget = 1): the direct light in the
//           ionised cells equals the photons the solution holds there (no
//           cutoff if the near field alone exceeds it). The cone part is
//           kept per half-octave distance bin, so lambda_x needs no second
//           march.
//
// h = P / tr P. Two corrections from the solution's own photons N:
//
//   otvet_split = 1 (thin or screen): per group g and cell, in the transport,
//     P_g = N_g I/3 + min(N_g, N_dir,g) (h - I/3), N_dir,g = q_g tr P/(4 pi)
//     (q_g the photons per unit L per substep injected into group g): the
//     direct light carries h, the rest of the photons is isotropic.
//   otvet_cone_blend = sigma > 0 (cone): h = h_s + G_sigma * (h_cone - h_s),
//     the small scales of the split thin tensor h_s (Gamma-weighted N at the
//     start of the hydro step) and the large scales of the cone; the
//     eigenvalues are then clamped to >= 0 (trace 1).
//
// Everything is deterministic and independent of the thread count.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "radray/constants.hpp"

namespace radray {

struct Config;

struct OtvetClosureParams {
  int closure = 0;  // 0 thin, 1 screen, 2 cone
  int split = 0;
  double screen_lambda = 0;  // RT cells; 0 = from the photon budget
  int cone_nside = 4;
  int cone_coarse = 4;
  double cone_rc = 8;
  double cone_taucap = 1;
  double cone_smax = 8;  // boxes
  int cone_budget = 1;
  double cone_blend = 4;  // RT cells; 0 = the cone alone
  // true if anything but the original thin tensor
  bool active() const { return closure != 0 || split != 0; }
};

// params.txt keys otvet_closure, otvet_split, otvet_screen_lambda,
// otvet_cone_{nside,coarse,rc,taucap,smax,budget,blend}. Returns false for
// any other key.
bool set_otvet_closure_param(Config& c, const std::string& key, const std::string& val);
void check_otvet_closure_param(const Config& c);
OtvetClosureParams otvet_closure_params(const Config& c);
std::string otvet_closure_summary(const Config& c);

struct OtvetSource {
  int i, j, k;  // 0-based RT cell
  double L;
};

struct OtvetClosureIn {
  int N = 0;
  std::vector<OtvetSource> src;
  // Gamma weighting: per group the photons injected per unit L per substep
  // (q_g, in the units of N) and sigma_HI(nu_g) in 1e-18 cm^2. Gamma-weighted
  // photons are sum_g sigma18_g N_g.
  std::array<double, kNfreq> q{}, sigma18{};
  // HI column per cell in 1e-18 cm^-2 (cone; may be null otherwise)
  const float* col = nullptr;
  // Gamma-weighted photons per cell of the solution (the budget of screen
  // and cone, and N of the blend's split); null: no budget (lambda and
  // lambda_x then need fixed values or are infinite) and h_s = h_thin.
  const double* nw = nullptr;
  // ionised fraction x_HII per cell (cone budget: only cells with x_HII > 0.5
  // count); null: every cell
  const float* xion = nullptr;
};

struct OtvetClosureOut {
  std::vector<float> h;     // 6 per cell
  std::vector<float> ndir;  // split: tr P / (4 pi) per cell (L units); else empty
  double lambda = 0;        // screen: lambda used (cells)
  double lam_x = 0;         // cone: lambda_x (cells; 0 = none)
  double ratio = 0;         // sum N_dir / sum N, Gamma-weighted (0 if no budget;
                            // cone: over the ionised cells)
  double seconds = 0;       // wall time
};

// The direct tensor h for the closure (and ndir if split). thin + split
// uses compute_otvet_tensor's h (method as otvet_tensor); screen and cone
// need FFTW3f.
void otvet_closure_tensor(const OtvetClosureParams& p, int tensor_method, const OtvetClosureIn& in,
                          OtvetClosureOut& out);

// Pieces, exposed for the tests and apps/otvet_tensor.cpp ----------------

// HEALPix ring-scheme pixel centres (unit vectors), in ring order.
std::vector<std::array<double, 3>> otvet_healpix_dirs(int nside);
// Generic FFT convolution with a kernel filled per component m by
// kfill(m, K) on the N^3 grid of minimum-image separations (K[i + N(j + N k)]
// for d = (mi(i), mi(j), mi(k))), single-precision FFTW; P[m][c] in double.
// fill must be a pure function of the separation (it runs in parallel).
using OtvetKernelFill = void (*)(int N, int m, const void* ctx, float* K);
void otvet_fft_convolve_kernel(int N, const std::vector<double>& rho, OtvetKernelFill fill,
                               const void* ctx, std::array<std::vector<double>, 6>& P);
// Smallest eigenvalue of the symmetric tensor (xx, yy, zz, xy, xz, yz).
double otvet_min_eigenvalue(const double h[6]);

}  // namespace radray
