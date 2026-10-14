// Moment-based radiative transfer: an alternative to the ray tracing
// (docs/rt_moments.md). Not in the original code.
//
// rt_method = m1 | otvet replaces the ray casting of every RT substep
// (ray_tracing_1/3) by an explicit two-moment solver for the photon number
// density N and flux F of each of the kNfreq frequency groups, on the RT
// grid, at the reduced speed of light of the run:
//
//   dN/dt + c_hat div F         = S - absorption
//   dF/dt + c_hat^2 div (D N)   =   - absorption
//
// with the Eddington tensor D from the M1 closure (Levermore 1984; local in
// (N, F)) or from the optically thin variable Eddington tensor of Gnedin &
// Abel 2001 (OTVET: the 1/r^2 weighted sum over the sources of the outer
// product of the unit vectors, recomputed every hydro step). Transport is
// the first-order Global Lax-Friedrichs (Rusanov) scheme with wave speed
// c_hat (as RAMSES-RT, Rosdahl et al. 2013, and ATON, Aubert & Teyssier
// 2008), moment_nsub subcycles per RT substep (an RT substep is one cell of
// light travel, so the Courant number is 1/moment_nsub). The photon-weighted
// group energy E = nu N is advected with N.
//
// The physics is the unchanged radiative_physics (photo_cell): the moment
// solver fills rt_grid%ngam/nu/dx as the ray casting does (N is the photon
// density in the cell during the substep, in units of <nH>, the quantity
// the rays deposit), and after photo_cell scales N, F and E by the photon
// fraction that survives (rt_grid%fgam), which is the absorption term above.
// The sources inject fResc*pop2_init*Luv*dtsec_rt photons per substep in
// their cell, the number ray_source emits.
//
// Everything here is deterministic and independent of the thread count:
// the transport is a Jacobi stencil, the sums are in a fixed order.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "radray/simulation.hpp"

namespace radray {

struct MomentState {
  int N = 0;
  std::int64_t ncell = 0;
  // Per group g and cell c at index g*ncell + c: photon density N in units
  // of <nH> (rt_grid%ngam), reduced flux F/c_hat (|F| <= N) and E = nu*N
  // (eV <nH>).
  std::vector<float> n, fx, fy, fz, e;
  // OTVET Eddington tensor (xx, yy, zz, xy, xz, yz per cell); empty for M1.
  std::vector<float> h;
  // Photon budget in photons per H atom of the box (sum over cells of N
  // divided by Nrt^3, the normalisation of rt_stat%ngam), cumulative:
  // emitted by the sources, absorbed in photo_cell (including the photons it
  // drops below ngam_min), lost in the redshift (groups below ngam_min).
  double injected = 0, absorbed = 0, lost = 0;
  // Diagnostics of the last OTVET tensor: source cells, method used
  // (1 direct, 2 FFT).
  std::int64_t nsrc = 0;
  int tensor_method = 0;
  // otvet_split = 1 (otvet_closure.hpp): the direct light tr P/(4 pi) per
  // cell in luminosity units (N_dir,g = q_g ndir); empty otherwise. Not
  // checkpointed (recomputed every hydro step, like h).
  std::vector<float> ndir;
  // Diagnostics of the last alternative tensor (otvet_closure.hpp)
  double closure_lambda = 0, closure_lam_x = 0, closure_ratio = 0, closure_seconds = 0;
  // Scratch of the transport (one group at a time; not checkpointed)
  std::vector<float> P, nn, nfx, nfy, nfz, ne;
};

// params.txt keys rt_method, moment_nsub, otvet_tensor. Returns false for
// any other key.
bool set_rt_method_param(Config& c, const std::string& key, const std::string& val);
void copy_rt_method(const Config& from, Config& to);
const char* rt_method_name(RTMethod m);
std::string rt_method_summary(const Config& c);
inline bool use_moments(const Config& c) { return c.rt_method != RTMethod::Rays; }

// Allocates s.moments (idempotent; called by prep_rt).
void moment_init(Simulation& s);
// Once per hydro step, after the sources are set (prep_rt): the OTVET tensor.
void moment_prep(Simulation& s);
// The moment version of rt_step_1: resets the radiation fields of the RT
// grid, updates tau, injects the source photons, transports N, F, E and
// deposits them in rt_grid%ngam/nu/dx for radiative_physics.
void moment_step_1(Simulation& s);
// After radiative_physics: N, F, E *= rt_grid%fgam.
void moment_step_3(Simulation& s);
// post_rt: redshift of the group spectra over the hydro step (the
// redshift_spectrum of the rays, per cell; F follows the photons).
void moment_redshift(Simulation& s);
// Sum of N over groups and cells (photons in flight, in units of <nH> per
// cell), in the fixed cell order.
double moment_photons_total(const Simulation& s);

// OTVET tensor for the sources of the RT grid (cells with sfr > sfr_min,
// weight Luv): h = P / tr P with P = sum_s L_s d d^T / |d|^4 over minimum-
// image separations d in cells (an isotropic L_s*4/3 in the source cell).
// method 1 = direct sum (double), 2 = FFT convolution (single precision
// FFTW; needs RADRAY_HAVE_FFTW), 0 = direct when Nsrc*Ncell is small.
// Returns the method used.
int compute_otvet_tensor(const RTGrid& G, const Units& u, int method, std::vector<float>& h,
                         std::int64_t* nsrc = nullptr);
bool otvet_fft_available();

// The convolution of the FFT method (src/rt_moments_fft.cpp): P[m][c] =
// sum over cells c' of rho[c'] K_m(c - c'), for the 6 kernel components.
void otvet_fft_convolve(int N, const std::vector<double>& rho,
                        std::array<std::vector<double>, 6>& P);
// The OTVET kernel K_m(d) (m = xx, yy, zz, xy, xz, yz) for a minimum-image
// separation d in (-N/2, N/2]^3. Off-diagonal terms with a component at
// exactly N/2 (two equidistant images of opposite sign) are 0.
void otvet_kernel(int N, int dx, int dy, int dz, double K[6]);

}  // namespace radray
