// Plane-parallel RT-domain mode (docs/plane_parallel.md): the parameters
// (Config::plane_*), the spectrum of the face sources, their normalisation
// and the photon bookkeeping shared by the serial RT (src/ray_tracing.cpp,
// src/rt_tools.cpp) and the AMReX backend (src/amrex/rt_level.cpp). New
// options, not in the original; everything here is inert unless
// plane_parallel = 1.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "radray/config.hpp"
#include "radray/constants.hpp"

namespace rrk {
struct RTP;  // include/radray/kernels/ray_kernels.hpp
}

namespace radray {

struct Simulation;

// The 5 frequency groups of the mode: J_nu ~ nu^-alpha between nu_min and
// nu_max (eV) cut into kNfreq bins of equal photon number, each
// represented by its photon-weighted mean frequency nu[n] (the effective
// frequency every ray carries, at which the cell physics evaluates the
// cross sections), with the photon fraction frac[n] = 1/kNfreq and the
// Verner et al. 1996 cross sections (rt_cell.hpp photoion_sigma) there.
struct PlaneSpectrum {
  double alpha = 1.5;
  double edge[kNfreq + 1] = {};        // bin edges, eV
  double nu[kNfreq] = {};              // effective frequency per bin, eV
  double frac[kNfreq] = {};            // photon fraction per bin (equal)
  double sig_HI[kNfreq] = {}, sig_HeI[kNfreq] = {}, sig_HeII[kNfreq] = {};  // cm^2
  // sum_n frac[n] sig_HI[n]: Gamma_HI = F * sig_eff for a photon flux F
  double sig_eff = 0;
};
PlaneSpectrum plane_spectrum(double alpha, double nu_min_eV, double nu_max_eV);
PlaneSpectrum plane_spectrum(const Config& c);  // 1 Ry (HI_eV) .. plane_nu_max_ry Ry

// The effective frequencies of the 5 groups the optical depths are built
// with (rt_step_1 / post_rt / the level's tau): the original's freq_init,
// or the mode's bins with plane_parallel = 1.
std::array<double, kNfreq> group_frequencies(const Config& c);

// Per-simulation state of the mode (Simulation::plane): the spectrum and
// the per-hydro-step photon bookkeeping of the conservation test
// (stats.txt plane_emitted / plane_absorbed / plane_escaped, in the rays'
// units: photons per <nH> per RT-cell volume, summed over the step's
// substeps; the escaped and absorbed sums are exact fixed-point sums, so
// they do not depend on the thread or rank count).
struct PlaneState {
  bool on = false;
  int L = 0;      // domain length in RT cells
  int faces = 0;  // Config::plane_faces
  int transport = 0;  // Config::plane_transport: 0 rays, 1 sweep
  PlaneSpectrum spec;
  double emitted = 0;
  __int128 absorbed_fx = 0, escaped_fx = 0;  // 2^kPlaneFxBits fixed point
  static constexpr int kFxBits = 40;
  double absorbed() const;
  double escaped() const;
  void reset_step();
  // plane_transport = sweep, serial RT: the packets in flight per lit face
  // (face_list order), kNfreq floats per cell in the RT grid's cell order
  // (the photons about to enter the cell along that face's axis; the
  // launch cells, 0 cells from the face, hold nothing: their packet is the
  // source packet), and their sum after the last advance (fixed point)
  std::vector<int> face_list;
  std::vector<std::vector<float>> field;
  __int128 inflight_fx = 0;
  double inflight() const;
};
__int128 plane_to_fx(double v);

// plane_init: fills Simulation::plane from cfg (after init_rt); a no-op
// without the mode. plane_check: the parameter rules (throws).
void plane_init(Simulation& s);
void check_plane_param(const Config& c);

// params.txt plane_* keys and the two step caps (dt_max_myr,
// rt_substeps_max): true if the key was handled.
bool set_plane_param(Config& c, const std::string& key, const std::string& val);
void copy_plane_param(const Config& from, Config& to);
std::string plane_summary(const Config& c);  // one line for the start-up echo
std::string plane_faces_string(int faces);   // "+x +y"

// Photons per ray this substep (RTP.plane.ngam), from the target Gamma_HI:
// N = Gamma * dtsec_rt / (nH_avg * len_rt * sig_eff) (per <nH> per RT-cell
// volume: a ray carries the photons that cross one cell face in one
// substep), split over the bins by frac. Independent of
// rt_cells_per_substep (k cells share a ray's photons and the rate uses
// the k-cell crossing speed).
void plane_fill_rtp(const Simulation& s, rrk::RTP& P);
// The photons emitted per substep: face cells x photons per ray
double plane_emitted_per_substep(const Simulation& s, const rrk::RTP& P);
// Whether the RT runs in this hydro step of apps/radhydro: the step's start
// is at or below plane_z_re (a_start >= a_re within 1e-9)
bool plane_rt_active(const Config& c, double a_start);
// <dir>/plane_spectrum.txt: the bins, the cross sections and the fluxes
void plane_write_spectrum(const std::string& dir, const Simulation& s);

// plane_transport = sweep (src/plane_sweep.cpp), the serial RT: the face
// list of a faces mask; the cast of every packet in flight (and of the
// source packet at the launch cells) into the RT grid's deposits in
// rt_step_1, with rrk::cast_1; the advance in rt_step_3, with rrk::cast_3:
// the transmitted packet moves rt_cells_per_substep cells along its axis,
// leaves at the domain's far face (escaped) or is dropped when faint.
// Thread-count independent (one column of cells per thread, fixed-point
// sums). The kernel helpers are shared with the backend.
std::vector<int> plane_face_list(int faces);
void plane_sweep_cast(Simulation& s);
void plane_sweep_advance(Simulation& s);

}  // namespace radray
