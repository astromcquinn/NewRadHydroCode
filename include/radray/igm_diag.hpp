// On-the-fly IGM diagnostics (diag_dz > 0, docs/plane_parallel.md "On-the-fly
// diagnostics"): the kinetic energy spectrum E(k) of the gas velocity and
// the temperature-density histogram, written to <outdir>/Diag every
// diag_dz in redshift instead of full snapshots. Read-only on the state:
// the run's evolution does not depend on them.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace radray {

// E(k) as tools/analysis/velocity_power.py: v (km/s, proper peculiar, lab
// frame) Fourier transformed on the N^3 grid, |v_hat|^2 summed over shells
// in nbins logarithmic bins of |k| from the fundamental to the Nyquist
// (comoving h/kpc), normalised so that sum E dk = <v^2>/2. Em: the same
// for sqrt(rho/<rho>) v (velocity_power.py --weight mass).
struct VelocitySpectrum {
  std::vector<double> edges, k, cnt, E, Em;
  double v2 = 0, v2m = 0;  // <v^2> and <rho v^2>/<rho>, (km/s)^2
};

// u: u_gas(5,N^3) (rho, rho v, E; i fastest; the lab frame, i.e. a zero grid
// velocity), vel_kms: km/s per code velocity, box: comoving Mpc/h.
// nthreads: FFTW threads (the transform is byte-identical for every count,
// docs/pm_refine.md); the shell sums are per k-plane and then in order.
void velocity_spectrum(const float* u, int N, double box, double vel_kms, int nbins, int nthreads,
                       VelocitySpectrum& out);
void write_velocity_spectrum(const std::string& path, double z, int N, double box,
                             const VelocitySpectrum& s);

// Volume and mass in bins of log10(Delta) (Delta = rho/<rho> of the gas)
// and log10(T [K]): kTDeltaND bins of 0.05 from -2, kTDeltaNT of 0.05 from
// 0; cells outside go to the edge bins.
constexpr int kTDeltaND = 120, kTDeltaNT = 160;
constexpr double kTDeltaD0 = -2, kTDeltaT0 = 0, kTDeltaBin = 0.05;
struct TDeltaHist {
  std::vector<std::int64_t> count;  // [iD * kTDeltaNT + iT]
  std::vector<double> mass;         // sum of Delta
};
// rho: u_gas(1,:) (stride 5), temp: the RT grid's T, ncell cells each
void tdelta_hist(const float* u, const float* temp, std::int64_t ncell, int N, TDeltaHist& out);
void write_tdelta(const std::string& path, double z, const TDeltaHist& h);

// Whether a diagnostic is due in a step from z_prev to z (z < z_prev):
// outputs at z0, z0 - dz, z0 - 2 dz, ... (none above z0).
bool diag_due(double z0, double dz, double z_prev, double z);

}  // namespace radray
