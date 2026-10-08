// Grid / ray-storage sizes. Mirrors the parameters of radhydro_global.f90,
// but with a runtime N (Ngas = Nrt = N per docs/validation_spec.md).
#pragma once

#include <cstdint>

namespace radray {

struct Config {
  // radhydro_global.f90: Ncpu. The reference is run with ONE thread, and the
  // port is serial; Ncpu enters ray-buffer sizes and the domain
  // decomposition, so we keep it explicit.
  int Ncpu = 1;

  int Ngas = 0;  // gas grid size per dimension
  int Nrt = 0;   // RT grid size per dimension (== Ngas here)
  int Hrt = 0;   // Nrt/2
  int Ndm = 0;   // DM particles per dimension (Npart = Ndm^3); default Ngas/2
  std::int64_t Npart = 0;  // Ndm^3

  // cosmo_parameters.f90: box (comoving Mpc/h). A compile-time parameter in
  // the original (5); a run-time value here (the Fortran harnesses set it
  // with -DBOX).
  double box = 5;
  std::int64_t Ncell = 0;  // Ngas^3
  std::int64_t Nrtc = 0;   // Nrt^3
  double ngasr = 0, nrtr = 0, hrtr = 0, nrtcr = 0, ncellr = 0;

  // N_rt_gas = Ngas/Nrt (gas cells per RT cell per dimension)
  int N_rt_gas = 1;
  int N_rt_gas_c = 1;
  double x_rt_gas = 1;

  double x_gas_box = 0;  // box/ngasr
  double x_dm_box = 0;   // box/ndmr
  double x_rt_box = 0;   // box/nrtr

  // Rays: Nray_max = nray_max_factor*Nrt^3 (64 in the original),
  // Nfree_max = Nray_max/(4*Ncpu)
  int nray_max_factor = 64;
  std::int64_t Nray_max = 0;
  std::int64_t Nfree_max = 0;

  // Runtime switch for the RT bookkeeping fixes B1-B4, B8, B9
  // (docs/port_notes.md section 5). true (params.txt legacy_rt = 1, the
  // default) reproduces the original Fortran bit-for-bit.
  bool legacy_rt = true;

  // Runtime switch for the atomic-physics fixes P1-P6 (docs/port_notes.md
  // section 5b). true (params.txt legacy_physics = 1, the default) keeps the
  // original rates and cooling sum bit-for-bit. The rate tables depend on it:
  // call calc_atomic_rates(rates, legacy_physics) after changing it.
  bool legacy_physics = true;

  // Reduced speed of light: c_hat = c_reduced * c. Enters only through
  // calc_timestep (the original's cscale): nstep_rt = nint(dl_c/dl_rt*cscale).
  // 1 (default) is the original.
  double c_reduced = 1.0;

  // Deterministic parallel ray tracing (docs/parallel.md). false (default):
  // the original sweep order (8 domains, one thread) and the Fortran ray
  // pool, bit-identical to the Fortran harness. true: ray_tracing_1/3 run
  // over 8^3-cell domains in colour-major order, in parallel within a
  // colour, with per-thread free lists; bit-identical for any thread count,
  // but not to the legacy order. Requires legacy_rt = false.
  bool det_parallel = false;

  // Ngas = Nrt = N; Ndm <= 0 means the original's Ngas/2; box in Mpc/h.
  static Config make(int N, int nray_max_factor = 64, int Ndm = -1, double box = 5);
};

}  // namespace radray
