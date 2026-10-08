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
  int Ndm = 0;   // Ngas/2 (only used for unit bookkeeping)
  std::int64_t Ncell = 0;  // Ngas^3
  std::int64_t Nrtc = 0;   // Nrt^3
  double ngasr = 0, nrtr = 0, hrtr = 0, nrtcr = 0, ncellr = 0;

  // N_rt_gas = Ngas/Nrt (gas cells per RT cell per dimension)
  int N_rt_gas = 1;
  int N_rt_gas_c = 1;
  double x_rt_gas = 1;

  double x_gas_box = 0;  // box/ngasr
  double x_rt_box = 0;   // box/nrtr

  // Rays: Nray_max = nray_max_factor*Nrt^3 (64 in the original),
  // Nfree_max = Nray_max/(4*Ncpu)
  int nray_max_factor = 64;
  std::int64_t Nray_max = 0;
  std::int64_t Nfree_max = 0;

  static Config make(int N, int nray_max_factor = 64);
};

}  // namespace radray
