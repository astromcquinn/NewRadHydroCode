// The complete state needed by the RT step, gathered in explicit structs
// (no hidden globals). This replaces the module variables of
// radhydro_global.f90 / rt_rates.f90 that the RT code uses.
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <ostream>
#include <vector>

#include "radray/config.hpp"
#include "radray/cosmology.hpp"
#include "radray/grid.hpp"
#include "radray/ray_quant.hpp"
#include "radray/rays.hpp"
#include "radray/rt_rates.hpp"

namespace radray {

struct RayDump;  // radray/ray_dump.hpp

// radhydro_global.f90: rt_stat_type
struct RTStat {
  std::int64_t nsfill = 0, ngfill = 0;
  double Luv = 0, sfr = 0, tau = 0, dtau = 0;
  std::array<double, 3> ngam{}, nion{}, ncol{}, nrec{}, nrab{};
  std::array<double, 3> gamma{}, lambda{};
};

// domain_tools.f90: RT domain decomposition (init_domain_rt).
// indx_dom_rt(2,3,n) -> box[n] = {i1,i2,j1,j2,k1,k2} (1-based, inclusive).
struct Domains {
  int Ndomain = 0;
  int Nx_dom = 0;
  std::vector<std::array<int, 6>> box;          // per domain
  std::vector<std::array<int, 27>> domain;      // domain_rt(27,n): status + 26 neighbours
};

struct Simulation {
  Config cfg;
  SimClock clk;
  Units units;
  RTGrid rt;
  GasGrid gas;
  RayState rays;
  AtomicRates rates;
  RTStat rt_stat;  // rt_stat(1)  (Ncpu = 1)
  RTStat rt_save;
  Domains dom;
  std::int64_t Nsource = 0;

  // ray.log (Fortran unit 201) and stdout-like summary; may be null.
  std::ostream* raylog = nullptr;
  std::ostream* out = nullptr;

  // Optional ray snapshots for the compression study (ray_dump.hpp); null
  // (the default) means off. Never read by the RT.
  std::shared_ptr<RayDump> raydump_state;
  // Lossy ray-precision experiment (ray_quant.hpp); off by default.
  RayQuantConfig ray_quant;
};

}  // namespace radray
