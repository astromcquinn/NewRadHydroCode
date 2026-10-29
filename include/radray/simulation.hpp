// The complete state needed by the RT step, gathered in explicit structs
// (no hidden globals). This replaces the module variables of
// radhydro_global.f90 / rt_rates.f90 that the RT code uses.
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <functional>
#include <ostream>
#include <vector>

#include "radray/config.hpp"
#include "radray/cosmology.hpp"
#include "radray/grid.hpp"
#include "radray/halos.hpp"
#include "radray/plane_parallel.hpp"
#include "radray/ray_quant.hpp"
#include "radray/rays.hpp"
#include "radray/rt_rates.hpp"

namespace radray {

struct RayDump;  // radray/ray_dump.hpp
struct MomentState;  // radray/rt_moments.hpp
struct RecombState;  // radray/recomb_rad.hpp
struct AmrexRT;
struct AmrexState;      // src/amrex/rt_level.cpp (rt_backend = amrex)
struct TensorDiag;   // radray/rt_tensor_diag.hpp

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

// deterministic_parallel ray tracing: cubic domains of d^3 cells (d = 8
// when possible; d >= 2k+2 for rt_cells_per_substep = k), n1 = Nrt/d per
// dimension (even), grouped into 8 colours by the parity of their domain
// indices. Domains of one colour are at least d cells apart, at least twice
// the reach of one cell's ray work (<= d/2 cells, checked at run time), so
// they can be swept concurrently.
struct DetDomains {
  int d = 0, n1 = 0;
  std::array<std::vector<std::array<int, 6>>, 8> colour;  // boxes per colour
};

// Wall-clock seconds per RT phase, accumulated over the run (diagnostics).
struct RTTimers {
  double prep = 0, calc_rt1 = 0, calc_ray1 = 0, trace1 = 0, physics = 0, trace3 = 0,
         sync = 0, post = 0;
  std::int64_t substeps = 0;
};
double wall_seconds();

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
  DetDomains ddom;
  RTTimers timers;
  std::int64_t Nsource = 0;
  // The sources of the current catalogue (apply_sources): position in RT
  // cells and luminosity, for the RT hierarchy, which places a source on
  // the finest level covering it (docs/stage_f.md section 3); every rank
  // holds the list
  struct RTSource {
    double x[3];
    float L, sfr;
  };
  std::vector<RTSource> rt_sources;

  // Called after every RT substep (after the clock update), e.g. for movie
  // frames; must not change the state.
  std::function<void(Simulation&)> substep_hook;

  // ray.log (Fortran unit 201) and stdout-like summary; may be null.
  std::ostream* raylog = nullptr;
  std::ostream* out = nullptr;

  // Optional ray snapshots for the compression study (ray_dump.hpp); null
  // (the default) means off. Never read by the RT.
  std::shared_ptr<RayDump> raydump_state;
  // Lossy ray-precision experiment (ray_quant.hpp); off by default.
  RayQuantConfig ray_quant;
  // Moment RT state (rt_method = m1 | otvet, rt_moments.hpp); null with the
  // ray tracing.
  std::shared_ptr<MomentState> moments;
  // Recombination radiation (recomb_rad = 1, recomb_rad.hpp); null otherwise.
  std::shared_ptr<RecombState> recomb;
  // The AMReX level of rt_backend = amrex (rt_backend.hpp); null otherwise.
  // Between hydro steps the rt and gas grids above hold its state; the rays
  // live only on the level (rays is then unused except for the counters
  // and thresholds). ngam_total_backend: sum of the live rays' photons after
  // post_rt (ngam_rays_total).
  std::shared_ptr<AmrexRT> amrex;
  double ngam_total_backend = 0;
  // amrex_state = 1 (amrex_state.hpp): the gas state on the levels; the
  // rt and gas grids above are then rank 0's mirror only.
  std::shared_ptr<AmrexState> state;
  // Eddington-tensor diagnostic (rt_tensor_diag.hpp); null (off) by default.
  // Never read by the RT.
  std::shared_ptr<TensorDiag> tdiag;
  // On-the-fly FoF cache (fof_dt_myr > 0, halos.hpp); unused otherwise.
  FofCache fof;
  // Plane-parallel RT-domain mode (plane_parallel.hpp, docs/plane_parallel.md):
  // off unless plane_parallel = 1 (plane_init).
  PlaneState plane;

  // The substep-averaged HI photoionization rate (fields / movie_fields
  // gammaHI_avg, docs/igm_diagnostics.md): per RT cell the sum over the
  // current hydro step's RT substeps of gamma(1) * dtsec_rt, and the
  // summed time, reset by radiative_transfer at the start of each hydro
  // step's RT (rt.hpp gamma_avg_*). Output only: never read by the RT,
  // not in the checkpoint (a restart begins at a step boundary; only a
  // step without RT right after a restart would write zeros where the
  // uninterrupted run writes the previous step's average). Serial RT;
  // the moment RT on rt_backend = amrex keeps acc on its level instead
  // (moment_level.cpp; acc stays empty, time and substeps are kept here).
  struct GammaAvg {
    bool on = false;
    std::vector<double> acc;  // per RT cell, s^-1 * s
    double time = 0;          // s
    std::int64_t substeps = 0;
  } gamma_avg;

  // X-ray background heating (xray_heating = 1, xray_background.hpp,
  // docs/xray_heating.md): the uniform rates per HI, HeI, HeII atom for the
  // current hydro step, set by xray_hook at the start of radiative_transfer
  // (the hook needs the source model, which this library does not link).
  // With on, radiative_physics processes every cell, also the neutral ones
  // without photons. Stateless (recomputed every step), so not in the
  // checkpoint. Serial RT only. Off by default: the original.
  struct Xray {
    bool on = false;
    double gamma[3] = {0, 0, 0};  // primary photoionisation rate, 1/s
    double phi[3] = {0, 0, 0};    // photoelectron power, erg/s
  } xray;
  std::function<void(Simulation&)> xray_hook;
};

}  // namespace radray
