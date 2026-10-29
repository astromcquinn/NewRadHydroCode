// Port of the RadHydro hydro solver: hydro_tools.f90 compiled with -DMF
// (moving-frame scheme), the gas domain decomposition of domain_tools.f90
// and the hydro parts of radhydro.f90 (calc_timestep, radhydro).
//
// State layout (as the Fortran arrays, component index fastest, cells in
// Fortran order i fastest; see grid.hpp):
//   GasGrid::u  = u_gas(5,N^3): rho, rho*(v - vg) (3), energy in the frame vg
//   GasGrid::w  = w_gas(4,N^3): rho_HI, rho_HeI, rho_HeII, rho_Z (passive)
//   HydroState::v = v_gas(3,N^3): grid (frame) velocity vg
// All in the supercomoving code units of cosmo_tools.f90: calc_sim_var.
//
// See docs/port_notes_hydro.md for the routine map and the quirks of the
// original that are reproduced here.
#pragma once

#include <array>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "radray/config.hpp"
#include "radray/cosmology.hpp"
#include "radray/fortran.hpp"
#include "radray/grid.hpp"
#include "radray/simulation.hpp"

namespace radray {

// radhydro_global.f90 time-stepping parameters
namespace hydrop {
using stepp::ra_max;         // cosmology.hpp
using stepp::dt_gas_factor;
using stepp::dt_pm_factor;
// hydro_tools.f90 module parameters
constexpr double rmf = 8;                // MF velocity smoothing scale [cells]
constexpr double cfl = dt_gas_factor;    // Courant factor of the sweeps
}  // namespace hydrop

// Diagnostic counters (not in the Fortran; they show which branches a
// problem exercises). Counted over the whole run.
struct HydroCounters {
  std::int64_t smooth_cells_1 = 0;   // calls on real(4) state
  std::int64_t smooth_cells_2 = 0;   // calls on real(8) half-step state
  std::int64_t climit_branch = 0;    // c >= c_limit branches (euler, smooth)
  std::int64_t freeze_smooth = 0;    // freeze: c > c_limit
  std::int64_t check_temp_fix = 0;   // check_temp: T < Tmin
  std::int64_t euler_pencils = 0, euler_subcycled = 0, euler_max_nt = 0;
  std::int64_t advect_pencils = 0, advect_subcycled = 0, advect_max_nt = 0;
};

// Counters of the unsplit scheme (hydro_unsplit.hpp), whole run. POD (it
// goes through the checkpoint).
struct UnsplitCounters {
  std::int64_t calls = 0;            // hydro_step_unsplit calls (2 per step)
  std::int64_t substeps = 0;         // CFL substeps
  std::int64_t max_substeps = 0;     // largest number of substeps of one call
  std::int64_t cells_from_e = 0;     // cell-substeps with the pressure from e
  std::int64_t cells_synced = 0;     // cell-substeps with e := E - KE
  std::int64_t cells_reset = 0;      // cell-substeps with E := KE + e
  std::int64_t floor_T = 0;          // temperature floor hits
  std::int64_t floor_rho = 0;        // density floor hits
  double dt_cfl_min = 1e300;         // smallest CFL step seen [code time]
  double frame_shift = 0;            // 1 if a moving-frame state was folded in
};

// Hydro variables of radhydro_global.f90 plus v_gas.
struct HydroState {
  int N = 0;
  std::vector<float> v;  // v_gas(3,N^3)
  double c_limit = 0;                               // set by hydro_step
  double dt_hy = 0;                                 // set by hydro_step
  std::array<double, 2> vmax_last{};  // hydro_step's local vmax (log only)
  HydroCounters counters;
  UnsplitCounters ucounters;  // hydro_scheme = unsplit only

  void allocate(int n) {
    N = n;
    v.assign(static_cast<std::size_t>(3) * n * n * n, 0.0f);
  }
  float* v_p(std::int64_t c) { return &v[c * 3]; }
};

// hydro_tools.f90: init_hydro (allocate and zero u_gas, v_gas, w_gas)
void init_hydro(const Config& cfg, GasGrid& gas, HydroState& hs);

enum class SweepOrder { Forward, Reverse };

// hydro_tools.f90: hydro_step(flag). One sweep along the first index for
// every pencil (j,k): smooth, euler, advect, freeze (forward) or smooth,
// advect, euler, freeze (reverse). dt_hy = dt_sim/2.
void hydro_step(const Config& cfg, const Units& un, const SimClock& clk,
                SweepOrder order, GasGrid& gas, HydroState& hs);

// hydro_tools.f90: transpose_ij, transpose_ik (in-place transposes of u_gas,
// v_gas and w_gas with the matching permutation of vector components).
void transpose_ij(GasGrid& gas, HydroState& hs);
void transpose_ik(GasGrid& gas, HydroState& hs);

// hydro_tools.f90: summarize_gas. Sets clk.cmax_gas, vmax_gas, vrms_gas and
// (side effects of the original) clamps u_gas(1) >= Dgas_min, w_gas >= 0.
struct GasSummary {
  // the unit-100 record "(50es13.5)": redshift, fHI v/m, fHeI v/m,
  // fHeII v/m, fe v/m, T v/m, Z v/m, clump(2:4,6:8,10:12,14:16,18:20)
  std::vector<double> record;
  // standard output lines: rho, vx, vy, vz, cs, T, HI, HeI, HeII, e, Z
  // (each: volume-weighted, "mass-weighted", min, max as printed)
  std::array<std::array<double, 4>, 11> table{};
};
GasSummary summarize_gas(const Config& cfg, const Units& un, SimClock& clk,
                         GasGrid& gas, HydroState& hs, std::ostream* out = nullptr);

// summarize_gas in parts (docs/rank0_streaming.md): the planes of the gas
// grid in k order, in runs of any length, then finish. summarize_gas runs
// it over the host arrays; with amrex_state = 1 rank 0 runs it over the
// slabs streamed from the level. add_planes takes nk consecutive planes,
// cells in the serial order: u (5 per cell) and w (4 per cell) are clamped
// in place as in summarize_gas; e (1 per cell) is the unsplit hydro's
// internal energy or null; v (3 per cell) the grid velocity or null (zero).
class GasSummer {
 public:
  GasSummer(const Config& cfg, const Units& un);
  void add_planes(int nk, float* u, float* w, const float* e, const float* v);
  GasSummary finish(SimClock& clk, std::ostream* out);

 private:
  struct Rec {
    double rho, vel[3], cs, tmp, fHI, fHeI, fHeII, fe, Z, nH, nHII, ne, alpha, x, cmax;
  };
  void acc(int m, double val, double rho) {
    x1[m] = x1[m] + val;
    x2[m] = x2[m] + val * rho;
    x3[m] = f::gmin(x3[m], val);
    x4[m] = f::gmax(x4[m], val);
  }
  const Config& cfg_;
  const Units& un_;
  double x1[11], x2[11], x3[11], x4[11];
  double cmax = 0, vrms = 0;
  double clump[21];  // 1-based
  std::vector<Rec> recs;
};

// Fortran "(50es13.5)" formatting of a record (gfortran style).
std::string format_es(double x, int w, int d);
std::string format_gas_record(const std::vector<double>& rec);

// radhydro.f90: calc_timestep is calc_timestep in cosmology.hpp (all drivers).

// radhydro.f90: radhydro without particle_mesh and radiative transfer:
// transposes + 3 forward sweeps, half clock step (calc_sim_var), 3 reverse
// sweeps + transposes, half clock step (calc_sim_var). With
// cfg.hydro_scheme = Unsplit each group of sweeps is one unsplit step over
// dt_sim/2 (hydro_half_step in hydro_unsplit.hpp).
void radhydro_hydro_step(const Config& cfg, SimClock& clk, Units& un, GasGrid& gas,
                         HydroState& hs);

// domain_tools.f90: init_domain_gas / set_domain_gas / end_domain_gas
// (same structure as the RT domains; Domains::box holds indx_dom_gas).
void init_domain_gas(const Config& cfg, Domains& dom, std::ostream* out = nullptr);
int set_domain_gas(Domains& dom);             // returns 0-based domain index
void end_domain_gas(Domains& dom, int idom);  // idom 0-based

}  // namespace radray
