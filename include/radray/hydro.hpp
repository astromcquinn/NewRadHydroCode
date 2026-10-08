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
#include "radray/grid.hpp"
#include "radray/simulation.hpp"

namespace radray {

// radhydro_global.f90 time-stepping parameters
namespace hydrop {
constexpr double ra_max = static_cast<double>(0.1f);  // real(8) from real(4) literal
constexpr double dt_gas_factor = 0.5;
constexpr double dt_pm_factor = 0.5;
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

// Hydro variables of radhydro_global.f90 plus v_gas.
struct HydroState {
  int N = 0;
  std::vector<float> v;  // v_gas(3,N^3)
  double cmax_gas = 0, vmax_gas = 0, vrms_gas = 0;  // set by summarize_gas
  double c_limit = 0;                               // set by hydro_step
  double dt_hy = 0;                                 // set by hydro_step
  std::array<double, 2> vmax_last{};  // hydro_step's local vmax (log only)
  HydroCounters counters;

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

// hydro_tools.f90: summarize_gas. Sets hs.cmax_gas, vmax_gas, vrms_gas and
// (side effects of the original) clamps u_gas(1) >= Dgas_min, w_gas >= 0.
struct GasSummary {
  // the unit-100 record "(50es13.5)": redshift, fHI v/m, fHeI v/m,
  // fHeII v/m, fe v/m, T v/m, Z v/m, clump(2:4,6:8,10:12,14:16,18:20)
  std::vector<double> record;
  // standard output lines: rho, vx, vy, vz, cs, T, HI, HeI, HeII, e, Z
  // (each: volume-weighted, "mass-weighted", min, max as printed)
  std::array<std::array<double, 4>, 11> table{};
};
GasSummary summarize_gas(const Config& cfg, const Units& un, const SimClock& clk,
                         GasGrid& gas, HydroState& hs, std::ostream* out = nullptr);

// Fortran "(50es13.5)" formatting of a record (gfortran style).
std::string format_es(double x, int w, int d);
std::string format_gas_record(const std::vector<double>& rec);

// radhydro.f90: calc_timestep, hydro part.
//   Fixed   : da_sim = da_fixed, dt from calc_dt_from_da (istep == 0 form)
//   Courant : the istep > 0 branch, dt_sim = min(2*dt_gas_factor/cmax_gas,
//             dt_pm_factor/sqrt(amax_dm)), da from calc_da_from_dt, limited
//             by ra_max, dt from calc_dt_from_da.
// Increments istep and converts vmax_gas, vrms_gas to km/s (as the
// original does, for printing). Returns dt_vg in years (log).
enum class TimestepMode { Fixed, Courant };
double calc_timestep_hydro(const Units& un, TimestepMode mode, double da_fixed,
                           double amax_dm, SimClock& clk, HydroState& hs,
                           std::ostream* out = nullptr);

// radhydro.f90: radhydro without particle_mesh and radiative transfer:
// transposes + 3 forward sweeps, half clock step (calc_sim_var), 3 reverse
// sweeps + transposes, half clock step (calc_sim_var).
void radhydro_hydro_step(const Config& cfg, SimClock& clk, Units& un, GasGrid& gas,
                         HydroState& hs);

// domain_tools.f90: init_domain_gas / set_domain_gas / end_domain_gas
// (same structure as the RT domains; Domains::box holds indx_dom_gas).
void init_domain_gas(const Config& cfg, Domains& dom, std::ostream* out = nullptr);
int set_domain_gas(Domains& dom);             // returns 0-based domain index
void end_domain_gas(Domains& dom, int idom);  // idom 0-based

}  // namespace radray
