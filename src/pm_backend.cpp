// The parameters of the PM backend switch (include/radray/pm_backend.hpp,
// docs/stage_e_pm.md): pm_backend = serial (the original's linked lists,
// float CIC and FFTW on the host, src/pm.cpp) | amrex (the particles in an
// AMReX ParticleContainer, fixed-point CIC deposits and the Poisson solve on
// a level of the PM mesh, src/amrex/pm_level.cpp, MPI + OpenMP). Also the
// parameters of the gravity on a refined gas hierarchy: pm_gas_levels (the
// uniform mesh fed by every level, docs/stage_f.md), Stage G's level Poisson
// solve pm_level_solve with its MLMG controls pm_mlmg_* (docs/stage_g.md)
// and the sub-cycled gravity of the refined levels pm_level_subcycle
// (docs/stage_g.md section 6), and amrex_pk, the power spectrum of the
// outputs on the levels (docs/rank0_streaming.md step 11). This file is
// AMReX-free and part of every build: set_pm_backend_param is called by the
// pm_driver and radhydro parameter parsers, the copy / summary / check
// functions by the drivers after the parameters are read. The meaning of
// every parameter is documented at its Config field
// (include/radray/config.hpp); the level's entry points live in
// src/amrex/pm_level.cpp or its stub.
#include "radray/pm_backend.hpp"

#include <stdexcept>

#include "radray/rt_backend.hpp"

namespace radray {

// One params.txt entry: true if key belongs to this group (the value is
// validated and stored, a bad value throws), false for any other key.
bool set_pm_backend_param(Config& c, const std::string& key, const std::string& val) {
  if (key == "pm_backend") {
    if (val == "serial")
      c.pm_backend = PMBackend::Serial;
    else if (val == "amrex")
      c.pm_backend = PMBackend::AMReX;
    else
      throw std::runtime_error("pm_backend must be serial or amrex");
    return true;
  }
  // box size of the PM mesh's decomposition; 0 = amrex_max_grid_size (with
  // amrex_state = 1 the mesh boxes are pm_refine times the gas boxes)
  if (key == "amrex_pm_max_grid_size") {
    c.amrex_pm_max_grid_size = std::stoi(val);
    if (c.amrex_pm_max_grid_size < 0) throw std::runtime_error("amrex_pm_max_grid_size must be >= 0");
    return true;
  }
  // the FFT of the Poisson solve on level 0: pencil (default; 1-D FFTW
  // transforms over the rank's pencils, docs/stage_e_pm.md) or AMReX's FFT
  if (key == "amrex_pm_fft") {
    if (val == "pencil")
      c.amrex_pm_fft = 0;
    else if (val == "amrex")
      c.amrex_pm_fft = 1;
    else
      throw std::runtime_error("amrex_pm_fft must be pencil or amrex");
    return true;
  }
  // the power spectrum of the data1 outputs with amrex_state = 1: on rank 0's
  // host PM mesh (host) or on the levels with exact bin sums (level)
  if (key == "amrex_pk") {
    if (val == "host")
      c.amrex_pk = 0;
    else if (val == "level")
      c.amrex_pk = 1;
    else
      throw std::runtime_error("amrex_pk must be host or level");
    return true;
  }
  // with a refined hierarchy (amrex_max_level > 0, amrex_state = 1): the gas
  // of every level (1) or of level 0 only (0) enters the uniform PM mesh
  if (key == "pm_gas_levels") {
    c.pm_gas_levels = std::stoi(val);
    if (c.pm_gas_levels != 0 && c.pm_gas_levels != 1) throw std::runtime_error("pm_gas_levels must be 0 or 1");
    return true;
  }
  // Stage G (docs/stage_g.md): the Poisson solve on the refined levels by
  // MLMG, with the parent level's potential as the boundary values and the
  // initial guess; pm_mlmg_iters V-cycles at most, stopping when the
  // residual is below pm_mlmg_tol times the initial one (0 = always the
  // full count); the bottom smoother runs at least pm_mlmg_bottom_sweeps
  // sweeps and, on a rectangular level, at least pm_mlmg_bottom_scale times
  // the square of its bottom grid's width
  if (key == "pm_level_solve") {
    c.pm_level_solve = std::stoi(val);
    if (c.pm_level_solve != 0 && c.pm_level_solve != 1) throw std::runtime_error("pm_level_solve must be 0 or 1");
    return true;
  }
  if (key == "pm_mlmg_iters") {
    c.pm_mlmg_iters = std::stoi(val);
    if (c.pm_mlmg_iters < 1) throw std::runtime_error("pm_mlmg_iters must be >= 1");
    return true;
  }
  if (key == "pm_mlmg_tol") {
    c.pm_mlmg_tol = std::stod(val);
    if (c.pm_mlmg_tol < 0) throw std::runtime_error("pm_mlmg_tol must be >= 0");
    return true;
  }
  if (key == "pm_mlmg_bottom_sweeps") {
    c.pm_mlmg_bottom_sweeps = std::stoi(val);
    if (c.pm_mlmg_bottom_sweeps < 1) throw std::runtime_error("pm_mlmg_bottom_sweeps must be >= 1");
    return true;
  }
  if (key == "pm_mlmg_bottom_scale") {
    c.pm_mlmg_bottom_scale = std::stod(val);
    if (!(c.pm_mlmg_bottom_scale >= 0)) throw std::runtime_error("pm_mlmg_bottom_scale must be >= 0");
    return true;
  }
  // docs/stage_g.md section 6: with pm_level_solve the refined levels' gravity
  // runs on their own (hydro) timescale inside the hydro calls
  if (key == "pm_level_subcycle") {
    c.pm_level_subcycle = std::stoi(val);
    if (c.pm_level_subcycle != 0 && c.pm_level_subcycle != 1) throw std::runtime_error("pm_level_subcycle must be 0 or 1");
    return true;
  }
  // 1 prints the MLMG residuals on rank 0; higher values print more
  if (key == "pm_mlmg_verbose") {
    c.pm_mlmg_verbose = std::stoi(val);
    return true;
  }
  return false;
}

// Every field of the group, from the parsed Config into the simulation's.
void copy_pm_backend_param(const Config& from, Config& to) {
  to.pm_backend = from.pm_backend;
  to.amrex_pm_max_grid_size = from.amrex_pm_max_grid_size;
  to.amrex_pm_fft = from.amrex_pm_fft;
  to.amrex_pk = from.amrex_pk;
  to.pm_gas_levels = from.pm_gas_levels;
  to.pm_level_solve = from.pm_level_solve;
  to.pm_mlmg_iters = from.pm_mlmg_iters;
  to.pm_mlmg_tol = from.pm_mlmg_tol;
  to.pm_mlmg_bottom_sweeps = from.pm_mlmg_bottom_sweeps;
  to.pm_mlmg_bottom_scale = from.pm_mlmg_bottom_scale;
  to.pm_mlmg_verbose = from.pm_mlmg_verbose;
  to.pm_level_subcycle = from.pm_level_subcycle;
}

// One line for the run header (pm_driver, radhydro); the hierarchy part
// only when there are refined levels, as the parameters only act then.
std::string pm_backend_summary(const Config& c) {
  if (!use_amrex_pm(c)) return "pm_backend = serial";
  const int mgs = c.amrex_pm_max_grid_size > 0 ? c.amrex_pm_max_grid_size : c.amrex_max_grid_size;
  return "pm_backend = amrex (max_grid_size = " + std::to_string(mgs) +
         ", fft = " + (c.amrex_pm_fft ? "amrex" : "pencil") +
         (c.amrex_max_level > 0
              ? (c.pm_level_solve ? std::string(", level solve (mlmg iters ") + std::to_string(c.pm_mlmg_iters) +
                                    (c.pm_level_subcycle ? ", sub-cycled" : "") + ")"
                                  : std::string(", gas levels ") + (c.pm_gas_levels ? "on" : "off (level 0)"))
              : "") +
         ")";
}

// Refuses what the backend does not implement; called once by the drivers
// after all parameters are set. The box checks keep every gas cell's
// pm_refine^3 mesh cells inside one box and the mesh a whole number of
// boxes.
void check_pm_backend_param(const Config& c) {
  if (!use_amrex_pm(c)) return;
  if (!amrex_available())
    throw std::runtime_error("pm_backend = amrex: this build has no AMReX (cmake -DRADRAY_AMREX=ON)");
  const int mgs = c.amrex_pm_max_grid_size > 0 ? c.amrex_pm_max_grid_size : c.amrex_max_grid_size;
  if (mgs < 4) throw std::runtime_error("pm_backend = amrex: the PM mesh box size must be >= 4");
  if (c.pm_refine > 1 && mgs % c.pm_refine != 0)
    throw std::runtime_error("pm_backend = amrex: amrex_pm_max_grid_size must be a multiple of pm_refine "
                             "(every gas cell's mesh cells lie in one box)");
  if (c.Ngas > 0 && (c.Ngas * c.pm_refine) % mgs != 0)
    throw std::runtime_error("pm_backend = amrex: the PM mesh size must be a multiple of the box size");
  if (c.pm_gas_levels && !c.pm_level_solve && c.amrex_max_level > 0 && c.amrex_state &&
      (c.pm_refine & (c.pm_refine - 1)) != 0)
    throw std::runtime_error("pm_gas_levels = 1 with refinement needs pm_refine a power of 2 (the level cells "
                             "must align with the mesh cells)");
  if (c.pm_level_solve && c.amrex_max_level > 0 && !c.amrex_state)
    throw std::runtime_error("pm_level_solve = 1 with refinement needs amrex_state = 1 (the particles and the "
                             "gas on the levels)");
  // the sub-cycle is driven from the hydro hierarchy's calls (pm_subcycle.hpp),
  // on the gas cells themselves
  if (c.pm_level_subcycle && c.pm_level_solve && c.amrex_max_level > 0) {
    if (c.pm_refine != 1) throw std::runtime_error("pm_level_subcycle = 1 needs pm_refine = 1");
    if (c.hydro_scheme != HydroScheme::Unsplit || !use_amrex_hydro(c))
      throw std::runtime_error("pm_level_subcycle = 1 needs hydro_scheme = unsplit and hydro_backend = amrex");
  }
}

}  // namespace radray
