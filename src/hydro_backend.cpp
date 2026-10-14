// hydro_backend parameters (hydro_backend.hpp); AMReX-free.
#include "radray/hydro_backend.hpp"

#include <stdexcept>

#include "radray/rt_backend.hpp"

namespace radray {

bool set_hydro_backend_param(Config& c, const std::string& key, const std::string& val) {
  if (key == "hydro_backend") {
    if (val == "serial")
      c.hydro_backend = HydroBackend::Serial;
    else if (val == "amrex")
      c.hydro_backend = HydroBackend::AMReX;
    else
      throw std::runtime_error("hydro_backend must be serial or amrex");
    return true;
  }
  if (key == "amrex_hydro_max_grid_size") {
    c.amrex_hydro_max_grid_size = std::stoi(val);
    if (c.amrex_hydro_max_grid_size != 0 && c.amrex_hydro_max_grid_size < 4)
      throw std::runtime_error("amrex_hydro_max_grid_size must be >= 4 (or 0: amrex_max_grid_size)");
    return true;
  }
  // Stage F: the hierarchy (docs/stage_f.md section 1)
  if (key == "amrex_max_level") {
    c.amrex_max_level = std::stoi(val);
    if (c.amrex_max_level < 0 || c.amrex_max_level > 8) throw std::runtime_error("amrex_max_level must be 0..8");
    return true;
  }
  if (key == "amrex_refine_mass") {
    c.amrex_refine_mass = std::stod(val);
    if (!(c.amrex_refine_mass > 0)) throw std::runtime_error("amrex_refine_mass must be > 0");
    return true;
  }
  if (key == "amrex_regrid_every") {
    c.amrex_regrid_every = std::stoi(val);
    if (c.amrex_regrid_every < 1) throw std::runtime_error("amrex_regrid_every must be >= 1");
    return true;
  }
  if (key == "amrex_fine_max_grid_size") {
    c.amrex_fine_max_grid_size = std::stoi(val);
    if (c.amrex_fine_max_grid_size < 4) throw std::runtime_error("amrex_fine_max_grid_size must be >= 4");
    return true;
  }
  if (key == "amrex_blocking_factor") {
    c.amrex_blocking_factor = std::stoi(val);
    if (c.amrex_blocking_factor < 2) throw std::runtime_error("amrex_blocking_factor must be >= 2");
    return true;
  }
  if (key == "amrex_n_error_buf") {
    c.amrex_n_error_buf = std::stoi(val);
    if (c.amrex_n_error_buf < 0) throw std::runtime_error("amrex_n_error_buf must be >= 0");
    return true;
  }
  if (key == "amrex_grid_eff") {
    c.amrex_grid_eff = std::stod(val);
    if (!(c.amrex_grid_eff > 0 && c.amrex_grid_eff <= 1)) throw std::runtime_error("amrex_grid_eff must be in (0, 1]");
    return true;
  }
  return false;
}

void copy_hydro_backend_param(const Config& from, Config& to) {
  to.hydro_backend = from.hydro_backend;
  to.amrex_hydro_max_grid_size = from.amrex_hydro_max_grid_size;
  to.amrex_max_level = from.amrex_max_level;
  to.amrex_refine_mass = from.amrex_refine_mass;
  to.amrex_regrid_every = from.amrex_regrid_every;
  to.amrex_fine_max_grid_size = from.amrex_fine_max_grid_size;
  to.amrex_blocking_factor = from.amrex_blocking_factor;
  to.amrex_n_error_buf = from.amrex_n_error_buf;
  to.amrex_grid_eff = from.amrex_grid_eff;
}

std::string hydro_backend_summary(const Config& c) {
  if (!use_amrex_hydro(c)) return "hydro_backend = serial";
  std::string s = "hydro_backend = amrex (max_grid_size = " +
                  std::to_string(c.amrex_hydro_max_grid_size ? c.amrex_hydro_max_grid_size : c.amrex_max_grid_size);
  if (c.amrex_max_level > 0)
    s += ", max_level = " + std::to_string(c.amrex_max_level) + ", refine_mass = " + std::to_string(c.amrex_refine_mass) +
         ", regrid_every = " + std::to_string(c.amrex_regrid_every) + ", fine boxes " +
         std::to_string(c.amrex_fine_max_grid_size) + "/" + std::to_string(c.amrex_blocking_factor);
  return s + ")";
}

void check_hydro_backend_param(const Config& c) {
  if (!use_amrex_hydro(c)) return;
  if (!amrex_available())
    throw std::runtime_error("hydro_backend = amrex: this build has no AMReX (cmake -DRADRAY_AMREX=ON)");
  if (c.hydro_scheme != HydroScheme::Unsplit)
    throw std::runtime_error("hydro_backend = amrex needs hydro_scheme = unsplit (the legacy split scheme "
                             "with its pencil FFT stays on the serial code)");
  if (c.amrex_max_level > 0) {
    if (c.amrex_fine_max_grid_size % c.amrex_blocking_factor != 0)
      throw std::runtime_error("amrex_fine_max_grid_size must be a multiple of amrex_blocking_factor");
    const int mgs0 = c.amrex_hydro_max_grid_size ? c.amrex_hydro_max_grid_size : c.amrex_max_grid_size;
    if (mgs0 % c.amrex_blocking_factor != 0)
      throw std::runtime_error("amrex_max_level > 0: the level-0 box size must be a multiple of amrex_blocking_factor");
  }
}

}  // namespace radray
