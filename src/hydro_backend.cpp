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
  return false;
}

void copy_hydro_backend_param(const Config& from, Config& to) {
  to.hydro_backend = from.hydro_backend;
  to.amrex_hydro_max_grid_size = from.amrex_hydro_max_grid_size;
}

std::string hydro_backend_summary(const Config& c) {
  if (!use_amrex_hydro(c)) return "hydro_backend = serial";
  return "hydro_backend = amrex (max_grid_size = " +
         std::to_string(c.amrex_hydro_max_grid_size ? c.amrex_hydro_max_grid_size : c.amrex_max_grid_size) + ")";
}

void check_hydro_backend_param(const Config& c) {
  if (!use_amrex_hydro(c)) return;
  if (!amrex_available())
    throw std::runtime_error("hydro_backend = amrex: this build has no AMReX (cmake -DRADRAY_AMREX=ON)");
  if (c.hydro_scheme != HydroScheme::Unsplit)
    throw std::runtime_error("hydro_backend = amrex needs hydro_scheme = unsplit (the legacy split scheme "
                             "with its pencil FFT stays on the serial code)");
}

}  // namespace radray
