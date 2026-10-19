// pm_backend parameters (pm_backend.hpp); AMReX-free.
#include "radray/pm_backend.hpp"

#include <stdexcept>

#include "radray/rt_backend.hpp"

namespace radray {

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
  if (key == "amrex_pm_max_grid_size") {
    c.amrex_pm_max_grid_size = std::stoi(val);
    if (c.amrex_pm_max_grid_size < 0) throw std::runtime_error("amrex_pm_max_grid_size must be >= 0");
    return true;
  }
  if (key == "amrex_pm_fft") {
    if (val == "pencil")
      c.amrex_pm_fft = 0;
    else if (val == "amrex")
      c.amrex_pm_fft = 1;
    else
      throw std::runtime_error("amrex_pm_fft must be pencil or amrex");
    return true;
  }
  if (key == "pm_gas_levels") {
    c.pm_gas_levels = std::stoi(val);
    if (c.pm_gas_levels != 0 && c.pm_gas_levels != 1) throw std::runtime_error("pm_gas_levels must be 0 or 1");
    return true;
  }
  return false;
}

void copy_pm_backend_param(const Config& from, Config& to) {
  to.pm_backend = from.pm_backend;
  to.amrex_pm_max_grid_size = from.amrex_pm_max_grid_size;
  to.amrex_pm_fft = from.amrex_pm_fft;
  to.pm_gas_levels = from.pm_gas_levels;
}

std::string pm_backend_summary(const Config& c) {
  if (!use_amrex_pm(c)) return "pm_backend = serial";
  const int mgs = c.amrex_pm_max_grid_size > 0 ? c.amrex_pm_max_grid_size : c.amrex_max_grid_size;
  return "pm_backend = amrex (max_grid_size = " + std::to_string(mgs) +
         ", fft = " + (c.amrex_pm_fft ? "amrex" : "pencil") +
         (c.amrex_max_level > 0 ? std::string(", gas levels ") + (c.pm_gas_levels ? "on" : "off (level 0)") : "") + ")";
}

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
  if (c.pm_gas_levels && c.amrex_max_level > 0 && c.amrex_state && (c.pm_refine & (c.pm_refine - 1)) != 0)
    throw std::runtime_error("pm_gas_levels = 1 with refinement needs pm_refine a power of 2 (the level cells "
                             "must align with the mesh cells)");
}

}  // namespace radray
