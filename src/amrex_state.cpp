// amrex_state parameters and checks (amrex_state.hpp); AMReX-free.
#include "radray/amrex_state.hpp"

#include <stdexcept>

#include "radray/hydro_backend.hpp"
#include "radray/pm_backend.hpp"
#include "radray/recomb_rad.hpp"
#include "radray/rt_backend.hpp"

namespace radray {

bool set_amrex_state_param(Config& c, const std::string& key, const std::string& val) {
  if (key != "amrex_state") return false;
  if (val == "1")
    c.amrex_state = true;
  else if (val == "0")
    c.amrex_state = false;
  else
    throw std::runtime_error("amrex_state must be 0 or 1");
  return true;
}

void copy_amrex_state_param(const Config& from, Config& to) { to.amrex_state = from.amrex_state; }

std::string amrex_state_summary(const Config& c) {
  return c.amrex_state ? "amrex_state = 1 (the gas on the levels, host mirror on rank 0)" : "amrex_state = 0";
}

void check_amrex_state_param(const Config& c) {
  if (!c.amrex_state) return;
  if (!amrex_available())
    throw std::runtime_error("amrex_state = 1: this build has no AMReX (cmake -DRADRAY_AMREX=ON)");
  if (!use_amrex(c) || !use_amrex_hydro(c) || !use_amrex_pm(c))
    throw std::runtime_error("amrex_state = 1 needs rt_backend = amrex, hydro_backend = amrex and pm_backend = amrex");
  if (c.hydro_scheme != HydroScheme::Unsplit)
    throw std::runtime_error("amrex_state = 1 needs hydro_scheme = unsplit");
  if (c.rt_method != RTMethod::Rays)
    throw std::runtime_error("amrex_state = 1 needs rt_method = rays");
  if (c.recomb_rad) throw std::runtime_error("amrex_state = 1 does not support recomb_rad = 1");
  if (c.N_rt_gas != 1) throw std::runtime_error("amrex_state = 1 needs N_rt_gas = 1 (Nrt = Ngas)");
  // one decomposition: the hydro and PM levels on the RT level's boxes
  if (c.amrex_hydro_max_grid_size != 0 && c.amrex_hydro_max_grid_size != c.amrex_max_grid_size)
    throw std::runtime_error("amrex_state = 1: amrex_hydro_max_grid_size must be 0 or equal to amrex_max_grid_size");
  if (c.amrex_pm_max_grid_size != 0 && c.amrex_pm_max_grid_size != c.pm_refine * c.amrex_max_grid_size)
    throw std::runtime_error("amrex_state = 1: amrex_pm_max_grid_size must be 0 or pm_refine * amrex_max_grid_size");
  if (c.Ngas % c.amrex_max_grid_size != 0)
    throw std::runtime_error("amrex_state = 1: amrex_max_grid_size must divide N");
}

}  // namespace radray
