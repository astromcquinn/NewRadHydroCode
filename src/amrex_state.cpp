// amrex_state parameters and checks (amrex_state.hpp); AMReX-free.
#include "radray/amrex_state.hpp"

#include <stdexcept>

#include "radray/hydro_backend.hpp"
#include "radray/pm_backend.hpp"
#include "radray/recomb_rad.hpp"
#include "radray/rt_backend.hpp"

namespace radray {

bool set_amrex_state_param(Config& c, const std::string& key, const std::string& val) {
  auto flag = [&](bool& f) {
    if (val == "1")
      f = true;
    else if (val == "0")
      f = false;
    else
      throw std::runtime_error(key + " must be 0 or 1");
  };
  if (key == "amrex_state") {
    flag(c.amrex_state);
    return true;
  }
  if (key == "amrex_balance") {
    flag(c.amrex_balance);
    return true;
  }
  if (key == "amrex_balance_ray_weight") {
    c.amrex_balance_ray_weight = std::stod(val);
    if (!(c.amrex_balance_ray_weight >= 0)) throw std::runtime_error("amrex_balance_ray_weight must be >= 0");
    return true;
  }
  if (key == "amrex_balance_min_gain") {
    c.amrex_balance_min_gain = std::stod(val);  // negative: remap at every RT step (tests)
    if (!(c.amrex_balance_min_gain >= -1)) throw std::runtime_error("amrex_balance_min_gain must be >= -1");
    return true;
  }
  if (key == "amrex_balance_method") {
    if (val == "sfc")
      c.amrex_balance_method = 0;
    else if (val == "knapsack")
      c.amrex_balance_method = 1;
    else
      throw std::runtime_error("amrex_balance_method must be sfc or knapsack");
    return true;
  }
  return false;
}

void copy_amrex_state_param(const Config& from, Config& to) {
  to.amrex_state = from.amrex_state;
  to.amrex_balance = from.amrex_balance;
  to.amrex_balance_ray_weight = from.amrex_balance_ray_weight;
  to.amrex_balance_min_gain = from.amrex_balance_min_gain;
  to.amrex_balance_method = from.amrex_balance_method;
}

std::string amrex_state_summary(const Config& c) {
  if (!c.amrex_state) return "amrex_state = 0";
  std::string s = "amrex_state = 1 (the gas on the levels, host mirror on rank 0)";
  if (c.amrex_balance)
    s += ", amrex_balance = 1 (" + std::string(c.amrex_balance_method ? "knapsack" : "sfc") +
         ", ray weight " + std::to_string(c.amrex_balance_ray_weight) + ", min gain " +
         std::to_string(c.amrex_balance_min_gain) + ")";
  return s;
}

void check_amrex_state_param(const Config& c) {
  if (c.amrex_balance && !c.amrex_state) throw std::runtime_error("amrex_balance = 1 needs amrex_state = 1");
  if (c.amrex_balance && c.amrex_slab)
    throw std::runtime_error("amrex_balance = 1 does not support amrex_slab = 1 (the slabs keep their z order)");
  if (!c.amrex_state) return;
  if (!amrex_available())
    throw std::runtime_error("amrex_state = 1: this build has no AMReX (cmake -DRADRAY_AMREX=ON)");
  if (!use_amrex(c) || !use_amrex_hydro(c) || !use_amrex_pm(c))
    throw std::runtime_error("amrex_state = 1 needs rt_backend = amrex, hydro_backend = amrex and pm_backend = amrex");
  if (c.hydro_scheme != HydroScheme::Unsplit)
    throw std::runtime_error("amrex_state = 1 needs hydro_scheme = unsplit");
  if (c.rt_method != RTMethod::Rays && !c.amrex_slab)
    throw std::runtime_error("amrex_state = 1 with rt_method = m1 | otvet needs amrex_slab = 1");
  if (c.recomb_rad && c.rt_method == RTMethod::Rays)
    throw std::runtime_error("amrex_state = 1 does not support recomb_rad = 1 with the rays");
  if (c.N_rt_gas != 1) throw std::runtime_error("amrex_state = 1 needs N_rt_gas = 1 (Nrt = Ngas)");
  // one decomposition: the hydro and PM levels on the RT level's boxes
  if (c.amrex_hydro_max_grid_size != 0 && c.amrex_hydro_max_grid_size != c.amrex_max_grid_size)
    throw std::runtime_error("amrex_state = 1: amrex_hydro_max_grid_size must be 0 or equal to amrex_max_grid_size");
  if (c.amrex_pm_max_grid_size != 0 && c.amrex_pm_max_grid_size != c.pm_refine * c.amrex_max_grid_size)
    throw std::runtime_error("amrex_state = 1: amrex_pm_max_grid_size must be 0 or pm_refine * amrex_max_grid_size");
  if (!c.amrex_slab && c.Ngas % c.amrex_max_grid_size != 0)
    throw std::runtime_error("amrex_state = 1: amrex_max_grid_size must divide N");
}

}  // namespace radray
