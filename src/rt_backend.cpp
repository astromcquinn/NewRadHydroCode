// rt_backend parameters (rt_backend.hpp); AMReX-free.
#include "radray/rt_backend.hpp"

#include <algorithm>
#include <stdexcept>

namespace radray {

bool set_rt_backend_param(Config& c, const std::string& key, const std::string& val) {
  if (key == "rt_backend") {
    if (val == "serial")
      c.rt_backend = RTBackend::Serial;
    else if (val == "amrex")
      c.rt_backend = RTBackend::AMReX;
    else
      throw std::runtime_error("rt_backend must be serial or amrex");
    return true;
  }
  if (key == "ray_encoding") {
    if (val == "float")
      c.ray_encoding = RayEncoding::Float;
    else if (val == "rstruct")
      c.ray_encoding = RayEncoding::RStruct;
    else
      throw std::runtime_error("ray_encoding must be float or rstruct");
    return true;
  }
  if (key == "amrex_max_grid_size") {
    c.amrex_max_grid_size = std::stoi(val);
    if (c.amrex_max_grid_size < 4) throw std::runtime_error("amrex_max_grid_size must be >= 4");
    return true;
  }
  if (key == "amrex_ghosts") {
    c.amrex_ghosts = std::stoi(val);
    if (c.amrex_ghosts < 2) throw std::runtime_error("amrex_ghosts must be >= 2 (the ray reach)");
    return true;
  }
  if (key == "amrex_slab") {
    if (val == "1")
      c.amrex_slab = true;
    else if (val == "0")
      c.amrex_slab = false;
    else
      throw std::runtime_error("amrex_slab must be 0 or 1");
    return true;
  }
  if (key == "amrex_moment_stats") {
    if (val == "exact")
      c.amrex_moment_stats = 0;
    else if (val == "fixed")
      c.amrex_moment_stats = 1;
    else
      throw std::runtime_error("amrex_moment_stats must be exact or fixed");
    return true;
  }
  if (key == "amrex_moment_fft_threads") {
    c.amrex_moment_fft_threads = std::stoi(val);
    if (c.amrex_moment_fft_threads < 0) throw std::runtime_error("amrex_moment_fft_threads must be >= 0");
    return true;
  }
  if (key == "amrex_slab_planes") {
    c.amrex_slab_planes = std::stoi(val);
    if (c.amrex_slab_planes < 0 || c.amrex_slab_planes % 2 != 0)
      throw std::runtime_error("amrex_slab_planes must be 0 (about N / ranks) or a positive even number");
    return true;
  }
  if (key == "amrex_split_on_entry") {
    c.amrex_split_on_entry = std::stoi(val);
    if (c.amrex_split_on_entry != 0 && c.amrex_split_on_entry != 1)
      throw std::runtime_error("amrex_split_on_entry must be 0 or 1");
    return true;
  }
  if (key == "amrex_ghost_rule") {
    if (val == "sqrt" || val == "0") c.amrex_ghost_rule = 0;
    else if (val == "parent" || val == "1") c.amrex_ghost_rule = 1;
    else if (val == "half" || val == "2") c.amrex_ghost_rule = 2;
    else if (val == "conserve" || val == "3") c.amrex_ghost_rule = 3;
    else throw std::runtime_error("amrex_ghost_rule must be sqrt, parent, half or conserve");
    return true;
  }
  if (key == "amrex_redist_local") {
    c.amrex_redist_local = std::stoi(val);
    if (c.amrex_redist_local < 0) throw std::runtime_error("amrex_redist_local must be >= 0");
    return true;
  }
  return false;
}

void copy_rt_backend_param(const Config& from, Config& to) {
  to.rt_backend = from.rt_backend;
  to.ray_encoding = from.ray_encoding;
  to.amrex_max_grid_size = from.amrex_max_grid_size;
  to.amrex_ghosts = from.amrex_ghosts;
  to.amrex_redist_local = from.amrex_redist_local;
  to.amrex_slab = from.amrex_slab;
  to.amrex_slab_planes = from.amrex_slab_planes;
  to.amrex_moment_stats = from.amrex_moment_stats;
  to.amrex_moment_fft_threads = from.amrex_moment_fft_threads;
  to.amrex_split_on_entry = from.amrex_split_on_entry;
  to.amrex_ghost_rule = from.amrex_ghost_rule;
}

std::string rt_backend_summary(const Config& c) {
  if (!use_amrex(c)) return "rt_backend = serial";
  // the values the level uses: at least the reach k+1 of one substep for
  // rt_cells_per_substep = k (src/amrex/rt_level.cpp)
  const int reach = c.rt_cells_per_substep + 1;
  return "rt_backend = amrex (ray_encoding = " +
         std::string(c.ray_encoding == RayEncoding::RStruct ? "rstruct" : "float") +
         (c.amrex_slab ? ", z-slabs of " +
                             (c.amrex_slab_planes ? std::to_string(c.amrex_slab_planes) : std::string("N / ranks")) +
                             " planes"
                       : ", max_grid_size = " + std::to_string(c.amrex_max_grid_size)) +
         (c.rt_method != RTMethod::Rays
              ? std::string(", moment stats = ") + (c.amrex_moment_stats ? "fixed" : "exact")
              : ", ghosts = " + std::to_string(std::max(c.amrex_ghosts, reach))) +
         ", redist_local = " +
         std::to_string(c.amrex_redist_local > 0 ? std::max(c.amrex_redist_local, reach) : 0) + ")";
}

void check_rt_backend_param(const Config& c) {
  if (!use_amrex(c)) {
    if (c.ray_encoding != RayEncoding::Float)
      throw std::runtime_error("ray_encoding = rstruct needs rt_backend = amrex");
    return;
  }
  if (!amrex_available())
    throw std::runtime_error("rt_backend = amrex: this build has no AMReX (cmake -DRADRAY_AMREX=ON)");
  if (c.legacy_rt)
    throw std::runtime_error(
        "rt_backend = amrex needs legacy_rt = 0: sorted ray arrays cannot reproduce the "
        "original's zombie-ray thresholds (docs/design_parallel_and_encoding.md section 0.1)");
  if (c.rt_method != RTMethod::Rays && !c.amrex_slab)
    throw std::runtime_error("rt_method = m1 | otvet with rt_backend = amrex needs amrex_slab = 1 (z-slabs)");
  if (c.recomb_rad && c.rt_method == RTMethod::Rays)
    throw std::runtime_error("rt_backend = amrex does not implement recomb_rad with the rays yet");
  if (c.det_parallel && c.rt_method == RTMethod::Rays)
    throw std::runtime_error("rt_backend = amrex: deterministic_parallel applies to the serial backend only");
  if (c.N_rt_gas != 1) throw std::runtime_error("rt_backend = amrex needs N_rt_gas = 1");
}

int amrex_slab_planes(const Config& c, int N, int nranks) {
  if (c.amrex_slab_planes > 0) return c.amrex_slab_planes;
  // about N / ranks, even, dividing N (N even); at least 2
  int p = std::max(2, (N / std::max(nranks, 1)) / 2 * 2);
  while (p > 2 && N % p != 0) p -= 2;
  return p;
}

}  // namespace radray
