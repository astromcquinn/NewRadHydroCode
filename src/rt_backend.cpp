// rt_backend parameters (rt_backend.hpp); AMReX-free.
#include "radray/rt_backend.hpp"

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
}

std::string rt_backend_summary(const Config& c) {
  if (!use_amrex(c)) return "rt_backend = serial";
  return "rt_backend = amrex (ray_encoding = " +
         std::string(c.ray_encoding == RayEncoding::RStruct ? "rstruct" : "float") +
         ", max_grid_size = " + std::to_string(c.amrex_max_grid_size) +
         ", ghosts = " + std::to_string(c.amrex_ghosts) +
         ", redist_local = " + std::to_string(c.amrex_redist_local) + ")";
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
  if (c.rt_method != RTMethod::Rays) throw std::runtime_error("rt_backend = amrex needs rt_method = rays");
  if (c.recomb_rad) throw std::runtime_error("rt_backend = amrex does not implement recomb_rad yet");
  if (c.det_parallel)
    throw std::runtime_error("rt_backend = amrex: deterministic_parallel applies to the serial backend only");
  if (c.N_rt_gas != 1) throw std::runtime_error("rt_backend = amrex needs N_rt_gas = 1");
}

}  // namespace radray
