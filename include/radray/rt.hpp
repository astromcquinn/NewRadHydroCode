// Port of rt_tools.f90 (+ domain_tools.f90 RT part, halo_tools.f90
// add_halo_to_rt).
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "radray/simulation.hpp"

namespace radray {

// A point source as given in sources.txt (1-based RT cell indices).
struct Source {
  int i = 0, j = 0, k = 0;
  float Luv = 0;  // source_type%L  (real(4))
  float sfr = 0;  // source_type%sfr (real(4))
};

// domain_tools.f90: init_domain_rt / set_domain_rt / end_domain_rt
void init_domain_rt(const Config& cfg, Domains& dom);
int set_domain_rt(Domains& dom);              // returns 0-based domain index
void end_domain_rt(Domains& dom, int idom);   // idom 0-based

// rt_tools.f90: init_rt (allocates rt grid, zeroes stats, atomic rates)
void init_rt(Simulation& s);

// Full initialisation of a Simulation for grid size N (init_variables,
// init_domain_rt, init_rt, init_ray_tracing, calc_sim_var).
void init_simulation(Simulation& s, int N, int nray_max_factor, double z_init);

// halo_tools.f90: calc_halo / add_halo_to_rt (metal injection skipped).
// Sources are summed per cell in the given order (= link_halo order).
void calc_halo(Simulation& s, const std::vector<Source>& sources);

// rt_tools.f90: radiative_transfer and its parts
void radiative_transfer(Simulation& s);
void prep_rt(Simulation& s);
void post_rt(Simulation& s);
void rt_step_1(Simulation& s);
void rt_step_2(Simulation& s);
void rt_step_3(Simulation& s);
void radiative_physics(Simulation& s, int krt);

// rt_tools.f90: summarize_rt (writes to s.out)
void summarize_rt(Simulation& s);

// Diagnostics for the validation outputs
std::int64_t count_live_rays(const Simulation& s);  // slots with nside > 0
double ngam_rays_total(const Simulation& s);         // sum ngam, live rays

}  // namespace radray
