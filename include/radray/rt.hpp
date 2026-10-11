// Port of rt_tools.f90 (+ domain_tools.f90 RT part, halo_tools.f90
// add_halo_to_rt).
#pragma once

#include <array>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "radray/constants.hpp"
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

// rt_cells_per_substep = k (Config, docs/rt_multicell.md): throws unless
// 1 <= k <= rrk::kMaxCellsPerSubstep, and for k > 1 unless legacy_rt = 0
// and rt_method = rays.
void check_rt_cells_param(const Config& c);

// rt_tools.f90: radiative_transfer and its parts
void radiative_transfer(Simulation& s);
void prep_rt(Simulation& s);
void post_rt(Simulation& s);
void rt_step_1(Simulation& s);
void rt_step_2(Simulation& s);
void rt_step_3(Simulation& s);
void radiative_physics(Simulation& s, int krt);
// The end of rt_step_3: the per-hydro-step sums of rt_stat / ray_stat
void accumulate_rt_stats(Simulation& s);

// rt_step_1 / calc_ray_1: the ngam_ray thresholds from the photon-count
// histogram ngam_pdf[1..N_pdf] (rrk::pdf_bin) of the live rays, and the F3
// ray-pool guard (legacy_rt = 0) from the histogram spl_pdf of the rays that
// may split. Shared with the AMReX backend; see rt_tools.cpp.
void calc_ngam_ray(const std::int64_t* ngam_pdf, int N_pdf, double ngam_avg, std::int64_t Nray,
                   std::int64_t Nray_max, bool legacy, std::array<double, 64>& ngam_ray,
                   std::ostream* log);
std::int64_t pool_split_candidates(const std::int64_t* spl_pdf, int N_pdf, double ngam_ray32);
// nray_src: rays per stellar emission (nray_init = 192; more with
// rt_cells_per_substep > 4, rrk::emit_nside).
void pool_guard(const std::int64_t* spl_pdf, int N_pdf, std::int64_t L, std::int64_t nsrc,
                std::int64_t nrec_rays, std::int64_t Nray_max, std::array<double, 64>& ngam_ray,
                std::int64_t& npool_limited, std::int64_t& npool_forced, std::ostream* log,
                int nray_src = rayp::nray_init);

// rt_tools.f90: summarize_rt (writes to s.out)
void summarize_rt(Simulation& s);

// One line with the accumulated RT phase timers (Simulation::timers)
std::string format_rt_timers(const RTTimers& t);

// Diagnostics for the validation outputs
std::int64_t count_live_rays(const Simulation& s);  // slots with nside > 0
double ngam_rays_total(const Simulation& s);         // sum ngam, live rays
double ngam_rays_total_serial(const Simulation& s);  // the serial chains' sum

}  // namespace radray
