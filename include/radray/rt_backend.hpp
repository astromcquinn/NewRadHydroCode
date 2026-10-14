// The RT backend switch (docs/stage_d.md): rt_backend = serial (the
// original's ray chains, src/ray_tracing.cpp) | amrex (the sorted-ray-array
// RT on an AMReX level, src/amrex/rt_level.cpp, MPI + OpenMP). This header
// is AMReX-free; with RADRAY_AMREX off the amrex_* functions are stubs that
// throw (src/amrex/rt_level_stub.cpp) and the process functions are no-ops.
#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "radray/config.hpp"
#include "radray/kernels/ray_kernels.hpp"

namespace radray {

struct Simulation;

// Parameters: rt_backend, ray_encoding, amrex_max_grid_size, amrex_ghosts,
// amrex_redist_local, amrex_slab, amrex_slab_planes; false for any other key.
bool set_rt_backend_param(Config& c, const std::string& key, const std::string& val);
void copy_rt_backend_param(const Config& from, Config& to);
std::string rt_backend_summary(const Config& c);
// Throws for a combination the backend does not implement (legacy_rt = 1,
// rt_method != rays, recomb_rad, deterministic_parallel).
void check_rt_backend_param(const Config& c);
// amrex_slab = 1: the k-planes per slab for an N^3 grid on nranks ranks
// (amrex_slab_planes, or for 0 about N / nranks: even, dividing N, >= 2)
int amrex_slab_planes(const Config& c, int N, int nranks);
inline bool use_amrex(const Config& c) { return c.rt_backend == RTBackend::AMReX; }

// Whether this build has the AMReX backend (RADRAY_AMREX).
bool amrex_available();

// Process set-up of the backend: amrex::Initialize (and MPI) / Finalize.
// No-ops without AMReX. Every driver calls them around its work; with MPI
// every rank runs the whole driver and only io_rank() writes files.
void backend_init(int& argc, char**& argv);
void backend_finalize();
bool io_rank();
int backend_rank();
int backend_nranks();
void backend_barrier();
// Min and max over the ranks of a double (identity without MPI)
void backend_reduce_minmax(double& lo, double& hi);
// io_rank()'s values of n doubles, or n bytes, onto every rank (identity
// without MPI)
void backend_bcast(double* v, int n);
void backend_bcast_bytes(void* p, std::size_t n);

// The RT of one hydro step on the level (called by rt_tools.cpp in place of
// prep_rt / rt_step_1-3 / post_rt when use_amrex(cfg)):
//   prep: the gas and sources of s (gas.u/w/e, rt.Luv/sfr) onto the level
//         and the RT grid gas state (calc_rt_i);
//   substep: rt_step_1 (thresholds, pool guard, emission/split/merge, the
//         casts), rt_step_2 (radiative_physics), rt_step_3 (advance,
//         migration); fills s.rays.stat/ngam_ray/Nray and s.rt_stat;
//   post: ray redshifting, calc_rt_f, and the level's fields back into
//         s.rt and s.gas on every rank;
//   sync_fields: the level's current fields into s.rt and s.gas (movie
//         frames inside a hydro step).
void amrex_rt_prep(Simulation& s);
void amrex_rt_substep(Simulation& s);
void amrex_rt_post(Simulation& s);
void amrex_rt_sync_fields(Simulation& s);
// Live rays on the level (all ranks) and the bytes of ray storage per ray
std::int64_t amrex_rt_nrays(const Simulation& s);
int amrex_rt_bytes_per_ray(const Simulation& s);
// Checkpoints: all rays gathered on io_rank() in the float layout
// (collective), and the reverse (every rank passes the full list and keeps
// the rays of its boxes).
void amrex_rt_gather_rays(const Simulation& s, std::vector<rrk::RayP>& out);
void amrex_rt_scatter_rays(Simulation& s, const std::vector<rrk::RayP>& in);
// The streaming forms (docs/stage_e_state.md: rank 0's transient is one
// chunk, not the pool): gather_chunks calls sink(chunk) on io_rank() for
// rank 0's rays and then every other rank's, in rank order, in chunks of
// at most 8 M rays (collective; the total is amrex_rt_nrays); scatter with
// first = true clears the level's rays, false appends, each call scattering
// one chunk held by io_rank() (collective, the same number of calls on
// every rank).
void amrex_rt_gather_rays_chunks(const Simulation& s,
                                 const std::function<void(const std::vector<rrk::RayP>&)>& sink);
void amrex_rt_scatter_rays(Simulation& s, const std::vector<rrk::RayP>& in, bool first);
// One line of the backend's phase timers
std::string amrex_rt_timers(const Simulation& s);
// rt_method = m1 | otvet on the level (checkpoint.cpp): the checkpoint's
// moment block (the five moment vectors and the budget) and recombination
// block (version >= 4) in the serial format, streamed slab by slab
// (collective): write calls sink(bytes) in file order (meaningful on
// io_rank() only); read calls rd(dst, bytes), which reads into dst, or
// skips the bytes for dst = null (the ranks other than io_rank() pass null
// for the data, real buffers for the headers and the budgets).
void amrex_moment_ckpt_write(const Simulation& s, const std::function<void(const void*, std::size_t)>& sink);
void amrex_moment_ckpt_read(Simulation& s, const std::function<void(void*, std::size_t)>& rd, int version);

}  // namespace radray
