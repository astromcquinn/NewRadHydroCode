// rt_backend = amrex without AMReX in the build (RADRAY_AMREX off): the
// process functions are no-ops and the level functions throw.
#include <stdexcept>

#include "radray/amrex_state.hpp"
#include "radray/rt_backend.hpp"

namespace radray {

namespace {
[[noreturn]] void no_amrex() {
  throw std::runtime_error("rt_backend = amrex: this build has no AMReX (cmake -DRADRAY_AMREX=ON)");
}
}  // namespace

bool amrex_available() { return false; }
void backend_init(int&, char**&) {}
void backend_finalize() {}
void backend_bcast(double*, int) {}
void backend_bcast_bytes(void*, std::size_t) {}
void amrex_state_gather_rt(Simulation&, bool) { no_amrex(); }
void amrex_state_scatter_rt(Simulation&) { no_amrex(); }
bool io_rank() { return true; }
int backend_rank() { return 0; }
int backend_nranks() { return 1; }
void backend_barrier() {}
void backend_reduce_minmax(double&, double&) {}
void amrex_rt_prep(Simulation&) { no_amrex(); }
void amrex_rt_substep(Simulation&) { no_amrex(); }
void amrex_rt_post(Simulation&) { no_amrex(); }
void amrex_rt_sync_fields(Simulation&) { no_amrex(); }
std::int64_t amrex_rt_nrays(const Simulation&) { return 0; }
int amrex_rt_bytes_per_ray(const Simulation&) { return 0; }
void amrex_rt_gather_rays(const Simulation&, std::vector<rrk::RayP>&) { no_amrex(); }
void amrex_rt_scatter_rays(Simulation&, const std::vector<rrk::RayP>&) { no_amrex(); }
std::string amrex_rt_timers(const Simulation&) { return ""; }

}  // namespace radray
