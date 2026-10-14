// hydro_backend = amrex without AMReX in the build (RADRAY_AMREX off): the
// level functions throw.
#include <stdexcept>

#include "radray/hydro_backend.hpp"

namespace radray {

void amrex_hydro_take_sink_events(std::vector<SinkEvent>& out) { out.clear(); }
void amrex_hydro_step(const Config&, const Units&, double, double, GasGrid&, HydroState&) {
  throw std::runtime_error("hydro_backend = amrex: this build has no AMReX (cmake -DRADRAY_AMREX=ON)");
}
std::string amrex_hydro_timers() { return ""; }
int amrex_hydro_finest_level() { return 0; }
void amrex_hydro_write_levels(const std::string&) {}
void amrex_hydro_read_levels(const Config&, int, const std::string&) { throw std::runtime_error("hydro_backend = amrex: this build has no AMReX"); }
void amrex_hydro_release() {}

}  // namespace radray
