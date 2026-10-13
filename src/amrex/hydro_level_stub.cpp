// hydro_backend = amrex without AMReX in the build (RADRAY_AMREX off): the
// level functions throw.
#include <stdexcept>

#include "radray/hydro_backend.hpp"

namespace radray {

void amrex_hydro_step(const Config&, const Units&, double, GasGrid&, HydroState&) {
  throw std::runtime_error("hydro_backend = amrex: this build has no AMReX (cmake -DRADRAY_AMREX=ON)");
}
std::string amrex_hydro_timers() { return ""; }
void amrex_hydro_release() {}

}  // namespace radray
