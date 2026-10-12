// pm_backend = amrex without AMReX in the build (RADRAY_AMREX off): the
// level functions throw (the process functions are in rt_level_stub.cpp).
#include <stdexcept>

#include "radray/pm_backend.hpp"

namespace radray {

namespace {
[[noreturn]] void no_amrex() {
  throw std::runtime_error("pm_backend = amrex: this build has no AMReX (cmake -DRADRAY_AMREX=ON)");
}
}  // namespace

void amrex_pm_step(pm::PMState&, float*, float*, SimClock&, const Units&) { no_amrex(); }
void amrex_pm_step_unsplit(pm::PMState&, float*, const float*, double, SimClock&, const Units&) {
  no_amrex();
}
void amrex_pm_fields(pm::PMState&, const float*, double) { no_amrex(); }
std::string amrex_pm_timers(const pm::PMState&) { return ""; }

}  // namespace radray
