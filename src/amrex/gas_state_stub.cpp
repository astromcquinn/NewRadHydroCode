// amrex_state = 1 without AMReX (RADRAY_AMREX off): check_amrex_state_param
// refuses it, so none of these is reached.
#include <stdexcept>

#include "radray/amrex_state.hpp"

namespace radray {

namespace {
[[noreturn]] void no_amrex() { throw std::runtime_error("amrex_state = 1: this build has no AMReX"); }
}  // namespace

void amrex_state_init(Simulation&) { no_amrex(); }
void amrex_state_scatter_gas(Simulation&) { no_amrex(); }
void amrex_state_gather_gas(Simulation&) { no_amrex(); }
void amrex_state_clamp_gas(Simulation&) { no_amrex(); }
void amrex_state_add_metals(Simulation&, const std::vector<std::int64_t>&, const std::vector<double>&) { no_amrex(); }
void amrex_state_release() {}

}  // namespace radray
