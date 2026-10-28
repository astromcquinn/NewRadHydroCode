// The thermal history of the gas before any heating, thermal_init = recfast
// (docs/thermal_history.md; quirk T1 of docs/port_notes.md 5b). New physics,
// off by default (thermal_init = floor, the original).
//
// The original starts the gas at the floor temp_min (an instantaneous
// decoupling at z ~ 216) and gives the neutral gas no electrons, so it never
// couples to the CMB: the unheated gas is 0.6-0.76x the RECFAST temperature
// at z <= 100. With thermal_init = recfast:
//  * the ICs start at T_RECFAST(z_init) (recfast_T);
//  * every hydro step, every gas cell couples to the CMB through the
//    residual electrons of recombination, x_e,res(z) (recfast_xe), with the
//    Compton coefficient of the cooling sum (coef_cool(14)):
//      dT/dt = (gamma-1) C14 n_e,add T_CMB^4 (T_CMB - T) / (k n_tot),
//    solved exactly over the step at the step's middle (T_CMB and x_e held):
//      T <- T_CMB + (T - T_CMB) exp(-dt/t_C).
//    n_e,add = max(x_e,res n_H - n_e, 0): only what the cell's own electrons
//    (ionised gas, X-ray electrons) lack, whose Compton term the RT physics
//    already applies. Thermal only: the ionisation state is unchanged.
// The table (src/recfast_table.inc, tools/thermal/make_recfast_table.py) is
// CAMB + RECFAST for the simulation cosmology of cosmo_parameters.f90.
#pragma once

#include "radray/config.hpp"
#include "radray/cosmology.hpp"
#include "radray/grid.hpp"

namespace radray {

// RECFAST gas temperature (K) and x_e = n_e/n_H at z, log-linear in ln(1+z);
// above the table (z > 2000) T = T_CMB and x_e at the table's top.
double recfast_T(double z);
double recfast_xe(double z);

// One hydro step of the residual-electron Compton coupling on every cell of
// gas (u(5), w; the dual-energy field e when allocated): dt_sec seconds at
// scale factor a (the middle of the step), with the units of that time and
// C14 = coef_cool(14) (AtomicRates). Thread-count independent (per cell).
// Returns the number of cells changed.
std::int64_t compton_residual_step(const Units& u, double a, double dt_sec, double C14,
                                   GasGrid& gas);

}  // namespace radray
