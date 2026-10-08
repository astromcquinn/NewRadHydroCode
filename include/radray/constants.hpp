// Physical constants and cosmological / RT parameters.
//
// Sources: physical_constants.f90, cosmo_parameters.f90,
//          ray_tracing_tools.f90 (HEALPix ray parameters).
//
// Many Fortran parameters are declared real(8) but initialised from
// SINGLE-precision literals, e.g. `real(8), parameter :: omegam = 0.27`.
// gfortran rounds the literal to real(4) first, then widens, so
// omegam == (double)0.27f != 0.27. We reproduce that with static_cast<double>
// of float literals.
//
// A few parameters are folded by gfortran at compile time with MPFR
// (correctly rounded transcendental functions / integer powers). For those
// we hardcode the exact bit pattern that gfortran produces (hex float
// literals); tests/fortran_consts/expected_consts.txt has the gfortran dump
// and the `fortran_consts` test checks every value below.
#pragma once

#include <array>
#include <cstdint>

namespace radray {

// Number of photon frequency groups (cosmo_parameters.f90: Nfreq).
constexpr int kNfreq = 5;

namespace pc {  // physical_constants.f90

constexpr double pi = 3.141592653589793238462643383279502884;  // acos(-1D0)

constexpr double yr2sec = 365 * 86400;  // integer product, exact
constexpr double sec2yr = 1.0 / yr2sec;

constexpr double c_cgs = 2.99792e+10;
constexpr double eV_cgs = 1.60218e-12;
constexpr double G_cgs = 6.67259e-08;
constexpr double h_cgs = 6.62608e-27;
constexpr double k_cgs = 1.38066e-16;
constexpr double mH_cgs = 1.67223e-24;
constexpr double mHe_cgs = 4 * mH_cgs;
constexpr double sT_cgs = 6.65246e-25;
constexpr double sS_cgs = 5.67051e-05;

constexpr double K2eV = k_cgs / eV_cgs;
constexpr double eV2K = 1.0 / K2eV;
constexpr double eV2erg = eV_cgs;

constexpr double Msun_cgs = 1.98900e+33;
constexpr double H0_cgs = 3.24086e-18;
constexpr double rhoc_cgs = 1.87890e-29;
constexpr double Mpc2cm = 3.08560e+24;

// real(8), parameter :: HI_eV = 13.60  (single literal!)
constexpr double HI_eV = static_cast<double>(13.60f);
constexpr double HeI_eV = static_cast<double>(24.59f);
constexpr double HeII_eV = static_cast<double>(54.42f);
constexpr double HI_K = HI_eV * eV2K;
constexpr double HeI_K = HeI_eV * eV2K;
constexpr double HeII_K = HeII_eV * eV2K;

}  // namespace pc

namespace cp {  // cosmo_parameters.f90

// box = 5 (Mpc/h, integer literal) in the original; a run-time value here
// (Config::box), this is only the default.
constexpr double box_default = 5;

constexpr double omegam = static_cast<double>(0.27f);
constexpr double omegal = 1 - omegam;
constexpr double omegab = static_cast<double>(0.045f);
constexpr double hubble0 = static_cast<double>(0.70f);
constexpr double wde = static_cast<double>(-1.0f);

constexpr double XHy = static_cast<double>(0.76f);
constexpr double YHe = static_cast<double>(0.24f);
// gamma_ideal = 5./3 : single-precision division
constexpr double gamma_ideal = static_cast<double>(5.0f / 3.0f);

constexpr double Tcmb0 = static_cast<double>(2.725f);
// omegar = 4.48E-7*(1+0.69)*Tcmb0**4/hubble0**2
// (4.48E-7*(1+0.69) is a real(4) product; Tcmb0**4 is folded by MPFR,
// correctly rounded.) gfortran value:
constexpr double omegar = 0x1.6559c34a8bc79p-14;

// RT parameters
constexpr std::array<double, kNfreq> freq_init = {
    static_cast<double>(15.71f), static_cast<double>(21.06f),
    static_cast<double>(29.64f), static_cast<double>(42.88f),
    static_cast<double>(74.12f)};
// pop2_init = 10**25.2 * (/0.4185,0.2942,0.2009,0.0811,0.0054/)
// evaluated entirely in real(4) at compile time (10**25.2 via MPFR).
constexpr std::array<double, kNfreq> pop2_init = {
    0x1.5f23100000000p+82, 0x1.edb0780000000p+81, 0x1.5120020000000p+81,
    0x1.102f000000000p+80, 0x1.21f89e0000000p+76};

constexpr double fResc = static_cast<double>(0.2f);
constexpr double rho_source = 50;
constexpr double eSF = static_cast<double>(0.05f);
constexpr double eSN = static_cast<double>(1e-5f);

// Gas density / temperature table parameters
constexpr double lgD_min = static_cast<double>(-4.0f);
constexpr double lgD_max = static_cast<double>(6.0f);
constexpr int Ntemp = 1000;
constexpr double lgT_min = 0;
constexpr double lgT_max = 10;
constexpr double lgT_del = static_cast<double>(0.01f);

}  // namespace cp

// Values gfortran folds at compile time with correctly rounded MPFR
// functions (used in cosmo_tools.f90 / ray_tracing_tools.f90 / rt_tools.f90).
namespace folded {
// cosmo_tools.f90 calc_sim_var: x = 1E3*(ob*h0**2)**0.4
constexpr double temp_x = 0x1.b2e76e0093428p+7;
// cosmo_tools.f90 calc_sim_var: 1D1**lgD_min, 1D1**lgD_max, 1D1**lgT_max
constexpr double ten_lgD_min = 0x1.a36e2eb1c432dp-14;
constexpr double ten_lgD_max = 1.0e6;
constexpr double ten_lgT_max = 1.0e10;
// ray_tracing_tools.f90 lside_from_radius: log(4.) is real(4)
constexpr double log4_r4 = static_cast<double>(0x1.62e43p+0f);
// rt_tools.f90 rt_step_1: lgn_del = log10(1.01) (real(4) log10)
constexpr double lgn_del = 0x1.1b348ep-8;
}  // namespace folded

namespace rayp {  // ray_tracing_tools.f90 module parameters
constexpr int lside_min = 2;
constexpr int nside_min = 1 << lside_min;
constexpr int lside_max = 13;
constexpr int nside_max = 1 << lside_max;
constexpr int lside_init = 2;
constexpr int nside_init = 1 << lside_init;
constexpr int nray_init = 12 * nside_init * nside_init;  // 192
constexpr double ray_per_cell = 1;
}  // namespace rayp

}  // namespace radray
