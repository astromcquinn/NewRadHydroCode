// A uniform X-ray background for preheating the IGM (docs/xray_heating.md),
// new physics behind xray_heating = 1 (default 0: the original, no X-rays).
//
//  1. Star formation: the comoving SFR density of all halos above the
//     atomic-cooling threshold, rho_SFR(z) = int dn/dlnM SFR(M, z) dlnM,
//     with the mass function the source model was calibrated with (the
//     Sheth-Tormen refit to the FoF halos, data/hmf/st_fof_40mpc_1024.txt)
//     and SFR = L_UV(M, z) / l_per_sfr of the abundance-matched source model
//     (sfr_model = abundance, docs/abundance_matching.md). Every halo in the
//     universe counts, not only those the box resolves.
//  2. Emission: L_X / SFR (log10, erg/s per Msun/yr) in the band
//     [e0, eband] keV, a power law L_E ~ E^-alpha from e0 to emax
//     (nothing below e0: absorbed in the host galaxy).
//  3. Transfer: the mean specific intensity J_E(z) of a homogeneous
//     background, redshifted and absorbed by the mean HI, HeI and HeII of
//     the IGM (Verner et al. 1996 cross-sections, as the RT):
//       J(E,z) = c/(4 pi) (1+z)^3 int_z^zmax dz'/(H(z')(1+z'))
//                eps(E (1+z')/(1+z), z') exp(-tau),
//     with eps the comoving emissivity per unit energy. The energy grid and
//     the ln(1+z) grid share one step, so a photon redshifts by exactly one
//     energy bin per redshift step (as in Mirocha 2014's ARES).
//  4. Rates per atom of species j (HI, HeI, HeII): the primary
//     photoionisation rate gamma_j = 4 pi int J/E sigma_j dE and the
//     photoelectron power phi_j = 4 pi int J/E sigma_j (E - I_j) dE. The
//     cell physics (rt_cell.hpp) splits the photoelectron energy into heat
//     and secondary ionisations of HI and HeI with the Shull & van
//     Steenberg (1985) fractions at the cell's ionised fraction.
//
// The defaults of the emission parameters are from memory and to be
// confirmed against the papers (docs/xray_heating.md section 5).
#pragma once

#include <string>
#include <vector>

#include "radray/abundance_matching.hpp"
#include "radray/mass_function.hpp"

namespace radray {

struct SourceModelParams;
struct Simulation;

struct XrayParams {
  bool on = false;            // xray_heating = 0 | 1
  double log_lx_sfr = 40.5;   // xray_log_lx_sfr: log10 L_X/SFR [erg/s/(Msun/yr)] in [e0, eband]
  double e0_keV = 0.5;        // xray_e0: lowest emitted energy (keV)
  double eband_keV = 2.0;     // xray_eband: upper edge of the normalisation band (keV)
  double emax_keV = 10.0;     // xray_emax: highest emitted energy (keV)
  double alpha = 1.0;         // xray_alpha: energy index, L_E ~ E^-alpha
  // Halos that emit: virial temperature (K) and mass (Msun) cuts; < 0
  // follows the source model's tvir_min, m_min
  double tvir_min = -1;       // xray_tvir_min
  double m_min = -1;          // xray_m_min
  double mu = 0.59;           // mean molecular weight of the virial temperature (halos.hpp)
  // Mass function (xray_st = A,a,p): the production refit
  // data/hmf/st_fof_40mpc_1024.txt (docs/abundance_matching.md)
  MassFunction::STParams st{0.345835, 0.616383, 0.3};
  double z_max = 40;          // xray_zmax: no emission above
  double dlog = 0.005;        // xray_dlog: step in ln E and in ln(1+z)
  // Absorbers (xray_opacity): box = the mass-weighted HI, HeI, HeII
  // fractions of the box at the start of the hydro step, applied at every
  // earlier redshift (stateless); neutral = a neutral IGM; none = no
  // absorption (tests)
  enum class Opacity { Box, Neutral, None } opacity = Opacity::Box;
};

// Sets one xray_* parameter from `key = value`; false if not an X-ray key.
bool set_xray_param(XrayParams& p, const std::string& key, const std::string& val);

// The uniform rates of one redshift (per atom: HI, HeI, HeII)
struct XrayRates {
  double z = 0;
  double gamma[3] = {0, 0, 0};  // primary photoionisation rate, 1/s
  double phi[3] = {0, 0, 0};    // photoelectron power, erg/s
  double xfrac[3] = {0, 0, 0};  // the absorber fractions used
  double sfrd = 0;              // Msun/yr/Mpc^3 comoving
  double eps = 0;               // erg/s/Mpc^3 comoving in [e0, eband]
};

class XrayBackground {
 public:
  // The source model must be sfr_model = abundance (its L_UV(M, z),
  // l_per_sfr, tvir_min and m_min). Tabulates sigma(M) and rho_SFR(z).
  // co: the simulation's cosmology (Config::cosmo): the expansion, the
  // mean densities and the mass function; the original's by default
  XrayBackground(const XrayParams& p, const SourceModelParams& src, const Cosmology& co = {});

  const XrayParams& params() const { return p_; }

  // Smallest emitting halo mass (Msun) at z: the virial-temperature cut
  // with the halo finder's virial radius (Bryan & Norman 1998 overdensity,
  // halos.cpp), or m_min if larger.
  double m_min_msun(double z) const;
  // rho_SFR (Msun/yr/Mpc^3 comoving), direct integration (slow: sigma table
  // but a fresh growth factor) and from the table (log-linear in ln(1+z))
  double sfrd_direct(double z) const;
  double sfrd(double z) const;
  // Comoving emissivity in [e0, eband], erg/s/Mpc^3
  double emissivity(double z) const;
  // Spectral shape per keV, normalised to 1 over [e0, eband]; 0 outside
  // [e0, emax]
  double shape(double E_keV) const;
  // The mean of shape over [E1, E2] in ln E (the bins of spectrum)
  double shape_avg(double E1_keV, double E2_keV) const;

  // J_E(E) at z0 (erg/s/cm^2/sr/keV) on the energy grid E_keV (ascending),
  // for the absorber fractions xfrac (multiples of the mean H and He
  // densities: HI/H, HeI/He, HeII/He)
  void spectrum(double z0, const double xfrac[3], std::vector<double>& E_keV,
                std::vector<double>& J) const;
  // The rates at z0 for the absorber fractions xfrac (ignored with
  // opacity = neutral (1, 1, 0) or none (0, 0, 0))
  XrayRates rates(double z0, const double xfrac[3]) const;

 private:
  XrayParams p_;
  Cosmology co_;
  LuvModel am_;
  double l_per_sfr_ = 8e27, tvir_min_ = 1e4, m_min_ = 0;
  MassFunction mf_;
  // sigma0 and dln(1/sigma)/dlnM on a ln M grid (M in Msun/h)
  double lnM0_ = 0, dlnM_ = 0;
  std::vector<double> sig0_, dsig_;
  // rho_SFR on the ln(1+z) grid 0 .. ln(1+z_max)
  double dlz_tab_ = 0;
  std::vector<double> lsfrd_;
};

// Shull & van Steenberg (1985) fits (as in Ricotti, Gnedin & Shull 2002) of
// the fractions of a fast photoelectron's energy that go into heat and into
// ionising HI and HeI, at ionised fraction x = n_HII/n_H (clamped to
// [1e-4, 1], the fits' range); in rt_cell.hpp (xray_secondary) for the cell
// physics.

// Mass-weighted HI/H, HeI/He, HeII/He of the gas grid (deterministic: per
// k-plane double sums added in plane order).
void box_ion_fractions(const Simulation& s, double x[3]);

}  // namespace radray
