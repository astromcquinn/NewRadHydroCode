// Ionising sources from halos: the source model (SFR, UV luminosity,
// stellar mass, gas fraction per halo), compatibility with the original
// external halo catalogs (io_tools.f90 read_halo, Halo/sfr.txt), and the
// halo_tools.f90 link_halo / add_halo_to_rt steps that bin a source list
// onto the RT (and gas) grid. See docs/halos_and_sources.md.
#pragma once

#include <array>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "radray/halos.hpp"
#include "radray/rt.hpp"
#include "radray/simulation.hpp"

namespace radray {

// radhydro_global.f90: source_type (56 bytes, no padding).
//   x(3)  comoving Mpc/h (link_halo bins with x/x_rt_box, x_rt_box=box/Nrt)
//   v(3)  peculiar km/s (not used by the original code; our convention)
//   N     number of particles
//   M     Msun; R comoving Mpc/h; T K   (not used by the original code)
//   L     UV luminosity, erg/s/Hz (ngam = fResc*pop2_init*L*dtsec_rt)
//   sfr   Msun/yr
//   mstar Msun; fgas (not used by the original code; our definitions)
struct SourceRecord {
  float x[3] = {0, 0, 0};
  float v[3] = {0, 0, 0};
  std::int32_t N = 0;
  float M = 0, R = 0, T = 0;
  float L = 0, sfr = 0, mstar = 0, fgas = 0;
};
static_assert(sizeof(SourceRecord) == 56, "source_type layout");

struct SourceCatalog {
  std::vector<SourceRecord> halo;  // Nhalo records, descending mass
  std::int64_t Nsource = 0;        // leading records with L > 0
  double redshift = 0;             // redshift of the catalog (z2 for legacy)
  double sfr_scale = 1;            // legacy rescaling applied to L and sfr
  std::string tag;                 // redshift tag, e.g. "z=08.0000"
};

// ---------------------------------------------------------------------------
// Legacy files
// ---------------------------------------------------------------------------

// Halo/sfr.txt rows: `read(un,*) z,dummy(8),sfr` (list-directed: a row
// may continue on following lines, the rest of the last line is skipped).
struct SfrHistoryRow {
  double z = 0;
  std::array<double, 8> dummy{};
  double sfr = 0;
};
std::vector<SfrHistoryRow> read_sfr_history(const std::string& path);

// io_tools.f90 read_halo redshift tag: z<10 -> 'z=0',f6.4; z<100 ->
// 'z=',f7.4 (Fortran F editing; gfortran rounds ties to even, as
// printf does). The original leaves the tag unset for z>=100; we throw.
std::string redshift_tag(double z);

// The sfr.txt scan of read_halo: z1/z2 = first row with z <= redshift and
// the row before it; catalog at z2; sfr_scale = sfr(redshift)/sfr2 by linear
// interpolation in z (extrapolation if the file ends first).
struct CatalogSelection {
  bool found = false;  // z2 > 0 .and. sfr2 > 0
  double z1 = 0, z2 = 0, sfr1 = 0, sfr2 = 0;
  double sfr = 0, sfr_scale = 0;
  std::string tag;
};
CatalogSelection select_halo_catalog(const std::vector<SfrHistoryRow>& hist, double redshift);

// Binary catalog: Nhalo (integer(8)) followed by Nhalo source_type records.
//  Stream:     no record markers. This is what the original reads: Intel
//              `form='binary'`, equivalent to gfortran access='stream'.
//  Sequential: Fortran unformatted sequential (4-byte record markers around
//              each record), what `form='unformatted'` would give.
//  Auto (read): decided from the file size.
// Little-endian, native IEEE floats.
enum class CatalogFormat { Stream, Sequential, Auto };

void write_halo_catalog(const std::string& path, const std::vector<SourceRecord>& halo,
                        CatalogFormat fmt = CatalogFormat::Stream);
std::vector<SourceRecord> read_halo_catalog(const std::string& path,
                                            CatalogFormat fmt = CatalogFormat::Auto);

// io_tools.f90 read_halo. dir_halo includes the trailing '/', as in
// radhydro_global.f90 (dir_halo = 'Halo/'). Leading records with L > 0 are
// the sources; their L and sfr are multiplied by sfr_scale in real(8) and
// stored back in real(4). No catalog (Nsource = 0) if the selection fails.
SourceCatalog read_halo(const std::string& dir_halo, double redshift, std::ostream* log = nullptr);

// ---------------------------------------------------------------------------
// Source model
// ---------------------------------------------------------------------------

enum class SfrModel {
  // Default, reconstruction of Trac & Cen (2007): SFR proportional to halo
  // mass for halos above an atomic-cooling threshold, normalised so that
  // the comoving SFR density of the box equals a global SFR history
  // (sfr.txt, or Hernquist & Springel 2003 analytic fit).
  GlobalMassProportional,
  // SFR = eSF * fb * M / t_dyn, t_dyn = sqrt(3 pi / (32 G rho_vir)),
  // the same form as the SN-heating star formation of rt_tools.f90.
  Dynamical,
};

struct SourceModelParams {
  SfrModel model = SfrModel::GlobalMassProportional;
  double Tvir_min = 1e4;  // K: halos below do not form stars
  double M_min = 0;       // Msun: optional extra mass cut

  // Global SFR history for GlobalMassProportional. If sfr_history_path is
  // set, column 10 of sfr.txt (Msun/yr/Mpc^3 comoving, assumed) is linearly
  // interpolated in z; otherwise the Hernquist & Springel (2003, eq. 46)
  // fit rho_sfr = rho_m chi^2/(1+alpha(chi-1)^3 exp(beta chi^(7/4))),
  // chi = (H(z)/H0)^(2/3).
  std::string sfr_history_path;
  double hs_rho_m = 0.013;  // Msun/yr/Mpc^3
  double hs_alpha = 0.012;
  double hs_beta = 0.041;

  // UV luminosity per unit SFR, erg/s/Hz per Msun/yr. 8.0e27 is the 1500 A
  // Salpeter value of Madau, Pozzetti & Dickinson (1998); with
  // pop2_init = 10**25.2 photons per erg/Hz it gives 1.3e53 ionising photons
  // per second per Msun/yr (before fResc), a standard Pop II value.
  double L_per_sfr = 8.0e27;

  double eSF = 0.05;  // Dynamical model efficiency (cosmo_parameters.f90)
  double z_max = 50;  // start of the integration for the stellar mass
};

// Global comoving SFR density (Msun/yr/Mpc^3) at redshift z: table if
// given, else the Hernquist & Springel fit.
double global_sfr_density(double z, const SourceModelParams& p,
                          const std::vector<SfrHistoryRow>* table);

// Comoving stellar mass density formed by redshift z (Msun/Mpc^3):
// integral of global_sfr_density over cosmic time from z_max.
double global_stellar_density(double z, const SourceModelParams& p,
                              const std::vector<SfrHistoryRow>* table);

// Cosmic time between z and z_max in years (flat LCDM of cosmo_parameters).
double cosmic_time_yr(double z, double z_max);

// Assign SFR, L, mstar and fgas to halos (already sorted by descending
// mass). Returns the source_type records in the legacy order: star-forming
// halos (L > 0) first, in the input order, then the others.
SourceCatalog assign_sources(const std::vector<Halo>& halos, const HaloUnits& hu,
                             const SourceModelParams& p,
                             const std::vector<SfrHistoryRow>* table = nullptr);

// ---------------------------------------------------------------------------
// halo_tools.f90: link_halo / calc_halo / add_halo_to_rt
// ---------------------------------------------------------------------------

// link_halo: RT cell of each of the first Nsource records,
// i = 1 + mod(int(x/x_rt_box), Nrt). The returned list is in record order,
// which is the order link_halo's chains visit the sources of a cell, so it
// can be passed to calc_halo (rt.hpp) unchanged.
std::vector<Source> link_halo(const std::vector<SourceRecord>& src, std::int64_t Nsource,
                              const Config& cfg);

// add_halo_to_rt metal injection: w_gas(4) += fZesc2*Zyield2*sfr*dt_yr*Msun
// /vol_gas/rho_unit at the gas cell 1+mod(int(x/x_gas_box),Ngas), in the
// original order (RT cells k,j,i; sources within a cell in record order).
void add_halo_metals(Simulation& s, const std::vector<SourceRecord>& src, std::int64_t Nsource);

// calc_halo for a source_type list: link_halo + rt.hpp calc_halo
// (+ metals if add_metals).
void calc_halo_records(Simulation& s, const std::vector<SourceRecord>& src,
                       std::int64_t Nsource, bool add_metals);

// ---------------------------------------------------------------------------
// Per-step driver
// ---------------------------------------------------------------------------

struct HaloSourceParams {
  enum class Mode { Legacy, OnTheFly };
  Mode mode = Mode::OnTheFly;
  std::string dir_halo = "Halo/";     // legacy input
  HaloFinderParams finder;
  SourceModelParams model;
  bool add_metals = true;             // the original always injects metals
  std::string write_catalog_dir;      // if set, on-the-fly catalogs are
                                      // written there in legacy format
};

// On-the-fly catalog: find_halos on xv_dm + assign_sources (+ optional
// legacy dump), for the redshift of clk; units from calc_sim_var. The same
// function serves the C++ driver and the C interface used by the Fortran
// full harness (radray_onthefly_sources, halo_capi.cpp).
SourceCatalog onthefly_catalog(const Config& cfg, const SimClock& clk, const Units& u,
                               const float* xv_dm, std::int64_t npart,
                               const HaloSourceParams& p);

// Sets one halo/source parameter from a `key = value` entry (params.txt of
// apps/radhydro and of the full Fortran harness). Returns false if the key
// is not a halo/source key. Keys: source_mode (legacy|onthefly), dir_halo,
// fof_b, fof_min_npart, fof_so (0|1), sfr_model (global|dynamical),
// tvir_min, m_min, l_per_sfr, sfr_history, write_catalogs (directory),
// add_metals (0|1).
bool set_halo_param(HaloSourceParams& p, const std::string& key, const std::string& val);

// What radhydro.f90 does once per hydro step before radiative_transfer:
// obtain the source list (read_halo, or FoF + source model on xv_dm) and,
// if Nsource > 0, calc_halo. Returns the catalog; the caller calls
// radiative_transfer when Nsource > 0. xv_dm may be null in legacy mode.
SourceCatalog update_sources(Simulation& s, const float* xv_dm, std::int64_t npart,
                             const HaloSourceParams& p, std::ostream* log = nullptr);

}  // namespace radray
