// C interface to the on-the-fly halo finder and source model, for the
// Fortran full harness (reference/harness_full/src/halo_fof_shim.f90). The
// Fortran run then gets exactly the source list the C++ driver computes
// from the same particles (onthefly_catalog in both).
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>

#include "radray/cosmology.hpp"
#include "radray/halo_sources.hpp"

namespace {
std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return std::string();
  return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
}  // namespace

extern "C" {

// xv:       xv_dm(6,npart), DM-grid units
// ngas,ndm: grid sizes; box: Mpc/h
// redshift, scalefactor: the clock at which the catalog is made (the end
//           of the hydro step, where read_halo is called)
// params:   NUL-terminated text of `key = value` lines; the halo/source keys
//           of set_halo_param are used, everything else is ignored
// nmax:     capacity of `records` (source_type records, 56 bytes each)
// nhalo, nsource: outputs
// Returns 0 on success, 1 if the catalog was truncated to nmax records,
// 2 on error (message on stderr).
int radray_onthefly_sources(const float* xv, std::int64_t npart, int ngas, int ndm, double box,
                            double redshift, double scalefactor, const char* params,
                            std::int64_t nmax, void* records, std::int64_t* nhalo,
                            std::int64_t* nsource) {
  using namespace radray;
  try {
    const Config cfg = Config::make(ngas, 64, ndm, box);
    SimClock clk;
    clk.redshift = redshift;
    clk.scalefactor = scalefactor;
    Units u;
    calc_sim_var(cfg, clk, u);

    HaloSourceParams p;
    std::istringstream in(params ? params : "");
    std::string line;
    while (std::getline(in, line)) {
      const auto h = line.find('#');
      if (h != std::string::npos) line.resize(h);
      const auto eq = line.find('=');
      if (eq == std::string::npos) continue;
      const std::string k = trim(line.substr(0, eq)), v = trim(line.substr(eq + 1));
      if (k != "source_mode" && k != "dir_halo") set_halo_param(p, k, v);
    }
    p.mode = HaloSourceParams::Mode::OnTheFly;

    const SourceCatalog cat = onthefly_catalog(cfg, clk, u, xv, npart, p);
    const std::int64_t n = static_cast<std::int64_t>(cat.halo.size());
    const std::int64_t ncopy = n < nmax ? n : nmax;
    if (ncopy > 0)
      std::memcpy(records, cat.halo.data(),
                  static_cast<std::size_t>(ncopy) * sizeof(SourceRecord));
    *nhalo = ncopy;
    *nsource = cat.Nsource < ncopy ? cat.Nsource : ncopy;
    return ncopy < n ? 1 : 0;
  } catch (const std::exception& e) {
    std::cerr << "radray_onthefly_sources: " << e.what() << "\n";
    return 2;
  }
}

}  // extern "C"
