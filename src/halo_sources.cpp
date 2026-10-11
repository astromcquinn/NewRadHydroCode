// Source model, legacy halo catalogs and link_halo/add_halo_to_rt
// (see include/radray/halo_sources.hpp and docs/halos_and_sources.md).
#include "radray/halo_sources.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include "radray/constants.hpp"
#include "radray/fortran.hpp"

namespace radray {

using namespace pc;
using namespace cp;

// ---------------------------------------------------------------------------
// sfr.txt
// ---------------------------------------------------------------------------

namespace {

// Split one line into list-directed value tokens (blanks and commas
// separate; a '/' ends the read; Fortran 'd' exponents accepted).
std::vector<std::string> list_tokens(const std::string& line) {
  std::vector<std::string> t;
  std::string cur;
  auto flush = [&] {
    if (!cur.empty()) t.push_back(cur);
    cur.clear();
  };
  for (char c : line) {
    if (c == '/') break;
    if (c == ' ' || c == '\t' || c == ',' || c == '\r') {
      flush();
    } else {
      cur.push_back((c == 'd' || c == 'D') ? 'e' : c);
    }
  }
  flush();
  return t;
}

double to_double(const std::string& s, const std::string& path) {
  char* end = nullptr;
  const double v = std::strtod(s.c_str(), &end);
  if (end == s.c_str() || *end != '\0')
    throw std::runtime_error(path + ": cannot parse value '" + s + "'");
  return v;
}

}  // namespace

std::vector<SfrHistoryRow> read_sfr_history(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<SfrHistoryRow> rows;
  std::string line;
  // Each READ consumes whole records until it has 10 values; values left on
  // its last record are discarded. A READ hitting EOF part way is an
  // end-of-file condition (the row is dropped), as with end=1.
  for (;;) {
    std::vector<double> v;
    bool eof = false;
    while (v.size() < 10) {
      if (!std::getline(in, line)) {
        eof = true;
        break;
      }
      for (const auto& tok : list_tokens(line)) {
        if (v.size() < 10) v.push_back(to_double(tok, path));
      }
    }
    if (eof) break;
    SfrHistoryRow r;
    r.z = v[0];
    for (int i = 0; i < 8; ++i) r.dummy[i] = v[1 + i];
    r.sfr = v[9];
    rows.push_back(r);
  }
  return rows;
}

// ---------------------------------------------------------------------------
// Redshift tag (Fortran F editing)
// ---------------------------------------------------------------------------

namespace {

// Fortran Fw.d of a non-negative value, right-justified in width w. gfortran
// rounds the exact binary value to nearest with ties to even, exactly like
// glibc printf (checked against gfortran in tests/test_halo_catalog.cpp).
std::string fortran_f(double x, int w, int d) {
  char buf[512];
  std::snprintf(buf, sizeof buf, "%.*f", d, x);
  std::string out(buf);
  // The leading zero of 0.xxxx is optional: dropped if the field is narrow
  if (static_cast<int>(out.size()) > w && out[0] == '0') out = out.substr(1);
  if (static_cast<int>(out.size()) > w) return std::string(w, '*');
  return std::string(w - out.size(), ' ') + out;
}

}  // namespace

std::string redshift_tag(double z) {
  if (!(z >= 0)) throw std::invalid_argument("redshift_tag: negative redshift");
  if (z < 10) return "z=0" + fortran_f(z, 6, 4);
  if (z < 100) return "z=" + fortran_f(z, 7, 4);
  throw std::invalid_argument("redshift_tag: z >= 100 (tag undefined in read_halo)");
}

// ---------------------------------------------------------------------------
// read_halo
// ---------------------------------------------------------------------------

CatalogSelection select_halo_catalog(const std::vector<SfrHistoryRow>& hist, double redshift) {
  CatalogSelection c;
  double z1 = 0, z2 = 0, sfr1 = 0, sfr2 = 0;
  for (const auto& r : hist) {
    z2 = z1;
    sfr2 = sfr1;
    z1 = r.z;
    sfr1 = r.sfr;
    if (r.z <= redshift) break;
  }
  c.z1 = z1;
  c.z2 = z2;
  c.sfr1 = sfr1;
  c.sfr2 = sfr2;
  if (z2 > 0 && sfr2 > 0) {
    c.found = true;
    c.tag = redshift_tag(z2);
    c.sfr = sfr2 + (sfr1 - sfr2) * (redshift - z2) / (z1 - z2);
    c.sfr_scale = c.sfr / sfr2;
  }
  return c;
}

void write_halo_catalog(const std::string& path, const std::vector<SourceRecord>& halo,
                        CatalogFormat fmt) {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot write " + path);
  const std::int64_t n = static_cast<std::int64_t>(halo.size());
  const bool seq = fmt == CatalogFormat::Sequential;
  auto marker = [&](std::int32_t m) {
    if (seq) out.write(reinterpret_cast<const char*>(&m), 4);
  };
  marker(8);
  out.write(reinterpret_cast<const char*>(&n), 8);
  marker(8);
  for (const auto& r : halo) {
    marker(56);
    out.write(reinterpret_cast<const char*>(&r), sizeof r);
    marker(56);
  }
  if (!out) throw std::runtime_error("write failed: " + path);
}

std::vector<SourceRecord> read_halo_catalog(const std::string& path, CatalogFormat fmt) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  in.seekg(0, std::ios::end);
  const std::int64_t size = static_cast<std::int64_t>(in.tellg());
  in.seekg(0);
  std::int64_t n = 0;
  if (fmt == CatalogFormat::Auto) {
    std::int64_t first = 0;
    in.read(reinterpret_cast<char*>(&first), 8);
    in.seekg(0);
    if (in && first >= 0 && size == 8 + 56 * first) {
      fmt = CatalogFormat::Stream;
    } else {
      std::int32_t m = 0;
      in.read(reinterpret_cast<char*>(&m), 4);
      std::int64_t nn = -1;
      in.read(reinterpret_cast<char*>(&nn), 8);
      in.seekg(0);
      if (m == 8 && nn >= 0 && size == 16 + 64 * nn)
        fmt = CatalogFormat::Sequential;
      else
        throw std::runtime_error(path + ": not a halo catalog (size does not match)");
    }
  }
  const bool seq = fmt == CatalogFormat::Sequential;
  auto marker = [&](std::int32_t expect) {
    if (!seq) return;
    std::int32_t m = 0;
    in.read(reinterpret_cast<char*>(&m), 4);
    if (m != expect) throw std::runtime_error(path + ": bad record marker");
  };
  marker(8);
  in.read(reinterpret_cast<char*>(&n), 8);
  marker(8);
  if (!in || n < 0) throw std::runtime_error(path + ": cannot read Nhalo");
  std::vector<SourceRecord> h(static_cast<std::size_t>(n));
  for (auto& r : h) {
    marker(56);
    in.read(reinterpret_cast<char*>(&r), sizeof r);
    marker(56);
  }
  if (!in) throw std::runtime_error(path + ": truncated");
  return h;
}

SourceCatalog read_halo(const std::string& dir_halo, double redshift, std::ostream* log) {
  SourceCatalog cat;
  const auto hist = read_sfr_history(dir_halo + "sfr.txt");
  const CatalogSelection sel = select_halo_catalog(hist, redshift);
  if (!sel.found) return cat;  // Nhalo = Nsource = 0
  cat.redshift = sel.z2;
  cat.sfr_scale = sel.sfr_scale;
  cat.tag = sel.tag;
  const std::string fn = dir_halo + "halo." + sel.tag;
  if (log) *log << " Reading " << fn << "\n";
  cat.halo = read_halo_catalog(fn);
  for (auto& r : cat.halo) {
    if (r.L > 0) {
      ++cat.Nsource;
      r.L = static_cast<float>(static_cast<double>(r.L) * sel.sfr_scale);
      r.sfr = static_cast<float>(static_cast<double>(r.sfr) * sel.sfr_scale);
    } else {
      break;
    }
  }
  if (log) *log << " " << cat.halo.size() << " " << cat.Nsource << "\n";
  return cat;
}

// ---------------------------------------------------------------------------
// Source model
// ---------------------------------------------------------------------------

namespace {

// H(z)/H0 for the cosmology of cosmo_parameters.f90 (flat, w=-1)
double E_of_z(double z) {
  const double a = 1 / (1 + z);
  return std::sqrt(omegam / (a * a * a) + omegal + omegar / (a * a * a * a));
}

// dt/dz in years (absolute value): 1/((1+z) H(z))
double dtdz_yr(double z) {
  const double H = H0_cgs * hubble0 * E_of_z(z);
  return sec2yr / ((1 + z) * H);
}

double table_sfr(double z, const std::vector<SfrHistoryRow>& t) {
  if (t.empty()) return 0;
  // Linear interpolation in z between neighbouring rows (any ordering);
  // zero outside the tabulated range.
  for (std::size_t i = 0; i + 1 < t.size(); ++i) {
    const double za = t[i].z, zb = t[i + 1].z;
    if ((z - za) * (z - zb) <= 0 && za != zb)
      return t[i].sfr + (t[i + 1].sfr - t[i].sfr) * (z - za) / (zb - za);
  }
  for (const auto& r : t)
    if (r.z == z) return r.sfr;
  return 0;
}

}  // namespace

double global_sfr_density(double z, const SourceModelParams& p,
                          const std::vector<SfrHistoryRow>* table) {
  if (table) return table_sfr(z, *table);
  const double chi = std::pow(E_of_z(z), 2.0 / 3.0);
  return p.hs_rho_m * chi * chi /
         (1 + p.hs_alpha * std::pow(chi - 1, 3) * std::exp(p.hs_beta * std::pow(chi, 1.75)));
}

double global_stellar_density(double z, const SourceModelParams& p,
                              const std::vector<SfrHistoryRow>* table) {
  if (z >= p.z_max) return 0;
  // Trapezoid in ln(1+z), 2000 steps
  const int n = 2000;
  const double u0 = std::log(1 + z), u1 = std::log(1 + p.z_max);
  const double du = (u1 - u0) / n;
  double s = 0;
  for (int i = 0; i <= n; ++i) {
    const double zz = std::exp(u0 + i * du) - 1;
    const double f = global_sfr_density(zz, p, table) * dtdz_yr(zz) * (1 + zz);
    s += (i == 0 || i == n) ? 0.5 * f : f;
  }
  return s * du;
}

double cosmic_time_yr(double z, double z_max) {
  if (z >= z_max) return 0;
  const int n = 2000;
  const double u0 = std::log(1 + z), u1 = std::log(1 + z_max);
  const double du = (u1 - u0) / n;
  double s = 0;
  for (int i = 0; i <= n; ++i) {
    const double zz = std::exp(u0 + i * du) - 1;
    const double f = dtdz_yr(zz) * (1 + zz);
    s += (i == 0 || i == n) ? 0.5 * f : f;
  }
  return s * du;
}

double f_esc_at(const SourceModelParams& p, double z) {
  if (p.f_esc_beta == 0) return p.f_esc;
  return std::min(p.f_esc_max, p.f_esc * std::pow((1 + z) / (1 + p.f_esc_z0), p.f_esc_beta));
}

SourceCatalog assign_sources(const std::vector<Halo>& halos, const HaloUnits& hu,
                             const SourceModelParams& p,
                             const std::vector<SfrHistoryRow>* table) {
  const double fb = omegab / omegam;
  const double z = hu.redshift;
  const std::size_t n = halos.size();
  std::vector<double> sfr(n, 0), mstar(n, 0), Lrt(n, 0);
  std::vector<char> active(n, 0);
  for (std::size_t i = 0; i < n; ++i)
    active[i] = halos[i].T >= p.Tvir_min && halos[i].M >= p.M_min && halos[i].M > 0;

  if (p.model == SfrModel::GlobalMassProportional) {
    double Msum = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (active[i]) Msum += halos[i].M;
    const double vbox = std::pow(hu.box / hubble0, 3);  // comoving Mpc^3
    const double sfr_tot = global_sfr_density(z, p, table) * vbox;
    const double mstar_tot = global_stellar_density(z, p, table) * vbox;
    for (std::size_t i = 0; i < n; ++i) {
      if (!active[i] || Msum <= 0) continue;
      const double w = halos[i].M / Msum;
      sfr[i] = w * sfr_tot;
      mstar[i] = std::min(w * mstar_tot, fb * halos[i].M);
    }
  } else if (p.model == SfrModel::AbundanceMatched) {
    double xi_pop2 = 0;
    for (double v : pop2_init) xi_pop2 += v;
    const double rt_scale = std::pow(10.0, p.log_xi_ion) / xi_pop2 * (f_esc_at(p, z) / fResc);
    const double tage = cosmic_time_yr(z, p.z_max);
    for (std::size_t i = 0; i < n; ++i) {
      if (!active[i]) continue;
      const double Luv = p.am.luv(halos[i].M, z);
      sfr[i] = Luv / p.L_per_sfr;
      Lrt[i] = Luv * rt_scale;
      mstar[i] = std::min(sfr[i] * tage, fb * halos[i].M);
    }
  } else {
    const double dvir = delta_vir_bryan_norman(hu.omegam_z);
    const double rho_vir = dvir * hu.rho_crit;
    const double tdyn_yr = std::sqrt(3 * pi / (32 * G_cgs * rho_vir)) * sec2yr;
    const double tage = cosmic_time_yr(z, p.z_max);
    for (std::size_t i = 0; i < n; ++i) {
      if (!active[i]) continue;
      sfr[i] = p.eSF * fb * halos[i].M / tdyn_yr;
      mstar[i] = std::min(sfr[i] * tage, fb * halos[i].M);
    }
  }

  SourceCatalog cat;
  cat.redshift = z;
  cat.tag = z < 100 ? redshift_tag(z) : "";
  cat.halo.reserve(n);
  auto make = [&](std::size_t i) {
    const Halo& h = halos[i];
    SourceRecord r;
    for (int d = 0; d < 3; ++d) {
      r.x[d] = static_cast<float>(h.x[d]);
      r.v[d] = static_cast<float>(h.v[d]);
    }
    r.N = static_cast<std::int32_t>(std::min<std::int64_t>(h.N, INT32_MAX));
    r.M = static_cast<float>(h.M);
    r.R = static_cast<float>(h.R);
    r.T = static_cast<float>(h.T);
    r.sfr = static_cast<float>(sfr[i]);
    r.L = static_cast<float>(p.model == SfrModel::AbundanceMatched ? Lrt[i]
                                                                   : p.L_per_sfr * sfr[i]);
    r.mstar = static_cast<float>(mstar[i]);
    r.fgas = h.M > 0 ? static_cast<float>(std::max(0.0, fb - mstar[i] / h.M)) : 0.0f;
    return r;
  };
  for (std::size_t i = 0; i < n; ++i)
    if (sfr[i] > 0) cat.halo.push_back(make(i));
  for (auto& r : cat.halo)
    if (r.L > 0) ++cat.Nsource;
  // Records whose L underflowed to 0 in real(4) are not sources either
  std::stable_partition(cat.halo.begin(), cat.halo.end(),
                        [](const SourceRecord& r) { return r.L > 0; });
  for (std::size_t i = 0; i < n; ++i)
    if (!(sfr[i] > 0)) cat.halo.push_back(make(i));
  return cat;
}

// ---------------------------------------------------------------------------
// link_halo / add_halo_to_rt
// ---------------------------------------------------------------------------

namespace {
int grid_index(float x, double cellsize, int N) {
  // 1 + mod(int(x/cellsize), N)   (real(4)/real(8) -> real(8); int truncates)
  const std::int32_t c = f::int4(static_cast<double>(x) / cellsize);
  const std::int32_t m = c % N;
  if (m < 0) throw std::runtime_error("source at negative position (out of bounds in the original)");
  return 1 + m;
}
}  // namespace

std::vector<Source> link_halo(const std::vector<SourceRecord>& src, std::int64_t Nsource,
                              const Config& cfg) {
  std::vector<Source> out;
  out.reserve(static_cast<std::size_t>(Nsource));
  for (std::int64_t ih = 0; ih < Nsource; ++ih) {
    const SourceRecord& r = src[ih];
    Source s;
    s.i = grid_index(r.x[0], cfg.x_rt_box, cfg.Nrt);
    s.j = grid_index(r.x[1], cfg.x_rt_box, cfg.Nrt);
    s.k = grid_index(r.x[2], cfg.x_rt_box, cfg.Nrt);
    s.Luv = r.L;
    s.sfr = r.sfr;
    out.push_back(s);
  }
  return out;
}

void add_halo_metals(Simulation& s, const std::vector<SourceRecord>& src, std::int64_t Nsource) {
  const Config& cfg = s.cfg;
  // fZesc2 = 0.10, Zyield2 = 0.02: real(8) parameters from real(4) literals
  const double fZesc2 = static_cast<double>(0.10f);
  const double Zyield2 = static_cast<double>(0.02f);
  const std::vector<Source> cells = link_halo(src, Nsource, cfg);
  // Visit order of add_halo_to_rt: RT cell (k slowest, i fastest), then the
  // cell's chain, which is in record order.
  std::vector<std::int64_t> order(cells.size());
  for (std::size_t n = 0; n < order.size(); ++n) order[n] = static_cast<std::int64_t>(n);
  auto key = [&](std::int64_t n) {
    const Source& c = cells[n];
    return (c.i - 1) + static_cast<std::int64_t>(cfg.Nrt) *
                           ((c.j - 1) + static_cast<std::int64_t>(cfg.Nrt) * (c.k - 1));
  };
  std::stable_sort(order.begin(), order.end(),
                   [&](std::int64_t a, std::int64_t b) { return key(a) < key(b); });
  for (std::int64_t n : order) {
    const SourceRecord& r = src[n];
    const int i = grid_index(r.x[0], cfg.x_gas_box, cfg.Ngas);
    const int j = grid_index(r.x[1], cfg.x_gas_box, cfg.Ngas);
    const int k = grid_index(r.x[2], cfg.x_gas_box, cfg.Ngas);
    const double dmstar = static_cast<double>(r.sfr) * s.clk.dt_yr * Msun_cgs;
    const double rhoZ = fZesc2 * Zyield2 * dmstar / s.units.vol_gas;
    float* w = s.gas.w_p(s.gas.cell(i, j, k));
    w[3] = static_cast<float>(static_cast<double>(w[3]) + rhoZ / s.units.rho_unit);
  }
}

void calc_halo_records(Simulation& s, const std::vector<SourceRecord>& src,
                       std::int64_t Nsource, bool add_metals) {
  calc_halo(s, link_halo(src, Nsource, s.cfg));
  if (add_metals && !s.gas.w.empty()) add_halo_metals(s, src, Nsource);
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

SourceCatalog onthefly_catalog(const Config& cfg, const SimClock& clk, const Units& u,
                               const float* xv_dm, std::int64_t npart,
                               const HaloSourceParams& p) {
  if (!xv_dm && npart > 0) throw std::invalid_argument("onthefly_catalog: xv_dm is null");
  const HaloUnits hu = make_halo_units(cfg, clk, u);
  const auto halos = find_halos(xv_dm, npart, hu, p.finder);
  std::vector<SfrHistoryRow> table;
  if (!p.model.sfr_history_path.empty()) table = read_sfr_history(p.model.sfr_history_path);
  SourceCatalog cat = assign_sources(halos, hu, p.model,
                                     p.model.sfr_history_path.empty() ? nullptr : &table);
  if (!p.write_catalog_dir.empty() && !cat.tag.empty()) {
    std::filesystem::create_directories(p.write_catalog_dir);
    write_halo_catalog(p.write_catalog_dir + "/halo." + cat.tag, cat.halo);
  }
  return cat;
}

SourceCatalog update_sources(Simulation& s, const float* xv_dm, std::int64_t npart,
                             const HaloSourceParams& p, std::ostream* log) {
  SourceCatalog cat;
  if (p.mode == HaloSourceParams::Mode::Legacy) {
    cat = read_halo(p.dir_halo, s.clk.redshift, log);
  } else {
    cat = onthefly_catalog(s.cfg, s.clk, s.units, xv_dm, npart, p);
    if (log) *log << " Halos " << cat.halo.size() << " " << cat.Nsource << "\n";
  }
  s.Nsource = cat.Nsource;
  if (cat.Nsource > 0) calc_halo_records(s, cat.halo, cat.Nsource, p.add_metals);
  return cat;
}

}  // namespace radray

namespace radray {

bool set_halo_param(HaloSourceParams& p, const std::string& key, const std::string& val) {
  auto num = [&]() {
    std::string v = val;
    for (auto& ch : v)
      if (ch == 'd' || ch == 'D') ch = 'e';
    return std::stod(v);
  };
  if (key == "source_mode") {
    if (val == "legacy")
      p.mode = HaloSourceParams::Mode::Legacy;
    else if (val == "onthefly")
      p.mode = HaloSourceParams::Mode::OnTheFly;
    else
      throw std::runtime_error("source_mode must be legacy or onthefly: " + val);
  } else if (key == "dir_halo") {
    p.dir_halo = val;
    if (!p.dir_halo.empty() && p.dir_halo.back() != '/') p.dir_halo += '/';
  } else if (key == "fof_b") {
    p.finder.b = num();
  } else if (key == "fof_min_npart") {
    p.finder.min_npart = static_cast<std::int64_t>(num());
  } else if (key == "fof_so") {
    p.finder.so_refine = num() != 0;
  } else if (key == "sfr_model") {
    if (val == "global")
      p.model.model = SfrModel::GlobalMassProportional;
    else if (val == "dynamical")
      p.model.model = SfrModel::Dynamical;
    else if (val == "abundance")
      p.model.model = SfrModel::AbundanceMatched;
    else
      throw std::runtime_error("sfr_model must be global, dynamical or abundance: " + val);
  } else if (key == "tvir_min") {
    p.model.Tvir_min = num();
  } else if (key == "m_min") {
    p.model.M_min = num();
  } else if (key == "l_per_sfr") {
    p.model.L_per_sfr = num();
  } else if (key == "sfr_history") {
    p.model.sfr_history_path = val;
  } else if (key == "write_catalogs") {
    p.write_catalog_dir = val;
  } else if (key == "add_metals") {
    p.add_metals = num() != 0;
  } else if (key == "am_z0") {
    p.model.am.z0 = num();
  } else if (key == "am_zlo") {
    p.model.am.z_lo = num();
  } else if (key == "am_zhi") {
    p.model.am.z_hi = num();
  } else if (key == "xi_ion") {
    p.model.log_xi_ion = num();
  } else if (key == "f_esc") {
    p.model.f_esc = num();
  } else if (key == "f_esc_beta") {
    p.model.f_esc_beta = num();
  } else if (key == "f_esc_z0") {
    p.model.f_esc_z0 = num();
  } else if (key == "f_esc_max") {
    p.model.f_esc_max = num();
  } else if (key == "am_norm_z" || key == "am_norm_dl") {
    std::vector<double> v;
    std::stringstream ss(val);
    std::string t;
    while (std::getline(ss, t, ',')) {
      for (auto& ch : t)
        if (ch == 'd' || ch == 'D') ch = 'e';
      if (t.find_first_not_of(" \t") != std::string::npos) v.push_back(std::stod(t));
    }
    (key == "am_norm_z" ? p.model.am.norm_z : p.model.am.norm_dlogL) = v;
  } else if (key.rfind("am_", 0) == 0) {
    for (int i = 0; i < LuvModel::kNpar; ++i) {
      if (key == LuvModel::name(i)) {
        p.model.am.p[i] = num();
        return true;
      }
    }
    return false;
  } else {
    return false;
  }
  return true;
}

}  // namespace radray
