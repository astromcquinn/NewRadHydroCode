// Plane-parallel RT-domain mode (plane_parallel.hpp, docs/plane_parallel.md).
#include "radray/plane_parallel.hpp"

#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>

#include "radray/kernels/ray_kernels.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/rt_backend.hpp"
#include "radray/rt_moments.hpp"
#include "radray/simulation.hpp"

namespace radray {

namespace {

// photon number per unit frequency of J_nu ~ nu^-alpha: nu^-(alpha+1)
double photons_between(double alpha, double a, double b) {
  if (std::fabs(alpha) < 1e-12) return std::log(b / a);
  return (std::pow(a, -alpha) - std::pow(b, -alpha)) / alpha;
}
// sum of nu * photons: nu^-alpha
double nu_photons_between(double alpha, double a, double b) {
  const double p = 1 - alpha;
  if (std::fabs(p) < 1e-12) return std::log(b / a);
  return (std::pow(b, p) - std::pow(a, p)) / p;
}
// the frequency nu with photons_between(nu_min, nu) = target
double invert_photons(double alpha, double nu_min, double target) {
  if (std::fabs(alpha) < 1e-12) return nu_min * std::exp(target);
  return std::pow(std::pow(nu_min, -alpha) - alpha * target, -1.0 / alpha);
}

}  // namespace

PlaneSpectrum plane_spectrum(double alpha, double nu_min_eV, double nu_max_eV) {
  if (!(nu_max_eV > nu_min_eV) || !(nu_min_eV > 0))
    throw std::runtime_error("plane_spectrum: need 0 < nu_min < nu_max");
  PlaneSpectrum S;
  S.alpha = alpha;
  const double total = photons_between(alpha, nu_min_eV, nu_max_eV);
  S.edge[0] = nu_min_eV;
  S.edge[kNfreq] = nu_max_eV;
  for (int n = 1; n < kNfreq; ++n) S.edge[n] = invert_photons(alpha, nu_min_eV, total * n / kNfreq);
  double nu[kNfreq];
  for (int n = 0; n < kNfreq; ++n) {
    const double a = S.edge[n], b = S.edge[n + 1];
    S.frac[n] = photons_between(alpha, a, b) / total;
    S.nu[n] = nu_photons_between(alpha, a, b) / photons_between(alpha, a, b);
    nu[n] = S.nu[n];
  }
  double sig[3][kNfreq];
  rrk::photoion_sigma(nu, sig);
  S.sig_eff = 0;
  for (int n = 0; n < kNfreq; ++n) {
    S.sig_HI[n] = sig[0][n];
    S.sig_HeI[n] = sig[1][n];
    S.sig_HeII[n] = sig[2][n];
    S.sig_eff += S.frac[n] * S.sig_HI[n];
  }
  if (!(S.sig_eff > 0)) throw std::runtime_error("plane_spectrum: no HI-ionising photons");
  return S;
}

PlaneSpectrum plane_spectrum(const Config& c) {
  return plane_spectrum(c.plane_alpha, pc::HI_eV, c.plane_nu_max_ry * pc::HI_eV);
}

std::array<double, kNfreq> group_frequencies(const Config& c) {
  std::array<double, kNfreq> nu = cp::freq_init;
  if (c.plane_parallel) {
    const PlaneSpectrum S = plane_spectrum(c);
    for (int n = 0; n < kNfreq; ++n) nu[n] = S.nu[n];
  }
  return nu;
}

__int128 plane_to_fx(double v) {
  return static_cast<__int128>(std::ldexp(v, PlaneState::kFxBits));
}
double PlaneState::absorbed() const { return std::ldexp(static_cast<double>(absorbed_fx), -kFxBits); }
double PlaneState::escaped() const { return std::ldexp(static_cast<double>(escaped_fx), -kFxBits); }
double PlaneState::inflight() const { return std::ldexp(static_cast<double>(inflight_fx), -kFxBits); }
std::vector<int> plane_face_list(int faces) {
  std::vector<int> v;
  for (int f = 0; f < 6; ++f)
    if ((faces >> f) & 1) v.push_back(f);
  return v;
}
void PlaneState::reset_step() {
  emitted = 0;
  absorbed_fx = 0;
  escaped_fx = 0;
}

std::string plane_faces_string(int faces) {
  std::string out;
  const char* names[6] = {"+x", "-x", "+y", "-y", "+z", "-z"};
  for (int f = 0; f < 6; ++f)
    if ((faces >> f) & 1) out += (out.empty() ? "" : " ") + std::string(names[f]);
  return out;
}

bool set_plane_param(Config& c, const std::string& key, const std::string& val) {
  if (key == "plane_parallel") {
    if (val == "0")
      c.plane_parallel = false;
    else if (val == "1")
      c.plane_parallel = true;
    else
      throw std::runtime_error("plane_parallel must be 0 or 1");
    return true;
  }
  if (key == "plane_transport") {
    if (val == "rays" || val == "0")
      c.plane_transport = 0;
    else if (val == "sweep" || val == "1")
      c.plane_transport = 1;
    else
      throw std::runtime_error("plane_transport must be rays or sweep");
    return true;
  }
  if (key == "plane_ndom") {
    c.plane_ndom = std::stoi(val);
    if (c.plane_ndom < 1) throw std::runtime_error("plane_ndom must be >= 1");
    return true;
  }
  if (key == "plane_faces") {
    // "+x +y", "x,-y", "+x+y": any of +x -x +y -y +z -z, a bare letter is +
    int faces = 0;
    std::string t = val;
    for (char& ch : t)
      if (ch == ',') ch = ' ';
    std::istringstream in(t);
    std::string w;
    while (in >> w) {
      std::size_t p = 0;
      while (p < w.size()) {
        int sign = 0;
        if (w[p] == '+') {
          sign = 0;
          ++p;
        } else if (w[p] == '-') {
          sign = 1;
          ++p;
        }
        if (p >= w.size() || (w[p] != 'x' && w[p] != 'y' && w[p] != 'z'))
          throw std::runtime_error("plane_faces: give faces among +x -x +y -y +z -z, e.g. \"+x +y\"");
        faces |= 1 << (2 * (w[p] - 'x') + sign);
        ++p;
      }
    }
    if (faces == 0) throw std::runtime_error("plane_faces: no face given");
    c.plane_faces = faces;
    return true;
  }
  if (key == "plane_gamma12") {
    c.plane_gamma12 = std::stod(val);
    if (!(c.plane_gamma12 > 0)) throw std::runtime_error("plane_gamma12 must be > 0");
    return true;
  }
  if (key == "plane_z_re") {
    c.plane_z_re = std::stod(val);
    return true;
  }
  if (key == "plane_alpha") {
    c.plane_alpha = std::stod(val);
    return true;
  }
  if (key == "plane_nu_max_ry") {
    c.plane_nu_max_ry = std::stod(val);
    if (!(c.plane_nu_max_ry > 1)) throw std::runtime_error("plane_nu_max_ry must be > 1");
    return true;
  }
  if (key == "dt_max_myr") {
    c.dt_max_myr = std::stod(val);
    if (!(c.dt_max_myr >= 0)) throw std::runtime_error("dt_max_myr must be >= 0 (0 = off)");
    return true;
  }
  if (key == "rt_substeps_max") {
    c.rt_substeps_max = std::stoi(val);
    if (c.rt_substeps_max < 0) throw std::runtime_error("rt_substeps_max must be >= 0 (0 = off)");
    return true;
  }
  if (key == "plane_turnon_cfl") {
    if (val == "0")
      c.plane_turnon_cfl = 0;
    else if (val == "1")
      c.plane_turnon_cfl = 1;
    else
      throw std::runtime_error("plane_turnon_cfl must be 0 or 1");
    return true;
  }
  if (key == "plane_turnon_temp") {
    c.plane_turnon_temp = std::stod(val);
    if (!(c.plane_turnon_temp > 0)) throw std::runtime_error("plane_turnon_temp must be > 0 (K)");
    return true;
  }
  if (key == "plane_turnon_courant") {
    c.plane_turnon_courant = std::stod(val);
    if (!(c.plane_turnon_courant > 0 && c.plane_turnon_courant <= 1))
      throw std::runtime_error("plane_turnon_courant must be in (0, 1]");
    return true;
  }
  if (key == "land_outputs") {
    if (val == "0")
      c.land_outputs = false;
    else if (val == "1")
      c.land_outputs = true;
    else
      throw std::runtime_error("land_outputs must be 0 or 1");
    return true;
  }
  if (key == "photon_consume_dense") {
    if (val == "0")
      c.photon_consume_dense = false;
    else if (val == "1")
      c.photon_consume_dense = true;
    else
      throw std::runtime_error("photon_consume_dense must be 0 or 1");
    return true;
  }
  return false;
}

void copy_plane_param(const Config& from, Config& to) {
  to.plane_parallel = from.plane_parallel;
  to.plane_transport = from.plane_transport;
  to.plane_ndom = from.plane_ndom;
  to.plane_faces = from.plane_faces;
  to.plane_gamma12 = from.plane_gamma12;
  to.plane_z_re = from.plane_z_re;
  to.plane_alpha = from.plane_alpha;
  to.plane_nu_max_ry = from.plane_nu_max_ry;
  to.dt_max_myr = from.dt_max_myr;
  to.rt_substeps_max = from.rt_substeps_max;
  to.photon_consume_dense = from.photon_consume_dense;
  to.land_outputs = from.land_outputs;
  to.plane_turnon_cfl = from.plane_turnon_cfl;
  to.plane_turnon_temp = from.plane_turnon_temp;
  to.plane_turnon_courant = from.plane_turnon_courant;
}

void check_plane_param(const Config& c) {
  if (c.dt_max_myr < 0) throw std::runtime_error("dt_max_myr must be >= 0");
  if (c.rt_substeps_max < 0) throw std::runtime_error("rt_substeps_max must be >= 0");
  if (!c.plane_parallel) return;
  if (c.legacy_rt) throw std::runtime_error("plane_parallel = 1 requires legacy_rt = 0");
  if (use_moments(c)) throw std::runtime_error("plane_parallel = 1 is implemented for rt_method = rays only");
  if (c.recomb_rad) throw std::runtime_error("plane_parallel = 1 does not allow recomb_rad (its rays would split)");
  if (c.ray_encoding != RayEncoding::Float)
    throw std::runtime_error("plane_parallel = 1 needs ray_encoding = float (rstruct assumes the Pop II bands)");
  if (c.amrex_max_level > 0) throw std::runtime_error("plane_parallel = 1 runs on the uniform grid only (amrex_max_level = 0)");
  if (c.plane_transport != 0 && c.plane_transport != 1) throw std::runtime_error("plane_transport must be rays or sweep");
  if (c.plane_ndom < 1 || c.Nrt % c.plane_ndom != 0)
    throw std::runtime_error("plane_ndom must divide N (" + std::to_string(c.Nrt) + ")");
  const int L = c.Nrt / c.plane_ndom;
  if (L % c.rt_cells_per_substep != 0)
    throw std::runtime_error("plane_parallel: the domain length " + std::to_string(L) +
                             " cells must be a multiple of rt_cells_per_substep");
  if (c.plane_faces == 0) throw std::runtime_error("plane_faces: no face given");
  if (!(c.plane_gamma12 > 0)) throw std::runtime_error("plane_gamma12 must be > 0");
  if (!(c.plane_nu_max_ry > 1)) throw std::runtime_error("plane_nu_max_ry must be > 1");
  if (c.plane_transport == 1) return;  // sweep: no ray pool
  // the pool: nfaces * N^3 / k rays in flight once the fronts have crossed
  int nfaces = 0;
  for (int f = 0; f < 6; ++f) nfaces += (c.plane_faces >> f) & 1;
  const double inflight = static_cast<double>(nfaces) * c.nrtcr / c.rt_cells_per_substep;  // nrtcr = Nrt^3
  if (static_cast<double>(c.Nray_max) < 1.25 * inflight)
    throw std::runtime_error("plane_parallel: up to " + std::to_string(static_cast<long long>(inflight)) +
                             " rays in flight (" + std::to_string(nfaces) + " faces); raise nray_max_factor to at least " +
                             std::to_string(static_cast<int>(std::ceil(1.25 * nfaces / c.rt_cells_per_substep)) + 1));
}

std::string plane_summary(const Config& c) {
  if (!c.plane_parallel) {
    std::string s;
    if (c.dt_max_myr > 0) s += "dt_max_myr = " + std::to_string(c.dt_max_myr);
    if (c.rt_substeps_max > 0) s += (s.empty() ? "" : ", ") + std::string("rt_substeps_max = ") + std::to_string(c.rt_substeps_max);
    if (c.photon_consume_dense) s += (s.empty() ? "" : ", ") + std::string("photon_consume_dense = 1");
    if (c.land_outputs) s += (s.empty() ? "" : ", ") + std::string("land_outputs = 1");
    return s;
  }
  char buf[256];
  std::snprintf(buf, sizeof buf, "plane_parallel = 1 (%s; %d^3 domains of %d cells, faces %s, Gamma_-12 = %g, z_re = %g, alpha = %g, nu_max = %g Ry)",
                c.plane_transport ? "sweep" : "rays", c.plane_ndom, c.Nrt / c.plane_ndom,
                plane_faces_string(c.plane_faces).c_str(), c.plane_gamma12, c.plane_z_re, c.plane_alpha, c.plane_nu_max_ry);
  std::string s = buf;
  if (c.dt_max_myr > 0) s += ", dt_max_myr = " + std::to_string(c.dt_max_myr);
  if (c.rt_substeps_max > 0) s += ", rt_substeps_max = " + std::to_string(c.rt_substeps_max);
  if (c.photon_consume_dense) s += ", photon_consume_dense = 1";
  if (c.land_outputs) s += ", land_outputs = 1";
  if (!c.plane_turnon_cfl)
    s += ", plane_turnon_cfl = 0";
  else if (c.plane_turnon_temp != 2e4 || c.plane_turnon_courant != 1)
    s += ", plane_turnon_temp = " + std::to_string(c.plane_turnon_temp) +
         ", plane_turnon_courant = " + std::to_string(c.plane_turnon_courant);
  return s;
}

void plane_init(Simulation& s) {
  PlaneState& P = s.plane;
  P.on = s.cfg.plane_parallel;
  if (!P.on) return;
  check_plane_param(s.cfg);
  P.L = s.cfg.Nrt / s.cfg.plane_ndom;
  P.faces = s.cfg.plane_faces;
  P.transport = s.cfg.plane_transport;
  P.spec = plane_spectrum(s.cfg);
  P.reset_step();
  P.face_list = plane_face_list(P.faces);
  P.field.clear();
  P.inflight_fx = 0;
  if (P.transport == 1 && !use_amrex(s.cfg)) {
    // the serial sweep's packets (the backend holds its own fields)
    P.field.assign(P.face_list.size(), std::vector<float>(static_cast<std::size_t>(kNfreq) * s.cfg.Nrtc, 0.0f));
  }
}

void plane_fill_rtp(const Simulation& s, rrk::RTP& P) {
  const PlaneState& ps = s.plane;
  P.plane.L = 0;
  if (!ps.on || !ps.active) return;  // no face sources (and no plane rules) before the turn-on
  P.plane.L = ps.L;
  // sweep: the kernels' per-cell emission (plane_rays) is off; the faces
  // are read from the PlaneState by the sweep itself
  P.plane.faces = ps.transport == 1 ? 0 : ps.faces;
  const double gamma = s.cfg.plane_gamma12 * 1e-12;
  const double ntot = gamma * s.clk.dtsec_rt / (s.units.nH_avg * s.units.len_rt * ps.spec.sig_eff);
  for (int n = 0; n < kNfreq; ++n) {
    P.plane.nu[n] = static_cast<float>(ps.spec.nu[n]);
    P.plane.ngam[n] = static_cast<float>(ntot * ps.spec.frac[n]);
  }
}

double plane_emitted_per_substep(const Simulation& s, const rrk::RTP& P) {
  if (!s.plane.on || !s.plane.active) return 0;
  int nfaces = 0;
  for (int f = 0; f < 6; ++f) nfaces += (s.plane.faces >> f) & 1;
  const double N = s.cfg.nrtr;
  double per_ray = 0;
  for (int n = 0; n < kNfreq; ++n) per_ray += static_cast<double>(P.plane.ngam[n]);
  return nfaces * N * N * s.cfg.plane_ndom * per_ray;
}

bool plane_rt_active(const Config& c, double a_start) {
  if (!c.plane_parallel) return false;
  if (c.plane_z_re < 0) return true;
  const double a_re = 1 / (1 + c.plane_z_re);
  return a_start >= a_re * (1 - 1e-9);
}

void plane_write_spectrum(const std::string& dir, const Simulation& s) {
  if (!s.plane.on) return;
  const PlaneSpectrum& S = s.plane.spec;
  std::FILE* fp = std::fopen((dir + "/plane_spectrum.txt").c_str(), "w");
  if (!fp) throw std::runtime_error("cannot write plane_spectrum.txt in " + dir);
  const double gamma = s.cfg.plane_gamma12 * 1e-12;
  const double F = gamma / S.sig_eff;  // photons cm^-2 s^-1 through a face
  std::fprintf(fp, "# plane-parallel RT-domain mode (docs/plane_parallel.md)\n");
  std::fprintf(fp, "alpha = %.17g\n", S.alpha);
  std::fprintf(fp, "ndom = %d\nL_cells = %d\nfaces = %s\n", s.cfg.plane_ndom, s.plane.L,
               plane_faces_string(s.plane.faces).c_str());
  std::fprintf(fp, "gamma_HI = %.17g\n", gamma);
  std::fprintf(fp, "sig_eff = %.17g\n", S.sig_eff);
  std::fprintf(fp, "flux = %.17g\n", F);
  std::fprintf(fp, "# bin  edge_lo_eV  edge_hi_eV  nu_eff_eV  frac  sig_HI  sig_HeI  sig_HeII  flux_bin\n");
  for (int n = 0; n < kNfreq; ++n)
    std::fprintf(fp, "bin %d = %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", n, S.edge[n], S.edge[n + 1],
                 S.nu[n], S.frac[n], S.sig_HI[n], S.sig_HeI[n], S.sig_HeII[n], F * S.frac[n]);
  std::fclose(fp);
}

}  // namespace radray
