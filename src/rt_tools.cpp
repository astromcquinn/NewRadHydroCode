// Port of rt_tools.f90, domain_tools.f90 (RT domains) and
// halo_tools.f90 (calc_halo / add_halo_to_rt).
#include "radray/rt.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

#include "radray/constants.hpp"
#include "radray/fortran.hpp"
#include "radray/parallel.hpp"
#include "radray/ray_tracing.hpp"

namespace radray {

using namespace pc;
using namespace cp;

// ---------------------------------------------------------------------------
// Domains (domain_tools.f90)
// ---------------------------------------------------------------------------

// domain_tools.f90: init_domain_rt
void init_domain_rt(const Config& cfg, Domains& dom) {
  const int Nrt = cfg.Nrt;
  // n = nint(Ncpu**(1./3)) : integer ** real(4)
  const int n0 = f::nint(static_cast<double>(
      std::pow(static_cast<float>(cfg.Ncpu), 1.0f / 3.0f)));
  int n1 = 0;
  for (int i = n0; i <= Nrt; ++i) {
    if (Nrt % (2 * i) == 0) {
      n1 = 2 * i;
      break;
    }
  }
  if (n1 == 0) throw std::runtime_error("init_domain_rt: no valid subdivision (N odd?)");
  const int n2 = n1, n3 = n1;
  dom.Nx_dom = Nrt / n1;
  dom.Ndomain = n1 * n2 * n3;
  dom.box.assign(dom.Ndomain, {});
  dom.domain.assign(dom.Ndomain, {});
  for (int n = 1; n <= dom.Ndomain; ++n) {
    const int i = 1 + (n - 1) % n1;
    const int j = 1 + ((n - 1) / n1) % n2;
    const int k = 1 + ((n - 1) / n1 / n2);
    auto& b = dom.box[n - 1];
    b[0] = 1 + (i - 1) * dom.Nx_dom;
    b[1] = i * dom.Nx_dom;
    b[2] = 1 + (j - 1) * dom.Nx_dom;
    b[3] = j * dom.Nx_dom;
    b[4] = 1 + (k - 1) * dom.Nx_dom;
    b[5] = k * dom.Nx_dom;
    auto& d = dom.domain[n - 1];
    int m = 1;
    d[0] = 0;
    for (int cc = -1; cc <= 1; ++cc) {
      const int k1 = 1 + (k + cc - 1 + n3) % n3;
      for (int bb = -1; bb <= 1; ++bb) {
        const int j1 = 1 + (j + bb - 1 + n2) % n2;
        for (int aa = -1; aa <= 1; ++aa) {
          const int i1 = 1 + (i + aa - 1 + n1) % n1;
          const int l = i1 + (j1 - 1) * n1 + (k1 - 1) * n1 * n2;
          if (l != n) {
            m = m + 1;
            d[m - 1] = l;
          }
        }
      }
    }
  }
}

// domain_tools.f90: set_domain_rt (single thread: picks the first domain
// with status 0 whose neighbours are not in progress)
int set_domain_rt(Domains& dom) {
  for (;;) {
    for (int n = 1; n <= dom.Ndomain; ++n) {
      auto& d = dom.domain[n - 1];
      if (d[0] == 0) {
        bool eval = true;
        for (int m = 2; m <= 27; ++m) {
          const int l = d[m - 1];
          if (dom.domain[l - 1][0] == 1) {
            eval = false;
            break;
          }
        }
        if (eval) {
          d[0] = 1;
          return n - 1;
        }
      }
    }
    // With one thread this cannot happen (it would spin forever in Fortran).
    throw std::runtime_error("set_domain_rt: no eligible domain");
  }
}

// domain_tools.f90: end_domain_rt
void end_domain_rt(Domains& dom, int idom) {
  dom.domain[idom][0] = 2;
  long sum = 0;
  for (auto& d : dom.domain) sum += d[0];
  if (sum == 2L * dom.Ndomain)
    for (auto& d : dom.domain) d[0] = 0;
}

// deterministic_parallel domains (simulation.hpp: DetDomains)
void init_det_domains(const Config& cfg, DetDomains& dd) {
  const int N = cfg.Nrt;
  dd = DetDomains{};
  for (int d = 8; d <= N / 2; ++d)
    if (N % d == 0 && (N / d) % 2 == 0) {
      dd.d = d;
      break;
    }
  if (dd.d == 0)
    throw std::runtime_error("deterministic_parallel needs Nrt = 2*m*d with d >= 8 (Nrt >= 16)");
  dd.n1 = N / dd.d;
  for (int k = 0; k < dd.n1; ++k)
    for (int j = 0; j < dd.n1; ++j)
      for (int i = 0; i < dd.n1; ++i) {
        const int col = (i % 2) + 2 * (j % 2) + 4 * (k % 2);
        dd.colour[col].push_back({1 + i * dd.d, (i + 1) * dd.d, 1 + j * dd.d, (j + 1) * dd.d,
                                  1 + k * dd.d, (k + 1) * dd.d});
      }
}

namespace {

// deterministic_parallel: sweep one ray-tracing pass over all domains,
// colour by colour; the domains of a colour are independent and run
// concurrently. The result is that of the serial colour-major sweep.
template <class F>
void det_sweep(Simulation& s, F&& sweep_domain) {
  for (int col = 0; col < 8; ++col) {
    const auto& boxes = s.ddom.colour[col];
    const int nb = static_cast<int>(boxes.size());
    std::string err;
#pragma omp parallel for schedule(dynamic, 1) if (nb > 1)
    for (int b = 0; b < nb; ++b) {
      try {
        sweep_domain(boxes[b]);
      } catch (const std::exception& e) {
#pragma omp critical(radray_det_err)
        err = e.what();
      }
    }
    if (!err.empty()) throw std::runtime_error(err);
  }
}

// deterministic_parallel: fold the per-thread ray counters into ray_stat(1)
// (integer sums and maxima: order independent) and update the live count.
void det_fold(RayState& R) {
  RayStat& st = R.stat;
  for (ThreadRays& t : R.tl) {
    st.ncst1 += t.ncst1;
    st.ncst2 = f::max(st.ncst2, t.ncst2);
    st.nmrg1 += t.nmrg1;
    st.nmrg2 = f::max(st.nmrg2, t.nmrg2);
    st.nspl1 += t.nspl1;
    st.nspl2 = f::max(st.nspl2, t.nspl2);
    st.ndel1 += t.ndel1;
    R.Nray += t.ncreated - t.ndel1;
    t.ncst1 = t.ncst2 = t.nmrg1 = t.nmrg2 = t.nspl1 = t.nspl2 = t.ndel1 = t.ncreated = 0;
  }
}

}  // namespace

double wall_seconds() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// ---------------------------------------------------------------------------
// init
// ---------------------------------------------------------------------------

void RTGrid::allocate(int n) {
  N = n;
  ncell = static_cast<std::int64_t>(n) * n * n;
  for (auto* v : {&rhog, &fHI, &fHeI, &fHeII, &temp, &Luv, &sfr})
    v->assign(ncell, 0.0f);
  for (auto* v : {&dx, &gamma, &lambda}) v->assign(3 * ncell, 0.0f);
  for (auto* v : {&nu, &ngam, &fgam, &tau}) v->assign(kNfreq * ncell, 0.0f);
  htoc.assign(4 * ncell, 0);
}

// rt_tools.f90: init_rt
void init_rt(Simulation& s) {
  s.rt.allocate(s.cfg.Nrt);
  s.rt_save = RTStat{};
  s.rt_stat = s.rt_save;
  calc_atomic_rates(s.rates, s.cfg.legacy_physics);
}

void init_simulation(Simulation& s, int N, int nray_max_factor, double z_init) {
  s.cfg = Config::make(N, nray_max_factor);
  init_variables(z_init, s.clk);
  init_domain_rt(s.cfg, s.dom);
  s.gas.allocate(s.cfg.Ngas);
  init_rt(s);
  init_ray_tracing(s);
  calc_sim_var(s.cfg, s.clk, s.units);
}

// ---------------------------------------------------------------------------
// Sources: halo_tools.f90 calc_halo / add_halo_to_rt
// ---------------------------------------------------------------------------

void calc_halo(Simulation& s, const std::vector<Source>& sources) {
  RTGrid& G = s.rt;
  s.rt_save.Luv = 0;
  s.rt_save.sfr = 0;
  s.rt_stat.Luv = 0;
  s.rt_stat.sfr = 0;
  std::fill(G.Luv.begin(), G.Luv.end(), 0.0f);
  std::fill(G.sfr.begin(), G.sfr.end(), 0.0f);

  // hoc_source/ll_source: per cell, sources in ascending input order
  // (link_halo reverses twice). Luv/sfr summed in real(8).
  std::vector<std::vector<int>> per_cell(G.ncell);
  for (std::size_t n = 0; n < sources.size(); ++n) {
    const Source& src = sources[n];
    if (src.i < 1 || src.i > G.N || src.j < 1 || src.j > G.N || src.k < 1 || src.k > G.N)
      throw std::runtime_error("source outside the grid");
    per_cell[G.cell(src.i, src.j, src.k)].push_back(static_cast<int>(n));
  }
  // add_halo_to_rt loops k, j, i; rt_stat accumulates in that order.
  for (int k = 1; k <= G.N; ++k)
    for (int j = 1; j <= G.N; ++j)
      for (int i = 1; i <= G.N; ++i) {
        const std::int64_t c = G.cell(i, j, k);
        double Luv = 0, sfr = 0;
        for (int n : per_cell[c]) {
          Luv = Luv + static_cast<double>(sources[n].Luv);
          sfr = sfr + static_cast<double>(sources[n].sfr);
        }
        G.Luv[c] = static_cast<float>(Luv);
        G.sfr[c] = static_cast<float>(sfr);
        s.rt_stat.Luv = s.rt_stat.Luv + Luv;
        s.rt_stat.sfr = s.rt_stat.sfr + sfr;
      }
  const double vbox = f::pow3(s.cfg.box / hubble0);
  s.rt_save.Luv = s.rt_stat.Luv / vbox;
  s.rt_save.sfr = s.rt_stat.sfr / vbox;
  s.Nsource = static_cast<std::int64_t>(sources.size());
}

// ---------------------------------------------------------------------------
// radiative_transfer
// ---------------------------------------------------------------------------

// rt_tools.f90: radiative_transfer
void radiative_transfer(Simulation& s) {
  SimClock& clk = s.clk;
  std::ostream* L = s.raylog;

  // Reset timing variables to the start of the hydro step
  clk.year = clk.year - clk.dt_yr;
  clk.scalefactor = clk.scalefactor - clk.da_sim;
  clk.redshift = 1 / clk.scalefactor - 1;
  calc_sim_var(s.cfg, clk, s.units);

  if (L) *L << " Hydro step : " << clk.istep << " " << static_cast<float>(clk.redshift) << "\n";

  prep_rt(s);

  for (clk.istep_rt = 1; clk.istep_rt <= clk.nstep_rt; ++clk.istep_rt) {
    if (L) *L << " RT step    : " << clk.istep_rt << " " << static_cast<float>(clk.redshift) << "\n";

    rt_step_1(s);
    rt_step_2(s);
    rt_step_3(s);

    const std::int64_t nr = s.rays.Nray - s.rays.save.nfree;
    const float fnr = static_cast<float>(static_cast<double>(nr) /
                                         static_cast<double>(s.rays.Nray_max));
    if (L) {
      const RayStat& rs = s.rays.save;
      *L << " Sources    : " << s.rt_save.nsfill << " "
         << static_cast<float>(s.rt_save.nsfill / s.cfg.nrtcr) << "\n";
      *L << " Photons    : " << s.rt_save.ngfill << " "
         << static_cast<float>(s.rt_save.ngfill / s.cfg.nrtcr) << "\n";
      *L << " Rays       : " << nr << " " << fnr << "\n";
      *L << " Buffer min : " << s.rays.stat.nfree << " " << s.rays.Nrbuf << "\n";
      *L << " Buffer max : " << s.rays.stat.nfree << "\n";
      *L << " Cast       : " << rs.ncst1 << " " << rs.ncst2 << "\n";
      *L << " Merge      : " << rs.nmrg1 << " " << rs.nmrg2 << "\n";
      *L << " Split      : " << rs.nspl1 << " " << rs.nspl2 << "\n";
      *L << " Delete     : " << rs.ndel1 << "\n";
      if (clk.istep_rt == clk.nstep_rt) *L << "\n";
    }

    // Update timing variables
    clk.year = clk.year + clk.dt_yr / clk.nstep_rt;
    clk.scalefactor = clk.scalefactor + clk.da_sim / clk.nstep_rt;
    clk.redshift = 1 / clk.scalefactor - 1;
    calc_sim_var(s.cfg, clk, s.units);
    // (make_data_2d / write_data_3: IO only, not ported)
    if (s.substep_hook) s.substep_hook(s);
  }
  // Leave istep_rt = nstep_rt+1 like the Fortran DO variable
  post_rt(s);
}

// ---------------------------------------------------------------------------
// prep_rt
// ---------------------------------------------------------------------------

namespace {

// Average of the Eulerian cells inside an RT cell + derived gas state.
// Shared by calc_rt_i (prep_rt), calc_rt_f (post_rt) and the end of
// photo_cell, which contain the same code.
struct GasState {
  double rhog, nH, nHe, nHI, nHeI, nHeII, nHII, nHeIII, ne, neff, tmp;
};

GasState gas_state_from_avg(const double us_in[5], const double ws[4],
                            const Units& u, bool clamp_density) {
  GasState g;
  double us1 = us_in[0];
  if (clamp_density) us1 = f::max(us1, u.Dgas_min);
  g.rhog = us1 * u.rho_unit;
  g.nH = g.rhog * XHy / mH_cgs;
  g.nHe = g.rhog * YHe / mHe_cgs;
  g.nHI = f::min(ws[0] * u.H_unit, g.nH);
  g.nHeI = f::min(ws[1] * u.He_unit, g.nHe);
  g.nHeII = f::min(ws[2] * u.He_unit, g.nHe - g.nHeI);
  g.nHII = g.nH - g.nHI;
  g.nHeIII = g.nHe - g.nHeI - g.nHeII;
  g.ne = g.nHII + g.nHeII + 2 * g.nHeIII;
  g.neff = g.nH + g.nHe + g.ne;
  // eth = (us(5)-sum(us(2:4)**2)/us(1)/2)*eden_unit   (real(8))
  const double ksum = ((0.0 + us_in[1] * us_in[1]) + us_in[2] * us_in[2]) + us_in[3] * us_in[3];
  const double eth = (us_in[4] - ksum / us1 / 2) * u.eden_unit;
  g.tmp = f::max((gamma_ideal - 1) * eth / k_cgs / g.neff, u.temp_min);
  return g;
}

void store_gas_state(RTGrid& G, std::int64_t c, const GasState& g, const Units& u) {
  G.rhog[c] = static_cast<float>(g.rhog / u.rhob_avg);
  G.fHI[c] = static_cast<float>(f::min(f::max(g.nHI / g.nH, 0.0), 1.0));
  G.fHeI[c] = static_cast<float>(f::min(f::max(g.nHeI / g.nHe, 0.0), 1.0));
  G.fHeII[c] = static_cast<float>(f::min(f::max(g.nHeII / g.nHe, 0.0), 1.0));
  G.temp[c] = static_cast<float>(g.tmp);
}

// Sum u_gas/w_gas over the Eulerian cells of RT cell (irt,jrt,krt) and
// normalise (us = us/N_rt_gas_c, us(1) = max(us(1),Dgas_min)).
void average_gas(const Simulation& s, int irt, int jrt, int krt, double us[5],
                 double ws[4]) {
  const int Nr = s.cfg.N_rt_gas;
  for (int m = 0; m < 5; ++m) us[m] = 0;
  for (int m = 0; m < 4; ++m) ws[m] = 0;
  for (int k = 1 + (krt - 1) * Nr; k <= krt * Nr; ++k)
    for (int j = 1 + (jrt - 1) * Nr; j <= jrt * Nr; ++j)
      for (int i = 1 + (irt - 1) * Nr; i <= irt * Nr; ++i) {
        const std::int64_t gc = s.gas.cell(i, j, k);
        const float* uu = &s.gas.u[gc * 5];
        const float* ww = &s.gas.w[gc * 4];
        for (int m = 0; m < 5; ++m) us[m] = us[m] + static_cast<double>(uu[m]);
        for (int m = 0; m < 4; ++m) ws[m] = ws[m] + static_cast<double>(ww[m]);
      }
  for (int m = 0; m < 5; ++m) us[m] = us[m] / s.cfg.N_rt_gas_c;
  for (int m = 0; m < 4; ++m) ws[m] = ws[m] / s.cfg.N_rt_gas_c;
  us[0] = f::max(us[0], s.units.Dgas_min);
}

}  // namespace

// rt_tools.f90: prep_rt (calc_rt_i + ray buffer)
void prep_rt(Simulation& s) {
  const double tw0 = wall_seconds();
  struct Acc {
    double& t;
    double t0;
    ~Acc() { t += wall_seconds() - t0; }
  } acc{s.timers.prep, tw0};
  s.rt_save.dtau = 0;
  s.rt_save.gamma = {};
  s.rt_save.lambda = {};
  s.rt_save.nsfill = 0;
  s.rt_save.ngfill = 0;
  // (random_seed: random numbers are never used in the RT, `.false.` branch)

  RTGrid& G = s.rt;
  const int Nrt = s.cfg.Nrt;
#pragma omp parallel for schedule(static) if (par::worth(s.cfg.Nrtc))
  for (int krt = 1; krt <= Nrt; ++krt)
    for (int jrt = 1; jrt <= Nrt; ++jrt)
      for (int irt = 1; irt <= Nrt; ++irt) {
        double us[5], ws[4];
        average_gas(s, irt, jrt, krt, us, ws);
        const GasState g = gas_state_from_avg(us, ws, s.units, false);
        store_gas_state(G, G.cell(irt, jrt, krt), g, s.units);
      }

  if (s.cfg.det_parallel) {
    // deterministic_parallel: per-thread free lists instead of the Fortran
    // buffers (Nrbuf, ray_buffer_*); Nray is the live count.
    if (s.cfg.legacy_rt) throw std::runtime_error("deterministic_parallel requires legacy_rt = 0");
    if (s.ddom.d == 0) init_det_domains(s.cfg, s.ddom);
    const int nt = par::max_threads();
    if (static_cast<int>(s.rays.tl.size()) < nt) s.rays.tl.resize(nt);
    s.rays.save.nfree = 0;
    s.rays.stat.nfree = 0;
    return;
  }

  // Nrbuf = min(max(int((0.98*Nray_max-Nray)/Ncpu,kind=8),0),Nfree_max/4)
  // 0.98*Nray_max is real(4) arithmetic.
  RayState& R = s.rays;
  float t = 0.98f * static_cast<float>(R.Nray_max);
  t = t - static_cast<float>(R.Nray);
  t = t / static_cast<float>(s.cfg.Ncpu);
  R.Nrbuf = f::min(f::max<std::int64_t>(f::int8(t), 0), R.Nfree_max / 4);
  ray_buffer_allocate(s);
}

// ---------------------------------------------------------------------------
// post_rt
// ---------------------------------------------------------------------------

// rt_tools.f90: post_rt
void post_rt(Simulation& s) {
  struct Acc {
    double& t;
    double t0;
    ~Acc() { t += wall_seconds() - t0; }
  } acc{s.timers.post, wall_seconds()};
  // deterministic_parallel: there are no Fortran buffers to return, and
  // freed slots are already zero.
  if (!s.cfg.det_parallel) ray_buffer_deallocate(s);
  ray_redshifting(s);

  RTGrid& G = s.rt;
  RayState& R = s.rays;
  const int Nrt = s.cfg.Nrt;
  const SigmaTable sig = photoion_sigma(freq_init);

  std::int64_t nray_pdf[24] = {0};
  // Fortran: dimension(-9:10), but 1E-9 is real(4) (9.99999972e-10) so empty
  // cells land in bin -10 -> out-of-bounds write in the original. Like the
  // harness (PATCHES.md 5) we use -10:10. Diagnostics only.
  std::int64_t ngam_pdf[21] = {0};  // index -10..10 -> 0..20
  std::int64_t nrad_pdf[21] = {0};

  // calc_rt_f
#pragma omp parallel for schedule(static) if (par::worth(s.cfg.Nrtc))
  for (int krt = 1; krt <= Nrt; ++krt)
    for (int jrt = 1; jrt <= Nrt; ++jrt)
      for (int irt = 1; irt <= Nrt; ++irt) {
        const std::int64_t c = G.cell(irt, jrt, krt);
        double us[5], ws[4];
        average_gas(s, irt, jrt, krt, us, ws);
        const GasState g = gas_state_from_avg(us, ws, s.units, false);
        store_gas_state(G, c, g, s.units);
        float* tau = G.tau_p(c);
        for (int i = 0; i < kNfreq; ++i) {
          const double sm = ((0.0 + g.nHI * sig[0][i]) + g.nHeI * sig[1][i]) + g.nHeII * sig[2][i];
          tau[i] = static_cast<float>(sm * s.units.len_rt);
        }
      }

  // calc_ray_f (diagnostic PDFs for ray.log)
  auto bin = [](double ngam) {
    ngam = f::min(f::max(ngam, static_cast<double>(1e-9f)), static_cast<double>(1e1f));
    int i;
    if (ngam < 0.5)
      i = f::floor4(std::log10(ngam));
    else
      i = f::floor4(ngam);
    return i + 10;
  };
  // integer histograms: per-thread partial counts, summed (order independent)
#pragma omp parallel if (par::worth(s.cfg.Nrtc))
  {
  std::int64_t nray_pdf_t[24] = {0}, ngam_pdf_t[21] = {0}, nrad_pdf_t[21] = {0};
  auto& nray_pdf_l = nray_pdf_t;
  auto& ngam_pdf_l = ngam_pdf_t;
  auto& nrad_pdf_l = nrad_pdf_t;
#pragma omp for schedule(static)
  for (int krt = 1; krt <= Nrt; ++krt)
    for (int jrt = 1; jrt <= Nrt; ++jrt)
      for (int irt = 1; irt <= Nrt; ++irt) {
        auto& ngam_pdf = ngam_pdf_l;
        auto& nrad_pdf = nrad_pdf_l;
        auto& nray_pdf = nray_pdf_l;
        const std::int64_t c = G.cell(irt, jrt, krt);
        ngam_pdf[bin(static_cast<double>(f::sum4<kNfreq>(G.ngam_p(c))))]++;
        std::int64_t nr = 0;
        std::int64_t ir = G.htoc_ray(1, c);
        while (ir > 0) {
          nr = nr + 1;
          nrad_pdf[bin(static_cast<double>(f::sum4<kNfreq>(R.at(ir).ngam)))]++;
          ir = R.at(ir).next;
        }
        if (static_cast<double>(G.sfr[c]) > s.units.sfr_min)
          nr = f::max<std::int64_t>(nr, rayp::nray_init);
        int b;
        if (nr <= 8) b = static_cast<int>(nr);
        else if (nr <= 12) b = 9;
        else if (nr <= 16) b = 10;
        else if (nr <= 24) b = 11;
        else if (nr <= 32) b = 12;
        else if (nr <= 48) b = 13;
        else if (nr <= 64) b = 14;
        else if (nr <= 96) b = 15;
        else if (nr <= 128) b = 16;
        else if (nr <= 192) b = 17;
        else if (nr <= 256) b = 18;
        else if (nr <= 384) b = 19;
        else if (nr <= 512) b = 20;
        else if (nr <= 768) b = 21;
        else if (nr <= 1024) b = 22;
        else b = 23;
        nray_pdf[b]++;
      }
#pragma omp critical(radray_calc_ray_f)
  {
    for (int b = 0; b < 24; ++b) nray_pdf[b] += nray_pdf_t[b];
    for (int b = 0; b < 21; ++b) ngam_pdf[b] += ngam_pdf_t[b];
    for (int b = 0; b < 21; ++b) nrad_pdf[b] += nrad_pdf_t[b];
  }
  }

  if (std::ostream* L = s.raylog) {
    *L << " PDF(nr)\n";
    for (int b = 0; b < 24; ++b) *L << "   " << b << " : " << nray_pdf[b] << "\n";
    *L << " PDF(ngam_RT)\n";
    for (int b = 1; b < 21; ++b) *L << "   " << (b - 10) << " : " << ngam_pdf[b] << "\n";
    *L << " PDF(ngam_ray)\n";
    for (int b = 1; b < 21; ++b) *L << "   " << (b - 10) << " : " << nrad_pdf[b] << "\n";
  }
}

// ---------------------------------------------------------------------------
// rt_step_1
// ---------------------------------------------------------------------------

// rt_tools.f90: rt_step_1
void rt_step_1(Simulation& s) {
  double tw = wall_seconds();
  auto lap = [&](double& t) {
    const double now = wall_seconds();
    t += now - tw;
    tw = now;
  };
  s.timers.substeps += 1;
  RTGrid& G = s.rt;
  RayState& R = s.rays;
  const Units& u = s.units;
  const int Nrt = s.cfg.Nrt;

  const double lgn_min = -6, lgn_max = +3;
  const double lgn_del = folded::lgn_del;
  const int N_pdf = f::nint((lgn_max - lgn_min) / lgn_del);

  // Initialize RT and ray variables (rt_stat(1), ray_stat(1))
  {
    RTStat& t = s.rt_stat;
    t.dtau = 0;
    t.ngam = {};
    t.nion = {};
    t.ncol = {};
    t.nrec = {};
    t.nrab = {};
    t.gamma = {};
    t.lambda = {};
    t.nsfill = 0;
    t.ngfill = 0;
    RayStat& r = R.stat;
    r.ncst1 = r.ncst2 = r.nmrg1 = r.nmrg2 = r.nspl1 = r.nspl2 = r.ndel1 = 0;
  }
  const bool det = s.cfg.det_parallel;

  // calc_rt_1: reset radiation fields, update optical depth of (partially)
  // ionised cells. Note: with legacy_rt, rt_grid%dx is NOT reset here
  // (original behaviour, B3) and tau is only updated where fHI or fHeI < 0.99
  // (B4).
  const bool legacy = s.cfg.legacy_rt;
  const SigmaTable sig = photoion_sigma(freq_init);
#pragma omp parallel for schedule(static) if (par::worth(s.cfg.Nrtc))
  for (int krt = 1; krt <= Nrt; ++krt)
    for (int jrt = 1; jrt <= Nrt; ++jrt)
      for (int irt = 1; irt <= Nrt; ++irt) {
        const std::int64_t c = G.cell(irt, jrt, krt);
        for (int n = 0; n < kNfreq; ++n) {
          G.nu[c * kNfreq + n] = 0;
          G.ngam[c * kNfreq + n] = 0;
          G.fgam[c * kNfreq + n] = 0;
        }
        for (int d = 0; d < 3; ++d) {
          G.gamma[c * 3 + d] = 0;
          G.lambda[c * 3 + d] = 0;
        }
        // B3 fix: also reset the photon-weighted direction field
        if (!legacy)
          for (int d = 0; d < 3; ++d) G.dx[c * 3 + d] = 0;
        // real(4) comparisons with real(4) literal 0.99.
        // B4 fix: without legacy_rt, update tau in every cell.
        if (G.fHI[c] < 0.99f || G.fHeI[c] < 0.99f || !legacy) {
          const double rhog = static_cast<double>(G.rhog[c]) * u.rhob_avg;
          const double nH = rhog * XHy / mH_cgs;
          const double nHe = rhog * YHe / mHe_cgs;
          const double nHI = nH * static_cast<double>(G.fHI[c]);
          const double nHeI = nHe * static_cast<double>(G.fHeI[c]);
          const double nHeII = nHe * static_cast<double>(G.fHeII[c]);
          float* tau = G.tau_p(c);
          for (int i = 0; i < kNfreq; ++i) {
            const double sm = ((0.0 + nHI * sig[0][i]) + nHeI * sig[1][i]) + nHeII * sig[2][i];
            tau[i] = static_cast<float>(sm * u.len_rt);
          }
        }
      }

  lap(s.timers.calc_rt1);
  // calc_ray_1: histogram of photon counts over ALL storage slots with
  // nside > 0 (includes freed rays that have not been zeroed yet).
  double ngam_avg = 0;
  std::vector<std::int64_t> ngam_pdf(N_pdf + 1, 0);  // 1-based
  // F3 guard (legacy_rt = 0): the same histogram for the rays that may split
  // this substep (dr < 0, nside < nside_max)
  std::vector<std::int64_t> spl_pdf(N_pdf + 1, 0);
  const std::int64_t top = R.hwm();
  if (!det) {
    for (std::int64_t ir = 1; ir <= top; ++ir) {
      const Ray& ry = R.at(ir);
      if (ry.nside > 0) {
        const double ntot = static_cast<double>(f::sum4<kNfreq>(ry.ngam));
        int n = 1 + f::int4((std::log10(ntot) - lgn_min) / lgn_del);
        n = f::min(f::max(n, 1), N_pdf);
        ngam_avg = ngam_avg + ntot;
        ngam_pdf[n] = ngam_pdf[n] + 1;
        if (ry.dr < 0 && ry.nside < rayp::nside_max) spl_pdf[n] = spl_pdf[n] + 1;
      }
    }
  } else {
    // deterministic_parallel: the slots of the live rays depend on thread
    // timing, so the histogram is summed per thread (integers) and ngam_avg
    // (ray.log only) in exact 2^-40 fixed point, both order independent.
    __int128 fix = 0;
#pragma omp parallel if (par::worth(top))
    {
      std::vector<std::int64_t> pdf_t(N_pdf + 1, 0), spl_t(N_pdf + 1, 0);
      __int128 fix_t = 0;
#pragma omp for schedule(static)
      for (std::int64_t ir = 1; ir <= top; ++ir) {
        const Ray& ry = R.at(ir);
        if (ry.nside > 0) {
          const double ntot = static_cast<double>(f::sum4<kNfreq>(ry.ngam));
          int n = 1 + f::int4((std::log10(ntot) - lgn_min) / lgn_del);
          n = f::min(f::max(n, 1), N_pdf);
          fix_t += static_cast<__int128>(std::ldexp(ntot, 40));
          pdf_t[n] = pdf_t[n] + 1;
          if (ry.dr < 0 && ry.nside < rayp::nside_max) spl_t[n] = spl_t[n] + 1;
        }
      }
#pragma omp critical(radray_calc_ray_1)
      {
        for (int n = 1; n <= N_pdf; ++n) ngam_pdf[n] += pdf_t[n];
        for (int n = 1; n <= N_pdf; ++n) spl_pdf[n] += spl_t[n];
        fix += fix_t;
      }
    }
    ngam_avg = std::ldexp(static_cast<double>(fix), -40);
  }

  // Calculate ngam_ray. With legacy_rt, thresholds that are not reached keep
  // their stale values (B2); otherwise they are recomputed from scratch.
  if (ngam_avg == 0) {
    R.ngam_ray.fill(0.0);
  } else {
    ngam_avg = ngam_avg / R.Nray;
    if (!legacy) R.ngam_ray.fill(0.0);  // B2 fix
    // legacy: int(Nray_max/64.) is a real(4) division (int64 product here;
    // the Fortran's int32 product overflows only for N >= ~323).
    // B9 fix: exact 64-bit integer Nray_max/64.
    const std::int64_t per =
        legacy ? f::int4(static_cast<double>(static_cast<float>(R.Nray_max) / 64.0f))
               : R.Nray_max / 64;
    for (int j = 1; j <= 64; ++j) {
      std::int64_t nr = 0;
      for (int i = N_pdf; i >= 1; --i) {
        nr = nr + ngam_pdf[i];
        if (nr > j * per) {
          R.ngam_ray[j - 1] = std::pow(1e1, (i - 1) * lgn_del + lgn_min);
          break;
        }
      }
      if (s.raylog && j % 4 == 0)
        *s.raylog << " Rays       : " << j << " " << static_cast<float>(R.ngam_ray[j - 1])
                  << " " << static_cast<float>(R.ngam_ray[j - 1] / ngam_avg) << "\n";
    }
  }

  // F3 (legacy_rt = 0 only; RT harness PATCHES.md 13): keep the ray pool from
  // overflowing. Rays are only created by splits (3 per split) and sources
  // (nray_init per source cell). From the live count L, the number of source
  // cells and the histogram of the rays that may split, limit the splits of
  // this substep so that even if every candidate splits at most 0.9 Nray_max
  // rays are live after it (the rest of the pool covers the free buffer,
  // <= Nray_max/16): ngam_ray(32) (the split threshold) is raised to the
  // lower edge of the brightest bins whose candidates fit in a third of the
  // remaining budget. If even the sources do not fit, splitting stops and
  // every ray except fresh source rays is merged (ngam_ray(16) =
  // ngam_ray(32) = 1e30). Integer arithmetic on the histogram only:
  // deterministic, the same in the Fortran harness, and inactive unless the
  // pool could overflow in this substep.
  if (!legacy) {
    std::int64_t L = 0;
    for (int n = 1; n <= N_pdf; ++n) L += ngam_pdf[n];
    std::int64_t nsrc = 0;
    for (std::int64_t c = 0; c < G.ncell; ++c)
      if (static_cast<double>(G.sfr[c]) > u.sfr_min) ++nsrc;
    const std::int64_t budget = (9 * R.Nray_max) / 10 - L - rayp::nray_init * nsrc;
    if (budget <= 0) {
      R.ngam_ray[16 - 1] = 1e30;
      R.ngam_ray[32 - 1] = 1e30;
      R.npool_forced += 1;
    } else {
      const std::int64_t smax = budget / 3;
      std::int64_t nr = 0;
      int i = N_pdf;
      for (; i >= 1; --i) {
        if (nr + spl_pdf[i] > smax) break;
        nr = nr + spl_pdf[i];
      }
      if (i >= 1) {
        const double thr = std::pow(1e1, i * lgn_del + lgn_min);
        if (thr > R.ngam_ray[32 - 1]) {
          R.ngam_ray[32 - 1] = thr;
          R.npool_limited += 1;
        }
      }
    }
    if (s.raylog && (R.npool_forced > 0 || R.npool_limited > 0))
      *s.raylog << " Pool guard : " << L << " " << nsrc << " " << budget << " "
                << R.npool_limited << " " << R.npool_forced << "\n";
  }

  lap(s.timers.calc_ray1);
  if (det) {
    // deterministic_parallel: colour-major sweep over the 8^3 domains
    det_sweep(s, [&](const std::array<int, 6>& box) { ray_tracing_1(s, box); });
    det_fold(R);
  } else {
    if (s.clk.istep_rt > 1) ray_buffer_supplement(s);

    // Ray casting, domain by domain (order 1..Ndomain with one thread)
    for (int n = 1; n <= s.dom.Ndomain; ++n) {
      const int i = set_domain_rt(s.dom);
      ray_tracing_1(s, s.dom.box[i]);
      end_domain_rt(s.dom, i);
    }
  }

  lap(s.timers.trace1);
  // Synchronise htoc_ray(1:2) = htoc_ray(3:4)
#pragma omp parallel for schedule(static) if (par::worth(s.cfg.Nrtc))
  for (std::int64_t c = 0; c < G.ncell; ++c) {
    G.htoc[c * 4 + 0] = G.htoc[c * 4 + 2];
    G.htoc[c * 4 + 1] = G.htoc[c * 4 + 3];
  }
  lap(s.timers.sync);
}

// ---------------------------------------------------------------------------
// rt_step_2 / radiative_physics
// ---------------------------------------------------------------------------

namespace {

// Accumulation of the rt_stat sums of radiative_physics. Direct mode adds
// each addend to rt_stat at once (the original); record mode stores the
// (accumulator, addend) sequence of one k-plane, which rt_step_2 replays in
// plane order after a parallel batch. Every accumulator then receives the
// same addends in the same order, so the double sums are bit-identical to
// the serial loop for any thread count.
struct RPSink {
  enum { kNgam = 0, kNion = 3, kNcol = 6, kNrec = 9, kNrab = 12, kGamma = 15, kLambda = 18,
         kDtau = 21, kNacc = 22 };
  double* direct[kNacc] = {};  // non-null: direct mode
  std::vector<std::uint8_t> id;
  std::vector<double> val;
  std::int64_t ngfill = 0, nsfill = 0;

  void bind(RTStat& t) {
    for (int sp = 0; sp < 3; ++sp) {
      direct[kNgam + sp] = &t.ngam[sp];
      direct[kNion + sp] = &t.nion[sp];
      direct[kNcol + sp] = &t.ncol[sp];
      direct[kNrec + sp] = &t.nrec[sp];
      direct[kNrab + sp] = &t.nrab[sp];
      direct[kGamma + sp] = &t.gamma[sp];
      direct[kLambda + sp] = &t.lambda[sp];
    }
    direct[kDtau] = &t.dtau;
  }
  void add(int a, double v) {
    if (direct[0]) {
      *direct[a] = *direct[a] + v;
    } else {
      id.push_back(static_cast<std::uint8_t>(a));
      val.push_back(v);
    }
  }
  void replay(RTStat& t) {
    RPSink d;
    d.bind(t);
    for (std::size_t n = 0; n < id.size(); ++n) *d.direct[id[n]] = *d.direct[id[n]] + val[n];
    t.ngfill += ngfill;
    t.nsfill += nsfill;
    id.clear();
    val.clear();
    ngfill = nsfill = 0;
  }
};

}  // namespace

void radiative_physics_plane(Simulation& s, int krt, RPSink& sink);

// rt_tools.f90: radiative_physics for one k-plane, accumulating directly
// into rt_stat(1)
void radiative_physics(Simulation& s, int krt) {
  RPSink sink;
  sink.bind(s.rt_stat);
  radiative_physics_plane(s, krt, sink);
  s.rt_stat.ngfill += sink.ngfill;
  s.rt_stat.nsfill += sink.nsfill;
}

// rt_tools.f90: rt_step_2. The k-planes are independent (each touches only
// its own RT and gas cells); with more than one thread they run in batches
// whose rt_stat addends are replayed in plane order (RPSink).
void rt_step_2(Simulation& s) {
  struct Acc {
    double& t;
    double t0;
    ~Acc() { t += wall_seconds() - t0; }
  } acc{s.timers.physics, wall_seconds()};
  const int N = s.cfg.Nrt;
  const int nt = par::max_threads();
  if (nt <= 1 || !par::worth(s.cfg.Nrtc)) {
    for (int k = 1; k <= N; ++k) radiative_physics(s, k);
    return;
  }
  const int batch = std::min(N, 4 * nt);
  std::vector<RPSink> sinks(batch);
  for (int k0 = 1; k0 <= N; k0 += batch) {
    const int nb = std::min(batch, N - k0 + 1);
    std::string err;
#pragma omp parallel for schedule(dynamic, 1)
    for (int b = 0; b < nb; ++b) {
      try {
        radiative_physics_plane(s, k0 + b, sinks[b]);
      } catch (const std::exception& e) {
#pragma omp critical(radray_rp_err)
        if (err.empty()) err = e.what();
      }
    }
    if (!err.empty()) throw std::runtime_error(err);
    for (int b = 0; b < nb; ++b) sinks[b].replay(s.rt_stat);
  }
}

namespace {

// rt_tools.f90: radiative_physics / make_index_gas_rt.
// Returns the Eulerian cells of the RT cell ordered by distance from the
// source (density peak) or from the upstream side of the radiation.
void make_index_gas_rt(const Simulation& s, int irt, int jrt, int krt,
                       std::vector<std::array<int, 3>>& out) {
  const int Nr = s.cfg.N_rt_gas;
  const RTGrid& G = s.rt;
  const std::int64_t c = G.cell(irt, jrt, krt);
  out.clear();
  if (Nr == 1) {
    out.push_back({irt, jrt, krt});
    return;
  }
  // General case (N_rt_gas > 1); not exercised by the validation problems.
  float x0[3] = {0, 0, 0};
  if (static_cast<double>(G.sfr[c]) > s.units.sfr_min) {
    const int i1 = 1 + (irt - 1) * Nr, j1 = 1 + (jrt - 1) * Nr, k1 = 1 + (krt - 1) * Nr;
    float rho = 0;
    for (int k = k1; k <= krt * Nr; ++k)
      for (int j = j1; j <= jrt * Nr; ++j)
        for (int i = i1; i <= irt * Nr; ++i) {
          const float uu = s.gas.u[s.gas.cell(i, j, k) * 5];
          if (uu > rho) {
            x0[0] = static_cast<float>(i - i1) + 0.5f;
            x0[1] = static_cast<float>(j - j1) + 0.5f;
            x0[2] = static_cast<float>(k - k1) + 0.5f;
            rho = uu;
          }
        }
  } else {
    for (int d = 0; d < 3; ++d)
      x0[d] = 0.5f * static_cast<float>(Nr) - G.dx[c * 3 + d] * static_cast<float>(Nr);
  }
  std::vector<float> rsq;
  std::vector<std::array<int, 3>> ijk;
  for (int k = 1; k <= Nr; ++k)
    for (int j = 1; j <= Nr; ++j)
      for (int i = 1; i <= Nr; ++i) {
        const float x1[3] = {static_cast<float>(i) - 0.5f, static_cast<float>(j) - 0.5f,
                             static_cast<float>(k) - 0.5f};
        float r2 = 0;
        for (int d = 0; d < 3; ++d) r2 += (x1[d] - x0[d]) * (x1[d] - x0[d]);
        rsq.push_back(r2);
        ijk.push_back({i, j, k});
      }
  std::vector<int> idx(rsq.size());
  for (std::size_t n = 0; n < idx.size(); ++n) idx[n] = static_cast<int>(n);
  // NOTE: Fortran uses the C library qsort (unstable); ties may be ordered
  // differently. Irrelevant for N_rt_gas = 1.
  std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return rsq[a] < rsq[b]; });
  for (int l : idx)
    out.push_back({(irt - 1) * Nr + ijk[l][0], (jrt - 1) * Nr + ijk[l][1],
                   (krt - 1) * Nr + ijk[l][2]});
}

}  // namespace

// rt_tools.f90: radiative_physics (with photo_source, photo_cell)
void radiative_physics_plane(Simulation& s, int krt, RPSink& sink) {
  RTGrid& G = s.rt;
  const Units& u = s.units;
  const SimClock& clk = s.clk;
  const AtomicRates& AR = s.rates;
  const int Nrt = s.cfg.Nrt;
  const double g = gamma_ideal;
  const double nu0[3] = {HI_eV, HeI_eV, HeII_eV};
  const double dtsec_rt = clk.dtsec_rt;
  const double Ncell = s.cfg.ncellr;
  const double Nc = s.cfg.N_rt_gas_c;
  std::vector<std::array<int, 3>> index_gas_rt;

  for (int jrt = 1; jrt <= Nrt; ++jrt) {
    for (int irt = 1; irt <= Nrt; ++irt) {
      const std::int64_t c = G.cell(irt, jrt, krt);
      float* gng = G.ngam_p(c);
      float* gnu = G.nu_p(c);
      float* gfg = G.fgam_p(c);
      float* gdx = G.dx_p(c);

      // Effective frequency
      for (int n = 0; n < kNfreq; ++n) {
        if (static_cast<double>(gng[n]) < u.ngam_min) {
          gnu[n] = static_cast<float>(freq_init[n]);
          gng[n] = 0;
        } else {
          gnu[n] = gnu[n] / gng[n];
        }
      }

      const double ngam_tot = static_cast<double>(f::sum4<kNfreq>(gng));

      // Effective unit direction vector
      if (ngam_tot < u.ngam_min) {
        gdx[0] = gdx[1] = gdx[2] = 0;
      } else {
        const float sng = f::sum4<kNfreq>(gng);
        double dx[3];
        for (int d = 0; d < 3; ++d) dx[d] = static_cast<double>(gdx[d] / sng);
        const double dl = std::sqrt(((0.0 + dx[0] * dx[0]) + dx[1] * dx[1]) + dx[2] * dx[2]);
        if (dl > static_cast<double>(1e-9f)) {
          for (int d = 0; d < 3; ++d) gdx[d] = static_cast<float>(dx[d] / dl);
        } else {
          gdx[0] = gdx[1] = gdx[2] = 0;
        }
      }

      if (!(ngam_tot > u.ngam_min || G.fHI[c] < 0.99f)) continue;

      if (ngam_tot > u.ngam_min) sink.ngfill += 1;

      make_index_gas_rt(s, irt, jrt, krt, index_gas_rt);

      // photo_source
      if (static_cast<double>(G.sfr[c]) > u.sfr_min) {
        for (int inu = 0; inu < kNfreq; ++inu) {
          const double nu = freq_init[inu];
          double ngam = fResc * pop2_init[inu] * (static_cast<double>(G.Luv[c]) * dtsec_rt);
          ngam = ngam / u.vol_rt / u.nH_avg;
          if (nu > HI_eV) sink.add(RPSink::kNgam + 0, ngam / s.cfg.nrtcr);
          if (nu > HeI_eV)
            sink.add(RPSink::kNgam + 1, ngam / s.cfg.nrtcr * (u.nH_avg / u.nHe_avg));
          if (nu > HeII_eV)
            sink.add(RPSink::kNgam + 2, ngam / s.cfg.nrtcr * (u.nH_avg / u.nHe_avg));
        }
        sink.nsfill += 1;
      }

      // ---- photo_cell ----
      double us[5] = {0, 0, 0, 0, 0}, ws[4] = {0, 0, 0, 0};
      std::array<double, kNfreq> nu;
      double nrad[kNfreq];
      for (int n = 0; n < kNfreq; ++n) {
        nu[n] = gnu[n];
        nrad[n] = static_cast<double>(gng[n]) * u.nH_avg * u.vol_rt;
      }
      const SigmaTable sigma_pht = photoion_sigma(nu);

      for (const auto& gi : index_gas_rt) {
        const std::int64_t gc = s.gas.cell(gi[0], gi[1], gi[2]);
        float* uu = &s.gas.u[gc * 5];
        float* ww = &s.gas.w[gc * 4];

        // Gas densities
        uu[0] = static_cast<float>(f::max(static_cast<double>(uu[0]), u.Dgas_min));
        const double rhog = static_cast<double>(uu[0]) * u.rho_unit;
        const double nH = rhog * XHy / mH_cgs;
        const double nHe = rhog * YHe / mHe_cgs;
        double nHI = f::min(static_cast<double>(ww[0]) * u.H_unit, nH);
        double nHeI = f::min(static_cast<double>(ww[1]) * u.He_unit, nHe);
        double nHeII = f::min(static_cast<double>(ww[2]) * u.He_unit, nHe - nHeI);
        double nHII = nH - nHI;
        double nHeIII = nHe - nHeI - nHeII;
        double ne = nHII + nHeII + 2 * nHeIII;
        double neff = nH + nHe + ne;
        double nff = nHII + nHeII + 4 * nHeIII;

        // Temperature. kin = sum(u(2:4)**2)/u(1)/2 is pure real(4).
        const float kin_f =
            (((0.0f + uu[1] * uu[1]) + uu[2] * uu[2]) + uu[3] * uu[3]) / uu[0] / 2.0f;
        const double kin = static_cast<double>(kin_f);
        double eth = (static_cast<double>(uu[4]) - kin) * u.eden_unit;
        // hydro_scheme = unsplit (docs/hydro_unsplit.md section 4): the
        // thermal energy is the dual-energy field, which gets the same
        // heating/cooling as u(5). Empty with the legacy scheme: no-op.
        float* ee = s.gas.e.empty() ? nullptr : &s.gas.e[gc];
        if (ee) eth = static_cast<double>(*ee) * u.eden_unit;
        double tmp = (g - 1) * eth / k_cgs / neff;
        if (tmp < u.temp_min) {
          tmp = u.temp_min;
          eth = neff * k_cgs * tmp / (g - 1);
        }

        // Temperature index for RT rates
        int iT = 1 + f::nint((std::log10(tmp) - lgT_min) / lgT_del);
        if (s.cfg.legacy_physics && (iT < 1 || iT > Ntemp))
          throw std::runtime_error("radiative_physics: temperature index out of table, iT=" +
                                   std::to_string(iT));
        // P9 (legacy_physics = 0; harness PATCHES.md 10): rates interpolated
        // in log T (rates_interp). With legacy_physics they are the table
        // values at the nearest node iT, as in the original.
        const RateSet rs = s.cfg.legacy_physics ? rates_node(AR, iT) : rates_interp(AR, tmp);

        double ngam[kNfreq];
        for (int n = 0; n < kNfreq; ++n) ngam[n] = static_cast<double>(gng[n]) * u.nH_avg;

        // Recombination radiation: disabled (`if (.false.)`)
        const double nrab[3] = {0, 0, 0};

        // Maximum photoionization rate in 1/s
        double gamma_pht[3][kNfreq], gamma_max[3][kNfreq];
        for (int inu = 0; inu < kNfreq; ++inu)
          for (int sp = 0; sp < 3; ++sp)
            gamma_pht[sp][inu] = ngam[inu] * sigma_pht[sp][inu] * (u.len_rt / dtsec_rt);
        for (int sp = 0; sp < 3; ++sp)
          for (int inu = 0; inu < kNfreq; ++inu) gamma_max[sp][inu] = gamma_pht[sp][inu];
        double gam_pht[3];
        auto sum_dim2 = [&](double out[3]) {
          for (int sp = 0; sp < 3; ++sp) {
            double t = 0;
            for (int inu = 0; inu < kNfreq; ++inu) t += gamma_pht[sp][inu];
            out[sp] = t;
          }
        };
        sum_dim2(gam_pht);

        // Implicit solution for ionisation balance
        int n = 0;
        double ne_prev = static_cast<double>(f::huge_r4);
        double nH1 = 0, nHe1 = 0, nHe2 = 0;
        while (std::fabs(ne / ne_prev - 1) > static_cast<double>(1e-3f) && n < 5) {
          n = n + 1;
          ne_prev = ne;
          const double AH2 = rs.alpha[3] * ne * dtsec_rt;
          const double AHe2 = rs.alpha[4] * ne * dtsec_rt;
          const double AHe3 = rs.alpha[5] * ne * dtsec_rt;
          const double BH1 = (gam_pht[0] + rs.gcol[0] * ne) * dtsec_rt;
          const double BHe1 = (gam_pht[1] + rs.gcol[1] * ne) * dtsec_rt;
          const double BHe2 = (gam_pht[2] + rs.gcol[2] * ne) * dtsec_rt;
          const double sub1 = 1 + AHe3 + AHe2 + BHe2;
          const double sub2 = AHe2;
          const double sub3 = BHe1 - AHe3;
          const double sub4 = 1 + BHe1;
          const double sub5 = AHe3 * nHe;
          const double sub6 = sub1 * sub4 - sub2 * sub3;
          nH1 = (nHI + AH2 * nH) / (1 + AH2 + BH1);
          nHe1 = (nHeI * sub1 + (nHeII + sub5) * sub2) / sub6;
          nHe2 = (nHeI * sub3 + (nHeII + sub5) * sub4) / sub6;
          ne = (nH - nH1) + nHe2 + 2 * (nHe - nHe1 - nHe2);

          // Photoionisations can not exceed the number of photons
          for (int inu = 0; inu < kNfreq; ++inu) {
            const double num1 = (((0.0 + gamma_pht[0][inu] * nH1) + gamma_pht[1][inu] * nHe1) +
                                 gamma_pht[2][inu] * nHe2) *
                                dtsec_rt * u.vol_gas;
            if (num1 > u.ngam_min) {
              const double x = nrad[inu] / num1;
              for (int sp = 0; sp < 3; ++sp)
                gamma_pht[sp][inu] = f::min(x * gamma_pht[sp][inu], gamma_max[sp][inu]);
            }
          }
          sum_dim2(gam_pht);
        }
        // (n == 0 cannot happen: the first test is ~1 > 1e-3)

        // Update H, He number densities
        nHI = nH1;
        nHII = nH - nHI;
        nHeI = nHe1;
        nHeII = nHe2;
        nHeIII = nHe - nHeI - nHeII;
        ne = nHII + nHeII + 2 * nHeIII;
        neff = nH + nHe + ne;
        nff = nHII + nHeII + 4 * nHeIII;

        // Count photoionisations, photon conservation
        double nion[3] = {0, 0, 0};
        for (int inu = 0; inu < kNfreq; ++inu) {
          const double d = rhog / u.rhob_avg;
          const double x = f::min(std::exp(-(d - rho_source) / rho_source), 1.0);
          const double num3[3] = {gamma_pht[0][inu] * nHI * dtsec_rt,
                                  gamma_pht[1][inu] * nHeI * dtsec_rt,
                                  gamma_pht[2][inu] * nHeII * dtsec_rt};
          for (int sp = 0; sp < 3; ++sp) nion[sp] = nion[sp] + num3[sp];
          const double snum3 = ((0.0 + num3[0]) + num3[1]) + num3[2];
          nrad[inu] = f::max(nrad[inu] - x * snum3 * u.vol_gas, 0.0);
        }

        // Collisional ionisations and recombinations
        const double nsp[3] = {nHI, nHeI, nHeII};
        const double nsp_rec[3] = {nHII, nHeII, nHeIII};
        double ncol[3], nrec[3];
        for (int sp = 0; sp < 3; ++sp) {
          ncol[sp] = rs.gcol[sp] * nsp[sp] * ne * dtsec_rt;
          nrec[sp] = rs.alpha[sp + 3] * nsp_rec[sp] * ne * dtsec_rt;
        }

        // Photoheating
        double lam_pht[3];
        for (int sp = 0; sp < 3; ++sp) {
          double t = 0;
          for (int inu = 0; inu < kNfreq; ++inu)
            t += gamma_pht[sp][inu] * (nu[inu] - nu0[sp]) * eV2erg;
          lam_pht[sp] = t;
        }
        const double heat_rate =
            ((0.0 + lam_pht[0] * nHI) + lam_pht[1] * nHeI) + lam_pht[2] * nHeII;

        // SN heating
        double sn_rate;
        if (tmp > 1e4 && rhog / u.rhob_avg > rho_source) {
          const double d = rhog / u.rhob_avg;
          const double x = f::max(1 - std::exp(-(d - rho_source) / rho_source), 0.0);
          const double tdyn = std::sqrt(3 * pi / (32 * G_cgs * rhog));
          sn_rate = x * eSN * eSF * rhog / tdyn * (c_cgs * c_cgs);
        } else {
          sn_rate = 0;
        }

        // Cooling
        const double Tcmb = u.Tcmb;
        double cool_rate;
        if (s.cfg.legacy_physics) {
          cool_rate = rs.cool[3] * ne * nHII + rs.cool[4] * ne * nHeII +
                      rs.cool[5] * ne * nHeIII + rs.cool[6] * ne * nHI +
                      rs.cool[7] * ne * nHeI + rs.cool[8] * ne * nHeII +
                      rs.cool[9] * ne * nHI + rs.cool[10] * ne * nHeI +
                      rs.cool[11] * ne * nHeII + rs.cool[12] * ne * nff +
                      rs.cool[13] * ne * (tmp - Tcmb) * f::pow4(Tcmb);
        } else {
          // P2/P3 (legacy_physics = 0): HeI collisional excitation (Cen 1992)
          // scales as ne^2*nHeII; slot 15 is HeII dielectronic recombination
          // cooling (Hui & Gnedin 1997), times ne*nHeII. Same term order as
          // the Fortran, left to right.
          cool_rate = rs.cool[3] * ne * nHII + rs.cool[4] * ne * nHeII +
                      rs.cool[5] * ne * nHeIII + rs.cool[6] * ne * nHI +
                      rs.cool[7] * ne * nHeI + rs.cool[8] * ne * nHeII +
                      rs.cool[9] * ne * nHI + rs.cool[10] * ne * ne * nHeII +
                      rs.cool[11] * ne * nHeII + rs.cool[12] * ne * nff +
                      rs.cool[13] * ne * (tmp - Tcmb) * f::pow4(Tcmb) +
                      rs.cool[14] * ne * nHeII + rs.cool[15] * ne * nHeI;
        }

        // Prevent overcooling/overheating
        {
          const double x = rhog / u.rhob_avg;
          if (x > 200) {
            cool_rate = cool_rate * f::pow2(200 / x);
            sn_rate = sn_rate * std::pow(200 / x, 1.5);
          }
        }

        const double net_rate = heat_rate + sn_rate - cool_rate;

        if (net_rate > 0) {
          eth = eth + net_rate * dtsec_rt;
        } else if (net_rate < 0) {
          eth = eth / (1 + std::log(1 - net_rate * dtsec_rt / eth));
        }

        tmp = (g - 1) * eth / k_cgs / neff;
        if (tmp < u.temp_min) {
          tmp = u.temp_min;
          eth = neff * k_cgs * tmp / (g - 1);
        }

        // Update u_gas and w_gas
        uu[4] = static_cast<float>(kin + eth / u.eden_unit);
        if (ee) *ee = static_cast<float>(eth / u.eden_unit);
        ww[0] = static_cast<float>(f::min(f::max(nHI, 0.0), nH) / u.H_unit);
        ww[1] = static_cast<float>(f::min(f::max(nHeI, 0.0), nHe) / u.He_unit);
        ww[2] = static_cast<float>(f::min(f::max(nHeII, 0.0), nHe) / u.He_unit);
        for (int m = 0; m < 5; ++m) us[m] = us[m] + static_cast<double>(uu[m]);
        for (int m = 0; m < 4; ++m) ws[m] = ws[m] + static_cast<double>(ww[m]);

        // Update radiation on RT grid
        float* gga = G.gamma_p(c);
        float* gla = G.lambda_p(c);
        for (int sp = 0; sp < 3; ++sp) {
          gga[sp] = static_cast<float>(static_cast<double>(gga[sp]) + gam_pht[sp] / Nc);
          gla[sp] = static_cast<float>(static_cast<double>(gla[sp]) + lam_pht[sp] / Nc);
        }

        // Update RT statistics
        const double navg[3] = {u.nH_avg, u.nHe_avg, u.nHe_avg};
        for (int sp = 0; sp < 3; ++sp) {
          sink.add(RPSink::kNion + sp, nion[sp] / navg[sp] / Ncell);
          sink.add(RPSink::kNcol + sp, ncol[sp] / navg[sp] / Ncell);
          sink.add(RPSink::kNrec + sp, nrec[sp] / navg[sp] / Ncell);
          sink.add(RPSink::kNrab + sp, nrab[sp] / navg[sp] / Ncell);
          sink.add(RPSink::kGamma + sp, gam_pht[sp] / Ncell);
          sink.add(RPSink::kLambda + sp, lam_pht[sp] / Ncell);
        }
        // Thomson optical depth: physical c (and physical dt), also with
        // c_reduced < 1
        sink.add(RPSink::kDtau, sT_cgs * ne * (c_cgs * dtsec_rt) / Ncell);
      }

      // Remaining photons -> fgam (numerator)
      for (int n = 0; n < kNfreq; ++n)
        gfg[n] = static_cast<float>(static_cast<double>(gfg[n]) + nrad[n] / u.vol_rt / u.nH_avg);

      // Update gas on RT grid (no Dgas_min clamp here, as in the original)
      for (int m = 0; m < 5; ++m) us[m] = us[m] / Nc;
      for (int m = 0; m < 4; ++m) ws[m] = ws[m] / Nc;
      const GasState gs = gas_state_from_avg(us, ws, u, false);
      store_gas_state(G, c, gs, u);
      // ---- end photo_cell ----

      // Update radiation on RT grid
      for (int m = 0; m < kNfreq; ++m) {
        if (static_cast<double>(gng[m]) < u.ngam_min) {
          gfg[m] = 0;
          gng[m] = 0;
        } else {
          gfg[m] = gfg[m] / gng[m];
          gng[m] = gng[m] * gfg[m];
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// rt_step_3
// ---------------------------------------------------------------------------

// rt_tools.f90: rt_step_3
void rt_step_3(Simulation& s) {
  RTGrid& G = s.rt;
  const double tw0 = wall_seconds();
  if (s.cfg.det_parallel) {
    det_sweep(s, [&](const std::array<int, 6>& box) { ray_tracing_3(s, box); });
    det_fold(s.rays);
  } else {
    for (int n = 1; n <= s.dom.Ndomain; ++n) {
      const int i = set_domain_rt(s.dom);
      ray_tracing_3(s, s.dom.box[i]);
      end_domain_rt(s.dom, i);
    }
  }
  const double tw1 = wall_seconds();
  s.timers.trace3 += tw1 - tw0;
  struct Acc {
    double& t;
    double t0;
    ~Acc() { t += wall_seconds() - t0; }
  } acc{s.timers.sync, tw1};
#pragma omp parallel for schedule(static) if (par::worth(s.cfg.Nrtc))
  for (std::int64_t c = 0; c < G.ncell; ++c) {
    G.htoc[c * 4 + 0] = G.htoc[c * 4 + 2];
    G.htoc[c * 4 + 1] = G.htoc[c * 4 + 3];
  }
  // (calc_gas_3: redshift of reionisation, zreion only -- not ported)

  RTStat& sv = s.rt_save;
  const RTStat& st = s.rt_stat;
  sv.nsfill = 0;
  sv.ngfill = 0;
  sv.gamma = {};
  sv.lambda = {};
  RayStat& rv = s.rays.save;
  const RayStat& rs = s.rays.stat;
  rv.nfree = 0;
  rv.ncst1 = rv.ncst2 = rv.nmrg1 = rv.nmrg2 = rv.nspl1 = rv.nspl2 = rv.ndel1 = 0;

  // icpu = 1
  sv.nsfill = sv.nsfill + st.nsfill;
  sv.ngfill = sv.ngfill + st.ngfill;
  sv.tau = sv.tau + st.dtau;
  sv.dtau = sv.dtau + st.dtau;
  for (int m = 0; m < 3; ++m) {
    sv.ngam[m] = sv.ngam[m] + st.ngam[m];
    sv.nion[m] = sv.nion[m] + st.nion[m];
    sv.ncol[m] = sv.ncol[m] + st.ncol[m];
    sv.nrec[m] = sv.nrec[m] + st.nrec[m];
    sv.nrab[m] = sv.nrab[m] + st.nrab[m];
    sv.gamma[m] = sv.gamma[m] + st.gamma[m];
    sv.lambda[m] = sv.lambda[m] + st.lambda[m];
  }
  rv.nfree = rv.nfree + rs.nfree;
  rv.ncst1 = rv.ncst1 + rs.ncst1;
  rv.ncst2 = f::max(rv.ncst2, rs.ncst2);
  rv.nmrg1 = rv.nmrg1 + rs.nmrg1;
  rv.nmrg2 = f::max(rv.nmrg2, rs.nmrg2);
  rv.nspl1 = rv.nspl1 + rs.nspl1;
  rv.nspl2 = f::max(rv.nspl2, rs.nspl2);
  rv.ndel1 = rv.ndel1 + rs.ndel1;
}

// ---------------------------------------------------------------------------
// summarize_rt
// ---------------------------------------------------------------------------

void summarize_rt(Simulation& s) {
  std::ostream* O = s.out;
  if (!O) return;
  const RTStat& v = s.rt_save;
  auto arr = [&](const char* name, const std::array<double, 3>& a) {
    *O << " " << name << static_cast<float>(a[0]) << " " << static_cast<float>(a[1]) << " "
       << static_cast<float>(a[2]) << "\n";
  };
  *O << " Source  : " << s.Nsource << " " << v.nsfill << " "
     << static_cast<float>(v.nsfill / s.cfg.nrtcr) << "\n";
  *O << " RT      : " << s.cfg.Nrtc << " " << v.ngfill << " "
     << static_cast<float>(v.ngfill / s.cfg.nrtcr) << "\n";
  *O << " Ray     : " << s.rays.Nray_max << " " << s.rays.Nray << " "
     << static_cast<float>(s.rays.Nray / static_cast<double>(s.rays.Nray_max)) << "\n";
  *O << " Luv     : " << static_cast<float>(v.Luv) << "\n";
  *O << " SFR     : " << static_cast<float>(v.sfr) << "\n";
  arr("ngam    : ", v.ngam);
  arr("nion    : ", v.nion);
  arr("ncol    : ", v.ncol);
  arr("nrec    : ", v.nrec);
  arr("nrab    : ", v.nrab);
  arr("Gamma   : ", v.gamma);
  arr("Lambda  : ", v.lambda);
  *O << " tau     : " << static_cast<float>(v.tau) << " " << static_cast<float>(v.dtau) << "\n";
}

std::string format_rt_timers(const RTTimers& t) {
  char buf[400];
  std::snprintf(buf, sizeof buf,
                "RT timers [s]: substeps %lld prep %.3f calc_rt_1 %.3f calc_ray_1 %.3f "
                "ray_tracing_1 %.3f radiative_physics %.3f ray_tracing_3 %.3f sync %.3f "
                "post_rt %.3f",
                static_cast<long long>(t.substeps), t.prep, t.calc_rt1, t.calc_ray1, t.trace1,
                t.physics, t.trace3, t.sync, t.post);
  return buf;
}

std::int64_t count_live_rays(const Simulation& s) {
  std::int64_t n = 0;
  for (std::int64_t ir = 1; ir <= s.rays.hwm(); ++ir)
    if (s.rays.at(ir).nside > 0) ++n;
  return n;
}

double ngam_rays_total(const Simulation& s) {
  // Live rays = those linked in the cell chains. Summed in storage order
  // over slots with nside > 0; valid after post_rt (dead slots zeroed).
  double t = 0;
  if (s.cfg.det_parallel) {
    // deterministic_parallel: storage slots depend on thread timing; sum in
    // the canonical order instead (cells k,j,i, each chain head to tail)
    for (std::int64_t c = 0; c < s.rt.ncell; ++c)
      for (std::int64_t ir = s.rt.htoc_ray(3, c); ir > 0; ir = s.rays.at(ir).next) {
        const Ray& r = s.rays.at(ir);
        double sr = 0;
        for (int n = 0; n < kNfreq; ++n) sr += static_cast<double>(r.ngam[n]);
        t = t + sr;
      }
    return t;
  }
  for (std::int64_t ir = 1; ir <= s.rays.hwm(); ++ir) {
    const Ray& r = s.rays.at(ir);
    if (r.nside > 0) {
      // ngam_tot + sum(real(ray(ir)%ngam,8)) of the harnesses: the ray's sum
      // first, then added (not five separate additions)
      double sr = 0;
      for (int n = 0; n < kNfreq; ++n) sr += static_cast<double>(r.ngam[n]);
      t = t + sr;
    }
  }
  return t;
}

}  // namespace radray
