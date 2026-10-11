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
#include "radray/kernels/ray_kernels.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/ray_dump.hpp"
#include "radray/ray_quant.hpp"
#include "radray/parallel.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/rt_backend.hpp"
#include "radray/rt_moments.hpp"
#include "radray/recomb_rad.hpp"

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

void check_rt_cells_param(const Config& c) {
  const int k = c.rt_cells_per_substep;
  if (k < 1 || k > rrk::kMaxCellsPerSubstep)
    throw std::runtime_error("rt_cells_per_substep must be 1.." +
                             std::to_string(rrk::kMaxCellsPerSubstep));
  if (k == 1) return;
  if (c.legacy_rt) throw std::runtime_error("rt_cells_per_substep > 1 requires legacy_rt = 0");
  if (use_moments(c))
    throw std::runtime_error("rt_cells_per_substep > 1 is implemented for rt_method = rays only");
}

// deterministic_parallel domains (simulation.hpp: DetDomains)
void init_det_domains(const Config& cfg, DetDomains& dd) {
  const int N = cfg.Nrt;
  dd = DetDomains{};
  // The reach of one cell's ray work is k+1 cells for k = rt_cells_per_substep
  // (a split or merge moves a ray by at most one cell, then a cast of k
  // cells; a ray moves k cells, plus one at a cell face, in ray_tracing_3),
  // and it must stay within d/2: d >= 2k+2. d = 8 (the original choice) for
  // k <= 3.
  const int dmin = std::max(8, 2 * cfg.rt_cells_per_substep + 2);
  for (int d = dmin; d <= N / 2; ++d)
    if (N % d == 0 && (N / d) % 2 == 0) {
      dd.d = d;
      break;
    }
  if (dd.d == 0)
    throw std::runtime_error("deterministic_parallel needs Nrt = 2*m*d with d >= " +
                             std::to_string(dmin) + " (Nrt >= " + std::to_string(2 * dmin) +
                             "; d >= 2k+2 for rt_cells_per_substep = k)");
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

    const bool moments = use_moments(s.cfg);
    if (use_amrex(s.cfg)) {
      // rt_backend = amrex (rt_backend.hpp): the whole substep on the level
      amrex_rt_substep(s);
    } else {
      if (s.raydump_state && !moments) {
        s.rays.prov_substep = ++s.raydump_state->substep;
      }
      rt_step_1(s);
      rt_step_2(s);
      rt_step_3(s);
      if (!moments) {
        if (s.ray_quant.enabled()) ray_quantise(s, s.ray_quant);  // research only
        if (s.raydump_state) ray_dump_after_substep(s);
      }
    }

    const std::int64_t nr = s.rays.Nray - s.rays.save.nfree;
    const float fnr = static_cast<float>(static_cast<double>(nr) /
                                         static_cast<double>(s.rays.Nray_max));
    if (L && moments) {
      // rt_method = m1 | otvet (rt_moments.hpp): no rays
      *L << " Sources    : " << s.rt_save.nsfill << " "
         << static_cast<float>(s.rt_save.nsfill / s.cfg.nrtcr) << "\n";
      *L << " Photons    : " << s.rt_save.ngfill << " "
         << static_cast<float>(s.rt_save.ngfill / s.cfg.nrtcr) << "\n";
      if (clk.istep_rt == clk.nstep_rt) *L << "\n";
    } else if (L) {
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

// The rt_grid fields of cell c as the kernels' pointer view
rrk::CellRef cell_ref(RTGrid& G, std::int64_t c) {
  rrk::CellRef r;
  r.rhog = &G.rhog[c];
  r.fHI = &G.fHI[c];
  r.fHeI = &G.fHeI[c];
  r.fHeII = &G.fHeII[c];
  r.temp = &G.temp[c];
  r.Luv = &G.Luv[c];
  r.sfr = &G.sfr[c];
  r.dx = G.dx_p(c);
  r.gamma = G.gamma_p(c);
  r.lambda = G.lambda_p(c);
  r.nu = G.nu_p(c);
  r.ngam = G.ngam_p(c);
  r.fgam = G.fgam_p(c);
  r.tau = G.tau_p(c);
  return r;
}

// The gas state of an RT cell from the averaged gas (calc_rt_i, calc_rt_f,
// the end of photo_cell) is the kernels' rrk::gas_state_from_avg /
// rrk::store_gas_state (include/radray/kernels/rt_cell.hpp).
using GasState = rrk::GasStateP;

// Sum u_gas/w_gas over the Eulerian cells of RT cell (irt,jrt,krt) and
// normalise (us = us/N_rt_gas_c, us(1) = max(us(1),Dgas_min)). Returns the
// averaged internal-energy density of hydro_scheme = unsplit, or -1.
double average_gas(const Simulation& s, int irt, int jrt, int krt, double us[5],
                   double ws[4]) {
  const int Nr = s.cfg.N_rt_gas;
  const bool has_e = !s.gas.e.empty();
  double es = 0;
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
        if (has_e) es = es + static_cast<double>(s.gas.e[gc]);
      }
  for (int m = 0; m < 5; ++m) us[m] = us[m] / s.cfg.N_rt_gas_c;
  for (int m = 0; m < 4; ++m) ws[m] = ws[m] / s.cfg.N_rt_gas_c;
  us[0] = f::max(us[0], s.units.Dgas_min);
  return has_e ? es / s.cfg.N_rt_gas_c : -1.0;
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
        const double es = average_gas(s, irt, jrt, krt, us, ws);
        const GasState g = rrk::gas_state_from_avg(us, ws, s.units, es);
        rrk::store_gas_state(cell_ref(G, G.cell(irt, jrt, krt)), g, s.units);
      }

  // recomb_rad = 1: buffers of the recombination radiation (recomb_rad.hpp)
  if (s.cfg.recomb_rad && !s.recomb) recomb_init(s);
  if (s.cfg.recomb_rad) check_recomb_param(s.cfg);

  // rt_backend = amrex: the level holds the cells and the rays
  // (rt_backend.hpp); the gas state above is recomputed there
  if (use_amrex(s.cfg)) {
    amrex_rt_prep(s);
    return;
  }

  // rt_method = m1 | otvet: no ray pool (rt_moments.hpp)
  if (use_moments(s.cfg)) {
    moment_prep(s);
    return;
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
  const bool moments = use_moments(s.cfg);
  if (use_amrex(s.cfg)) {
    amrex_rt_post(s);  // redshift, calc_rt_f, fields back to s (rt_backend.hpp)
    return;
  }
  if (moments) {
    moment_redshift(s);  // rt_method = m1 | otvet (rt_moments.hpp)
  } else {
    // deterministic_parallel: there are no Fortran buffers to return, and
    // freed slots are already zero.
    if (!s.cfg.det_parallel) ray_buffer_deallocate(s);
    ray_redshifting(s);
  }

  RTGrid& G = s.rt;
  RayState& R = s.rays;
  const int Nrt = s.cfg.Nrt;
  double sig[3][kNfreq];
  rrk::photoion_sigma(freq_init.data(), sig);

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
        const double es = average_gas(s, irt, jrt, krt, us, ws);
        const GasState g = rrk::gas_state_from_avg(us, ws, s.units, es);
        rrk::store_gas_state(cell_ref(G, c), g, s.units);
        rrk::tau_from_state(G.tau_p(c), g, sig, s.units.len_rt);
      }

  if (moments) return;  // (the ray PDFs below)

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
          nr = f::max<std::int64_t>(
              nr, rrk::nray_from_nside(rrk::emit_nside(s.cfg.rt_cells_per_substep)));
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
// ngam_ray thresholds and the pool guard (rt_step_1 / calc_ray_1), shared
// with the AMReX backend (src/amrex/rt_level.cpp)
// ---------------------------------------------------------------------------

void calc_ngam_ray(const std::int64_t* ngam_pdf, int N_pdf, double ngam_avg, std::int64_t Nray,
                   std::int64_t Nray_max, bool legacy, std::array<double, 64>& ngam_ray,
                   std::ostream* log) {
  const double lgn_min = rrk::kLgnMin;
  const double lgn_del = folded::lgn_del;
  // With legacy_rt, thresholds that are not reached keep their stale values
  // (B2); otherwise they are recomputed from scratch.
  if (ngam_avg == 0) {
    ngam_ray.fill(0.0);
    return;
  }
  ngam_avg = ngam_avg / Nray;
  if (!legacy) ngam_ray.fill(0.0);  // B2 fix
  // legacy: int(Nray_max/64.) is a real(4) division (int64 product here;
  // the Fortran's int32 product overflows only for N >= ~323).
  // B9 fix: exact 64-bit integer Nray_max/64.
  const std::int64_t per =
      legacy ? f::int4(static_cast<double>(static_cast<float>(Nray_max) / 64.0f)) : Nray_max / 64;
  for (int j = 1; j <= 64; ++j) {
    std::int64_t nr = 0;
    for (int i = N_pdf; i >= 1; --i) {
      nr = nr + ngam_pdf[i];
      if (nr > j * per) {
        ngam_ray[j - 1] = std::pow(1e1, (i - 1) * lgn_del + lgn_min);
        break;
      }
    }
    if (log && j % 4 == 0)
      *log << " Rays       : " << j << " " << static_cast<float>(ngam_ray[j - 1]) << " "
           << static_cast<float>(ngam_ray[j - 1] / ngam_avg) << "\n";
  }
}

std::int64_t pool_split_candidates(const std::int64_t* spl_pdf, int N_pdf, double ngam_ray32) {
  std::int64_t nspl = 0;
  for (int n = 1; n <= N_pdf; ++n)
    if (std::pow(1e1, n * folded::lgn_del + rrk::kLgnMin) > ngam_ray32) nspl += spl_pdf[n];
  return nspl;
}

// F3 (legacy_rt = 0 only; RT harness PATCHES.md 13): keep the ray pool from
// overflowing. Rays are only created by splits (3 per split) and sources
// (nray_init per source cell, nray_src with rt_cells_per_substep > 4:
// rrk::emit_nside). From the live count L, the number of source
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
void pool_guard(const std::int64_t* spl_pdf, int N_pdf, std::int64_t L, std::int64_t nsrc,
                std::int64_t nrec_rays, std::int64_t Nray_max, std::array<double, 64>& ngam_ray,
                std::int64_t& npool_limited, std::int64_t& npool_forced, std::ostream* log,
                int nray_src) {
  const double lgn_min = rrk::kLgnMin;
  const double lgn_del = folded::lgn_del;
  const std::int64_t budget = (9 * Nray_max) / 10 - L - nray_src * nsrc - nrec_rays;
  if (budget <= 0) {
    ngam_ray[16 - 1] = 1e30;
    ngam_ray[32 - 1] = 1e30;
    npool_forced += 1;
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
      if (thr > ngam_ray[32 - 1]) {
        ngam_ray[32 - 1] = thr;
        npool_limited += 1;
      }
    }
  }
  if (log && (npool_forced > 0 || npool_limited > 0))
    *log << " Pool guard : " << L << " " << nsrc << " " << budget << " " << npool_limited << " "
         << npool_forced << "\n";
}

// ---------------------------------------------------------------------------
// rt_step_1
// ---------------------------------------------------------------------------

// rt_tools.f90: rt_step_1
void rt_step_1(Simulation& s) {
  if (use_moments(s.cfg)) {
    moment_step_1(s);  // rt_method = m1 | otvet (rt_moments.hpp)
    return;
  }
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
  double sig[3][kNfreq];
  rrk::photoion_sigma(freq_init.data(), sig);
#pragma omp parallel for schedule(static) if (par::worth(s.cfg.Nrtc))
  for (int krt = 1; krt <= Nrt; ++krt)
    for (int jrt = 1; jrt <= Nrt; ++jrt)
      for (int irt = 1; irt <= Nrt; ++irt)
        rrk::calc_rt_1_cell(cell_ref(G, G.cell(irt, jrt, krt)), u, sig, legacy);

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

  // Calculate ngam_ray (calc_ngam_ray below)
  calc_ngam_ray(ngam_pdf.data(), N_pdf, ngam_avg, R.Nray, R.Nray_max, legacy, R.ngam_ray,
                s.raylog);

  // F3 (legacy_rt = 0 only): the ray pool guard (pool_guard below)
  if (!legacy) {
    std::int64_t L = 0;
    for (int n = 1; n <= N_pdf; ++n) L += ngam_pdf[n];
    std::int64_t nsrc = 0;
    for (std::int64_t c = 0; c < G.ncell; ++c)
      if (static_cast<double>(G.sfr[c]) > u.sfr_min) ++nsrc;
    // recomb_rad = 1 (recomb_rad.hpp): choose the blocks that launch
    // recombination rays this substep and absorb the rest on the spot; their
    // rays count against the budget like the stellar source rays
    std::int64_t nrec_rays = 0;
    if (s.recomb) {
      // stellar splits expected this substep: candidates above the split
      // threshold ngam_ray(32), 3 new rays each; the recombination rays only
      // get the room left after them
      const std::int64_t nspl = pool_split_candidates(spl_pdf.data(), N_pdf, R.ngam_ray[32 - 1]);
      nrec_rays = recomb_prepare(s, L, nsrc, 3 * nspl);
    }
    pool_guard(spl_pdf.data(), N_pdf, L, nsrc, nrec_rays, R.Nray_max, R.ngam_ray,
               R.npool_limited, R.npool_forced, s.raylog,
               rrk::nray_from_nside(rrk::emit_nside(s.cfg.rt_cells_per_substep)));
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
  const int Nrt = s.cfg.Nrt;
  const double dtsec_rt = clk.dtsec_rt;
  const double Ncell = s.cfg.ncellr;
  const double Nc = s.cfg.N_rt_gas_c;
  std::vector<std::array<int, 3>> index_gas_rt;
  // The per-gas-cell physics is the kernels' photo_gas_cell
  // (include/radray/kernels/rt_cell.hpp)
  rrk::PhysCtx X;
  X.u = u;
  X.rates.cool = s.rates.coef_cool_.data();
  X.rates.alpha = s.rates.alpha_rec_.data();
  X.rates.gcol = s.rates.gamma_col_.data();
  X.rates.legacy_physics = s.cfg.legacy_physics;
  X.dtsec_rt = dtsec_rt;
  X.cells = s.cfg.rt_cells_per_substep;
  X.photon_conserve = s.cfg.photon_conserve;
  X.want_nrab = static_cast<bool>(s.recomb);

  for (int jrt = 1; jrt <= Nrt; ++jrt) {
    for (int irt = 1; irt <= Nrt; ++irt) {
      const std::int64_t c = G.cell(irt, jrt, krt);
      const rrk::CellRef cr = cell_ref(G, c);

      // Effective frequency and unit direction vector; skip neutral cells
      // without photons
      double ngam_tot = 0;
      if (!rrk::cell_radiation_prep(cr, u, ngam_tot)) continue;

      if (ngam_tot > u.ngam_min) sink.ngfill += 1;

      make_index_gas_rt(s, irt, jrt, krt, index_gas_rt);

      // photo_source
      if (static_cast<double>(G.sfr[c]) > u.sfr_min) {
        double ngam_sp[3][kNfreq];
        rrk::source_stats(G.Luv[c], dtsec_rt, u, ngam_sp);
        for (int inu = 0; inu < kNfreq; ++inu)
          for (int sp = 0; sp < 3; ++sp)
            if (ngam_sp[sp][inu] != 0) sink.add(RPSink::kNgam + sp, ngam_sp[sp][inu] / s.cfg.nrtcr);
        sink.nsfill += 1;
      }

      // ---- photo_cell ----
      double us[5] = {0, 0, 0, 0, 0}, ws[4] = {0, 0, 0, 0}, es = 0;
      double nu[kNfreq], nrad[kNfreq];
      for (int n = 0; n < kNfreq; ++n) {
        nu[n] = cr.nu[n];
        nrad[n] = static_cast<double>(cr.ngam[n]) * u.nH_avg * u.vol_rt;
      }
      double sigma_pht[3][kNfreq];
      rrk::photoion_sigma(nu, sigma_pht);

      // recomb_case = A: case A rates and cooling (alpha_rec(1:3),
      // coef_cool(1:3)); B (default, the original): alpha_rec(4:6),
      // coef_cool(4:6). With recomb_rad = 1, cells whose block is on the
      // spot this substep use case B.
      const bool onspot = s.recomb && recomb_onspot(*s.recomb, c);
      X.ia = (s.cfg.recomb_case_A && !onspot) ? 0 : 3;

      for (const auto& gi : index_gas_rt) {
        const std::int64_t gc = s.gas.cell(gi[0], gi[1], gi[2]);
        float* uu = &s.gas.u[gc * 5];
        float* ww = &s.gas.w[gc * 4];
        // hydro_scheme = unsplit (docs/hydro_unsplit.md section 4): the
        // thermal energy is the dual-energy field, which gets the same
        // heating/cooling as u(5). Empty with the legacy scheme: no-op.
        float* ee = s.gas.e.empty() ? nullptr : &s.gas.e[gc];

        rrk::PhotoStats st;
        const int err = rrk::photo_gas_cell(uu, ww, ee, cr.ngam, nu, sigma_pht, nrad, X, st);
        if (err == rrk::kPhotoBadTemperatureIndex)
          throw std::runtime_error("radiative_physics: temperature index out of table");
        if (s.recomb)
          recomb_count(s, c, st.nrab, s.cfg.recomb_rad_lines ? st.nrB : nullptr, st.T_in, u.vol_gas);

        if (ee) es = es + static_cast<double>(*ee);
        for (int m = 0; m < 5; ++m) us[m] = us[m] + static_cast<double>(uu[m]);
        for (int m = 0; m < 4; ++m) ws[m] = ws[m] + static_cast<double>(ww[m]);

        // Update radiation on RT grid
        for (int sp = 0; sp < 3; ++sp) {
          cr.gamma[sp] = static_cast<float>(static_cast<double>(cr.gamma[sp]) + st.gam[sp] / Nc);
          cr.lambda[sp] = static_cast<float>(static_cast<double>(cr.lambda[sp]) + st.lam[sp] / Nc);
        }

        // Update RT statistics
        const double navg[3] = {u.nH_avg, u.nHe_avg, u.nHe_avg};
        for (int sp = 0; sp < 3; ++sp) {
          sink.add(RPSink::kNion + sp, st.nion[sp] / navg[sp] / Ncell);
          sink.add(RPSink::kNcol + sp, st.ncol[sp] / navg[sp] / Ncell);
          sink.add(RPSink::kNrec + sp, st.nrec[sp] / navg[sp] / Ncell);
          sink.add(RPSink::kNrab + sp, st.nrab[sp] / navg[sp] / Ncell);
          sink.add(RPSink::kGamma + sp, st.gam[sp] / Ncell);
          sink.add(RPSink::kLambda + sp, st.lam[sp] / Ncell);
        }
        // Thomson optical depth: physical c (and physical dt), also with
        // c_reduced < 1
        sink.add(RPSink::kDtau, st.dtau / Ncell);
      }

      // Remaining photons -> fgam, the gas on the RT grid, the transmitted
      // fraction (photo_cell_finish)
      rrk::photo_cell_finish(cr, nrad, us, ws, s.gas.e.empty() ? -1.0 : es, Nc, u);
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
  if (use_moments(s.cfg)) {
    moment_step_3(s);  // rt_method = m1 | otvet (rt_moments.hpp)
  } else if (s.cfg.det_parallel) {
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
  if (!use_moments(s.cfg)) s.timers.trace3 += tw1 - tw0;
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
  accumulate_rt_stats(s);
}

// rt_step_3: the per-hydro-step sums of the substep's rt_stat(1) and
// ray_stat(1) (also called by the AMReX backend after its substep)
void accumulate_rt_stats(Simulation& s) {
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
  // rt_method = m1 | otvet: the photons of the moment solver
  if (use_moments(s.cfg)) return moment_photons_total(s);
  // rt_backend = amrex: the fixed-point sum over the level's rays
  if (use_amrex(s.cfg)) return s.ngam_total_backend;
  return ngam_rays_total_serial(s);
}

double ngam_rays_total_serial(const Simulation& s) {
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
