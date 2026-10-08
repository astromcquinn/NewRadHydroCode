// Host-side driver logic shared by the spike's serial reference, AMReX and
// Parthenon prototypes: problem input, the clock of rt_driver /
// radiative_transfer, ngam_ray thresholds from the global histogram, the
// redshift table, and output in the format of docs/validation_spec.md (so
// tools/compare/compare.py works unchanged). Uses the radray library for
// the clock, units and rate tables.
#pragma once

#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "radray/config.hpp"
#include "radray/cosmology.hpp"
#include "radray/grid.hpp"
#include "radray/io.hpp"
#include "radray/rt.hpp"
#include "radray/rt_rates.hpp"
#include "ray_kernels.hpp"
#include "rt_cell.hpp"

namespace rrhost {

struct SourceCell {
  int i, j, k;  // 0-based
  float Luv, sfr;
};

struct StepStats {
  long long nray_active = 0, ncst1 = 0, nmrg1 = 0, nspl1 = 0, ndel1 = 0;
  long long ngam_tot_fx = 0;  // fixed point, kScaleTot
};

class RTRun {
 public:
  radray::ProblemParams P;
  radray::Config cfg;
  radray::SimClock clk;
  radray::Units units;
  radray::AtomicRates rates;
  std::vector<SourceCell> sources;
  double sig_init[3][5];
  std::array<double, 64> ngam_ray{};
  int N_pdf = 0;

  void init(const std::string& pdir, const std::vector<std::string>& sets) {
    P = radray::read_params(pdir + "/params.txt");
    for (const std::string& kv : sets) {
      const auto eq = kv.find('=');
      if (eq == std::string::npos) throw std::runtime_error("--set needs key=value: " + kv);
      radray::set_param(P, kv.substr(0, eq), kv.substr(eq + 1));
    }
    if (P.legacy_rt != 0)
      throw std::runtime_error("the spike implements legacy_rt = 0 only (sorted ray arrays)");
    cfg = radray::Config::make(P.N, P.nray_max_factor);
    cfg.legacy_rt = false;
    cfg.legacy_physics = (P.legacy_physics == 1);
    radray::init_variables(P.z_init, clk);
    radray::calc_atomic_rates(rates, cfg.legacy_physics);
    radray::calc_sim_var(cfg, clk, units);
    const std::vector<radray::Source> src = radray::read_sources(pdir + "/sources.txt");
    // add_halo_to_rt: per cell, sum Luv/sfr in input order in real(8)
    std::map<long long, std::pair<double, double>> acc;
    std::map<long long, int> order;
    for (const auto& s : src) {
      const long long c = (s.i - 1) + static_cast<long long>(P.N) * ((s.j - 1) + static_cast<long long>(P.N) * (s.k - 1));
      auto& a = acc[c];
      a.first = a.first + static_cast<double>(s.Luv);
      a.second = a.second + static_cast<double>(s.sfr);
    }
    for (const auto& kv : acc) {
      const long long c = kv.first;
      const int i = static_cast<int>(c % P.N), j = static_cast<int>((c / P.N) % P.N),
                k = static_cast<int>(c / (static_cast<long long>(P.N) * P.N));
      sources.push_back({i, j, k, static_cast<float>(kv.second.first), static_cast<float>(kv.second.second)});
    }
    const std::array<double, 5> fi = radray::cp::freq_init;
    double nu[5];
    for (int n = 0; n < 5; ++n) nu[n] = fi[n];
    rrk::photoion_sigma(nu, sig_init);
    N_pdf = rrk::pdf_nbins();
    if (N_pdf + 2 > rrk::kNpdfMax) throw std::runtime_error("kNpdfMax too small");
  }

  // Gas input (9,N,N,N) as GasGrid (the caller copies its part)
  void read_gas(const std::string& pdir, radray::GasGrid& g) const {
    g.allocate(P.N);
    radray::read_gas(pdir + "/gas.bin", g);
  }

  // rt_driver: start of a hydro step, then the rewind of radiative_transfer.
  // Returns nstep_rt.
  int begin_step() {
    radray::calc_sim_var(cfg, clk, units);
    radray::calc_timestep(cfg, P.da_sim, clk);
    radray::advance_clock_hydro_step(cfg, clk, units);
    clk.year = clk.year - clk.dt_yr;
    clk.scalefactor = clk.scalefactor - clk.da_sim;
    clk.redshift = 1 / clk.scalefactor - 1;
    radray::calc_sim_var(cfg, clk, units);
    return clk.nstep_rt;
  }
  // radiative_transfer: clock update after each substep
  void end_substep() {
    clk.year = clk.year + clk.dt_yr / clk.nstep_rt;
    clk.scalefactor = clk.scalefactor + clk.da_sim / clk.nstep_rt;
    clk.redshift = 1 / clk.scalefactor - 1;
    radray::calc_sim_var(cfg, clk, units);
  }

  rrk::PhysCtx phys() const {
    rrk::PhysCtx X;
    X.u = units;
    X.rates.cool = rates.coef_cool_.data();
    X.rates.alpha = rates.alpha_rec_.data();
    X.rates.gcol = rates.gamma_col_.data();
    X.rates.legacy_physics = cfg.legacy_physics;
    X.dtsec_rt = clk.dtsec_rt;
    return X;
  }
  rrk::RTP rtp() const {
    rrk::RTP p;
    p.N = cfg.Nrt;
    p.nrtr = cfg.nrtr;
    p.hrtr = cfg.hrtr;
    p.seg_min = units.seg_min;
    p.ngam_min = units.ngam_min;
    p.fgam_min = units.fgam_min;
    p.sfr_min = units.sfr_min;
    p.vol_rt = units.vol_rt;
    p.nH_avg = units.nH_avg;
    p.dtsec_rt = clk.dtsec_rt;
    p.nseg = 4;
    p.ngam_ray16 = ngam_ray[15];
    p.ngam_ray32 = ngam_ray[31];
    return p;
  }

  // rt_step_1/calc_ray_1 with legacy_rt = 0. hist[1..N_pdf]: global counts
  // of live rays per bin; npos: number of live rays with sum(ngam) > 0.
  void thresholds(const long long* hist, long long npos) {
    ngam_ray.fill(0.0);
    if (npos == 0) return;
    const long long per = cfg.Nray_max / 64;
    const double lgn_del = radray::folded::lgn_del;
    for (int j = 1; j <= 64; ++j) {
      long long nr = 0;
      for (int i = N_pdf; i >= 1; --i) {
        nr = nr + hist[i];
        if (nr > j * per) {
          ngam_ray[j - 1] = std::pow(1e1, (i - 1) * lgn_del + -6.0);
          break;
        }
      }
    }
  }

  // ray_redshifting table (after the last substep of a hydro step)
  rrk::RedshiftTable redshift_table() const {
    using namespace radray;
    rrk::RedshiftTable T;
    T.bmin = -6;
    T.bmax = 6;
    T.bdel = static_cast<double>(0.1f);
    const int Nbin = 1 + static_cast<int>(std::round((T.bmax - T.bmin) / T.bdel));
    if (Nbin != rrk::kZbin) throw std::runtime_error("redshift table size");
    T.fscale = (clk.scalefactor - clk.da_sim) / clk.scalefactor;
    T.ngam_min = units.ngam_min;
    double* freq = T.freq;
    freq[0] = pc::HI_eV;
    freq[1] = std::sqrt(pc::HI_eV * pc::HeI_eV);
    freq[2] = pc::HeI_eV;
    freq[3] = std::sqrt(pc::HeI_eV * pc::HeII_eV);
    freq[4] = pc::HeII_eV;
    freq[5] = static_cast<double>(1e2f);
    double rmax = 0., rmin = static_cast<double>(3.40282346638528859812e+38F);
    for (int j = 1; j <= 4; ++j) {
      const double fj = freq[j - 1], fj1 = freq[j], fj2 = freq[j + 1];
      for (int i = 1; i <= Nbin; ++i) {
        const double b = (i - 1) * T.bdel + T.bmin;
        const double c = b + 1;
        double v;
        if (std::fabs(c) < T.bdel / 1e2)
          v = std::log(fj2 / fj1) / std::log(fj1 / fj);
        else
          v = (std::pow(fj2, c) - std::pow(fj1, c)) / (std::pow(fj1, c) - std::pow(fj, c));
        T.tab[j - 1][2 * (i - 1) + 1] = b;
        T.tab[j - 1][2 * (i - 1) + 0] = v;
        rmax = rrf::max(rmax, v);
        rmin = rrf::min(rmin, v);
      }
    }
    T.rmin = rmin;
    T.rmax = rmax;
    return T;
  }

  // Output of one hydro step: cells in global order (i fastest)
  void write_step(const std::string& odir, const std::vector<rrk::CellData>& c,
                  const StepStats& st) const {
    char name[32];
    std::snprintf(name, sizeof name, "step_%04d", clk.istep);
    const std::string dir = odir + "/" + name;
    std::filesystem::create_directories(dir);
    auto put = [&](const char* f, int nc, auto get) {
      std::vector<float> v(c.size() * nc);
      for (std::size_t n = 0; n < c.size(); ++n)
        for (int m = 0; m < nc; ++m) v[n * nc + m] = get(c[n], m);
      radray::write_raw_f32(dir + "/" + f + ".bin", v);
    };
    put("fHI", 1, [](const rrk::CellData& x, int) { return x.fHI; });
    put("fHeI", 1, [](const rrk::CellData& x, int) { return x.fHeI; });
    put("fHeII", 1, [](const rrk::CellData& x, int) { return x.fHeII; });
    put("temp", 1, [](const rrk::CellData& x, int) { return x.temp; });
    put("rhog", 1, [](const rrk::CellData& x, int) { return x.rhog; });
    put("ngam", 5, [](const rrk::CellData& x, int m) { return x.ngam[m]; });
    put("nu", 5, [](const rrk::CellData& x, int m) { return x.nu[m]; });
    put("fgam", 5, [](const rrk::CellData& x, int m) { return x.fgam[m]; });
    put("tau", 5, [](const rrk::CellData& x, int m) { return x.tau[m]; });
    put("gamma", 3, [](const rrk::CellData& x, int m) { return x.gamma[m]; });
    put("lambda", 3, [](const rrk::CellData& x, int m) { return x.lambda[m]; });
    put("dx", 3, [](const rrk::CellData& x, int m) { return x.dx[m]; });
    put("u_gas", 5, [](const rrk::CellData& x, int m) { return x.u[m]; });
    put("w_gas", 4, [](const rrk::CellData& x, int m) { return x.w[m]; });
    std::FILE* fp = std::fopen((dir + "/stats.txt").c_str(), "w");
    if (!fp) throw std::runtime_error("cannot write stats.txt in " + dir);
    std::fprintf(fp, "istep = %d\n", clk.istep);
    std::fprintf(fp, "redshift = %.17g\n", clk.redshift);
    std::fprintf(fp, "scalefactor = %.17g\n", clk.scalefactor);
    std::fprintf(fp, "nstep_rt = %d\n", clk.nstep_rt);
    std::fprintf(fp, "nray_active = %lld\n", st.nray_active);
    std::fprintf(fp, "ncst1 = %lld\n", st.ncst1);
    std::fprintf(fp, "nmrg1 = %lld\n", st.nmrg1);
    std::fprintf(fp, "nspl1 = %lld\n", st.nspl1);
    std::fprintf(fp, "ndel1 = %lld\n", st.ndel1);
    std::fprintf(fp, "ngam_rays_total = %.17g\n", rrk::from_fx(st.ngam_tot_fx, rrk::kScaleTot));
    std::fclose(fp);
  }
};

// Initial CellData of global cell c from the gas input
inline rrk::CellData initial_cell(const radray::GasGrid& g, long long c) {
  rrk::CellData d{};
  for (int m = 0; m < 5; ++m) d.u[m] = g.u[c * 5 + m];
  for (int m = 0; m < 4; ++m) d.w[m] = g.w[c * 4 + m];
  return d;
}

}  // namespace rrhost
