// The gas state on the level (amrex_state = 1; amrex_state.hpp,
// docs/stage_e_state.md): one FabArray of floats (u, e, w) on the
// decomposition shared by the RT, hydro and PM levels, and the slab
// transfers to and from rank 0's host mirror.
#include "gas_state.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <AMReX.H>
#include <AMReX_MFIter.H>

#include <AMReX_DistributionMapping.H>

#include "radray/amrex_state.hpp"
#include "radray/fortran.hpp"
#include "radray/hydro_backend.hpp"
#include "radray/rt_backend.hpp"
#include "radray/simulation.hpp"

namespace radray {

namespace {
gs::GasState* g_active = nullptr;

gs::GasState& state_of(Simulation& s) {
  if (!s.state) throw std::runtime_error("amrex_state = 1: amrex_state_init was not called");
  return s.state->gas;
}

// Rank 0's host arrays must hold the grid (the mirror)
void check_mirror(const Simulation& s) {
  if (!io_rank()) return;
  const GasGrid& g = s.gas;
  if (static_cast<std::int64_t>(g.u.size()) != 5 * g.ncell || static_cast<std::int64_t>(g.w.size()) != 4 * g.ncell ||
      static_cast<std::int64_t>(g.e.size()) != g.ncell)
    throw std::runtime_error("amrex_state = 1: rank 0's host gas arrays are not allocated (u, w, e)");
}
}  // namespace

gs::GasState* gs::active() { return g_active; }

AmrexState::~AmrexState() {
  if (g_active == &gas) g_active = nullptr;
}

void amrex_state_release() { g_active = nullptr; }

void amrex_state_init(Simulation& s) {
  namespace ax = amrex;
  check_amrex_state_param(s.cfg);
  if (s.state) return;
  if (!amrex::Initialized()) throw std::runtime_error("amrex_state = 1: backend_init was not called");
  auto st = std::make_shared<AmrexState>();
  gs::GasState& g = st->gas;
  g.N = s.cfg.Ngas;
  g.domain = ax::Box(ax::IntVect(0), ax::IntVect(g.N - 1));
  const ax::RealBox rb({0., 0., 0.}, {double(g.N), double(g.N), double(g.N)});  // dx = 1
  const ax::Array<int, 3> per{1, 1, 1};
  g.geom = ax::Geometry(g.domain, rb, ax::CoordSys::cartesian, per);
  g.period = g.geom.periodicity();
  g.ba = ax::BoxArray(g.domain);
  g.ba.maxSize(s.cfg.amrex_max_grid_size);
  g.dm = ax::DistributionMapping(g.ba);
  g.G.define(g.ba, g.dm, gs::kNC, 0);
  g.G.setVal(0.0f);
  s.state = st;
  g_active = &st->gas;
}

void amrex_state_scatter_gas(Simulation& s) {
  gs::GasState& g = state_of(s);
  check_mirror(s);
  const int N = g.N;
  const float *u = s.gas.u.data(), *w = s.gas.w.data(), *e = s.gas.e.data();
  gs::scatter_slabs(g.G, N, gs::kNC, [&](int k0, int nk, amrex::Array4<float> const& a) {
#pragma omp parallel for schedule(static) collapse(2)
    for (int k = k0; k < k0 + nk; ++k)
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
          const std::int64_t c = i + static_cast<std::int64_t>(N) * (j + static_cast<std::int64_t>(N) * k);
          for (int m = 0; m < 5; ++m) a(i, j, k, gs::kU + m) = u[c * 5 + m];
          a(i, j, k, gs::kE) = e[c];
          for (int m = 0; m < 4; ++m) a(i, j, k, gs::kW + m) = w[c * 4 + m];
        }
  });
}

void amrex_state_gather_gas(Simulation& s) {
  gs::GasState& g = state_of(s);
  check_mirror(s);
  const int N = g.N;
  float *u = s.gas.u.data(), *w = s.gas.w.data(), *e = s.gas.e.data();
  gs::gather_slabs(g.G, N, gs::kNC, [&](int k0, int nk, amrex::Array4<float const> const& a) {
#pragma omp parallel for schedule(static) collapse(2)
    for (int k = k0; k < k0 + nk; ++k)
      for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
          const std::int64_t c = i + static_cast<std::int64_t>(N) * (j + static_cast<std::int64_t>(N) * k);
          for (int m = 0; m < 5; ++m) u[c * 5 + m] = a(i, j, k, gs::kU + m);
          e[c] = a(i, j, k, gs::kE);
          for (int m = 0; m < 4; ++m) w[c * 4 + m] = a(i, j, k, gs::kW + m);
        }
  });
}

void amrex_state_clamp_gas(Simulation& s) {
  namespace ax = amrex;
  gs::GasState& g = state_of(s);
  const double rho_min = s.units.Dgas_min;
  for (ax::MFIter mfi(g.G); mfi.isValid(); ++mfi) {
    auto a = g.G.array(mfi);
    ax::ParallelForOMP(mfi.validbox(), [=](int i, int j, int k) {
      // summarize_gas: u(1) = max(u(1), Dgas_min) in double, w = max(w, 0) in float
      a(i, j, k, gs::kU) = static_cast<float>(f::gmax(static_cast<double>(a(i, j, k, gs::kU)), rho_min));
      for (int m = 0; m < 4; ++m) a(i, j, k, gs::kW + m) = f::gmax(a(i, j, k, gs::kW + m), 0.0f);
    });
  }
}

// ---- load balancing (amrex_balance = 1) --------------------------------------

namespace {
// mean over max of the per-rank cost sums under a mapping
double efficiency(const amrex::Vector<amrex::Real>& cost, const amrex::DistributionMapping& dm, int nranks) {
  std::vector<double> per(static_cast<std::size_t>(nranks), 0.0);
  for (int b = 0; b < static_cast<int>(cost.size()); ++b) per[static_cast<std::size_t>(dm[b])] += cost[b];
  double mx = 0, sum = 0;
  for (double c : per) {
    mx = std::max(mx, c);
    sum += c;
  }
  return mx > 0 ? sum / nranks / mx : 1.0;
}
}  // namespace

bool amrex_balance_step(Simulation& s, pm::PMState& pm, double& eff_old, double& eff_new) {
  namespace ax = amrex;
  gs::GasState& g = state_of(s);
  const int nranks = ax::ParallelDescriptor::NProcs();
  const int nbox = static_cast<int>(g.ba.size());
  // the cost of every box: its cells plus the weighted live rays
  ax::Vector<ax::Real> cost(nbox, 0.0);
  rt_box_ray_counts(s, cost);
  for (int b = 0; b < nbox; ++b)
    cost[b] = static_cast<ax::Real>(g.ba[b].numPts()) + s.cfg.amrex_balance_ray_weight * cost[b];
  eff_old = efficiency(cost, g.dm, nranks);
  ax::DistributionMapping dm2 = s.cfg.amrex_balance_method == 1 ? ax::DistributionMapping::makeKnapSack(cost)
                                                                 : ax::DistributionMapping::makeSFC(cost, g.ba);
  eff_new = efficiency(cost, dm2, nranks);
  if (nranks == 1 || !(eff_new > eff_old + s.cfg.amrex_balance_min_gain)) {
    eff_new = eff_old;
    return false;
  }
  // the gas onto the new mapping, then the levels
  gs::FFab G2(g.ba, dm2, gs::kNC, 0);
  G2.ParallelCopy(g.G, 0, 0, gs::kNC);
  g.G = std::move(G2);
  g.dm = dm2;
  rt_remap(s, dm2);
  pm_remap(pm, dm2);
  amrex_hydro_remap(dm2);  // level 0 on the new mapping, the refined levels kept
  return true;
}

void amrex_state_add_metals(Simulation& s, const std::vector<std::int64_t>& cells, const std::vector<double>& dz) {
  namespace ax = amrex;
  gs::GasState& g = state_of(s);
  if (cells.size() != dz.size()) throw std::invalid_argument("amrex_state_add_metals: sizes differ");
  const int N = g.N;
  for (ax::MFIter mfi(g.G); mfi.isValid(); ++mfi) {
    const ax::Box bx = mfi.validbox();
    auto a = g.G.array(mfi);
    for (std::size_t n = 0; n < cells.size(); ++n) {
      const std::int64_t c = cells[n];
      const int i = static_cast<int>(c % N), j = static_cast<int>((c / N) % N), k = static_cast<int>(c / (static_cast<std::int64_t>(N) * N));
      if (!bx.contains(ax::IntVect(i, j, k))) continue;
      float& w3 = a(i, j, k, gs::kW + 3);
      w3 = static_cast<float>(static_cast<double>(w3) + dz[n]);
    }
  }
}

}  // namespace radray
