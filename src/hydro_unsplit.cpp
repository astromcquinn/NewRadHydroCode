// Host driver of the unsplit MUSCL-Hancock + HLLC scheme with the
// dual-energy formulation (see hydro_unsplit.hpp, docs/hydro_unsplit.md).
//
// Determinism (docs/parallel.md): every pass computes each cell (or face)
// from a fixed stencil of the previous pass with the same function, in
// per-cell or per-slab buffers; the only reductions are a max (the CFL
// speed) and integer counters, so the result is the same bits for any
// thread count.
#include "radray/hydro_unsplit.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

#include "radray/hydro_backend.hpp"
#include "radray/kernels/hydro_mh.hpp"
#include "radray/parallel.hpp"

namespace radray {

using rrh::NG;
using rrh::NQ;
using rrh::NS;
using rrh::NV;

// --------------------------------------------------------------------------
// Parameters
// --------------------------------------------------------------------------
bool set_hydro_scheme_param(Config& cfg, const std::string& key, const std::string& val) {
  auto bad = [&](const std::string& what) {
    throw std::runtime_error(key + ": " + what + " (got '" + val + "')");
  };
  if (key == "hydro_scheme") {
    if (val == "legacy") cfg.hydro_scheme = HydroScheme::Legacy;
    else if (val == "unsplit") cfg.hydro_scheme = HydroScheme::Unsplit;
    else bad("must be legacy or unsplit");
  } else if (key == "unsplit_cfl") {
    cfg.unsplit_cfl = std::stod(val);
    if (!(cfg.unsplit_cfl > 0 && cfg.unsplit_cfl <= 1)) bad("must be in (0, 1]");
  } else if (key == "dual_energy") {
    const int v = std::stoi(val);
    if (v != 0 && v != 1) bad("must be 0 or 1");
    cfg.dual_energy = v == 1;
  } else if (key == "dual_eta1") {
    cfg.dual_eta1 = std::stod(val);
    if (!(cfg.dual_eta1 > 0 && cfg.dual_eta1 < 1)) bad("must be in (0, 1)");
  } else if (key == "dual_eta2") {
    cfg.dual_eta2 = std::stod(val);
    if (!(cfg.dual_eta2 > 0 && cfg.dual_eta2 < 1)) bad("must be in (0, 1)");
  } else if (key == "dual_nb_radius") {
    cfg.dual_nb_radius = std::stoi(val);
    if (cfg.dual_nb_radius < 1 || cfg.dual_nb_radius > 2) bad("must be 1 or 2");
  } else if (key == "dual_reset_rule") {
    cfg.dual_reset_rule = std::stoi(val);
    if (cfg.dual_reset_rule < 0 || cfg.dual_reset_rule > 4) bad("must be 0..4");
  } else if (key == "dual_eta2_local") {
    const int v = std::stoi(val);
    if (v != 0 && v != 1) bad("must be 0 or 1");
    cfg.dual_eta2_local = v == 1;
  } else if (key == "dual_p_old") {
    const int v = std::stoi(val);
    if (v != 0 && v != 1) bad("must be 0 or 1");
    cfg.dual_p_old = v == 1;
  } else if (key == "dual_eta1_local") {
    const int v = std::stoi(val);
    if (v != 0 && v != 1) bad("must be 0 or 1");
    cfg.dual_eta1_local = v == 1;
  } else {
    return false;
  }
  return true;
}

void copy_hydro_scheme(const Config& from, Config& to) {
  to.hydro_scheme = from.hydro_scheme;
  to.unsplit_cfl = from.unsplit_cfl;
  to.dual_energy = from.dual_energy;
  to.dual_eta1 = from.dual_eta1;
  to.dual_eta2 = from.dual_eta2;
  to.dual_eta1_local = from.dual_eta1_local;
  to.dual_nb_radius = from.dual_nb_radius;
  to.dual_reset_rule = from.dual_reset_rule;
  to.dual_eta2_local = from.dual_eta2_local;
  to.dual_p_old = from.dual_p_old;
}

std::string hydro_scheme_summary(const Config& cfg) {
  if (cfg.hydro_scheme == HydroScheme::Legacy) return "hydro_scheme = legacy";
  return "hydro_scheme = unsplit, unsplit_cfl = " + std::to_string(cfg.unsplit_cfl) +
         ", dual_energy = " + (cfg.dual_energy ? "1" : "0") +
         ", dual_eta1 = " + std::to_string(cfg.dual_eta1) +
         ", dual_eta2 = " + std::to_string(cfg.dual_eta2) +
         ", dual_eta1_local = " + (cfg.dual_eta1_local ? "1" : "0") +
         ", dual_nb_radius = " + std::to_string(cfg.dual_nb_radius) +
         ", dual_reset_rule = " + std::to_string(cfg.dual_reset_rule) +
         ", dual_eta2_local = " + (cfg.dual_eta2_local ? "1" : "0") +
         ", dual_p_old = " + (cfg.dual_p_old ? "1" : "0");
}

namespace {

// Work arrays of one step, kept between calls (one simulation per process).
struct Work {
  int n = 0, m = 0;
  std::vector<double> U, S;    // padded (m^3): conserved state, scalars
  std::vector<double> Es;      // padded (m^3): selected thermal energy density
  std::vector<double> Q;       // (n+2)^3: predictor states
  std::vector<double> Un, Sn;  // interior (n^3): updated state, scalars

  void resize(int nn) {
    if (nn == n) return;
    n = nn;
    m = n + 2 * NG;
    const std::int64_t mc = static_cast<std::int64_t>(m) * m * m;
    const std::int64_t qc = static_cast<std::int64_t>(n + 2) * (n + 2) * (n + 2);
    const std::int64_t nc = static_cast<std::int64_t>(n) * n * n;
    U.assign(mc * NV, 0.0);
    S.assign(mc * NS, 0.0);
    Es.assign(mc, 0.0);
    Q.assign(qc * NQ, 0.0);
    Un.assign(nc * NV, 0.0);
    Sn.assign(nc * NS, 0.0);
  }
  // padded index, i,j,k in [-NG, n+NG)
  std::int64_t at(int i, int j, int k) const {
    return (i + NG) + static_cast<std::int64_t>(m) * ((j + NG) + static_cast<std::int64_t>(m) * (k + NG));
  }
  // predictor index, i,j,k in [-1, n+1)
  std::int64_t aq(int i, int j, int k) const {
    const std::int64_t mq = n + 2;
    return (i + 1) + mq * ((j + 1) + mq * (k + 1));
  }
  // interior index, i,j,k in [0, n)
  std::int64_t ai(int i, int j, int k) const {
    return i + static_cast<std::int64_t>(n) * (j + static_cast<std::int64_t>(n) * k);
  }
};

Work& work(int n) {
  static std::unique_ptr<Work> w;
  if (!w) w = std::make_unique<Work>();
  w->resize(n);
  return *w;
}

// Periodic ghost cells of U and S (state = true) or of Es
void fill_ghosts(Work& W, bool state) {
  const int n = W.n;
  const bool par_ok = par::worth(static_cast<long long>(n) * n * n);
#pragma omp parallel for schedule(static) if (par_ok)
  for (int k = -NG; k < n + NG; ++k)
    for (int j = -NG; j < n + NG; ++j)
      for (int i = -NG; i < n + NG; ++i) {
        if (i >= 0 && i < n && j >= 0 && j < n && k >= 0 && k < n) continue;
        const int ii = (i + n) % n, jj = (j + n) % n, kk = (k + n) % n;
        const std::int64_t d = W.at(i, j, k), s = W.at(ii, jj, kk);
        if (state) {
          for (int c = 0; c < NV; ++c) W.U[d * NV + c] = W.U[s * NV + c];
          for (int c = 0; c < NS; ++c) W.S[d * NS + c] = W.S[s * NS + c];
        } else {
          W.Es[d] = W.Es[s];
        }
      }
}

// The selected thermal energy density of every interior cell (the eta1
// rule against the cell's E or the maximum over its face neighbours),
// from the ghost-filled U and S; then its ghosts.
void select_pass(Work& W, const rrh::DualEnergy& de) {
  const int n = W.n;
  const bool par_ok = par::worth(static_cast<long long>(n) * n * n);
  auto Ea = [&W](int i, int j, int k) { return W.U[W.at(i, j, k) * NV + 4]; };
#pragma omp parallel for schedule(static) if (par_ok)
  for (int k = 0; k < n; ++k)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        const std::int64_t c = W.at(i, j, k);
        const double* U = &W.U[c * NV];
        double Eref = U[4], e1;
        if (!de.eta1_local) Eref = rrh::face_neighbour_max(Ea, i, j, k, de.nb_radius, e1);
        W.Es[c] = rrh::select_thermal(U, W.S[c * NS], Eref, de.eta1);
      }
  fill_ghosts(W, false);
}

using Face = rrh::FaceData;

// The face between cell a and cell b = a + e_d, direction d (rrh::face_all
// on the work arrays)
inline void compute_face(const Work& W, int ia, int ja, int ka, int d, double dtdx, Face& f) {
  auto Qa = [&W](int i, int j, int k) { return &W.Q[W.aq(i, j, k) * NQ]; };
  auto Ua = [&W](int i, int j, int k, int c) { return W.U[W.at(i, j, k) * NV + c]; };
  auto Sa = [&W](int i, int j, int k, int s) { return W.S[W.at(i, j, k) * NS + s]; };
  rrh::face_all(Qa, Ua, Sa, ia, ja, ka, d, dtdx, f);
}

// Thread-local face buffers of one slab
struct SlabFaces {
  std::vector<Face> fx, fy, fz0, fz1;
  void resize(int n) {
    fx.resize(static_cast<std::size_t>(n + 1) * n);
    fy.resize(static_cast<std::size_t>(n) * (n + 1));
    fz0.resize(static_cast<std::size_t>(n) * n);
    fz1.resize(static_cast<std::size_t>(n) * n);
  }
};

struct SubstepCounts {
  std::int64_t from_e = 0, synced = 0, floor_T = 0, floor_rho = 0, reset = 0;
};

// One substep of length dt on the padded state (ghosts filled)
SubstepCounts substep(Work& W, double dt, const rrh::DualEnergy& de, double rho_min,
                      double Tmin, bool p_old) {
  const int n = W.n;
  const double dtdx = dt;  // dx = 1 cell
  const bool par_ok = par::worth(static_cast<long long>(n) * n * n);
  // 1. predictor on the interior plus one ghost layer
  auto Ua = [&W](int i, int j, int k, int c) { return W.U[W.at(i, j, k) * NV + c]; };
  auto Ea = [&W](int i, int j, int k) { return W.Es[W.at(i, j, k)]; };
#pragma omp parallel for schedule(static) if (par_ok)
  for (int k = -1; k < n + 1; ++k)
    for (int j = -1; j < n + 1; ++j)
      for (int i = -1; i < n + 1; ++i)
        rrh::predict_cell(Ua, Ea, i, j, k, dtdx, &W.Q[W.aq(i, j, k) * NQ]);

  // 2. fluxes and update, slab by slab (the z faces of a slab are computed
  // by that slab, so no data is shared between iterations)
  const int nth = par_ok ? par::max_threads() : 1;
  std::vector<SlabFaces> faces(nth);
  for (auto& f : faces) f.resize(n);
#pragma omp parallel for schedule(static) if (par_ok)
  for (int k = 0; k < n; ++k) {
    SlabFaces& sf = faces[par::thread_id()];
    for (int j = 0; j < n; ++j)
      for (int i = 0; i <= n; ++i) compute_face(W, i - 1, j, k, 0, dtdx, sf.fx[i + (n + 1) * j]);
    for (int j = 0; j <= n; ++j)
      for (int i = 0; i < n; ++i) compute_face(W, i, j - 1, k, 1, dtdx, sf.fy[i + n * j]);
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        compute_face(W, i, j, k - 1, 2, dtdx, sf.fz0[i + n * j]);
        compute_face(W, i, j, k, 2, dtdx, sf.fz1[i + n * j]);
      }
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        const std::int64_t c = W.at(i, j, k), cn = W.ai(i, j, k);
        const Face& x0 = sf.fx[i + (n + 1) * j];
        const Face& x1 = sf.fx[i + 1 + (n + 1) * j];
        const Face& y0 = sf.fy[i + n * j];
        const Face& y1 = sf.fy[i + n * (j + 1)];
        const Face& z0 = sf.fz0[i + n * j];
        const Face& z1 = sf.fz1[i + n * j];
        // the cell pressure of the work term: time-centred (Hancock) or at
        // the old time (Bryan et al. 2014 eq. 43)
        const double p_half = p_old ? (rrh::kGamma - 1) * rrf::max(W.Es[c], 0.0)
                                    : W.Q[W.aq(i, j, k) * NQ + 4];
        rrh::update_from_faces(&W.U[c * NV], &W.S[c * NS], x0, x1, y0, y1, z0, z1, p_half, dtdx,
                               &W.Un[cn * NV], &W.Sn[cn * NS]);
      }
  }

  // 3. floors, dual-energy selection and synchronisation; back into U, S
  std::int64_t n_from_e = 0, n_synced = 0, n_floor_T = 0, n_floor_rho = 0, n_reset = 0;
#pragma omp parallel for schedule(static) if (par_ok) \
    reduction(+ : n_from_e, n_synced, n_floor_T, n_floor_rho, n_reset)
  for (int k = 0; k < n; ++k)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        const std::int64_t cn = W.ai(i, j, k), c = W.at(i, j, k);
        auto Ena = [&W, n](int ii, int jj, int kk) {
          return W.Un[W.ai((ii + n) % n, (jj + n) % n, (kk + n) % n) * NV + 4];
        };
        double Enb1;
        const double Enb = rrh::face_neighbour_max(Ena, i, j, k, de.nb_radius, Enb1);
        double* U = &W.U[c * NV];
        for (int m = 0; m < NV; ++m) U[m] = W.Un[cn * NV + m];
        double eD = W.Sn[cn * NS];
        const rrh::CellFlags f = rrh::finalize_cell(U, eD, Enb, Enb1, de, rho_min, Tmin);
        W.S[c * NS] = eD;
        for (int s = 1; s < NS; ++s) W.S[c * NS + s] = W.Sn[cn * NS + s];
        n_from_e += f.from_e;
        n_synced += f.synced;
        n_floor_T += f.floor_T;
        n_floor_rho += f.floor_rho;
        n_reset += f.reset;
      }
  return SubstepCounts{n_from_e, n_synced, n_floor_T, n_floor_rho, n_reset};
}

// Largest signal speed over the interior (after select_pass)
double max_signal(const Work& W) {
  const int n = W.n;
  const bool par_ok = par::worth(static_cast<long long>(n) * n * n);
  double smax = 0;
#pragma omp parallel for schedule(static) if (par_ok) reduction(max : smax)
  for (int k = 0; k < n; ++k)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        const std::int64_t c = W.at(i, j, k);
        smax = rrf::max(smax, rrh::signal_sum(&W.U[c * NV], W.Es[c]));
      }
  return smax;
}

// Folds a non-zero grid velocity into u (lab frame) and zeroes it. Returns
// whether anything was changed.
bool fold_frame(GasGrid& gas, HydroState& hs) {
  bool any = false;
  for (float x : hs.v)
    if (x != 0.0f) {
      any = true;
      break;
    }
  if (!any) return false;
  for (std::int64_t c = 0; c < gas.ncell; ++c) {
    float* u = gas.u_p(c);
    float* v = hs.v_p(c);
    double U[NV];
    for (int m = 0; m < NV; ++m) U[m] = u[m];
    const double rho = U[0];
    const double ke_frame = rrh::kinetic_energy(U);
    for (int m = 0; m < 3; ++m) U[1 + m] = U[1 + m] + rho * static_cast<double>(v[m]);
    U[4] = U[4] - ke_frame + rrh::kinetic_energy(U);
    for (int m = 1; m < NV; ++m) u[m] = static_cast<float>(U[m]);
    v[0] = v[1] = v[2] = 0.0f;
  }
  return true;
}

double eta1_of(const Config& cfg) { return cfg.dual_energy ? cfg.dual_eta1 : 0.0; }

}  // namespace

// --------------------------------------------------------------------------
void init_unsplit_state(const Config& cfg, const Units& un, GasGrid& gas, HydroState& hs) {
  if (cfg.hydro_scheme != HydroScheme::Unsplit) return;
  if (cfg.amrex_state && gas.u.empty()) {
    // amrex_state = 1 on a rank other than 0: the gas (and e) is on the level
    gas.e_on = true;
    return;
  }
  if (fold_frame(gas, hs)) hs.ucounters.frame_shift = 1;
  if (gas.has_e()) return;
  const double Tmin = un.Thydro_min;  // thermal_floor (1 K by default)
  gas.e.assign(static_cast<std::size_t>(gas.ncell), 0.0f);
  for (std::int64_t c = 0; c < gas.ncell; ++c) {
    const float* u = gas.u_p(c);
    double U[NV];
    for (int m = 0; m < NV; ++m) U[m] = u[m];
    const double rho = rrf::max(U[0], un.Dgas_min);
    double e = U[4] - rrh::kinetic_energy(U);
    if (!(e >= rho * Tmin)) e = rho * Tmin;
    gas.e[c] = static_cast<float>(e);
  }
}

void hydro_step_unsplit(const Config& cfg, const Units& un, double dt, GasGrid& gas,
                        HydroState& hs, double scalefactor) {
  if (!gas.has_e()) throw std::runtime_error("hydro_step_unsplit: call init_unsplit_state first");
  if (!(cfg.unsplit_cfl > 0)) throw std::runtime_error("unsplit_cfl must be > 0");
  const int n = gas.N;
  // amrex_state = 1: the gas is on the level and v is identically zero
  // (the unsplit scheme never shifts the frame), nothing to fold
  if (!cfg.amrex_state && fold_frame(gas, hs)) hs.ucounters.frame_shift = 1;
  hs.dt_hy = dt;
  if (use_amrex_hydro(cfg)) {  // the same step on the AMReX level (hydro_backend.hpp)
    amrex_hydro_step(cfg, un, dt, scalefactor, gas, hs);
    return;
  }
  Work& W = work(n);
  UnsplitCounters& C = hs.ucounters;

  // float -> double
  const bool par_ok = par::worth(static_cast<long long>(n) * n * n);
#pragma omp parallel for schedule(static) if (par_ok)
  for (int k = 0; k < n; ++k)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        const std::int64_t cg = gas.cell(i + 1, j + 1, k + 1), c = W.at(i, j, k);
        const float* u = &gas.u[cg * 5];
        const float* w = &gas.w[cg * 4];
        for (int m = 0; m < NV; ++m) W.U[c * NV + m] = u[m];
        W.S[c * NS] = gas.e[cg];
        for (int s = 0; s < 4; ++s) W.S[c * NS + 1 + s] = w[s];
      }

  const rrh::DualEnergy de{eta1_of(cfg), cfg.dual_eta2, cfg.dual_eta1_local, cfg.dual_nb_radius,
                           cfg.dual_reset_rule, cfg.dual_eta2_local};
  const double rho_min = un.Dgas_min;
  const double Tmin = un.Thydro_min;  // thermal_floor (1 K by default)

  double remaining = dt;
  std::int64_t nsub_call = 0;
  while (remaining > 0) {
    fill_ghosts(W, true);
    select_pass(W, de);
    const double smax = max_signal(W);
    if (!(smax > 0) || !std::isfinite(smax))
      throw std::runtime_error("hydro_step_unsplit: non-finite or zero signal speed");
    const double dt_cfl = cfg.unsplit_cfl / smax;
    C.dt_cfl_min = std::min(C.dt_cfl_min, dt_cfl);
    int nsub = static_cast<int>(std::ceil(remaining / dt_cfl - 1e-9));
    if (nsub < 1) nsub = 1;
    const double dt_sub = nsub == 1 ? remaining : remaining / nsub;
    const SubstepCounts sc = substep(W, dt_sub, de, rho_min, Tmin, cfg.dual_p_old);
    C.cells_from_e += sc.from_e;
    C.cells_synced += sc.synced;
    C.floor_T += sc.floor_T;
    C.floor_rho += sc.floor_rho;
    C.cells_reset += sc.reset;
    ++nsub_call;
    remaining = nsub == 1 ? 0.0 : remaining - dt_sub;
  }
  C.calls += 1;
  C.substeps += nsub_call;
  C.max_substeps = std::max(C.max_substeps, nsub_call);

  // double -> float
#pragma omp parallel for schedule(static) if (par_ok)
  for (int k = 0; k < n; ++k)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        const std::int64_t cg = gas.cell(i + 1, j + 1, k + 1), c = W.at(i, j, k);
        float* u = &gas.u[cg * 5];
        float* w = &gas.w[cg * 4];
        for (int m = 0; m < NV; ++m) u[m] = static_cast<float>(W.U[c * NV + m]);
        gas.e[cg] = static_cast<float>(W.S[c * NS]);
        for (int s = 0; s < 4; ++s) w[s] = static_cast<float>(W.S[c * NS + 1 + s]);
      }
}

void hydro_half_step(const Config& cfg, const Units& un, const SimClock& clk, HydroHalf half,
                     GasGrid& gas, HydroState& hs) {
  if (cfg.hydro_scheme == HydroScheme::Unsplit) {
    hydro_step_unsplit(cfg, un, clk.dt_sim / 2, gas, hs, clk.scalefactor);
    return;
  }
  if (half == HydroHalf::First) {
    // radhydro.f90: forward hydro step with transposes in between
    transpose_ij(gas, hs);  // x,y,z -> y,x,z
    transpose_ik(gas, hs);  // y,x,z -> z,x,y
    hydro_step(cfg, un, clk, SweepOrder::Forward, gas, hs);
    transpose_ik(gas, hs);  // z,x,y -> y,x,z
    hydro_step(cfg, un, clk, SweepOrder::Forward, gas, hs);
    transpose_ij(gas, hs);  // y,x,z -> x,y,z
    hydro_step(cfg, un, clk, SweepOrder::Forward, gas, hs);
  } else {
    // Reverse hydro step with transposes in between
    hydro_step(cfg, un, clk, SweepOrder::Reverse, gas, hs);
    transpose_ij(gas, hs);  // x,y,z -> y,x,z
    hydro_step(cfg, un, clk, SweepOrder::Reverse, gas, hs);
    transpose_ik(gas, hs);  // y,x,z -> z,x,y
    hydro_step(cfg, un, clk, SweepOrder::Reverse, gas, hs);
    transpose_ik(gas, hs);  // z,x,y -> y,x,z
    transpose_ij(gas, hs);  // y,x,z -> x,y,z
  }
}

}  // namespace radray
