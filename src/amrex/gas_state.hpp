// The gas state on the level (amrex_state = 1, docs/stage_e_state.md) and
// the slab transfers between a level and rank 0's host arrays. Internal to
// the AMReX backends (src/amrex/*.cpp); the AMReX-free interface is
// include/radray/amrex_state.hpp.
#pragma once
#include <cstdint>
#include <memory>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace radray {

struct Config;

namespace gs {

namespace ax = amrex;
using FFab = ax::FabArray<ax::BaseFab<float>>;

// The components of the gas MultiFab: u(1:5), e, w(1:4) (the order of the
// hydro level's U)
constexpr int kU = 0, kE = 5, kW = 6, kNC = 10;

struct GasState {
  int N = 0;
  ax::Box domain;
  ax::Geometry geom;
  ax::Periodicity period;
  ax::BoxArray ba;
  ax::DistributionMapping dm;
  FFab G;
  // Stage F3a (docs/stage_f.md section 3): the refined levels of the
  // hydro hierarchy (amrex_max_level > 0), float like G, on the
  // hierarchy's BoxArray and mapping of each level; the hydro copies them
  // in and out at every call and (re)defines them at a regrid, so the gas
  // of every level lives here between steps (the RT of F3b reads them).
  // fine[l - 1] is level l; finest is the deepest level in use.
  struct FineLevel {
    ax::BoxArray ba;
    ax::DistributionMapping dm;
    FFab G;
  };
  std::vector<FineLevel> fine;
  int finest = 0;
};

// The active gas state (null unless amrex_state = 1 and amrex_state_init
// has run); the levels take their BoxArray and DistributionMapping from it.
GasState* active();

// The decomposition of an N^3 level (cfg: amrex_slab, amrex_slab_planes,
// amrex_max_grid_size): boxes of amrex_max_grid_size with AMReX's default
// mapping, or for amrex_slab = 1 z-slabs (whole i, j planes) mapped to the
// ranks in z order, slab s to rank s * nranks / nslab (the exact chains of
// docs/stage_e_moments.md need that order). Same on every rank.
void make_decomposition(const Config& cfg, int N, ax::BoxArray& ba, ax::DistributionMapping& dm);

// k-planes per slab of the transfers (a one-box FabArray of the slab on
// rank 0, ParallelCopy): rank 0's transient is one slab, not the grid.
constexpr int kSlabPlanes = 8;

// Rank 0 visits the whole level slab by slab: fn(k0, nk, a) on rank 0 with
// a = Array4<float const> of the slab (indices i, j in 0..N-1, k in
// k0..k0+nk-1, components 0..ncomp-1). Collective.
template <class F>
void gather_slabs(const FFab& src, int N, int ncomp, F fn) {
  const ax::DistributionMapping dm0(ax::Vector<int>{0});
  for (int k0 = 0; k0 < N; k0 += kSlabPlanes) {
    const int nk = std::min(kSlabPlanes, N - k0);
    const ax::Box slab(ax::IntVect(0, 0, k0), ax::IntVect(N - 1, N - 1, k0 + nk - 1));
    FFab s0(ax::BoxArray(slab), dm0, ncomp, 0);
    s0.ParallelCopy(src, 0, 0, ncomp);
    if (ax::ParallelDescriptor::IOProcessor()) fn(k0, nk, s0.const_array(0));
  }
}

// The reverse: rank 0 fills each slab (fill(k0, nk, a) with a =
// Array4<float>) and it is copied onto the level. Collective.
template <class F>
void scatter_slabs(FFab& dst, int N, int ncomp, F fill) {
  const ax::DistributionMapping dm0(ax::Vector<int>{0});
  for (int k0 = 0; k0 < N; k0 += kSlabPlanes) {
    const int nk = std::min(kSlabPlanes, N - k0);
    const ax::Box slab(ax::IntVect(0, 0, k0), ax::IntVect(N - 1, N - 1, k0 + nk - 1));
    FFab s0(ax::BoxArray(slab), dm0, ncomp, 0);
    if (ax::ParallelDescriptor::IOProcessor()) fill(k0, nk, s0.array(0));
    dst.ParallelCopy(s0, 0, 0, ncomp);
  }
}

}  // namespace gs

// The owner, held by Simulation::state
struct AmrexState {
  gs::GasState gas;
  ~AmrexState();
};

// Load balancing (amrex_balance = 1): implemented by the levels, called by
// amrex_balance_step in gas_state.cpp. rt_box_costs adds the live ray count
// of every box (global box index) into cost; rt_remap and pm_remap move a
// level's persistent data onto a new mapping of the same BoxArray; the
// hydro level is rebuilt lazily (amrex_hydro_release).
struct Simulation;
namespace pm {
struct PMState;
}
void rt_box_ray_counts(Simulation& s, amrex::Vector<amrex::Real>& count);
void rt_remap(Simulation& s, const amrex::DistributionMapping& dm);
void pm_remap(pm::PMState& p, const amrex::DistributionMapping& dm);
// The hydro hierarchy's level 0 onto the new mapping (its data is copied
// in from the gas state at every step; the refined levels keep their own
// mapping); nothing without a hierarchy.
void amrex_hydro_remap(const amrex::DistributionMapping& dm);
// The hydro's level-0 reference for its push-down (docs/stage_f.md section
// 1) set to the state's current level 0: called by the RT hierarchy after
// its step, which changed the refined levels itself and averaged them down
// (so that change is not pushed down again); the PM's later changes still
// are.
void amrex_hydro_level0_is_current();

}  // namespace radray
