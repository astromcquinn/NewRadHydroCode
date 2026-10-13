// The gas state on the level (amrex_state = 1, docs/stage_e_state.md) and
// the slab transfers between a level and rank 0's host arrays. Internal to
// the AMReX backends (src/amrex/*.cpp); the AMReX-free interface is
// include/radray/amrex_state.hpp.
#pragma once
#include <cstdint>
#include <memory>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>

namespace radray {

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
};

// The active gas state (null unless amrex_state = 1 and amrex_state_init
// has run); the levels take their BoxArray and DistributionMapping from it.
GasState* active();

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

}  // namespace radray
