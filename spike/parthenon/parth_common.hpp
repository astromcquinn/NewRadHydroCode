// Small helpers shared by the two Parthenon prototypes: ghost exchange as a
// task collection, block geometry, rank-0 gather of interior data.
#pragma once

#include <mpi.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>

#include "bvals/comms/bvals_in_one.hpp"

namespace sp {

using parthenon::BoundaryType;
using parthenon::Mesh;
using parthenon::MeshBlock;
using parthenon::MeshData;
using parthenon::Real;
using parthenon::TaskCollection;
using parthenon::TaskID;
using parthenon::TaskRegion;
using parthenon::TaskStatus;

inline double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

[[noreturn]] inline void die(const char* msg) {
  std::fprintf(stderr, "fatal: %s\n", msg);
  MPI_Abort(MPI_COMM_WORLD, 3);
  std::abort();
}

// Ghost exchange (owner -> ghosts, periodic) of the FillGhost variables in
// `fields` (empty: all of them), with Parthenon's boundary-communication
// tasks: one task list per default block partition.
inline void exchange_ghosts(Mesh* pmesh, const std::string& label,
                            const std::vector<std::string>& fields) {
  TaskCollection tc;
  TaskID none(0);
  auto partitions = pmesh->GetDefaultBlockPartitions();
  TaskRegion& region = tc.AddRegion(partitions.size());
  for (std::size_t p = 0; p < partitions.size(); ++p) {
    auto& tl = region[p];
    auto& base = pmesh->mesh_data.Add("base", partitions[p]);
    std::shared_ptr<MeshData<Real>> md =
        fields.empty() ? base : pmesh->mesh_data.AddShallow(label, base, fields);
    auto start = tl.AddTask(none, parthenon::StartReceiveBoundBufs<BoundaryType::any>, md);
    parthenon::AddBoundaryExchangeTasks(start, tl, md, pmesh->multilevel);
  }
  if (tc.Execute() != parthenon::TaskListStatus::complete) die("ghost exchange failed");
}

// Interior geometry of a block of a uniform mesh with nglobal^3 cells:
// global index of the first interior cell and the interior size.
struct BlockGeom {
  int off[3];  // global 0-based cell index of the first interior cell
  int n[3];    // interior cells
  int is, js, ks;  // first interior index in Parthenon's (ghost-including) arrays
};
inline BlockGeom block_geom(MeshBlock* pmb, int nglobal, double xmin, double xmax) {
  BlockGeom g;
  using parthenon::IndexDomain;
  const auto& cb = pmb->cellbounds;
  g.n[0] = cb.ncellsi(IndexDomain::interior);
  g.n[1] = cb.ncellsj(IndexDomain::interior);
  g.n[2] = cb.ncellsk(IndexDomain::interior);
  g.is = cb.is(IndexDomain::interior);
  g.js = cb.js(IndexDomain::interior);
  g.ks = cb.ks(IndexDomain::interior);
  const double h = (xmax - xmin) / nglobal;
  g.off[0] = static_cast<int>(std::lround((pmb->block_size.xmin(parthenon::X1DIR) - xmin) / h));
  g.off[1] = static_cast<int>(std::lround((pmb->block_size.xmin(parthenon::X2DIR) - xmin) / h));
  g.off[2] = static_cast<int>(std::lround((pmb->block_size.xmin(parthenon::X3DIR) - xmin) / h));
  return g;
}

// Rank-0 gather of per-cell records: every rank passes (global cell index,
// record bytes) pairs; rank 0 gets a dense array of `ncell` records.
inline std::vector<char> gather_cells(const std::vector<long long>& gidx, const std::vector<char>& rec,
                                      std::size_t rec_size, long long ncell, MPI_Comm comm) {
  int rank, nr;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nr);
  const std::size_t per = sizeof(long long) + rec_size;
  std::vector<char> send(gidx.size() * per);
  for (std::size_t c = 0; c < gidx.size(); ++c) {
    std::memcpy(&send[c * per], &gidx[c], sizeof(long long));
    std::memcpy(&send[c * per + sizeof(long long)], &rec[c * rec_size], rec_size);
  }
  int mine = static_cast<int>(send.size());
  std::vector<int> cnt(nr), dsp(nr);
  MPI_Gather(&mine, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, comm);
  std::vector<char> recv;
  if (rank == 0) {
    long long tot = 0;
    for (int r = 0; r < nr; ++r) {
      dsp[r] = static_cast<int>(tot);
      tot += cnt[r];
    }
    if (tot != ncell * static_cast<long long>(per)) die("gather: wrong number of cells");
    recv.resize(tot);
  }
  MPI_Gatherv(send.data(), mine, MPI_BYTE, recv.data(), cnt.data(), dsp.data(), MPI_BYTE, 0, comm);
  std::vector<char> out;
  if (rank == 0) {
    out.assign(ncell * rec_size, 0);
    std::vector<char> seen(ncell, 0);
    for (long long c = 0; c < ncell; ++c) {
      long long g;
      std::memcpy(&g, &recv[c * per], sizeof(long long));
      if (g < 0 || g >= ncell || seen[g]) die("gather: bad cell index");
      seen[g] = 1;
      std::memcpy(&out[g * rec_size], &recv[c * per + sizeof(long long)], rec_size);
    }
  }
  return out;
}

inline int num_threads() { return Kokkos::DefaultExecutionSpace().concurrency(); }

}  // namespace sp
