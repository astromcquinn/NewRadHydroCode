// Port of ray_tracing_tools.f90 (default code path: RAYFAST and FREQ1
// undefined). See include/radray/ray_tracing.hpp for conventions.
#include "radray/ray_tracing.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "radray/constants.hpp"
#include "radray/fortran.hpp"
#include "radray/healpix.hpp"
#include "radray/kernels/ray_kernels.hpp"
#include "radray/parallel.hpp"
#include "radray/recomb_rad.hpp"
#include "radray/rt_backend.hpp"
#include "radray/rt_moments.hpp"
#include "radray/rt_tensor_diag.hpp"

namespace radray {

using namespace rayp;

// ---------------------------------------------------------------------------
// RayState (storage)
// ---------------------------------------------------------------------------

void RayBuffer::reset(std::int64_t n) {
  std::free(p);
  p = nullptr;
  cap = n;
  // calloc: lazily mapped zero pages; all-zero bits == Ray{}
  p = static_cast<Ray*>(std::calloc(static_cast<std::size_t>(n), sizeof(Ray)));
  if (!p) throw std::runtime_error("cannot reserve ray storage (" + std::to_string(n) + " rays)");
}

void RayState::init(std::int64_t nray_max, std::int64_t nfree_max) {
  if (nray_max >= (std::int64_t(1) << 31) - 1)
    throw std::runtime_error("Nray_max >= 2^31: ray links are int32 (see rays.hpp)");
  Nray_max = nray_max;
  Nfree_max = nfree_max;
  Nray = 0;
  Nrbuf = 0;
  ray.reset(Nray_max + 1);
  hwm_ = 0;
  ray_index_.clear();
  ray_index_.shrink_to_fit();
  ray_index_.reserve(static_cast<std::size_t>(Nray_max) + 1);
  ray_index_.resize(1);  // position 0 unused
  prov.reset(track_prov ? Nray_max + 1 : 0);
  ray_index_[0] = 0;
  stat = RayStat{};
  save = RayStat{};
  ngam_ray.fill(0.0);
  npool_limited = npool_forced = 0;
  pool_free.clear();
  next_fresh = 1;
  tl.clear();
}

void RayState::ensure(std::int64_t i) {
  if (i < 1 || i > Nray_max)
    throw std::runtime_error("ray storage index out of range: " + std::to_string(i) +
                             " (ray pool overflow, Nray_max = " + std::to_string(Nray_max) + ")");
  // ifree is sized on first use (Nfree_max = Nray_max/4 entries)
  if (stat.ifree.empty()) stat.ifree.assign(static_cast<std::size_t>(Nfree_max) + 1, 0);
  const std::int64_t old = static_cast<std::int64_t>(ray_index_.size());
  if (i >= old) {
    ray_index_.resize(static_cast<std::size_t>(i) + 1);
    for (std::int64_t p = old; p <= i; ++p) ray_index_[p] = static_cast<std::int32_t>(p);
  }
  if (i > hwm_) hwm_ = i;
}

// deterministic_parallel: take nr free slots for thread t. Slots come from
// the thread's list, refilled in chunks from the shared pool (freed slots
// first, then never-used ones). Which slots a ray gets depends on thread
// timing, but no result depends on storage indices with legacy_rt = 0
// (docs/parallel.md).
void RayState::det_take(ThreadRays& t, int nr, std::int64_t* indx) {
  constexpr std::int64_t kChunk = 4096;
  if (static_cast<std::int64_t>(t.free.size()) < nr) {
    std::lock_guard<std::mutex> lk(pool_mtx);
    while (static_cast<std::int64_t>(t.free.size()) < kChunk && !pool_free.empty()) {
      t.free.push_back(pool_free.back());
      pool_free.pop_back();
    }
    const std::int64_t want = kChunk - static_cast<std::int64_t>(t.free.size());
    const std::int64_t avail = Nray_max + 1 - next_fresh;
    const std::int64_t n = want < avail ? want : avail;
    for (std::int64_t m = 0; m < n; ++m) t.free.push_back(static_cast<std::int32_t>(next_fresh + m));
    next_fresh += n > 0 ? n : 0;
    if (next_fresh - 1 > hwm_) hwm_ = next_fresh - 1;
    if (static_cast<std::int64_t>(t.free.size()) < nr)
      throw std::runtime_error("ray pool overflow (Nray_max = " + std::to_string(Nray_max) +
                               "); raise nray_max_factor");
  }
  for (int m = 0; m < nr; ++m) {
    indx[m] = t.free.back();
    t.free.pop_back();
  }
  t.ncreated += nr;
}

void RayState::det_give(ThreadRays& t, std::int64_t ir) {
  constexpr std::size_t kSpill = 16384, kKeep = 8192;
  t.free.push_back(static_cast<std::int32_t>(ir));
  t.ndel1 += 1;
  if (t.free.size() > kSpill) {
    std::lock_guard<std::mutex> lk(pool_mtx);
    while (t.free.size() > kKeep) {
      pool_free.push_back(t.free.back());
      t.free.pop_back();
    }
  }
}

// ---------------------------------------------------------------------------
// init_ray_tracing
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: init_ray_tracing
// (The Fortran also writes Log/nside_of_r.dat and Log/r_of_nside.dat; these
// diagnostics are not reproduced.)
void init_ray_tracing(Simulation& s) {
  // the moment RT on the AMReX level holds no rays: no ray storage (Nray_max
  // N^3 rays would be calloc'ed on every rank; s.cfg.Nray_max, which the
  // checkpoint records, is unchanged)
  const bool no_rays = use_amrex(s.cfg) && use_moments(s.cfg);
  s.rays.init(no_rays ? 0 : s.cfg.Nray_max, s.cfg.Nfree_max);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

namespace {

// 1 + mod(int(x),Nrt): RT cell index (1-based) of a position in RT units.
inline int cell_of(double x, int Nrt) { return 1 + f::mod(f::int4(x), Nrt); }

// deterministic_parallel: the cell being processed by this thread, and the
// reach check. Every chain insertion must stay within dmax = d/2 - cast
// cells (periodic Chebyshev distance) of it, where cast is the reach of a
// cast that follows the insertion (k = rt_cells_per_substep cells in
// ray_tracing_1, 1 in the original; none in ray_tracing_3), so that all
// the ray work of a cell stays within d/2 and concurrent domains of one
// colour (>= d cells apart) never touch the same cell or ray chain.
thread_local int det_cell[3] = {0, 0, 0};
thread_local int det_cast = 1;
inline void det_set_cell(int i, int j, int k, int cast) {
  det_cell[0] = i;
  det_cell[1] = j;
  det_cell[2] = k;
  det_cast = cast;
}
void det_check_reach(const Simulation& s, int i, int j, int k) {
  const int N = s.cfg.Nrt, dmax = s.ddom.d / 2 - det_cast;
  const int ijk[3] = {i, j, k};
  for (int m = 0; m < 3; ++m) {
    int dd = std::abs(ijk[m] - det_cell[m]);
    dd = std::min(dd, N - dd);
    if (dd > dmax)
      throw std::runtime_error("deterministic_parallel: ray moved " + std::to_string(dd) +
                               " cells from the cell being processed (limit " +
                               std::to_string(dmax) + ")");
  }
}

inline double sum_ngam_r4(const Ray& r) {
  // sum(ray(ir)%ngam): real(4) sum, then widened
  return rrk::sum_ngam(r);
}

// 1-based cell of a float position (the ray's stored x)
inline void home_cell(const Ray& ry, int Nrt, int& i, int& j, int& k) {
  i = cell_of(ry.x[0], Nrt);
  j = cell_of(ry.x[1], Nrt);
  k = cell_of(ry.x[2], Nrt);
}

}  // namespace

rrk::RTP rt_params(const Simulation& s) {
  rrk::RTP p;
  p.N = s.cfg.Nrt;
  p.nrtr = s.cfg.nrtr;
  p.hrtr = s.cfg.hrtr;
  p.seg_min = s.units.seg_min;
  p.ngam_min = s.units.ngam_min;
  p.fgam_min = s.units.fgam_min;
  p.sfr_min = s.units.sfr_min;
  p.vol_rt = s.units.vol_rt;
  p.nH_avg = s.units.nH_avg;
  p.dtsec_rt = s.clk.dtsec_rt;
  // B8 fix: without legacy_rt, up to 4 segments (the original uses 3); 3k+1
  // for a path of k = rt_cells_per_substep cells (legacy_rt = 0 only)
  p.cells = s.cfg.rt_cells_per_substep;
  p.nseg = rrk::nseg_for(s.cfg.legacy_rt, p.cells);
  p.deposit = s.cfg.ray_deposit;  // share (0, the original) or exposure (1)
  p.nside_merge = s.rays.nside_merge;
  p.ngam_ray16 = s.rays.ngam_ray[16 - 1];
  p.ngam_ray32 = s.rays.ngam_ray[32 - 1];
  // B16 fix ray_rotate: the keys of the frame rotation (rrk::frame_rotation)
  p.ray_rotate = s.cfg.ray_rotate;
  p.rot_seed = s.cfg.ray_rotate_seed;
  static_assert(sizeof(p.rot_abits) == sizeof(s.clk.scalefactor), "double");
  std::memcpy(&p.rot_abits, &s.clk.scalefactor, sizeof(p.rot_abits));
  p.rot_istep = s.clk.istep_rt;
  return p;
}

// ---------------------------------------------------------------------------
// ray_tracing_1
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_tracing_1
void ray_tracing_1(Simulation& s, const std::array<int, 6>& box) {
  RTGrid& G = s.rt;
  RayState& R = s.rays;
  RayStat& rys = R.stat;
  const Units& u = s.units;

  // merge1(0:191): Fortran capacity 10000 rays per cone
  constexpr int kMergeCap = 10000;
  static thread_local std::vector<std::vector<std::int64_t>> merge1_ir(192);
  int merge1_nr[192];
  const bool det = s.cfg.det_parallel;
  const rrk::RTP P = rt_params(s);
  // recomb_rad = 1 (recomb_rad.hpp): cells that launch recombination rays
  const std::uint8_t* rec_emit = s.recomb ? s.recomb->emit.data() : nullptr;

  for (int krt = box[4]; krt <= box[5]; ++krt) {
    for (int jrt = box[2]; jrt <= box[3]; ++jrt) {
      for (int irt = box[0]; irt <= box[1]; ++irt) {
        const std::int64_t c = G.cell(irt, jrt, krt);
        if (!(G.htoc_ray(1, c) > 0 ||
              static_cast<double>(G.sfr[c]) > u.sfr_min || (rec_emit && rec_emit[c])))
          continue;
        if (det) det_set_cell(irt, jrt, krt, P.cells);

        std::int64_t ncst_tot = 0, nspl_tot = 0, nmrg_tot = 0;
        // 2 in the original (both branches of its choice give 2, B12)
        const int nside_merge = P.nside_merge;

        // Source?
        if (static_cast<double>(G.sfr[c]) > u.sfr_min) {
          ray_source(s, irt, jrt, krt);
        }
        // Recombination radiation (recomb_rad = 1)
        if (rec_emit && rec_emit[c]) ray_source_recomb(s, irt, jrt, krt);

        // Use htoc_ray(1:2) instead of htoc_ray(3:4) for head & tail
        const std::int64_t head = G.htoc_ray(1, c);
        const std::int64_t tail = G.htoc_ray(2, c);
        std::int64_t ir = head;

        bool loop;
        int npix_merge;
        if (ir > 0) {
          loop = true;
          npix_merge = 12 * nside_merge * nside_merge;
          for (int ic = 0; ic < npix_merge; ++ic) merge1_nr[ic] = 0;
        } else {
          loop = false;
          npix_merge = 0;
        }

        while (loop) {
          if (ir <= 0)
            throw std::runtime_error("ray_tracing_1: broken ray chain");
          // Exit loop after reaching original tail ray
          if (ir == tail) loop = false;

          const std::int64_t next = R.at(ir).next;

          const Ray& ry = R.at(ir);
          // classify: cast / split / merge / remove (last match wins)
          switch (rrk::classify(ry, P)) {
            case rrk::kCast:
              ray_cast_1(s, ir, irt, jrt, krt);
              ncst_tot += 1;
              break;
            case rrk::kSplit:
              ray_cast_split(s, ir, irt, jrt, krt);
              ncst_tot += 4;
              nspl_tot += 4;
              break;
            case rrk::kMerge: {
              const int ic = rrk::merge_cone(ry, nside_merge);
              const int nr = merge1_nr[ic] + 1;
              if (nr > kMergeCap)
                throw std::runtime_error(
                    "ray_tracing_1: merge1 overflow (Fortran bound 10000)");
              merge1_nr[ic] = nr;
              if (static_cast<int>(merge1_ir[ic].size()) < nr)
                merge1_ir[ic].resize(nr);
              merge1_ir[ic][nr - 1] = ir;
              break;
            }
            default:  // rrk::kRemove
              ray_remove(s, ir, irt, jrt, krt, true);
              break;
          }

          ir = next;
        }

        // Ray merging?
        int nsum = 0;
        for (int ic = 0; ic < npix_merge; ++ic) nsum += merge1_nr[ic];
        if (nsum > 0) {
          for (int ic = 0; ic < npix_merge; ++ic) {
            const int nr = merge1_nr[ic];
            if (nr == 1) {
              ray_cast_1(s, merge1_ir[ic][0], irt, jrt, krt);
              ncst_tot += 1;
            } else if (nr > 1) {
              ray_cast_merge(s, nr, merge1_ir[ic].data(), irt, jrt, krt);
              ncst_tot += 1;
              nmrg_tot += (nr - 1);
            }
          }
        }

        // Update ray stat (deterministic_parallel: the thread's counters)
        if (det) {
          ThreadRays& t = R.tl[par::thread_id()];
          t.ncst1 += ncst_tot;
          t.ncst2 = f::max(t.ncst2, ncst_tot);
          t.nspl1 += nspl_tot;
          t.nspl2 = f::max(t.nspl2, nspl_tot);
          t.nmrg1 += nmrg_tot;
          t.nmrg2 = f::max(t.nmrg2, nmrg_tot);
        } else {
          rys.ncst1 = rys.ncst1 + ncst_tot;
          rys.ncst2 = f::max(rys.ncst2, ncst_tot);
          rys.nspl1 = rys.nspl1 + nspl_tot;
          rys.nspl2 = f::max(rys.nspl2, nspl_tot);
          rys.nmrg1 = rys.nmrg1 + nmrg_tot;
          rys.nmrg2 = f::max(rys.nmrg2, nmrg_tot);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ray_source
// ---------------------------------------------------------------------------

// B16 fix: the frame rotation of ray_source (ray_tracing.hpp); the
// arithmetic is rrk::frame_rotation, shared with the AMReX backend
void ray_frame_rotation(const Simulation& s, int irt, int jrt, int krt, double rot[3][3],
                        std::uint64_t stream) {
  rrk::frame_rotation(rt_params(s), irt, jrt, krt, stream, rot);
}

// ray_tracing_tools.f90: ray_source
void ray_source(Simulation& s, int irt, int jrt, int krt) {
  RTGrid& G = s.rt;
  RayState& R = s.rays;
  const std::int64_t c = G.cell(irt, jrt, krt);
  const rrk::RTP P = rt_params(s);

  if (rrk::source_is_thick(G.tau_p(c))) {
    // Optically thick source cell: add photons directly to the cell
    double v[rrk::kNdep];
    rrk::source_thick_deposit(G.Luv[c], P, v);
    float* gng = G.ngam_p(c);
    float* gnu = G.nu_p(c);
    for (int n = 0; n < kNfreq; ++n) {
      gng[n] = static_cast<float>(static_cast<double>(gng[n]) + v[n]);
      gnu[n] = static_cast<float>(static_cast<double>(gnu[n]) + v[5 + n]);
    }
    if (s.tdiag && s.tdiag->acc) tensor_diag_iso(*s.tdiag, c, v);  // diagnostic
    return;
  }

  // Generate new rays: nray_init, more with rt_cells_per_substep > 4
  // (rrk::emit_nside)
  rrk::SourceSpec spec = rrk::source_spec(irt - 1, jrt - 1, krt - 1, G.Luv[c], P);
  const int nray = rrk::nray_from_nside(spec.nside);
  static thread_local std::vector<std::int64_t> indx_buf;
  if (static_cast<int>(indx_buf.size()) < nray) indx_buf.resize(nray);
  std::int64_t* indx = indx_buf.data();
  ray_create(s, nray, indx);
  // B16 fix (ray_rotate != 0): rotated HEALPix frame; the kernel's source_ray
  // rotates the direction and makes ipix the grid-frame pixel containing it,
  // so that ray_cast_split keeps the deviation dx - centre(ipix), as for a
  // merged ray
  const bool rotate = s.cfg.ray_rotate != 0;
  double rot[3][3];
  if (rotate) ray_frame_rotation(s, irt, jrt, krt, rot);
  rrk::set_rotation(spec, rotate ? rot : nullptr);

  for (int i = 1; i <= nray; ++i) {
    const int ipix = i - 1;
    const std::int64_t ir = indx[i - 1];
    Ray& ry = R.at(ir);
    static_cast<rrk::RayP&>(ry) = rrk::source_ray(spec, ipix);
    if (R.track_prov) {
      RayProv& pv = R.prov[ir];
      pv = RayProv{};
      pv.src = static_cast<std::int32_t>(c);
      pv.birth = static_cast<std::int32_t>(R.prov_substep);
    }
    // Insert rays into head of chain
    ray_insert_head(s, ir, irt, jrt, krt);
  }

  // Update head in htoc_ray(1) since source rays added at head
  G.htoc_ray(1, c) = indx[nray - 1];
  // Only update tail in htoc_ray(2) if initially empty
  if (G.htoc_ray(2, c) == 0) G.htoc_ray(2, c) = indx[0];
}

// ---------------------------------------------------------------------------
// ray_cast_split
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_cast_split
void ray_cast_split(Simulation& s, std::int64_t ir, int irt, int jrt, int krt) {
  RayState& R = s.rays;
  const int Nrt = s.cfg.Nrt;
  const rrk::RTP P = rt_params(s);

  // The 4 daughters (rrk::split); the first keeps the parent's slot
  const rrk::RayP parent = R.at(ir);
  rrk::RayP child[4];
  bool clamped = false;
  rrk::split(parent, P, child, &clamped);

  // Ray-dump provenance (not used by the RT)
  RayProv pv_child;
  if (R.track_prov) {
    pv_child = R.prov[ir];
    if (clamped) {
      pv_child.flags |= RayProv::kClamped;
      if (pv_child.nclamp < 255) ++pv_child.nclamp;
    }
  }

  // Storage indices: parent + 3 new rays
  std::int64_t indx[4];
  indx[0] = ir;
  ray_create(s, 3, indx + 1);

  for (int n = 1; n <= 4; ++n) {
    const std::int64_t irn = indx[n - 1];
    if (R.track_prov) R.prov[irn] = pv_child;
    Ray& ry = R.at(irn);
    static_cast<rrk::RayP&>(ry) = child[n - 1];
    // Safer to get x from ray(ir)%x
    int i, j, k;
    home_cell(ry, Nrt, i, j, k);

    if (n == 1) {
      // First daughter keeps the parent's storage index; relink only if it
      // moved into another cell.
      if (i != irt || j != jrt || k != krt) {
        ray_remove(s, irn, irt, jrt, krt, false);
        ray_insert_tail(s, irn, i, j, k);
      }
    } else {
      ray_insert_tail(s, irn, i, j, k);
    }

    ray_cast_1(s, irn, i, j, k);
  }
}

// ---------------------------------------------------------------------------
// ray_cast_merge
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_cast_merge
void ray_cast_merge(Simulation& s, int nrc, const std::int64_t* irc, int irt,
                    int jrt, int krt) {
  RayState& R = s.rays;
  const int Nrt = s.cfg.Nrt;
  const rrk::RTP P = rt_params(s);

  // Ray-dump provenance of the merged ray (not used by the RT)
  RayProv pv_merged;
  if (R.track_prov) {
    pv_merged = R.prov[irc[0]];
    int nm = pv_merged.nmerge;
    for (int i = 1; i < nrc; ++i) {
      const RayProv& q = R.prov[irc[i]];
      if (q.src != pv_merged.src) pv_merged.src = -1;
      if (q.birth < pv_merged.birth) pv_merged.birth = q.birth;
      pv_merged.flags |= q.flags;
      if (q.nclamp > pv_merged.nclamp) pv_merged.nclamp = q.nclamp;
      nm = std::max(nm, static_cast<int>(q.nmerge));
    }
    pv_merged.nmerge = static_cast<std::uint16_t>(std::min(nm + 1, 65535));
    pv_merged.flags |= RayProv::kMerged;
    if (pv_merged.src < 0) pv_merged.flags |= RayProv::kMixedSrc;
  }

  // Photon-weighted sums over the rays of the cone, in chain order
  // (rrk::merge_add); the other rays are removed entirely
  float cc[3];
  rrk::cell_centre(irt - 1, jrt - 1, krt - 1, cc);
  rrk::MergeAcc acc;
  for (int i = 1; i <= nrc; ++i) {
    const std::int64_t ir = irc[i - 1];
    rrk::merge_add(acc, R.at(ir), cc, P);
    if (i > 1) ray_remove(s, ir, irt, jrt, krt, true);
  }

  // Merged ray properties (rrk::merge_finish)
  const std::int64_t ir = irc[0];
  if (R.track_prov) R.prov[ir] = pv_merged;
  Ray& ry = R.at(ir);
  static_cast<rrk::RayP&>(ry) = rrk::merge_finish(acc, cc, P);
  int i, j, k;
  home_cell(ry, Nrt, i, j, k);

  if (i != irt || j != jrt || k != krt) {
    ray_remove(s, ir, irt, jrt, krt, false);
    ray_insert_tail(s, ir, i, j, k);
  }

  ray_cast_1(s, ir, i, j, k);
}

// ---------------------------------------------------------------------------
// ray_tracing_3
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_tracing_3
void ray_tracing_3(Simulation& s, const std::array<int, 6>& box) {
  RTGrid& G = s.rt;
  RayState& R = s.rays;
  const Units& u = s.units;

  for (int krt = box[4]; krt <= box[5]; ++krt) {
    for (int jrt = box[2]; jrt <= box[3]; ++jrt) {
      for (int irt = box[0]; irt <= box[1]; ++irt) {
        const std::int64_t c = G.cell(irt, jrt, krt);
        if (!(G.htoc_ray(1, c) > 0 ||
              static_cast<double>(G.sfr[c]) > u.sfr_min))
          continue;
        if (s.cfg.det_parallel) det_set_cell(irt, jrt, krt, 0);

        const std::int64_t head = G.htoc_ray(1, c);
        const std::int64_t tail = G.htoc_ray(2, c);
        std::int64_t ir = head;
        bool loop = ir > 0;

        while (loop) {
          if (ir <= 0)
            throw std::runtime_error("ray_tracing_3: broken ray chain");
          if (ir == tail) loop = false;
          const std::int64_t next = R.at(ir).next;
          ray_cast_3(s, ir, irt, jrt, krt);
          ir = next;
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ray_cast_1 / ray_cast_3 (non-RAYFAST)
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_cast_1 (#else RAYFAST branch): rrk::cast_1
// with the deposits added to the float32 RT grid in segment order
void ray_cast_1(Simulation& s, std::int64_t ir, int /*irt*/, int /*jrt*/,
                int /*krt*/) {
  RTGrid& G = s.rt;
  const rrk::RTP P = rt_params(s);
  auto tau = [&](int i, int j, int k) -> const float* { return G.tau_p(G.cell(i + 1, j + 1, k + 1)); };
  TensorDiag* TD = (s.tdiag && s.tdiag->acc) ? s.tdiag.get() : nullptr;  // diagnostic
  auto dep = [&](int i, int j, int k, const double* v) {
    const std::int64_t c = G.cell(i + 1, j + 1, k + 1);
    float* gng = G.ngam_p(c);
    float* gnu = G.nu_p(c);
    float* gdx = G.dx_p(c);
    for (int q = 0; q < kNfreq; ++q)
      gng[q] = static_cast<float>(static_cast<double>(gng[q]) + v[q]);
    for (int q = 0; q < kNfreq; ++q)
      gnu[q] = static_cast<float>(static_cast<double>(gnu[q]) + v[5 + q]);
    for (int d = 0; d < 3; ++d)
      gdx[d] = static_cast<float>(static_cast<double>(gdx[d]) + v[10 + d]);
    if (TD) tensor_diag_ray(*TD, c, v);
  };
  rrk::cast_1(s.rays.at(ir), P, tau, dep);
}

// ray_tracing_tools.f90: ray_cast_3 (#else RAYFAST branch): rrk::cast_3,
// then the removal or the relink to the new home cell
void ray_cast_3(Simulation& s, std::int64_t& ir, int irt, int jrt, int krt) {
  RTGrid& G = s.rt;
  const int Nrt = s.cfg.Nrt;
  const rrk::RTP P = rt_params(s);
  auto tf = [&](int i, int j, int k, const float*& tau, const float*& fgam) {
    const std::int64_t c = G.cell(i + 1, j + 1, k + 1);
    tau = G.tau_p(c);
    fgam = G.fgam_p(c);
  };
  Ray& ry = s.rays.at(ir);
  if (!rrk::cast_3(ry, P, tf)) {
    ray_remove(s, ir, irt, jrt, krt, true);
    ir = 0;
  } else {
    int i, j, k;
    home_cell(ry, Nrt, i, j, k);
    if (i != irt || j != jrt || k != krt) {
      ray_remove(s, ir, irt, jrt, krt, false);
      ray_insert_tail(s, ir, i, j, k);
    }
  }
}

// ray_tracing_tools.f90: ray_segment (rrk::ray_segment)
std::array<double, 4> ray_segment(const double x0[3], const double dx[3],
                                  double seg_min, int nseg) {
  if (nseg < 1 || nseg > 4) throw std::invalid_argument("ray_segment: nseg must be 1..4");
  std::array<double, 4> seg{};
  rrk::ray_segment(x0, dx, seg_min, nseg, seg.data());
  return seg;
}

std::vector<double> ray_segment_path(const double x0[3], const double dx[3], double seg_min,
                                     int cells) {
  if (cells < 1 || cells > rrk::kMaxCellsPerSubstep)
    throw std::invalid_argument("ray_segment_path: cells must be 1.." +
                                std::to_string(rrk::kMaxCellsPerSubstep));
  std::vector<double> seg(static_cast<std::size_t>(rrk::nseg_for(false, cells)), 0.0);
  rrk::ray_segment(x0, dx, seg_min, static_cast<int>(seg.size()), seg.data(), cells);
  return seg;
}

// ---------------------------------------------------------------------------
// Storage / chain operations
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_create
void ray_create(Simulation& s, int nr, std::int64_t* indx) {
  RayState& R = s.rays;
  if (s.cfg.det_parallel) {
    R.det_take(R.tl[par::thread_id()], nr, indx);
    return;
  }
  RayStat& st = R.stat;

  if (st.nfree < nr) {
    const std::int64_t n = R.Nrbuf / 4;
    const std::int64_t if1 = st.nfree + 1;
    const std::int64_t if2 = st.nfree + n;
    if (if2 > R.Nfree_max)
      throw std::runtime_error("ray_create: free buffer overflow (Fortran bound)");
    st.nfree = if2;
    const std::int64_t ir1 = R.Nray + 1;
    const std::int64_t ir2 = R.Nray + n;
    R.Nray = ir2;
    for (std::int64_t m = 0; m < n; ++m) st.ifree[if1 + m] = R.ray_index(ir1 + m);
  }
  if (st.nfree < nr)
    throw std::runtime_error("ray_create: not enough free rays (Nrbuf too small)");

  // indx = ifree(nfree-nr+1:nfree)
  for (int m = 0; m < nr; ++m) indx[m] = st.ifree[st.nfree - nr + 1 + m];
  st.nfree = st.nfree - nr;
  for (int m = 0; m < nr; ++m) R.ensure(indx[m]);
}

// ray_tracing_tools.f90: ray_remove
void ray_remove(Simulation& s, std::int64_t ir, int irt, int jrt, int krt,
                bool flag) {
  RayState& R = s.rays;
  RTGrid& G = s.rt;
  RayStat& st = R.stat;
  const std::int64_t c = G.cell(irt, jrt, krt);

  // Remove link and reconnect chain
  const std::int64_t id = R.at(ir).next;
  const std::int64_t iu = R.at(ir).prev;
  if (id != 0) R.at(id).prev = iu;
  if (iu != 0) R.at(iu).next = id;

  // Head or tail of chain?
  if (ir == G.htoc_ray(3, c)) G.htoc_ray(3, c) = id;
  if (ir == G.htoc_ray(4, c)) G.htoc_ray(4, c) = iu;

  if (flag && s.cfg.det_parallel) {
    // freed slots are always cleared (legacy_rt = 0 is required)
    R.at(ir) = Ray{};
    R.det_give(R.tl[par::thread_id()], ir);
  } else if (flag) {
    // Free buffer full: return a quarter of Nrbuf to ray_index
    if (st.nfree == R.Nfree_max) {
      const std::int64_t n = R.Nrbuf / 4;
      const std::int64_t if2 = st.nfree;
      const std::int64_t if1 = if2 - (n - 1);
      st.nfree = st.nfree - n;
      const std::int64_t ir2 = R.Nray;
      const std::int64_t ir1 = ir2 - (n - 1);
      R.Nray = R.Nray - n;
      for (std::int64_t m = 0; m < n; ++m) R.set_ray_index(ir1 + m, st.ifree[if1 + m]);
    }
    st.nfree = st.nfree + 1;
    st.ifree[st.nfree] = ir;
    st.ndel1 = st.ndel1 + 1;

    // B1 fix: without legacy_rt, clear the freed slot (as
    // ray_buffer_deallocate does), so calc_ray_1 (nside > 0) ignores it.
    if (!s.cfg.legacy_rt) R.at(ir) = Ray{};
  }
}

// ray_tracing_tools.f90: ray_insert_head
void ray_insert_head(Simulation& s, std::int64_t ir, int irt, int jrt, int krt) {
  RayState& R = s.rays;
  RTGrid& G = s.rt;
  const std::int64_t c = G.cell(irt, jrt, krt);
  const std::int64_t head = G.htoc_ray(3, c);
  if (head == 0) {
    R.at(ir).next = 0;
    R.at(ir).prev = 0;
    G.htoc_ray(3, c) = ir;
    G.htoc_ray(4, c) = ir;
  } else {
    R.at(head).prev = ir;
    R.at(ir).next = head;
    R.at(ir).prev = 0;
    G.htoc_ray(3, c) = ir;
  }
}

// ray_tracing_tools.f90: ray_insert_tail
void ray_insert_tail(Simulation& s, std::int64_t ir, int irt, int jrt, int krt) {
  RayState& R = s.rays;
  RTGrid& G = s.rt;
  if (s.cfg.det_parallel) det_check_reach(s, irt, jrt, krt);
  const std::int64_t c = G.cell(irt, jrt, krt);
  const std::int64_t tail = G.htoc_ray(4, c);
  if (tail == 0) {
    R.at(ir).next = 0;
    R.at(ir).prev = 0;
    G.htoc_ray(3, c) = ir;
    G.htoc_ray(4, c) = ir;
  } else {
    R.at(tail).next = ir;
    R.at(ir).next = 0;
    R.at(ir).prev = tail;
    G.htoc_ray(4, c) = ir;
  }
}

// ---------------------------------------------------------------------------
// nside / radius relations
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: lside_from_radius etc. (the kernels' functions)
int lside_from_radius(double r) { return rrk::lside_from_radius(r); }
int nside_from_radius(double r) { return rrk::nside_from_radius(r); }
int nray_from_nside(int nside) { return rrk::nray_from_nside(nside); }
double rmin_from_nside(int nside) { return rrk::rmin_from_nside(nside); }
double rmax_from_nside(int nside) { return rrk::rmax_from_nside(nside); }

// ---------------------------------------------------------------------------
// interpolate (general_tools.f90), 'lin' branch
// ---------------------------------------------------------------------------

double interpolate_lin(double x, const double* fa, int n, int a, int b) {
  return rrk::interp_lin(x, fa, n, a, b);
}

// ---------------------------------------------------------------------------
// ray_redshifting (non-FREQ1)
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: the tables of ray_redshifting (Nfreq = 5 case)
RedshiftTable make_redshift_table() { return rrk::build_redshift_table(); }

// ray_tracing_tools.f90: ray_redshifting / redshift_ray (#else FREQ1)
void ray_redshifting(Simulation& s) {
  RayState& R = s.rays;
  if (R.Nray <= 0) return;
  const double ngam_min = s.units.ngam_min;
  // (a - da)/a, or per distance with redshift_per_distance (cosmology.hpp)
  const double fscale = redshift_fscale(s.cfg, s.clk);
  const RedshiftTable T = make_redshift_table();

  // every ray is redshifted independently
  const std::int64_t top = R.hwm();
#pragma omp parallel for schedule(static) if (par::worth(top))
  for (std::int64_t ir = 1; ir <= top; ++ir) {
    Ray& ry = R.at(ir);
    if (!(ry.nside > 0)) continue;
    rrk::redshift_ray(ry, T, fscale, ngam_min);
  }
}

// ---------------------------------------------------------------------------
// Ray buffer management (Ncpu = 1)
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_buffer_allocate
void ray_buffer_allocate(Simulation& s) {
  RayState& R = s.rays;
  RayStat& st = R.stat;
  if (st.ifree.empty()) st.ifree.assign(static_cast<std::size_t>(R.Nfree_max) + 1, 0);
  const std::int64_t i1 = R.Nray + 1;  // + (icpu-1)*nrbuf, icpu = 1
  st.nfree = R.Nrbuf;
  for (std::int64_t m = 0; m < R.Nrbuf; ++m) st.ifree[1 + m] = R.ray_index(i1 + m);
  R.Nray = R.Nray + R.Nrbuf * s.cfg.Ncpu;
  R.save.nfree = st.nfree;
}

// ray_tracing_tools.f90: ray_buffer_supplement
void ray_buffer_supplement(Simulation& s) {
  RayState& R = s.rays;
  RayStat& st = R.stat;
  const int Ncpu = s.cfg.Ncpu;

  // int(0.98*Nray_max,kind=8): real(4) product
  const std::int64_t nlim = f::max<std::int64_t>(
      f::int8(static_cast<double>(0.98f * static_cast<float>(R.Nray_max))) -
          (R.Nray - st.nfree),
      0);
  const std::int64_t navg = f::min(R.Nrbuf, nlim / Ncpu);
  const std::int64_t ndel = navg / 4;

  // First reduce buffer where needed
  {
    const std::int64_t n = f::max<std::int64_t>(st.nfree - navg, 0);
    const std::int64_t nr = n;  // nr(icpu); maxval(nr) == n with one thread
    if (n > 0 && nr >= ndel) {
      const std::int64_t i1 = R.Nray + 1 - nr;  // sum(nr(icpu:Ncpu))
      for (std::int64_t m = 0; m < n; ++m)
        R.set_ray_index(i1 + m, st.ifree[st.nfree - n + 1 + m]);
      st.nfree = st.nfree - n;
    }
    if (nr >= ndel) R.Nray = R.Nray - nr;
  }

  // Second increase buffer where needed
  {
    const std::int64_t n = f::min(f::max<std::int64_t>(navg - st.nfree, 0), ndel);
    const std::int64_t nr = n;
    if (n > 0 && nr >= ndel) {
      const std::int64_t i2 = R.Nray + nr;  // sum(nr(1:icpu))
      const std::int64_t i1 = i2 - n + 1;
      if (st.nfree + n > R.Nfree_max)
        throw std::runtime_error("ray_buffer_supplement: free buffer overflow");
      for (std::int64_t m = 0; m < n; ++m)
        st.ifree[st.nfree + 1 + m] = R.ray_index(i1 + m);
      st.nfree = st.nfree + n;
    }
    if (nr >= ndel) R.Nray = R.Nray + nr;
  }

  R.save.nfree = st.nfree;
}

// ray_tracing_tools.f90: ray_buffer_deallocate
void ray_buffer_deallocate(Simulation& s) {
  RayState& R = s.rays;
  RayStat& st = R.stat;

  if (st.nfree > 0) {
    const std::int64_t i1 = R.Nray + 1 - st.nfree;  // sum(ray_stat(icpu:Ncpu)%nfree)
    for (std::int64_t m = 0; m < st.nfree; ++m) R.set_ray_index(i1 + m, st.ifree[1 + m]);
  }
  R.Nray = R.Nray - st.nfree;
  st.nfree = 0;
  R.save.nfree = 0;

  // Zero inactive ray storage: ray_index(i) for i = max(1,Nray+1)..Nray_max.
  // Positions beyond the materialised range hold untouched (zero) rays.
  const std::int64_t i1 = f::max<std::int64_t>(1, R.Nray + 1);
  const std::int64_t i2 = static_cast<std::int64_t>(R.ray_index_.size()) - 1;
  for (std::int64_t i = i1; i <= i2; ++i) {
    const std::int64_t ir = R.ray_index_[i];
    R.ray[ir] = Ray{};
  }
}

}  // namespace radray
