// Port of ray_tracing_tools.f90 (default code path: RAYFAST and FREQ1
// undefined). See include/radray/ray_tracing.hpp for conventions.
#include "radray/ray_tracing.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "radray/constants.hpp"
#include "radray/fortran.hpp"
#include "radray/healpix.hpp"

namespace radray {

using namespace rayp;

// ---------------------------------------------------------------------------
// RayState (storage)
// ---------------------------------------------------------------------------

void RayState::init(std::int64_t nray_max, std::int64_t nfree_max) {
  Nray_max = nray_max;
  Nfree_max = nfree_max;
  Nray = 0;
  Nrbuf = 0;
  ray.clear();
  ray_index_.clear();
  ray.reserve(static_cast<std::size_t>(Nray_max) + 1);
  ray_index_.reserve(static_cast<std::size_t>(Nray_max) + 1);
  ray.resize(1);         // slot 0: dummy
  ray_index_.resize(1);  // position 0 unused
  prov.clear();
  if (track_prov) prov.resize(1);
  ray_index_[0] = 0;
  stat = RayStat{};
  stat.ifree.assign(static_cast<std::size_t>(Nfree_max) + 1, 0);
  save = RayStat{};
  ngam_ray.fill(0.0);
}

void RayState::ensure(std::int64_t i) {
  if (i < 1 || i > Nray_max)
    throw std::runtime_error("ray storage index out of range: " + std::to_string(i));
  const std::int64_t old = static_cast<std::int64_t>(ray_index_.size());
  if (i >= old) {
    ray_index_.resize(static_cast<std::size_t>(i) + 1);
    for (std::int64_t p = old; p <= i; ++p) ray_index_[p] = p;
  }
  if (i >= static_cast<std::int64_t>(ray.size()))
    ray.resize(static_cast<std::size_t>(i) + 1);
  if (track_prov && i >= static_cast<std::int64_t>(prov.size()))
    prov.resize(static_cast<std::size_t>(i) + 1);
}

// ---------------------------------------------------------------------------
// init_ray_tracing
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: init_ray_tracing
// (The Fortran also writes Log/nside_of_r.dat and Log/r_of_nside.dat; these
// diagnostics are not reproduced.)
void init_ray_tracing(Simulation& s) {
  s.rays.init(s.cfg.Nray_max, s.cfg.Nfree_max);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

namespace {

// 1 + mod(int(x),Nrt): RT cell index (1-based) of a position in RT units.
inline int cell_of(double x, int Nrt) { return 1 + f::mod(f::int4(x), Nrt); }

inline double sum_ngam_r4(const Ray& r) {
  // sum(ray(ir)%ngam): real(4) sum, then widened
  return static_cast<double>(f::sum4<kNfreq>(r.ngam));
}

}  // namespace

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

  for (int krt = box[4]; krt <= box[5]; ++krt) {
    for (int jrt = box[2]; jrt <= box[3]; ++jrt) {
      for (int irt = box[0]; irt <= box[1]; ++irt) {
        const std::int64_t c = G.cell(irt, jrt, krt);
        if (!(G.htoc_ray(1, c) > 0 ||
              static_cast<double>(G.sfr[c]) > u.sfr_min))
          continue;

        std::int64_t ncst_tot = 0, nspl_tot = 0, nmrg_tot = 0;
        const int nside_merge = 2;

        // Source?
        if (static_cast<double>(G.sfr[c]) > u.sfr_min) {
          ray_source(s, irt, jrt, krt);
        }

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

          enum Case { kMerge, kCast, kSplit, kRemove };
          Case ray_case = kMerge;
          const Ray& ry = R.at(ir);
          const double ntot = sum_ngam_r4(ry);

          if (ntot > R.ngam_ray[16 - 1]) ray_case = kCast;
          if (ry.r == 0) ray_case = kCast;
          if (ry.dr < 0 && ntot > R.ngam_ray[32 - 1]) ray_case = kSplit;
          if (ntot < u.ngam_min) ray_case = kRemove;
          if (ry.nside == nside_max) ray_case = kRemove;

          switch (ray_case) {
            case kCast:
              ray_cast_1(s, ir, irt, jrt, krt);
              ncst_tot += 1;
              break;
            case kSplit:
              ray_cast_split(s, ir, irt, jrt, krt);
              ncst_tot += 4;
              nspl_tot += 4;
              break;
            case kMerge: {
              double dxv[3] = {static_cast<double>(ry.dx[0]),
                               static_cast<double>(ry.dx[1]),
                               static_cast<double>(ry.dx[2])};
              const int ic = healpix::vec2pix_nest(nside_merge, dxv);
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
            case kRemove:
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

        // Update ray stat
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

// ---------------------------------------------------------------------------
// ray_source
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_source
void ray_source(Simulation& s, int irt, int jrt, int krt) {
  RTGrid& G = s.rt;
  RayState& R = s.rays;
  const Units& u = s.units;
  const std::int64_t c = G.cell(irt, jrt, krt);
  float* tau = G.tau_p(c);

  // minval(rtg%tau): real(4)
  float tmin = tau[0];
  for (int n = 1; n < kNfreq; ++n)
    if (tau[n] < tmin) tmin = tau[n];

  // Photons per frequency group, normalised by <nH>:
  // ngam = fResc*pop2_init*(rtg%Luv*dtsec_rt)
  const double Ldt = static_cast<double>(G.Luv[c]) * s.clk.dtsec_rt;

  // exp of a real(4) argument is real(4) (expf); 1E-2 is real(4)
  if (std::exp(-tmin) < 1e-2f) {
    // Optically thick source cell: add photons directly to the cell
    float* gng = G.ngam_p(c);
    float* gnu = G.nu_p(c);
    for (int n = 0; n < kNfreq; ++n) {
      const double nu = cp::freq_init[n];
      double ngam = cp::fResc * cp::pop2_init[n] * Ldt;
      ngam = ngam / u.vol_rt / u.nH_avg;
      gng[n] = static_cast<float>(static_cast<double>(gng[n]) + ngam);
      gnu[n] = static_cast<float>(static_cast<double>(gnu[n]) + ngam * nu);
    }
    return;
  }

  // Generate new rays
  const double r = 0;
  const double dr = rmax_from_nside(nside_init);
  // x = (/irt,jrt,krt/) - 0.5 : real(4) arithmetic
  const float xs[3] = {static_cast<float>(irt) - 0.5f,
                       static_cast<float>(jrt) - 0.5f,
                       static_cast<float>(krt) - 0.5f};

  std::int64_t indx[nray_init];
  ray_create(s, nray_init, indx);

  double nu[kNfreq], ngam[kNfreq];
  for (int n = 0; n < kNfreq; ++n) {
    nu[n] = cp::freq_init[n];
    ngam[n] = cp::fResc * cp::pop2_init[n] * Ldt;
    ngam[n] = ngam[n] / u.vol_rt / u.nH_avg / nray_init;
  }

  for (int i = 1; i <= nray_init; ++i) {
    const int ipix = i - 1;
    double dx[3];
    healpix::pix2vec_nest(nside_init, ipix, dx);

    const std::int64_t ir = indx[i - 1];
    Ray& ry = R.at(ir);
    ry.nside = nside_init;
    ry.ipix = ipix;
    ry.r = static_cast<float>(r);
    ry.dr = static_cast<float>(dr);
    for (int d = 0; d < 3; ++d) {
      ry.x[d] = xs[d];
      ry.dx[d] = static_cast<float>(dx[d]);
    }
    for (int n = 0; n < kNfreq; ++n) {
      ry.nu[n] = static_cast<float>(nu[n]);
      ry.ngam[n] = static_cast<float>(ngam[n]);
    }
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
  G.htoc_ray(1, c) = indx[nray_init - 1];
  // Only update tail in htoc_ray(2) if initially empty
  if (G.htoc_ray(2, c) == 0) G.htoc_ray(2, c) = indx[0];
}

// ---------------------------------------------------------------------------
// ray_cast_split
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_cast_split
void ray_cast_split(Simulation& s, std::int64_t ir, int irt, int jrt, int krt) {
  RayState& R = s.rays;
  const double nrtr = s.cfg.nrtr;
  const int Nrt = s.cfg.Nrt;

  // Parent ray
  int nside, ipix;
  double r, dr, x[3], dx[3], nu[kNfreq], ngam[kNfreq];
  {
    const Ray& p = R.at(ir);
    nside = p.nside;
    ipix = p.ipix;
    r = p.r;
    dr = p.dr;
    for (int d = 0; d < 3; ++d) {
      x[d] = p.x[d];
      dx[d] = p.dx[d];
    }
    for (int n = 0; n < kNfreq; ++n) {
      nu[n] = p.nu[n];
      ngam[n] = p.ngam[n];
    }
  }

  // Current indices: deviation of the ray direction from its pixel centre
  const int ipixc = ipix;
  double dc[3];
  healpix::pix2vec_nest(nside, ipixc, dc);
  for (int d = 0; d < 3; ++d) dc[d] = dx[d] - dc[d];

  // Daughter rays
  nside = f::min(2 * nside, nside_max);
  const double rmin = rmin_from_nside(nside);
  const double rmax = rmax_from_nside(nside);
  const double r_parent = r;
  r = f::min(f::max(r, rmin), rmin + 1);
  dr = rmax - r;

  // Ray-dump provenance (not used by the RT)
  RayProv pv_child;
  if (R.track_prov) {
    pv_child = R.prov[ir];
    if (r != r_parent) {
      pv_child.flags |= RayProv::kClamped;
      if (pv_child.nclamp < 255) ++pv_child.nclamp;
    }
  }

  // Effective origin, position corresponding to r = 0
  double x0[3];
  for (int d = 0; d < 3; ++d)
    x0[d] = f::mod(x[d] - f::mod(r * dx[d], nrtr) + nrtr, nrtr);

  // Divide photon count
  for (int n = 0; n < kNfreq; ++n) ngam[n] = ngam[n] / 4;

  // Storage indices: parent + 3 new rays
  std::int64_t indx[4];
  indx[0] = ir;
  ray_create(s, 3, indx + 1);

  for (int n = 1; n <= 4; ++n) {
    ipix = 4 * ipixc + (n - 1);

    // New direction: default direction of new pixel + parent deviation
    healpix::pix2vec_nest(nside, ipix, dx);
    for (int d = 0; d < 3; ++d) dx[d] = dx[d] + dc[d];
    const double nrm = std::sqrt(f::sum8<3>(
        std::array<double, 3>{dx[0] * dx[0], dx[1] * dx[1], dx[2] * dx[2]}.data()));
    for (int d = 0; d < 3; ++d) dx[d] = dx[d] / nrm;
    ipix = healpix::vec2pix_nest(nside, dx);

    // New position
    for (int d = 0; d < 3; ++d)
      x[d] = f::mod(x0[d] + f::mod(r * dx[d], nrtr) + nrtr, nrtr);

    const std::int64_t irn = indx[n - 1];
    if (R.track_prov) R.prov[irn] = pv_child;
    {
      Ray& ry = R.at(irn);
      ry.nside = nside;
      ry.ipix = ipix;
      ry.r = static_cast<float>(r);
      ry.dr = static_cast<float>(dr);
      for (int d = 0; d < 3; ++d) {
        ry.x[d] = static_cast<float>(x[d]);
        ry.dx[d] = static_cast<float>(dx[d]);
      }
      for (int k = 0; k < kNfreq; ++k) {
        ry.nu[k] = static_cast<float>(nu[k]);
        ry.ngam[k] = static_cast<float>(ngam[k]);
      }
      // Safer to get x from ray(ir)%x
      for (int d = 0; d < 3; ++d) x[d] = ry.x[d];
    }
    const int i = cell_of(x[0], Nrt);
    const int j = cell_of(x[1], Nrt);
    const int k = cell_of(x[2], Nrt);

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
  const double nrtr = s.cfg.nrtr;
  const double hrt = s.cfg.hrtr;
  const int Nrt = s.cfg.Nrt;
  const Units& u = s.units;

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

  double r = 0, dr = 0, ntot = 0;
  double x[3] = {0, 0, 0}, dx[3] = {0, 0, 0}, xd[3];
  double nu[kNfreq] = {0, 0, 0, 0, 0}, ngam[kNfreq] = {0, 0, 0, 0, 0};
  // (/irt,jrt,krt/)-0.5 is real(4)
  const float cc[3] = {static_cast<float>(irt) - 0.5f,
                       static_cast<float>(jrt) - 0.5f,
                       static_cast<float>(krt) - 0.5f};

  for (int i = 1; i <= nrc; ++i) {
    const std::int64_t ir = irc[i - 1];
    {
      const Ray& ry = R.at(ir);
      // xd = ray(ir)%x - ((/irt,jrt,krt/)-0.5): real(4) subtraction
      for (int d = 0; d < 3; ++d) xd[d] = static_cast<double>(ry.x[d] - cc[d]);
      for (int d = 0; d < 3; ++d) {
        if (std::fabs(xd[d]) > hrt) {
          if (xd[d] > 0)
            xd[d] = xd[d] - nrtr;
          else
            xd[d] = xd[d] + nrtr;
        }
      }
      const double w = sum_ngam_r4(ry);
      ntot = ntot + w;
      for (int d = 0; d < 3; ++d) x[d] = x[d] + w * xd[d];
      // ray(ir)%r**2 is a real(4) product
      r = r + w * static_cast<double>(ry.r * ry.r);
      dr = dr + w * static_cast<double>(ry.dr);
      for (int d = 0; d < 3; ++d) dx[d] = dx[d] + w * static_cast<double>(ry.dx[d]);
      // ray(ir)%ngam*ray(ir)%nu is a real(4) product
      for (int n = 0; n < kNfreq; ++n)
        nu[n] = nu[n] + static_cast<double>(ry.ngam[n] * ry.nu[n]);
      for (int n = 0; n < kNfreq; ++n)
        ngam[n] = ngam[n] + static_cast<double>(ry.ngam[n]);
    }
    // Remove other rays entirely
    if (i > 1) ray_remove(s, ir, irt, jrt, krt, true);
  }

  // Merged ray properties
  const std::int64_t ir = irc[0];
  r = r / ntot;
  dr = dr / ntot;
  for (int d = 0; d < 3; ++d) xd[d] = x[d] / ntot;
  for (int d = 0; d < 3; ++d) dx[d] = dx[d] / ntot;
  r = std::sqrt(r);
  {
    const double nrm = std::sqrt((0.0 + dx[0] * dx[0]) + dx[1] * dx[1] + dx[2] * dx[2]);
    for (int d = 0; d < 3; ++d) dx[d] = dx[d] / nrm;
  }
  for (int d = 0; d < 3; ++d)
    x[d] = f::mod(static_cast<double>(cc[d]) + xd[d] + nrtr, nrtr);

  // Conserve photon number; conserve energy via effective frequency
  for (int n = 0; n < kNfreq; ++n) {
    if (ngam[n] < u.ngam_min) {
      nu[n] = cp::freq_init[n];
      ngam[n] = 0;
    } else {
      nu[n] = nu[n] / ngam[n];
    }
  }

  // nside from mean radius, then distance from splitting, final radius
  const int nside = nside_from_radius(r);
  const double rmin = rmin_from_nside(nside);
  const double rmax = rmax_from_nside(nside);
  r = f::max(rmax - dr, rmin);
  dr = rmax - r;

  const int ipix = healpix::vec2pix_nest(nside, dx);

  if (R.track_prov) R.prov[ir] = pv_merged;
  {
    Ray& ry = R.at(ir);
    ry.nside = nside;
    ry.ipix = ipix;
    ry.r = static_cast<float>(r);
    ry.dr = static_cast<float>(dr);
    for (int d = 0; d < 3; ++d) {
      ry.x[d] = static_cast<float>(x[d]);
      ry.dx[d] = static_cast<float>(dx[d]);
    }
    for (int n = 0; n < kNfreq; ++n) {
      ry.nu[n] = static_cast<float>(nu[n]);
      ry.ngam[n] = static_cast<float>(ngam[n]);
    }
    for (int d = 0; d < 3; ++d) x[d] = ry.x[d];
  }
  const int i = cell_of(x[0], Nrt);
  const int j = cell_of(x[1], Nrt);
  const int k = cell_of(x[2], Nrt);

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

// ray_tracing_tools.f90: ray_cast_1 (#else RAYFAST branch)
void ray_cast_1(Simulation& s, std::int64_t ir, int /*irt*/, int /*jrt*/,
                int /*krt*/) {
  RTGrid& G = s.rt;
  const double nrtr = s.cfg.nrtr;
  const int Nrt = s.cfg.Nrt;

  double x[3], dx[3], nu[kNfreq], ngam[kNfreq];
  {
    const Ray& ry = s.rays.at(ir);
    for (int d = 0; d < 3; ++d) {
      x[d] = ry.x[d];
      dx[d] = ry.dx[d];
    }
    for (int n = 0; n < kNfreq; ++n) {
      nu[n] = ry.nu[n];
      ngam[n] = ry.ngam[n];
    }
  }
  // B8 fix: without legacy_rt, up to 4 segments (the original uses 3)
  const int nseg = s.cfg.legacy_rt ? 3 : 4;
  const std::array<double, 4> dl = ray_segment(x, dx, s.units.seg_min, nseg);

  double l = 0;
  for (int n = 0; n < nseg; ++n) {
    l = l + dl[n];
    if (dl[n] >= s.units.seg_min) {
      double x1[3];
      for (int d = 0; d < 3; ++d)
        x1[d] = f::mod(x[d] + (l - dl[n] / 2) * dx[d] + nrtr, nrtr);
      const int i = cell_of(x1[0], Nrt);
      const int j = cell_of(x1[1], Nrt);
      const int k = cell_of(x1[2], Nrt);
      const std::int64_t c = G.cell(i, j, k);

      // fdl = dl(n)/sum(dl(n:))
      double sdl = 0;
      for (int m = n; m < nseg; ++m) sdl += dl[m];
      const double fdl = dl[n] / sdl;

      const float* tau = G.tau_p(c);
      double dng[kNfreq];
      for (int q = 0; q < kNfreq; ++q) {
        const double fabs_ = 1 - std::exp(static_cast<double>(-tau[q]) * dl[n]);
        if (fabs_ < fdl) {
          dng[q] = ngam[q] * fdl;
          ngam[q] = ngam[q] - dng[q];
        } else {
          dng[q] = ngam[q] * fabs_;
          ngam[q] = ngam[q] - dng[q];
        }
      }

      float* gng = G.ngam_p(c);
      float* gnu = G.nu_p(c);
      float* gdx = G.dx_p(c);
      for (int q = 0; q < kNfreq; ++q)
        gng[q] = static_cast<float>(static_cast<double>(gng[q]) + dng[q]);
      for (int q = 0; q < kNfreq; ++q)
        gnu[q] = static_cast<float>(static_cast<double>(gnu[q]) + dng[q] * nu[q]);
      const double sdng = f::sum8<kNfreq>(dng);
      for (int d = 0; d < 3; ++d)
        gdx[d] = static_cast<float>(static_cast<double>(gdx[d]) + dx[d] * sdng);
    }
  }
}

// ray_tracing_tools.f90: ray_cast_3 (#else RAYFAST branch)
void ray_cast_3(Simulation& s, std::int64_t& ir, int irt, int jrt, int krt) {
  RTGrid& G = s.rt;
  const double nrtr = s.cfg.nrtr;
  const int Nrt = s.cfg.Nrt;
  const Units& u = s.units;

  double r, dr, x[3], dx[3], ngam[kNfreq];
  {
    const Ray& ry = s.rays.at(ir);
    r = ry.r;
    dr = ry.dr;
    for (int d = 0; d < 3; ++d) {
      x[d] = ry.x[d];
      dx[d] = ry.dx[d];
    }
    for (int n = 0; n < kNfreq; ++n) ngam[n] = ry.ngam[n];
  }
  // B8 fix: without legacy_rt, up to 4 segments (the original uses 3)
  const int nseg = s.cfg.legacy_rt ? 3 : 4;
  const std::array<double, 4> dl = ray_segment(x, dx, u.seg_min, nseg);

  // Save original
  const double r0 = r, dr0 = dr;
  const double x0[3] = {x[0], x[1], x[2]};
  double ng1[kNfreq], ng2[kNfreq];
  for (int n = 0; n < kNfreq; ++n) {
    ng1[n] = ngam[n];
    ng2[n] = ngam[n];
  }
  const double ngam_cut = f::max(u.ngam_min, f::sum8<kNfreq>(ngam) * u.fgam_min);

  double l = 0;
  for (int n = 0; n < kNfreq; ++n) ngam[n] = 0;

  for (int n = 0; n < nseg; ++n) {
    l = l + dl[n];
    if (dl[n] >= u.seg_min && f::sum8<kNfreq>(ng2) > ngam_cut) {
      double x1[3];
      for (int d = 0; d < 3; ++d)
        x1[d] = f::mod(x0[d] + (l - dl[n] / 2) * dx[d] + nrtr, nrtr);
      const int i = cell_of(x1[0], Nrt);
      const int j = cell_of(x1[1], Nrt);
      const int k = cell_of(x1[2], Nrt);
      const std::int64_t c = G.cell(i, j, k);

      double sdl = 0;
      for (int m = n; m < nseg; ++m) sdl += dl[m];
      const double fdl = dl[n] / sdl;

      const float* tau = G.tau_p(c);
      const float* fgam = G.fgam_p(c);
      double dng[kNfreq];
      for (int q = 0; q < kNfreq; ++q) {
        const double fabs_ = 1 - std::exp(static_cast<double>(-tau[q]) * dl[n]);
        if (fabs_ < fdl) {
          dng[q] = ng1[q] * fdl;
          ng1[q] = ng1[q] - dng[q];
        } else {
          dng[q] = ng1[q] * fabs_;
          ng1[q] = ng1[q] - dng[q];
        }
      }

      // Consume photons along ray segment: (1-fgam) is real(4) arithmetic
      for (int q = 0; q < kNfreq; ++q)
        ng2[q] = ng2[q] - dng[q] * static_cast<double>(1.0f - fgam[q]);

      if (f::sum8<kNfreq>(ng2) > ngam_cut) {
        // Advance ray to end of segment (and slightly beyond); 1E-6 real(4)
        for (int d = 0; d < 3; ++d)
          x[d] = f::mod(x0[d] + (l + static_cast<double>(1e-6f)) * dx[d] + nrtr, nrtr);
        r = r0 + l;
        dr = dr0 - l;
        for (int q = 0; q < kNfreq; ++q)
          ngam[q] = ngam[q] + dng[q] * static_cast<double>(fgam[q]);
      }
    }
  }

  if (f::sum8<kNfreq>(ngam) < ngam_cut) {
    ray_remove(s, ir, irt, jrt, krt, true);
    ir = 0;
  } else {
    {
      Ray& ry = s.rays.at(ir);
      ry.r = static_cast<float>(r);
      ry.dr = static_cast<float>(dr);
      for (int d = 0; d < 3; ++d) ry.x[d] = static_cast<float>(x[d]);
      for (int q = 0; q < kNfreq; ++q) ry.ngam[q] = static_cast<float>(ngam[q]);
      for (int d = 0; d < 3; ++d) x[d] = ry.x[d];
    }
    const int i = cell_of(x[0], Nrt);
    const int j = cell_of(x[1], Nrt);
    const int k = cell_of(x[2], Nrt);
    if (i != irt || j != jrt || k != krt) {
      ray_remove(s, ir, irt, jrt, krt, false);
      ray_insert_tail(s, ir, i, j, k);
    }
  }
}

// ray_tracing_tools.f90: ray_segment
std::array<double, 4> ray_segment(const double x0[3], const double dx[3],
                                  double seg_min, int nseg) {
  if (nseg < 1 || nseg > 4) throw std::invalid_argument("ray_segment: nseg must be 1..4");
  std::array<double, 4> seg{};
  double l = 0;
  double x[3] = {x0[0], x0[1], x0[2]};
  const double eps = static_cast<double>(1e-6f);  // 1E-6 is real(4)

  for (int n = 0; n < nseg; ++n) {
    double dl = 1 - l;
    if (dl > seg_min) {
      for (int i = 0; i < 3; ++i) {
        if (std::fabs(dx[i]) > eps) {
          double ds;
          if (dx[i] < 0) {
            ds = f::floor4(x[i]) - x[i];
            if (ds > -eps) ds = -1;
          } else {
            ds = f::ceiling4(x[i]) - x[i];
            if (ds < +eps) ds = +1;
          }
          dl = f::min(dl, ds / dx[i]);
        }
      }
    }
    l = l + dl;
    for (int i = 0; i < 3; ++i) x[i] = x[i] + dl * dx[i];
    seg[n] = dl;
  }
  return seg;
}

// ---------------------------------------------------------------------------
// Storage / chain operations
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_create
void ray_create(Simulation& s, int nr, std::int64_t* indx) {
  RayState& R = s.rays;
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

  if (flag) {
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

// ray_tracing_tools.f90: lside_from_radius
int lside_from_radius(double r) {
  if (r < 1) return lside_min;
  // log(pi*ray_per_cell/3.*r**2)/log(4.) ; log(4.) is real(4)
  const int lside = f::ceiling4(
      std::log(pc::pi * ray_per_cell / 3.0 * f::pow2(r)) / folded::log4_r4);
  return f::min(f::max(lside, lside_min), lside_max);
}

// ray_tracing_tools.f90: nside_from_radius
int nside_from_radius(double r) { return 1 << lside_from_radius(r); }

// ray_tracing_tools.f90: nray_from_nside
int nray_from_nside(int nside) {
  if (nside < nside_min) return nray_init;
  return 12 * nside * nside;
}

// ray_tracing_tools.f90: rmin_from_nside
double rmin_from_nside(int nside) {
  if (nside == nside_min) return 0;
  return (nside / 2) * std::sqrt(3 / (pc::pi * ray_per_cell));
}

// ray_tracing_tools.f90: rmax_from_nside
double rmax_from_nside(int nside) {
  return nside * std::sqrt(3 / (pc::pi * ray_per_cell));
}

// ---------------------------------------------------------------------------
// interpolate (general_tools.f90), 'lin' branch
// ---------------------------------------------------------------------------

double interpolate_lin(double x, const double* fa, int n, int a, int b) {
  auto F = [&](int row, int i) { return fa[2 * (i - 1) + (row - 1)]; };
  const double dxmin = static_cast<double>(f::epsilon_r4);
  int i1 = 1, i2 = n;
  if (F(a, 1) < F(a, n)) {
    while (i2 - i1 > 1) {
      const int i = (i1 + i2) / 2;
      if (x > F(a, i))
        i1 = i;
      else
        i2 = i;
    }
  } else {
    while (i2 - i1 > 1) {
      const int i = (i1 + i2) / 2;
      if (x < F(a, i))
        i1 = i;
      else
        i2 = i;
    }
  }
  const double x1 = F(a, i1), y1 = F(b, i1), x2 = F(a, i2), y2 = F(b, i2);
  const double ddx = x2 - x1;
  if (std::fabs(ddx) < dxmin) return (y1 + y2) / 2;
  return y1 + (y2 - y1) * ((x - x1) / (x2 - x1));
}

// ---------------------------------------------------------------------------
// ray_redshifting (non-FREQ1)
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_redshifting / redshift_ray (#else FREQ1)
void ray_redshifting(Simulation& s) {
  RayState& R = s.rays;
  if (R.Nray <= 0) return;
  const Units& u = s.units;

  const double bmin = -6, bmax = 6;
  const double bdel = static_cast<double>(0.1f);
  const int Nbin = 1 + f::nint((bmax - bmin) / bdel);  // 121
  constexpr int NF = kNfreq;

  const double fscale = (s.clk.scalefactor - s.clk.da_sim) / s.clk.scalefactor;

  // Frequency edges (Nfreq = 5 case)
  double freq[NF + 1];
  freq[0] = pc::HI_eV;
  freq[1] = std::sqrt(pc::HI_eV * pc::HeI_eV);
  freq[2] = pc::HeI_eV;
  freq[3] = std::sqrt(pc::HeI_eV * pc::HeII_eV);
  freq[4] = pc::HeII_eV;
  freq[5] = static_cast<double>(1e2f);

  // Lookup table indx(2,Nbin,Nfreq-1) stored as tab[j][2*(i-1)+(row-1)]
  std::vector<std::vector<double>> tab(NF - 1, std::vector<double>(2 * Nbin));
  double rmax = 0.;
  double rmin = static_cast<double>(f::huge_r4);
  for (int j = 1; j <= NF - 1; ++j) {
    // freq(j), freq(j1), freq(j2) with j1=j+1, j2=j+2 (1-based)
    const double fj = freq[j - 1], fj1 = freq[j], fj2 = freq[j + 1];
    for (int i = 1; i <= Nbin; ++i) {
      const double b = (i - 1) * bdel + bmin;
      const double c = b + 1;
      double v;
      if (std::fabs(c) < bdel / 1e2) {
        v = std::log(fj2 / fj1) / std::log(fj1 / fj);
      } else {
        v = (std::pow(fj2, c) - std::pow(fj1, c)) / (std::pow(fj1, c) - std::pow(fj, c));
      }
      tab[j - 1][2 * (i - 1) + 1] = b;
      tab[j - 1][2 * (i - 1) + 0] = v;
      rmax = f::max(rmax, v);
      rmin = f::min(rmin, v);
    }
  }

  const std::int64_t top = R.hwm();
  for (std::int64_t ir = 1; ir <= top; ++ir) {
    Ray& ry = R.at(ir);
    if (!(ry.nside > 0)) continue;

    double nu[NF], ngam[NF], egam[NF], ngam1[NF], egam1[NF];
    for (int i = 0; i < NF; ++i) {
      nu[i] = ry.nu[i];
      ngam[i] = ry.ngam[i];
      egam[i] = nu[i] * ngam[i];
      ngam1[i] = 0;
      egam1[i] = 0;
    }

    for (int i = 1; i <= NF; ++i) {
      const int im = i - 1, ip = i + 1;
      if (ngam[i - 1] > u.ngam_min) {
        const double fq = freq[i - 1] / fscale;
        double rN, rE, b;
        if (i < NF) {
          rN = ngam[ip - 1] / ngam[i - 1];
          rN = f::min(f::max(rN, rmin), rmax);
          b = interpolate_lin(rN, tab[i - 1].data(), Nbin, 1, 2);
          b = f::min(f::max(b, bmin), bmax);
        } else {
          rN = ngam[i - 1] / f::max(ngam[im - 1], u.ngam_min);
          rN = f::min(f::max(rN, rmin), rmax);
          b = interpolate_lin(rN, tab[im - 1].data(), Nbin, 1, 2);
          b = f::min(f::max(b, bmin), bmax);
        }
        const double c = b + 1;
        const double d = b + 2;
        const double fi = freq[i - 1], fip = freq[ip - 1];

        if (std::fabs(c) < bdel / 1e2)
          rN = std::log(fq / fi) / std::log(fip / fi);
        else
          rN = (std::pow(fq, c) - std::pow(fi, c)) / (std::pow(fip, c) - std::pow(fi, c));

        if (std::fabs(d) < bdel / 1e2)
          rE = std::log(fq / fi) / std::log(fip / fi);
        else
          rE = (std::pow(fq, d) - std::pow(fi, d)) / (std::pow(fip, d) - std::pow(fi, d));

        rN = f::min(f::max(rN, 0.0), 1.0);
        rE = f::min(f::max(rE, 0.0), 1.0);
        ngam1[i - 1] = (1 - rN) * ngam[i - 1];
        egam1[i - 1] = (1 - rE) * egam[i - 1];

        if (i > 1) {
          ngam1[im - 1] = ngam1[im - 1] + rN * ngam[i - 1];
          egam1[im - 1] = egam1[im - 1] + rE * egam[i - 1];
        }
      }
    }

    for (int i = 1; i <= NF; ++i) {
      if (ngam1[i - 1] < u.ngam_min) {
        nu[i - 1] = cp::freq_init[i - 1];
        ngam[i - 1] = 0;
      } else {
        nu[i - 1] = egam1[i - 1] / ngam1[i - 1] * fscale;
        nu[i - 1] = f::min(f::max(nu[i - 1], freq[i - 1]), freq[i]);
        ngam[i - 1] = ngam1[i - 1];
      }
    }
    for (int i = 0; i < NF; ++i) {
      ry.nu[i] = static_cast<float>(nu[i]);
      ry.ngam[i] = static_cast<float>(ngam[i]);
    }

    // Enforce minimum for stability (real(4) sum)
    if (static_cast<double>(f::sum4<NF>(ry.ngam)) < u.ngam_min)
      ry.ngam[0] = static_cast<float>(u.ngam_min);
  }
}

// ---------------------------------------------------------------------------
// Ray buffer management (Ncpu = 1)
// ---------------------------------------------------------------------------

// ray_tracing_tools.f90: ray_buffer_allocate
void ray_buffer_allocate(Simulation& s) {
  RayState& R = s.rays;
  RayStat& st = R.stat;
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
