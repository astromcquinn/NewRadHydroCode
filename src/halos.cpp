// On-the-fly FoF halo finder (see include/radray/halos.hpp and
// docs/halos_and_sources.md).
#include "radray/halos.hpp"
#include "radray/fortran.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "radray/constants.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace radray {

namespace {

using i64 = std::int64_t;

// Minimum-image separation along one axis (any |d|, any wrap count).
inline double min_image(double d, double L) { return d - L * std::nearbyint(d / L); }

inline double wrap(double x, double L) {
  if (x >= 0 && x < L) return x;
  double y = x - L * std::floor(x / L);
  if (y >= L || y < 0) y = 0;  // rounding at the upper edge
  return y;
}

// ---------------------------------------------------------------------------
// Lock-free union-find. Links always go from the larger root to the smaller
// one, so parent[x] <= x, there are no cycles and the final root of every
// component is its smallest element: the result is independent of the order
// in which unions are performed (and therefore of the thread count).
// ---------------------------------------------------------------------------
struct UnionFind {
  std::vector<i64> parent;
  explicit UnionFind(i64 n) : parent(n) { std::iota(parent.begin(), parent.end(), i64{0}); }

  i64 load(i64 x) const { return __atomic_load_n(&parent[x], __ATOMIC_RELAXED); }

  i64 find(i64 x) {
    for (;;) {
      const i64 p = load(x);
      if (p == x) return x;
      const i64 gp = load(p);
      if (gp != p) {
        i64 expected = p;  // path halving; failure is harmless
        __atomic_compare_exchange_n(&parent[x], &expected, gp, false, __ATOMIC_RELAXED,
                                    __ATOMIC_RELAXED);
      }
      x = gp;
    }
  }

  void unite(i64 a, i64 b) {
    for (;;) {
      a = find(a);
      b = find(b);
      if (a == b) return;
      if (a < b) std::swap(a, b);  // a > b: hang a below b
      i64 expected = a;
      if (__atomic_compare_exchange_n(&parent[a], &expected, b, false, __ATOMIC_ACQ_REL,
                                      __ATOMIC_RELAXED))
        return;
    }
  }
};

// Linked-cell grid: particles sorted by cell (counting sort, stable in
// particle index), cell size >= the search radius it was built for.
struct CellGrid {
  int nc = 1;
  double L = 1;
  double inv_cell = 1;
  std::vector<i64> start;  // nc^3 + 1
  std::vector<i64> idx;    // particle indices, by cell then ascending index

  int cell1(double x) const {
    int c = static_cast<int>(wrap(x, L) * inv_cell);
    return c >= nc ? nc - 1 : c;
  }
  i64 cell_id(int i, int j, int k) const {
    return i + static_cast<i64>(nc) * (j + static_cast<i64>(nc) * k);
  }
};

CellGrid build_grid(const float* pos, i64 n, int stride, double L, double rmin, int ncell_max) {
  CellGrid g;
  g.L = L;
  int nc = static_cast<int>(std::floor(L / rmin));
  if (ncell_max <= 0) {
    // ~8 cells per particle on average at most, never more than 1024^3.
    const double cap = 2.0 * std::cbrt(static_cast<double>(std::max<i64>(n, 1)));
    ncell_max = static_cast<int>(std::min(1024.0, std::max(1.0, cap)));
  }
  nc = std::max(1, std::min(nc, ncell_max));
  g.nc = nc;
  g.inv_cell = nc / L;
  const i64 ncell = static_cast<i64>(nc) * nc * nc;
  std::vector<i64> cell(n);
#pragma omp parallel for schedule(static)
  for (i64 p = 0; p < n; ++p) {
    const float* x = pos + p * stride;
    cell[p] = g.cell_id(g.cell1(x[0]), g.cell1(x[1]), g.cell1(x[2]));
  }
  g.start.assign(ncell + 1, 0);
  for (i64 p = 0; p < n; ++p) ++g.start[cell[p] + 1];
  for (i64 c = 0; c < ncell; ++c) g.start[c + 1] += g.start[c];
  g.idx.resize(n);
  std::vector<i64> fill(g.start.begin(), g.start.end() - 1);
  for (i64 p = 0; p < n; ++p) g.idx[fill[cell[p]]++] = p;
  return g;
}

// Unique (sorted) periodic neighbour cells of cell (i,j,k) at offsets
// -r..r (deduplicated for small nc).
void neighbour_cells(const CellGrid& g, int i, int j, int k, int r, std::vector<i64>& out) {
  out.clear();
  const int nc = g.nc;
  for (int dk = -r; dk <= r; ++dk)
    for (int dj = -r; dj <= r; ++dj)
      for (int di = -r; di <= r; ++di) {
        const int a = ((i + di) % nc + nc) % nc;
        const int b = ((j + dj) % nc + nc) % nc;
        const int c = ((k + dk) % nc + nc) % nc;
        out.push_back(g.cell_id(a, b, c));
      }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
}

inline double dist2(const float* a, const float* b, double L) {
  const double dx = min_image(static_cast<double>(a[0]) - b[0], L);
  const double dy = min_image(static_cast<double>(a[1]) - b[1], L);
  const double dz = min_image(static_cast<double>(a[2]) - b[2], L);
  return dx * dx + dy * dy + dz * dz;
}

FofGroups collect_groups(UnionFind& uf, i64 n, i64 nmin) {
  FofGroups G;
  G.label.resize(n);
#pragma omp parallel for schedule(static)
  for (i64 p = 0; p < n; ++p) G.label[p] = uf.find(p);
  // Count members per root (reuse uf.parent as counter storage)
  std::vector<i64>& cnt = uf.parent;
  std::fill(cnt.begin(), cnt.end(), 0);
  for (i64 p = 0; p < n; ++p) ++cnt[G.label[p]];
  // Group index per root, -1 for small groups
  i64 ng = 0;
  for (i64 p = 0; p < n; ++p) {
    if (G.label[p] == p && cnt[p] >= nmin) {
      G.group_label.push_back(p);
      cnt[p] = ng++;
    } else if (G.label[p] == p) {
      cnt[p] = -1;
    }
  }
  G.offset.assign(ng + 1, 0);
  for (i64 p = 0; p < n; ++p) {
    const i64 g = cnt[G.label[p]];
    if (g >= 0) ++G.offset[g + 1];
  }
  for (i64 g = 0; g < ng; ++g) G.offset[g + 1] += G.offset[g];
  G.members.resize(G.offset[ng]);
  std::vector<i64> fill(G.offset.begin(), G.offset.end() - 1);
  for (i64 p = 0; p < n; ++p) {
    const i64 g = cnt[G.label[p]];
    if (g >= 0) G.members[fill[g]++] = p;
  }
  return G;
}

}  // namespace

// ---------------------------------------------------------------------------
// FoF
// ---------------------------------------------------------------------------

FofGroups fof_partition(const float* pos, i64 npart, int stride, double L, double ll,
                        i64 nmin, int ncell_max) {
  if (!(L > 0) || !(ll > 0)) throw std::invalid_argument("fof_partition: L and ll must be > 0");
  if (ll >= 0.5 * L) throw std::invalid_argument("fof_partition: ll must be < L/2");
  UnionFind uf(npart);
  const CellGrid g = build_grid(pos, npart, stride, L, ll, ncell_max);
  const int nc = g.nc;
  const i64 ncell = static_cast<i64>(nc) * nc * nc;
  const double ll2 = ll * ll;

#pragma omp parallel
  {
    std::vector<i64> nb;
#pragma omp for schedule(dynamic, 64)
    for (i64 c = 0; c < ncell; ++c) {
      const i64 b0 = g.start[c], b1 = g.start[c + 1];
      if (b0 == b1) continue;
      const int i = static_cast<int>(c % nc);
      const int j = static_cast<int>((c / nc) % nc);
      const int k = static_cast<int>(c / (static_cast<i64>(nc) * nc));
      neighbour_cells(g, i, j, k, 1, nb);
      for (i64 c2 : nb) {
        if (c2 < c) continue;  // each unordered cell pair once
        const i64 e0 = g.start[c2], e1 = g.start[c2 + 1];
        for (i64 a = b0; a < b1; ++a) {
          const i64 p = g.idx[a];
          const float* xp = pos + p * stride;
          for (i64 b = (c2 == c ? a + 1 : e0); b < e1; ++b) {
            const i64 q = g.idx[b];
            if (dist2(xp, pos + q * stride, L) <= ll2) uf.unite(p, q);
          }
        }
      }
    }
  }
  return collect_groups(uf, npart, nmin);
}

FofGroups fof_partition_bruteforce(const float* pos, i64 npart, int stride, double L, double ll,
                                   i64 nmin) {
  UnionFind uf(npart);
  const double ll2 = ll * ll;
  for (i64 p = 0; p < npart; ++p)
    for (i64 q = p + 1; q < npart; ++q)
      if (dist2(pos + p * stride, pos + q * stride, L) <= ll2) uf.unite(p, q);
  return collect_groups(uf, npart, nmin);
}

// ---------------------------------------------------------------------------
// Units and virial quantities
// ---------------------------------------------------------------------------

double particle_mass_msun(i64 npart_total, double box) {
  using namespace pc;
  using namespace cp;
  // cosmo_tools.f90: mp_unit = rhoc_cgs*om*(box*Mpc2cm)**3/h0/Npart
  // (same expression as Units::mp_unit in calc_sim_var)
  const double mp_unit = rhoc_cgs * omegam * f::pow3(box * Mpc2cm) / hubble0 /
                         static_cast<double>(npart_total);
  return mp_unit / Msun_cgs;
}

HaloUnits make_halo_units(const Config& cfg, const SimClock& clk, const Units& u) {
  HaloUnits h;
  const i64 Npart = static_cast<i64>(cfg.Ndm) * cfg.Ndm * cfg.Ndm;
  h.L = cfg.Ndm;
  h.box = cfg.box;
  h.pos_to_mpch = cfg.x_dm_box;                       // box/ndmr
  const double x_dm_gas = static_cast<double>(cfg.Ngas / cfg.Ndm);
  h.vel_to_kms = x_dm_gas * u.vel_unit / 1e5;         // pm_tools.f90 summarize_dm
  h.mpart_msun = u.Msun_unit;                         // calc_sim_var
  (void)Npart;
  h.redshift = clk.redshift;
  h.hubble0 = cp::hubble0;
  h.rho_crit = u.rho_crit;
  h.rhom_mean = u.rho_unit;                           // rhoc*om*h0^2/a^3
  h.omegam_z = h.rhom_mean / h.rho_crit;
  return h;
}

double delta_vir_bryan_norman(double omegam_z) {
  const double x = omegam_z - 1;
  return 18 * pc::pi * pc::pi + 82 * x - 39 * x * x;
}

double virial_radius_mpch(double M, double delta_crit, const HaloUnits& hu) {
  // M = 4pi/3 delta rho_crit R_phys^3
  const double Mg = M * pc::Msun_cgs;
  const double Rphys = std::cbrt(3 * Mg / (4 * pc::pi * delta_crit * hu.rho_crit));
  return Rphys * (1 + hu.redshift) / pc::Mpc2cm * hu.hubble0;
}

double virial_temperature(double M, double R_mpch, double mu, const HaloUnits& hu) {
  const double Rphys = R_mpch / hu.hubble0 * pc::Mpc2cm / (1 + hu.redshift);
  return mu * pc::mH_cgs * pc::G_cgs * M * pc::Msun_cgs / (2 * pc::k_cgs * Rphys);
}

// ---------------------------------------------------------------------------
// Halo properties
// ---------------------------------------------------------------------------

namespace {

// Centre of mass of a member list, unwrapped around the first member.
void centre_of_mass(const float* xv, const i64* mem, i64 n, double L, double c[3]) {
  const float* x0 = xv + 6 * mem[0];
  double s[3] = {0, 0, 0};
  for (i64 a = 0; a < n; ++a) {
    const float* x = xv + 6 * mem[a];
    for (int d = 0; d < 3; ++d) s[d] += min_image(static_cast<double>(x[d]) - x0[d], L);
  }
  for (int d = 0; d < 3; ++d) c[d] = wrap(x0[d] + s[d] / n, L);
}

void potential_minimum(const float* xv, const i64* mem, i64 n, double L, double eps,
                       i64 nmax, double c[3]) {
  const i64 step = (n + nmax - 1) / nmax;
  const double eps2 = eps * eps;
  double best = 0;
  i64 ibest = 0;
  for (i64 a = 0; a < n; ++a) {
    const float* xa = xv + 6 * mem[a];
    double phi = 0;
    for (i64 b = 0; b < n; b += step) {
      if (b == a) continue;
      phi -= 1 / std::sqrt(dist2(xa, xv + 6 * mem[b], L) + eps2);
    }
    if (a == 0 || phi < best) {  // strict: ties keep the smaller index
      best = phi;
      ibest = a;
    }
  }
  const float* x = xv + 6 * mem[ibest];
  for (int d = 0; d < 3; ++d) c[d] = wrap(x[d], L);
}

void mean_velocity(const float* xv, const i64* mem, i64 n, double v[3]) {
  double s[3] = {0, 0, 0};
  for (i64 a = 0; a < n; ++a)
    for (int d = 0; d < 3; ++d) s[d] += xv[6 * mem[a] + 3 + d];
  for (int d = 0; d < 3; ++d) v[d] = s[d] / n;
}

// All particles within r of c (particle units), sorted by (distance, index).
void gather_sphere(const float* xv, const CellGrid& g, const double c[3], double r,
                   std::vector<std::pair<double, i64>>& out) {
  out.clear();
  const int nc = g.nc;
  const double cs = g.L / nc;
  const int nr = std::min(nc / 2, static_cast<int>(std::ceil(r / cs)));
  std::vector<i64> nb;
  neighbour_cells(g, g.cell1(c[0]), g.cell1(c[1]), g.cell1(c[2]), nr, nb);
  const double r2 = r * r;
  for (i64 cell : nb)
    for (i64 a = g.start[cell]; a < g.start[cell + 1]; ++a) {
      const i64 p = g.idx[a];
      const float* x = xv + 6 * p;
      const double dx = min_image(x[0] - c[0], g.L);
      const double dy = min_image(x[1] - c[1], g.L);
      const double dz = min_image(x[2] - c[2], g.L);
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 <= r2) out.emplace_back(d2, p);
    }
  std::sort(out.begin(), out.end());
}

}  // namespace

std::vector<Halo> find_halos(const float* xv, i64 npart, const HaloUnits& hu,
                             const HaloFinderParams& p) {
  if (npart <= 0) return {};
  const double L = hu.L;
  const double spacing = L / std::cbrt(static_cast<double>(npart));
  const double ll = p.b * spacing;
  const FofGroups G = fof_partition(xv, npart, 6, L, ll, std::max<i64>(p.min_npart, 1),
                                    p.ncell_max);
  const i64 ng = G.ngroups();

  const double dvir = p.so_delta_crit > 0 ? p.so_delta_crit : delta_vir_bryan_norman(hu.omegam_z);

  // SO threshold in particles per (particle unit)^3:
  // rho_com [Msun/(Mpc/h)^3] = delta*rho_crit*a^3 * (Mpc2cm/h)^3 / Msun
  const double a = 1 / (1 + hu.redshift);
  const double rho_thr_com = dvir * hu.rho_crit * a * a * a *
                             std::pow(pc::Mpc2cm / hu.hubble0, 3) / pc::Msun_cgs;
  const double n_thr = rho_thr_com * std::pow(hu.pos_to_mpch, 3) / hu.mpart_msun;

  CellGrid grid;
  if (p.so_refine) grid = build_grid(xv, npart, 6, L, spacing, p.ncell_max);

  std::vector<Halo> out(ng);
  std::vector<char> keep(ng, 1);

#pragma omp parallel
  {
    std::vector<std::pair<double, i64>> sph;
    std::vector<i64> inner;
#pragma omp for schedule(dynamic, 1)
    for (i64 gi = 0; gi < ng; ++gi) {
      const i64* mem = &G.members[G.offset[gi]];
      const i64 n = G.size(gi);
      Halo& h = out[gi];
      h.id = G.group_label[gi];
      h.N_fof = n;
      double c[3];
      if (p.centre == HaloCentre::PotentialMinimum)
        potential_minimum(xv, mem, n, L, 0.1 * ll, std::max<i64>(p.potmin_nmax, 1), c);
      else
        centre_of_mass(xv, mem, n, L, c);
      double v[3];
      mean_velocity(xv, mem, n, v);
      i64 N = n;
      double R_pu = 0;  // radius in particle units (SO only)

      if (p.so_refine) {
        // Search radius: enclose all FoF members and twice the virial
        // radius of the FoF mass, grown until the outermost crossing is
        // inside the search sphere.
        double rmax2 = 0;
        for (i64 m = 0; m < n; ++m) {
          const float* x = xv + 6 * mem[m];
          const double dx = min_image(x[0] - c[0], L), dy = min_image(x[1] - c[1], L),
                       dz = min_image(x[2] - c[2], L);
          rmax2 = std::max(rmax2, dx * dx + dy * dy + dz * dz);
        }
        const double rvir_pu = virial_radius_mpch(n * hu.mpart_msun, dvir, hu) / hu.pos_to_mpch;
        double rs = std::max(std::sqrt(rmax2), 2 * rvir_pu) + spacing;
        i64 kso = 0;
        double rso = 0;
        for (int it = 0; it <= p.so_iterations; ++it) {
          for (;;) {
            gather_sphere(xv, grid, c, rs, sph);
            kso = 0;
            rso = 0;
            // Outermost radius where the enclosed mean density >= threshold
            for (std::size_t m = 0; m < sph.size(); ++m) {
              const double r = std::sqrt(sph[m].first);
              const double k = static_cast<double>(m + 1);
              if (r > 0 && k >= n_thr * (4.0 / 3.0) * pc::pi * r * r * r) {
                kso = static_cast<i64>(m + 1);
                rso = r;
              }
            }
            if (kso < static_cast<i64>(sph.size()) || rs >= 0.45 * L) break;
            rs = std::min(2 * rs, 0.45 * L);
          }
          if (kso == 0) break;
          inner.clear();
          for (i64 m = 0; m < kso; ++m) inner.push_back(sph[m].second);
          std::sort(inner.begin(), inner.end());
          if (it == p.so_iterations) break;
          centre_of_mass(xv, inner.data(), kso, L, c);
        }
        N = kso;
        R_pu = rso;
        if (kso > 0) mean_velocity(xv, inner.data(), kso, v);
        if (N < p.min_npart) keep[gi] = 0;
      }

      h.N = N;
      h.M = N * hu.mpart_msun;
      for (int d = 0; d < 3; ++d) {
        h.x[d] = c[d] * hu.pos_to_mpch;
        h.v[d] = v[d] * hu.vel_to_kms;
      }
      h.R = p.so_refine ? R_pu * hu.pos_to_mpch : virial_radius_mpch(h.M, dvir, hu);
      h.T = h.M > 0 ? virial_temperature(h.M, h.R > 0 ? h.R : 1, p.mu, hu) : 0;
    }
  }

  std::vector<Halo> res;
  res.reserve(ng);
  for (i64 gi = 0; gi < ng; ++gi)
    if (keep[gi]) res.push_back(out[gi]);
  std::sort(res.begin(), res.end(), [](const Halo& A, const Halo& B) {
    if (A.M != B.M) return A.M > B.M;
    if (A.N != B.N) return A.N > B.N;
    return A.id < B.id;
  });
  return res;
}

}  // namespace radray
