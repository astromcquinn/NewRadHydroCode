// Case A recombination and recombination radiation (recomb_rad.hpp,
// docs/recomb_radiation.md). New code, not in the original.
#include "radray/recomb_rad.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "radray/constants.hpp"
#include "radray/fortran.hpp"
#include "radray/healpix.hpp"
#include "radray/parallel.hpp"
#include "radray/ray_tracing.hpp"
#include "radray/simulation.hpp"

namespace radray {

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

bool set_recomb_param(Config& c, const std::string& key, const std::string& val) {
  if (key == "recomb_case") {
    if (val == "A" || val == "a")
      c.recomb_case_A = true;
    else if (val == "B" || val == "b")
      c.recomb_case_A = false;
    else
      throw std::runtime_error("recomb_case must be A or B");
    return true;
  }
  if (key == "recomb_rad") {
    c.recomb_rad = std::stoi(val) != 0;
    return true;
  }
  if (key == "recomb_rad_nside") {
    const int n = std::stoi(val);
    if (n != 1 && n != 2 && n != 4)
      throw std::runtime_error("recomb_rad_nside must be 1, 2 or 4");
    c.recomb_rad_nside = n;
    return true;
  }
  if (key == "recomb_rad_block") {
    c.recomb_rad_block = std::stoi(val);
    if (c.recomb_rad_block < 1) throw std::runtime_error("recomb_rad_block must be >= 1");
    return true;
  }
  if (key == "recomb_rad_pool") {
    c.recomb_rad_pool = std::stod(val);
    if (!(c.recomb_rad_pool > 0 && c.recomb_rad_pool < 0.9))
      throw std::runtime_error("recomb_rad_pool must be in (0, 0.9)");
    return true;
  }
  return false;
}

void copy_recomb_param(const Config& from, Config& to) {
  to.recomb_case_A = from.recomb_case_A;
  to.recomb_rad = from.recomb_rad;
  to.recomb_rad_nside = from.recomb_rad_nside;
  to.recomb_rad_block = from.recomb_rad_block;
  to.recomb_rad_pool = from.recomb_rad_pool;
}

void check_recomb_param(const Config& c) {
  if (!c.recomb_rad) return;
  if (!c.recomb_case_A)
    throw std::runtime_error(
        "recomb_rad = 1 needs recomb_case = A: with case B the ground-state "
        "recombinations are already absorbed on the spot");
  if (c.rt_method != RTMethod::Rays) throw std::runtime_error("recomb_rad = 1 needs rt_method = rays");
  if (c.legacy_rt) throw std::runtime_error("recomb_rad = 1 needs legacy_rt = 0 (the pool guard)");
  if (c.Nrt % c.recomb_rad_block != 0)
    throw std::runtime_error("recomb_rad_block must divide N");
}

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

void recomb_init(Simulation& s) {
  if (!s.cfg.recomb_rad) return;
  check_recomb_param(s.cfg);
  if (!s.recomb) s.recomb = std::make_shared<RecombState>();
  RecombState& Rc = *s.recomb;
  const std::size_t n = static_cast<std::size_t>(s.cfg.Nrtc) * kNfreq;
  if (Rc.ngam.size() != n) {
    Rc.ngam.assign(n, 0.0f);
    Rc.egam.assign(n, 0.0f);
  }
  Rc.emit_ngam.assign(n, 0.0f);
  Rc.emit_egam.assign(n, 0.0f);
  Rc.emit.assign(static_cast<std::size_t>(s.cfg.Nrtc), 0);
}

// photo_cell. Mean photon energy of a ground-state recombination: the
// threshold plus the mean kinetic energy of the recombining electrons, kT
// (sigma_rec ~ 1/E against a Maxwellian), kept inside the group.
void recomb_count(Simulation& s, std::int64_t c, const double nrab[3], double T,
                  double vol_gas) {
  RecombState& Rc = *s.recomb;
  const Units& u = s.units;
  static constexpr int grp[3] = {0, 2, 4};
  const double nu0[3] = {pc::HI_eV, pc::HeI_eV, pc::HeII_eV};
  const double nu_hi[3] = {std::sqrt(pc::HI_eV * pc::HeI_eV),
                           std::sqrt(pc::HeI_eV * pc::HeII_eV), 1e2};
  const double kT = pc::K2eV * T;
  for (int sp = 0; sp < 3; ++sp) {
    if (!(nrab[sp] > 0)) continue;
    const double ng = nrab[sp] * vol_gas / (u.nH_avg * u.vol_rt);
    const double nu = std::min(nu0[sp] + kT, nu_hi[sp]);
    const std::size_t k = static_cast<std::size_t>(c) * kNfreq + grp[sp];
    Rc.ngam[k] = static_cast<float>(static_cast<double>(Rc.ngam[k]) + ng);
    Rc.egam[k] = static_cast<float>(static_cast<double>(Rc.egam[k]) + ng * nu);
  }
}

// ---------------------------------------------------------------------------
// recomb_prepare (rt_step_1)
// ---------------------------------------------------------------------------

std::int64_t recomb_prepare(Simulation& s, std::int64_t L, std::int64_t nsrc, std::int64_t nsplit) {
  RecombState& Rc = *s.recomb;
  RTGrid& G = s.rt;
  const Units& u = s.units;
  const int Nrt = s.cfg.Nrt;
  const int b = s.cfg.recomb_rad_block;
  const int nb = Nrt / b;
  const int ns = s.cfg.recomb_rad_nside;
  const std::int64_t nray = 12 * ns * ns;
  const std::int64_t nblk = static_cast<std::int64_t>(nb) * nb * nb;
  // The emitting cell of every block rotates through the block
  const int b3 = b * b * b;
  const int rot = static_cast<int>((31LL * s.clk.istep + s.clk.istep_rt) % b3);
  const int oi = rot % b, oj = (rot / b) % b, ok = rot / (b * b);

  // Pass 1: per block, the photons per ray and the class
  // (0 none, 1 too faint or optically thick: on the spot, 2 candidate)
  std::vector<float> pray(static_cast<std::size_t>(nblk), 0.0f);
  std::vector<std::uint8_t> cls(static_cast<std::size_t>(nblk), 0);
  const double lmin = -12, ldel = 0.01;
  const int NH = 1600;  // log10(photons per ray) from 1e-12 to 1e4
  std::vector<std::int64_t> hist(NH + 1, 0);
#pragma omp parallel if (par::worth(s.cfg.Nrtc))
  {
    std::vector<std::int64_t> hist_t(NH + 1, 0);
#pragma omp for schedule(static)
    for (int bk = 0; bk < nb; ++bk)
      for (int bj = 0; bj < nb; ++bj)
        for (int bi = 0; bi < nb; ++bi) {
          const std::int64_t ib = (static_cast<std::int64_t>(bk) * nb + bj) * nb + bi;
          double tot = 0;
          for (int k = 0; k < b; ++k)
            for (int j = 0; j < b; ++j)
              for (int i = 0; i < b; ++i) {
                const std::int64_t c = G.cell(bi * b + i + 1, bj * b + j + 1, bk * b + k + 1);
                for (int n = 0; n < kNfreq; ++n)
                  tot += static_cast<double>(Rc.ngam[static_cast<std::size_t>(c) * kNfreq + n]);
              }
          if (!(tot > 0)) continue;
          const double p = tot / static_cast<double>(nray);
          pray[ib] = static_cast<float>(p);
          const std::int64_t ce = G.cell(bi * b + oi + 1, bj * b + oj + 1, bk * b + ok + 1);
          const float* tau = G.tau_p(ce);
          float tmin = tau[0];
          for (int n = 1; n < kNfreq; ++n) tmin = std::min(tmin, tau[n]);
          // ray_source's test for an optically thick source cell; rays below
          // ngam_min would be removed at once
          if (std::exp(-tmin) < 1e-2f || p < u.ngam_min) {
            cls[ib] = 1;
          } else {
            cls[ib] = 2;
            int h = 1 + static_cast<int>(std::floor((std::log10(p) - lmin) / ldel));
            h = std::min(std::max(h, 1), NH);
            hist_t[h] += 1;
          }
        }
#pragma omp critical(radray_recomb_hist)
    for (int h = 1; h <= NH; ++h) hist[h] += hist_t[h];
  }

  // Pool budget: launch only while L + stellar source rays + the rays of the
  // expected stellar splits + recombination rays stay below
  // recomb_rad_pool * Nray_max; the brightest blocks first
  // (whole histogram bins, so the choice does not depend on the order).
  const RayState& R = s.rays;
  const std::int64_t room = static_cast<std::int64_t>(s.cfg.recomb_rad_pool * static_cast<double>(R.Nray_max)) -
                            L - rayp::nray_init * nsrc - nsplit;
  const std::int64_t nmax = room > 0 ? room / nray : 0;
  std::int64_t ncand = 0;
  for (int h = 1; h <= NH; ++h) ncand += hist[h];
  int hcut = 1;  // candidates in bins >= hcut launch
  if (ncand > nmax) {
    std::int64_t nr = 0;
    hcut = NH + 1;
    for (int h = NH; h >= 1; --h) {
      if (nr + hist[h] > nmax) break;
      nr += hist[h];
      hcut = h;
    }
  }
  Rc.nblock_last = ncand;
  // Coarser merging of the faint rays (12 instead of 48 direction bins per
  // cell) while the pool is crowded, so that the live count can fall
  s.rays.nside_merge = (L > static_cast<std::int64_t>(s.cfg.recomb_rad_pool * static_cast<double>(R.Nray_max))) ? 1 : 2;

  // Pass 2: launch or absorb on the spot; per-plane sums, added in order
  std::vector<double> sum_prod(nb, 0), sum_launch(nb, 0), sum_thick(nb, 0), sum_pool(nb, 0);
  std::vector<std::int64_t> nemit(nb, 0);
#pragma omp parallel for schedule(static) if (par::worth(s.cfg.Nrtc))
  for (int bk = 0; bk < nb; ++bk)
    for (int bj = 0; bj < nb; ++bj)
      for (int bi = 0; bi < nb; ++bi) {
        const std::int64_t ib = (static_cast<std::int64_t>(bk) * nb + bj) * nb + bi;
        if (cls[ib] == 0) continue;
        bool launch = false;
        if (cls[ib] == 2) {
          int h = 1 + static_cast<int>(std::floor((std::log10(static_cast<double>(pray[ib])) - lmin) / ldel));
          h = std::min(std::max(h, 1), NH);
          launch = h >= hcut;
        }
        const std::int64_t ce = G.cell(bi * b + oi + 1, bj * b + oj + 1, bk * b + ok + 1);
        double tot = 0;
        for (int k = 0; k < b; ++k)
          for (int j = 0; j < b; ++j)
            for (int i = 0; i < b; ++i) {
              const std::int64_t c = G.cell(bi * b + i + 1, bj * b + j + 1, bk * b + k + 1);
              float* ng = &Rc.ngam[static_cast<std::size_t>(c) * kNfreq];
              float* eg = &Rc.egam[static_cast<std::size_t>(c) * kNfreq];
              for (int n = 0; n < kNfreq; ++n) {
                if (!(ng[n] > 0)) continue;
                tot += static_cast<double>(ng[n]);
                if (launch) {
                  const std::size_t ke = static_cast<std::size_t>(ce) * kNfreq + n;
                  Rc.emit_ngam[ke] += ng[n];
                  Rc.emit_egam[ke] += eg[n];
                } else {
                  // on the spot: into the cell's photon field before
                  // photo_cell (gnu holds ngam*nu until photo_cell)
                  float* gng = G.ngam_p(c);
                  float* gnu = G.nu_p(c);
                  gng[n] = static_cast<float>(static_cast<double>(gng[n]) + static_cast<double>(ng[n]));
                  gnu[n] = static_cast<float>(static_cast<double>(gnu[n]) + static_cast<double>(eg[n]));
                }
                ng[n] = 0;
                eg[n] = 0;
              }
            }
        sum_prod[bk] += tot;
        if (launch) {
          Rc.emit[ce] = 1;
          sum_launch[bk] += tot;
          nemit[bk] += 1;
        } else if (cls[ib] == 1) {
          sum_thick[bk] += tot;
        } else {
          sum_pool[bk] += tot;
        }
      }
  const double norm = 1.0 / static_cast<double>(s.cfg.Nrtc);
  std::int64_t ne = 0;
  for (int bk = 0; bk < nb; ++bk) {
    Rc.produced += sum_prod[bk] * norm;
    Rc.launched += sum_launch[bk] * norm;
    Rc.local_thick += sum_thick[bk] * norm;
    Rc.local_pool += sum_pool[bk] * norm;
    ne += nemit[bk];
  }
  Rc.nemit_last = ne;
  if (s.raylog)
    *s.raylog << " Recomb     : produced " << static_cast<float>(Rc.produced) << " launched "
              << static_cast<float>(Rc.launched) << " thick " << static_cast<float>(Rc.local_thick)
              << " pool " << static_cast<float>(Rc.local_pool) << " blocks " << ncand << " emit "
              << ne << " merge_nside " << s.rays.nside_merge << "\n";
  return ne * nray;
}

// ---------------------------------------------------------------------------
// ray_source_recomb (ray_tracing_1)
// ---------------------------------------------------------------------------

void ray_source_recomb(Simulation& s, int irt, int jrt, int krt) {
  RTGrid& G = s.rt;
  RayState& R = s.rays;
  RecombState& Rc = *s.recomb;
  const std::int64_t c = G.cell(irt, jrt, krt);
  Rc.emit[c] = 0;
  const int ns = s.cfg.recomb_rad_nside;
  const int nray = 12 * ns * ns;
  const double dr = rmax_from_nside(ns);
  const float xs[3] = {static_cast<float>(irt) - 0.5f, static_cast<float>(jrt) - 0.5f,
                       static_cast<float>(krt) - 0.5f};
  double nu[kNfreq], ngam[kNfreq];
  for (int n = 0; n < kNfreq; ++n) {
    const std::size_t k = static_cast<std::size_t>(c) * kNfreq + n;
    const double ng = Rc.emit_ngam[k];
    nu[n] = ng > 0 ? static_cast<double>(Rc.emit_egam[k]) / ng : cp::freq_init[n];
    ngam[n] = ng / nray;
    Rc.emit_ngam[k] = 0;
    Rc.emit_egam[k] = 0;
  }
  std::int64_t indx[rayp::nray_init];
  ray_create(s, nray, indx);
  for (int i = 0; i < nray; ++i) {
    double dx[3];
    healpix::pix2vec_nest(ns, i, dx);
    const std::int64_t ir = indx[i];
    Ray& ry = R.at(ir);
    ry.nside = ns;
    ry.ipix = i;
    ry.r = 0;
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
    ray_insert_head(s, ir, irt, jrt, krt);
  }
  // as ray_source: new rays at the head; the tail only if the chain was empty
  G.htoc_ray(1, c) = indx[nray - 1];
  if (G.htoc_ray(2, c) == 0) G.htoc_ray(2, c) = indx[0];
}

}  // namespace radray
