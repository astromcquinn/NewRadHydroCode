// Case A recombination and recombination radiation (recomb_rad.hpp,
// docs/recomb_radiation.md). New code, not in the original.
#include "radray/recomb_rad.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "radray/constants.hpp"
#include "radray/fortran.hpp"
#include "radray/kernels/ray_kernels.hpp"
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
  const std::size_t nc = static_cast<std::size_t>(s.cfg.Nrtc);
  const std::size_t n = nc * kNfreq;
  if (Rc.ngam.size() != n) {
    Rc.ngam.assign(n, 0.0f);
    Rc.egam.assign(n, 0.0f);
    Rc.ots.assign(nc, 0.0f);
    Rc.mode.assign(nc, 1);  // launch until a block has shown its production
  }
  Rc.emit_ngam.assign(n, 0.0f);
  Rc.emit_egam.assign(n, 0.0f);
  Rc.emit.assign(nc, 0);
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
  const bool onspot = recomb_onspot(Rc, c);
  for (int sp = 0; sp < 3; ++sp) {
    if (!(nrab[sp] > 0)) continue;
    const double ng = nrab[sp] * vol_gas / (u.nH_avg * u.vol_rt);
    if (onspot) {
      Rc.ots[c] = static_cast<float>(static_cast<double>(Rc.ots[c]) + ng);
      continue;
    }
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
  auto cell_in = [&](int bi, int bj, int bk, int i, int j, int k) {
    return G.cell(bi * b + i + 1, bj * b + j + 1, bk * b + k + 1);
  };

  // Pass 1, per block (blocks are independent; per-plane sums added in
  // order): launch the buffer of a launch-mode block from its emitting cell,
  // tally the on-the-spot count, and keep the block's production of the
  // last substep (photons per ray) and class for the choice of modes
  // (0 no production, 1 on the spot by optical depth or faintness,
  // 2 candidate for launch).
  std::vector<float> pray(static_cast<std::size_t>(nblk), 0.0f);
  std::vector<std::uint8_t> cls(static_cast<std::size_t>(nblk), 0);
  std::vector<double> sum_launch(nb, 0), sum_drop(nb, 0), sum_ots(nb, 0);
  std::vector<std::int64_t> nemit(nb, 0);
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
          const std::int64_t ce = cell_in(bi, bj, bk, oi, oj, ok);
          double tot = 0, tots = 0;
          for (int k = 0; k < b; ++k)
            for (int j = 0; j < b; ++j)
              for (int i = 0; i < b; ++i) {
                const std::int64_t c = cell_in(bi, bj, bk, i, j, k);
                for (int n = 0; n < kNfreq; ++n)
                  tot += static_cast<double>(Rc.ngam[static_cast<std::size_t>(c) * kNfreq + n]);
                tots += static_cast<double>(Rc.ots[c]);
                Rc.ots[c] = 0;
              }
          sum_ots[bk] += tots;
          if (tot > 0) {
            const bool launch = tot / static_cast<double>(nray) >= u.ngam_min;
            for (int k = 0; k < b; ++k)
              for (int j = 0; j < b; ++j)
                for (int i = 0; i < b; ++i) {
                  const std::int64_t c = cell_in(bi, bj, bk, i, j, k);
                  float* ng = &Rc.ngam[static_cast<std::size_t>(c) * kNfreq];
                  float* eg = &Rc.egam[static_cast<std::size_t>(c) * kNfreq];
                  for (int n = 0; n < kNfreq; ++n) {
                    if (launch && ng[n] > 0) {
                      const std::size_t ke = static_cast<std::size_t>(ce) * kNfreq + n;
                      Rc.emit_ngam[ke] += ng[n];
                      Rc.emit_egam[ke] += eg[n];
                    }
                    ng[n] = 0;
                    eg[n] = 0;
                  }
                }
            if (launch) {
              Rc.emit[ce] = 1;
              sum_launch[bk] += tot;
              nemit[bk] += 1;
            } else {
              sum_drop[bk] += tot;
            }
          }
          // the block's production in the last substep, in either mode
          const double p = (tot + tots) / static_cast<double>(nray);
          if (!(p > 0)) continue;
          pray[ib] = static_cast<float>(p);
          const float* tau = G.tau_p(ce);
          float tmin = tau[0];
          for (int n = 1; n < kNfreq; ++n) tmin = std::min(tmin, tau[n]);
          // ray_source's test for an optically thick source cell; rays
          // below ngam_min would be removed at once
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
  std::int64_t ne = 0;
  for (int bk = 0; bk < nb; ++bk) ne += nemit[bk];

  // Pool budget for the blocks in launch mode this substep (their rays are
  // created at the next substep): L + stellar source rays + the rays of the
  // expected stellar splits + this substep's recombination rays + 12 nside^2
  // per launch block stay below recomb_rad_pool * Nray_max. The brightest
  // blocks first, in whole histogram bins, so the choice does not depend on
  // the order. Blocks without production stay in launch mode at no cost.
  const RayState& R = s.rays;
  const std::int64_t cap = static_cast<std::int64_t>(s.cfg.recomb_rad_pool * static_cast<double>(R.Nray_max));
  const std::int64_t room = cap - L - rayp::nray_init * nsrc - nsplit - ne * nray;
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
  // Coarser merging of the faint rays (12 instead of 48 direction bins per
  // cell) while the pool is crowded, so that the live count can fall
  s.rays.nside_merge = (L > cap) ? 1 : 2;

  // Pass 2: the modes of this substep
  std::vector<std::int64_t> nlaunch(nb, 0);
#pragma omp parallel for schedule(static) if (par::worth(s.cfg.Nrtc))
  for (int bk = 0; bk < nb; ++bk)
    for (int bj = 0; bj < nb; ++bj)
      for (int bi = 0; bi < nb; ++bi) {
        const std::int64_t ib = (static_cast<std::int64_t>(bk) * nb + bj) * nb + bi;
        bool launch = cls[ib] == 0;
        if (cls[ib] == 2) {
          int h = 1 + static_cast<int>(std::floor((std::log10(static_cast<double>(pray[ib])) - lmin) / ldel));
          h = std::min(std::max(h, 1), NH);
          launch = h >= hcut;
        }
        if (launch && cls[ib] == 2) nlaunch[bk] += 1;
        const std::uint8_t m = launch ? 1 : 0;
        for (int k = 0; k < b; ++k)
          for (int j = 0; j < b; ++j)
            for (int i = 0; i < b; ++i) Rc.mode[cell_in(bi, bj, bk, i, j, k)] = m;
      }

  const double norm = 1.0 / static_cast<double>(s.cfg.Nrtc);
  std::int64_t nl = 0;
  for (int bk = 0; bk < nb; ++bk) {
    Rc.launched += sum_launch[bk] * norm;
    Rc.dropped += sum_drop[bk] * norm;
    Rc.onspot += sum_ots[bk] * norm;
    nl += nlaunch[bk];
  }
  Rc.nemit_last = ne;
  Rc.nlaunch_last = nl;
  if (s.raylog)
    *s.raylog << " Recomb     : launched " << static_cast<float>(Rc.launched) << " dropped "
              << static_cast<float>(Rc.dropped) << " onspot " << static_cast<float>(Rc.onspot)
              << " emit " << ne << " launch_mode " << nl << " of " << ncand << " merge_nside "
              << s.rays.nside_merge << "\n";
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
  double nu[kNfreq], ngam[kNfreq];
  for (int n = 0; n < kNfreq; ++n) {
    const std::size_t k = static_cast<std::size_t>(c) * kNfreq + n;
    const double ng = Rc.emit_ngam[k];
    nu[n] = ng > 0 ? static_cast<double>(Rc.emit_egam[k]) / ng : cp::freq_init[n];
    ngam[n] = ng / nray;
    Rc.emit_ngam[k] = 0;
    Rc.emit_egam[k] = 0;
  }
  // 12 nside^2 rays from the cell centre with this spectrum (the kernels'
  // emission, as ray_source)
  const rrk::SourceSpec spec = rrk::source_spec_spectrum(irt - 1, jrt - 1, krt - 1, ns, nu, ngam);
  std::int64_t indx[rayp::nray_init];
  ray_create(s, nray, indx);
  for (int i = 0; i < nray; ++i) {
    const std::int64_t ir = indx[i];
    Ray& ry = R.at(ir);
    static_cast<rrk::RayP&>(ry) = rrk::source_ray(spec, i);
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
