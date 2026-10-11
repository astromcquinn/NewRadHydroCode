// Case A recombination and recombination radiation (recomb_rad.hpp,
// docs/recomb_radiation.md). New code, not in the original.
#include "radray/recomb_rad.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "radray/constants.hpp"
#include "radray/fortran.hpp"
#include "radray/healpix.hpp"
#include "radray/rt_rates.hpp"
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
    const int v = std::stoi(val);
    if (v < 0 || v > 1) throw std::runtime_error("recomb_rad must be 0 or 1");
    c.recomb_rad = v != 0;
    return true;
  }
  if (key == "recomb_rad_lines") {
    const int v = std::stoi(val);
    if (v < 0 || v > 1) throw std::runtime_error("recomb_rad_lines must be 0 or 1");
    c.recomb_rad_lines = v != 0;
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
  to.recomb_rad_lines = from.recomb_rad_lines;
  to.recomb_rad_nside = from.recomb_rad_nside;
  to.recomb_rad_block = from.recomb_rad_block;
  to.recomb_rad_pool = from.recomb_rad_pool;
}

void check_recomb_param(const Config& c) {
  if (c.recomb_rad_lines && !c.recomb_rad)
    throw std::runtime_error("recomb_rad_lines = 1 needs recomb_rad = 1");
  if (!c.recomb_rad) return;
  if (!c.recomb_case_A)
    throw std::runtime_error(
        "recomb_rad = 1 needs recomb_case = A: with case B the ground-state "
        "recombinations are already absorbed on the spot");
  // rt_method = m1 | otvet: every cell emits its own recombinations into the
  // moments (moment_step_1); the ray-only settings below are not used
  if (c.rt_method != RTMethod::Rays) return;
  if (c.legacy_rt) throw std::runtime_error("recomb_rad = 1 needs legacy_rt = 0 (the pool guard)");
  if (c.Nrt % c.recomb_rad_block != 0)
    throw std::runtime_error("recomb_rad_block must divide N");
}

// ---------------------------------------------------------------------------
// Photon energies
// ---------------------------------------------------------------------------

namespace {

// Group edges [eV] (the original's: thresholds and their geometric means)
const double kEdge[kNfreq + 1] = {pc::HI_eV, std::sqrt(pc::HI_eV * pc::HeI_eV), pc::HeI_eV,
                                  std::sqrt(pc::HeI_eV * pc::HeII_eV), pc::HeII_eV, 1e2};

int group_of(double e) {
  for (int g = 0; g < kNfreq; ++g)
    if (e < kEdge[g + 1]) return g;
  return kNfreq - 1;
}

// x e^x E1(x), E1 the exponential integral (asymptotic series for large x)
double xexE1(double x) {
  if (x > 50) return 1 - 1 / x + 2 / (x * x) - 6 / (x * x * x) + 24 / (x * x * x * x);
  return x * std::exp(x) * (-std::expint(-x));
}

// Two-photon decay 2s -> 1s of a hydrogenic ion: photons per unit y = E/E12,
// the fit of Nussbaumer & Schmutz (1984, A&A 138, 495) for HI, unnormalised
// (coefficients from memory, to be checked against the paper). Used for HeII
// and, approximately, for HeI 2^1S.
double two_photon_shape(double y) {
  const double w = y * (1 - y);
  const double q = std::pow(4 * w, 0.8);
  return w * (1 - q) + 0.88 * std::pow(w, 1.53) * q;
}

// Adds n photons per recombination of a two-photon decay of total energy
// E12 [eV] to the groups (only the part above 13.6 eV ionises).
void add_two_photon(double n, double E12, double* ln, double* le) {
  const int M = 4000;
  double tot = 0;
  double gn[kNfreq] = {}, ge[kNfreq] = {};
  for (int i = 0; i < M; ++i) {
    const double y = (i + 0.5) / M;
    const double a = two_photon_shape(y);
    tot += a;
    const double e = y * E12;
    if (e < pc::HI_eV) continue;
    const int g = group_of(e);
    gn[g] += a;
    ge[g] += a * e;
  }
  for (int g = 0; g < kNfreq; ++g) {  // 2 photons per decay
    ln[g] += 2 * n * gn[g] / tot;
    le[g] += 2 * n * ge[g] / tot;
  }
}

void add_line(double n, double e, double* ln, double* le) {
  const int g = group_of(e);
  ln[g] += n;
  le[g] += n * e;
}

// Fraction of the case-B recombinations HeIII -> HeII that capture directly
// to n = 2 (Balmer continuum, E >= 13.6 eV), with Kramers' partial
// recombination coefficients alpha_n ~ (1/n) x_n e^x_n E1(x_n),
// x_n = Z^2 chi_H / (n^2 kT) (no Gaunt factors: an estimate to 10-20%).
double heII_fBC(double T) {
  const double x1 = pc::HeII_eV / (pc::K2eV * T);
  double sum = 0, a2 = 0;
  for (int n = 2; n <= 2000; ++n) {
    const double a = xexE1(x1 / (static_cast<double>(n) * n)) / n;
    if (n == 2) a2 = a;
    sum += a;
  }
  return a2 / sum;
}

}  // namespace

double recomb_mean_excess(int sp, double T) {
  const double chi = sp == 0 ? pc::HI_eV : sp == 1 ? pc::HeI_eV : sp == 2 ? pc::HeII_eV : 0.25 * pc::HeII_eV;
  const double kT = pc::K2eV * T;
  // photons per dE ~ E^2 sigma(E) exp(-x), x = (E - chi)/kT, on x in [0, 60]
  const int M = 3000;
  const double dx = 60.0 / M;
  double num = 0, den = 0;
  for (int i0 = 0; i0 < M; i0 += kNfreq) {
    std::array<double, kNfreq> e{};
    for (int j = 0; j < kNfreq; ++j) e[j] = chi + (i0 + j + 0.5) * dx * kT;
    SigmaTable sig{};
    if (sp < 3) sig = photoion_sigma(e);
    for (int j = 0; j < kNfreq; ++j) {
      const double x = (i0 + j + 0.5) * dx;
      const double sg = sp < 3 ? sig[sp][j] : std::pow(chi / e[j], 3);
      const double w = e[j] * e[j] * sg * std::exp(-x);
      num += w * x;
      den += w;
    }
  }
  return den > 0 ? kT * num / den : kT;
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
  if (Rc.dE.empty()) {
    const int NT = RecombState::kNT;
    Rc.dE.assign(4 * NT, 0.0);
    Rc.fBC.assign(NT, 0.0);
    for (int it = 0; it < NT; ++it) {
      const double T = std::pow(10.0, RecombState::kLogT0 + it * RecombState::kDLogT);
      for (int sp = 0; sp < 4; ++sp) Rc.dE[sp * NT + it] = recomb_mean_excess(sp, T);
      Rc.fBC[it] = heII_fBC(T);
    }
    // recomb_rad_lines: H-ionising photons per case-B recombination. The
    // branching fractions are the usual low-density estimates (Osterbrock &
    // Ferland 2006, ch. 4; from memory, to be checked):
    // HeII -> HeI: 3/4 to triplets, all ending in 2^3S -> 1^1S (19.82 eV);
    //   1/4 to singlets, 2/3 of them via 2^1P -> 1^1S (21.22 eV) and 1/3
    //   via the 2^1S two-photon decay (20.62 eV in all).
    // HeIII -> HeII: 2/3 via 2P -> 1S (Lya, 40.81 eV), 1/3 via the 2S
    //   two-photon decay (40.81 eV in all); plus the Balmer continuum
    //   (fBC, added in recomb_count).
    for (auto* a : {&Rc.line_n[0][0], &Rc.line_e[0][0], &Rc.line_n[1][0], &Rc.line_e[1][0]})
      for (int g = 0; g < kNfreq; ++g) a[g] = 0;
    add_line(0.75, 19.82, Rc.line_n[0], Rc.line_e[0]);
    add_line(0.25 * 2.0 / 3.0, 21.22, Rc.line_n[0], Rc.line_e[0]);
    add_two_photon(0.25 / 3.0, 20.62, Rc.line_n[0], Rc.line_e[0]);
    add_line(2.0 / 3.0, 40.81, Rc.line_n[1], Rc.line_e[1]);
    add_two_photon(1.0 / 3.0, 40.81, Rc.line_n[1], Rc.line_e[1]);
  }
}

// photo_cell. Mean photon energy of a ground-state recombination: the
// threshold plus recomb_mean_excess(T) (about kT; Milne relation with the
// code's cross-sections), kept inside the group.
void recomb_count(Simulation& s, std::int64_t c, const double nrab[3], const double nrB[3],
                  double T, double vol_gas) {
  RecombState& Rc = *s.recomb;
  const Units& u = s.units;
  static constexpr int grp[3] = {0, 2, 4};
  const double nu0[3] = {pc::HI_eV, pc::HeI_eV, pc::HeII_eV};
  const bool onspot = recomb_onspot(Rc, c);
  // log T interpolation in the energy tables
  const int NT = RecombState::kNT;
  double ft = (std::log10(std::max(T, 1.0)) - RecombState::kLogT0) / RecombState::kDLogT;
  ft = std::min(std::max(ft, 0.0), static_cast<double>(NT - 1) - 1e-9);
  const int it = static_cast<int>(ft);
  const double wt = ft - it;
  auto tab = [&](const double* t) { return (1 - wt) * t[it] + wt * t[it + 1]; };
  const double conv = vol_gas / (u.nH_avg * u.vol_rt);
  auto add = [&](int g, double ng, double nu) {
    const std::size_t k = static_cast<std::size_t>(c) * kNfreq + g;
    Rc.ngam[k] = static_cast<float>(static_cast<double>(Rc.ngam[k]) + ng);
    Rc.egam[k] = static_cast<float>(static_cast<double>(Rc.egam[k]) + ng * nu);
  };
  for (int sp = 0; sp < 3; ++sp) {
    if (!(nrab[sp] > 0)) continue;
    const double ng = nrab[sp] * conv;
    if (onspot) {
      Rc.ots[c] = static_cast<float>(static_cast<double>(Rc.ots[c]) + ng);
      continue;
    }
    const double nu = std::min(nu0[sp] + tab(&Rc.dE[sp * NT]), kEdge[grp[sp] + 1]);
    add(grp[sp], ng, nu);
  }
  // recomb_rad_lines: the H-ionising helium photons of the case-B
  // recombinations (in the case-B cells they are counted on the spot: not
  // transported, as in the original)
  if (!nrB) return;
  for (int h = 0; h < 2; ++h) {
    const int sp = h + 1;
    if (!(nrB[sp] > 0)) continue;
    const double nb = nrB[sp] * conv;
    double tot = 0;
    for (int g = 0; g < kNfreq; ++g) tot += Rc.line_n[h][g];
    if (sp == 2) tot += tab(Rc.fBC.data());
    if (onspot) {
      Rc.ots[c] = static_cast<float>(static_cast<double>(Rc.ots[c]) + nb * tot);
      continue;
    }
    for (int g = 0; g < kNfreq; ++g)
      if (Rc.line_n[h][g] > 0) add(g, nb * Rc.line_n[h][g], Rc.line_e[h][g] / Rc.line_n[h][g]);
    if (sp == 2) {
      const double nu = std::min(0.25 * pc::HeII_eV + tab(&Rc.dE[3 * NT]), kEdge[1]);
      add(0, nb * tab(Rc.fBC.data()), nu);
    }
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
  // B16 fix (ray_rotate != 0): rotated HEALPix frame, as in ray_source
  const bool rotate = s.cfg.ray_rotate != 0;
  double rot[3][3];
  if (rotate) ray_frame_rotation(s, irt, jrt, krt, rot, 1);
  for (int i = 0; i < nray; ++i) {
    double dx[3];
    healpix::pix2vec_nest(ns, i, dx);
    int ipix = i;
    if (rotate) {
      const double v[3] = {dx[0], dx[1], dx[2]};
      for (int d = 0; d < 3; ++d) dx[d] = rot[d][0] * v[0] + rot[d][1] * v[1] + rot[d][2] * v[2];
      ipix = healpix::vec2pix_nest(ns, dx);
    }
    const std::int64_t ir = indx[i];
    Ray& ry = R.at(ir);
    ry.nside = ns;
    ry.ipix = ipix;
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
