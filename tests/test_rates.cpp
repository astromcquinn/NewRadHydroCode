// Checks the legacy_physics = 0 atomic rates (docs/port_notes.md section 5b)
// against independent evaluations of the published fits, and the
// interpolated rate lookup (P9).
//
// usage: test_rates [--csv FILE]
//
// For every coefficient changed by P1-P8 and a set of table temperatures it
// prints the table value, the published fit evaluated here from the paper's
// own constants, and their ratio, and checks:
//   (a) the transcription: the published fit evaluated with the code's
//       threshold temperatures (HI_K, HeI_K, HeII_K: 13.6/24.59/54.42 eV) and
//       k_cgs agrees with the table to 1e-12 (1e-6 where the table keeps the
//       original's single-precision exponent literals);
//   (b) the published constants: the same fit with the paper's threshold
//       temperatures (Hui & Gnedin 1997: 157807, 285335, 631515 K; CODATA k)
//       agrees within the tolerance listed per row, which comes only from
//       the threshold temperatures differing by up to 1.4e-4 relative;
//   (c) the dielectronic terms against Cen (1992)'s form (470000 K, 94000 K,
//       1.24e-13) at T >= 1e5 K, where it matters, to 5%;
//   (d) legacy_physics = 1 tables are unchanged by the switch, and every
//       coefficient that P1-P8 do not touch is identical in both modes;
//   (e) atomic_rates_at(T) (the formulas at any T) reproduces every table
//       entry exactly at the nodes, in both modes;
//   (f) the P9 interpolation: at a node (w = 0) it returns the node values
//       exactly; at mid-interval, for every slot and every interval of the
//       table where the rate is positive, it matches the direct formula
//       evaluation at that temperature to within the expected error of
//       linear interpolation of ln R in log T,
//         ln(interp/formula) = +h^2/8 d^2 lnR/du^2  (u = log10 T, h = lgT_del),
//       estimated from the table's second differences; it is far more
//       accurate than the nearest-node lookup; zeros, the edges, and the
//       T -> (i, w) mapping are checked separately.
// --csv writes logT and all coef_cool, alpha_rec, gamma_col entries in both
// modes (used by tools/plot_cooling.py for docs/figures/cooling_rates.png).
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

#include "check.hpp"
#include "radray/constants.hpp"
#include "radray/rt_rates.hpp"

using namespace radray;

namespace {

// Threshold temperatures and k for evaluating a fit.
struct Consts {
  double THI, THeI, THeII, k;
};
// Hui & Gnedin (1997, MNRAS 292, 27), Appendix A; k from CODATA 2018.
const Consts kPaper{157807., 285335., 631515., 1.380649e-16};
// The code's values (physical_constants.f90).
const Consts kCode{pc::HI_K, pc::HeI_K, pc::HeII_K, pc::k_cgs};

// ---- Published fits, written out as printed in the papers ----
// lambda_X = 2 T_X / T.
// HG97: case B HII recombination cooling [erg cm^3 s^-1]
double hg97_recool_B_HII(double T, const Consts& c) {
  const double l = 2 * c.THI / T;
  return 3.435e-30 * T * std::pow(l, 1.970) / std::pow(1 + std::pow(l / 2.250, 0.376), 3.720);
}
// HG97: case B HeII recombination cooling = k T alpha_HeII^B
double hg97_recool_B_HeII(double T, const Consts& c) {
  const double l = 2 * c.THeI / T;
  return c.k * T * 1.26e-14 * std::pow(l, 0.750);
}
// HG97: case B HeIII recombination cooling = 8 x HII fit at lambda_HeII
double hg97_recool_B_HeIII(double T, const Consts& c) {
  const double l = 2 * c.THeII / T;
  return 8 * 3.435e-30 * T * std::pow(l, 1.970) / std::pow(1 + std::pow(l / 2.250, 0.376), 3.720);
}
// HG97: case B HeIII recombination = 2 x alpha_HII^B(lambda_HeII) [cm^3 s^-1]
double hg97_alpha_B_HeIII(double T, const Consts& c) {
  const double l = 2 * c.THeII / T;
  return 2 * 2.753e-14 * std::pow(l, 1.500) / std::pow(1 + std::pow(l / 2.740, 0.407), 2.242);
}
// HG97: collisional ionisation rates (Lotz 1967 fits) [cm^3 s^-1]
double hg97_ci(double T, double TX, double A, double a, double b, double c, double d) {
  const double l = 2 * TX / T;
  return A * std::pow(T, -1.5) * std::exp(-l / 2) * std::pow(l, -a) /
         std::pow(1 + std::pow(l / b, c), d);
}
double hg97_ci_HI(double T, const Consts& c) {
  return hg97_ci(T, c.THI, 21.11, 1.089, 0.354, 0.874, 1.101);
}
double hg97_ci_HeI(double T, const Consts& c) {
  return hg97_ci(T, c.THeI, 32.38, 1.146, 0.416, 0.987, 1.056);
}
double hg97_ci_HeII(double T, const Consts& c) {
  return hg97_ci(T, c.THeII, 19.95, 1.089, 0.553, 0.735, 1.275);
}
// HG97: collisional ionisation cooling = k T_X Gamma_X
double hg97_cicool_HI(double T, const Consts& c) { return c.k * c.THI * hg97_ci_HI(T, c); }
double hg97_cicool_HeI(double T, const Consts& c) { return c.k * c.THeI * hg97_ci_HeI(T, c); }
double hg97_cicool_HeII(double T, const Consts& c) { return c.k * c.THeII * hg97_ci_HeII(T, c); }
// HG97 (Aldrovandi & Pequignot 1973): HeII dielectronic recombination
double hg97_alpha_di(double T, const Consts& c) {
  const double l = 2 * c.THeII / T;
  return 1.90e-3 * std::pow(T, -1.5) * std::exp(-0.75 * l / 2) *
         (1 + 0.3 * std::exp(-0.15 * l / 2));
}
double hg97_dicool(double T, const Consts& c) { return 0.75 * c.k * c.THeII * hg97_alpha_di(T, c); }
// Glover & Jappsen (2007, ApJ 666, 1, Table 7): He(1^1S) collisional
// excitation cooling [erg cm^3 s^-1], times ne*nHeI (P10)
double gj07_ceHeI_ground(double T) {
  return 1.1e-19 * std::pow(T, 0.082) * std::exp(-230000. / T);
}
// Cen (1992): the same two processes in Cen's form
double cen92_alpha_di(double T) {
  return 1.9e-3 * std::pow(T, -1.5) * std::exp(-470000. / T) * (1 + 0.3 * std::exp(-94000. / T));
}
double cen92_dicool(double T) {
  return 1.24e-13 * std::pow(T, -1.5) * std::exp(-470000. / T) * (1 + 0.3 * std::exp(-94000. / T));
}
// Cen (1992): He(2^3S) collisional excitation cooling [erg cm^6 s^-1],
// multiplies ne^2 nHeII
double cen92_ceHeI(double T) {
  return 9.10e-27 * std::pow(T, -0.1687) / (1 + std::sqrt(T / 1e5)) * std::exp(-13179. / T);
}
// HG97: HeII radiative recombination (case A/B, Burgess & Seaton fits)
double hg97_alpha_B_HeII(double T, const Consts& c) {
  return 1.26e-14 * std::pow(2 * c.THeI / T, 0.750);
}
double hg97_alpha_A_HeII(double T, const Consts& c) {
  return 3.0e-14 * std::pow(2 * c.THeI / T, 0.654);
}
// HG97: HeII case A recombination cooling = k T alpha_HeII^A
double hg97_recool_A_HeII(double T, const Consts& c) { return c.k * T * hg97_alpha_A_HeII(T, c); }
// Katz, Weinberg & Hernquist (1996), Table 1: bremsstrahlung with
// g_ff = 1.1 + 0.34 exp(-(5.5 - log10 T)^2/3), per ne*(nHII+nHeII+4nHeIII)
double kwh96_gff(double T) {
  return 1.1 + 0.34 * std::exp(-std::pow(5.5 - std::log10(T), 2) / 3);
}
double kwh96_brems(double T, const Consts&) { return 1.42e-27 * std::sqrt(T) * kwh96_gff(T); }

double table_T(int i) { return std::pow(1e1, (i - 1) * cp::lgT_del + cp::lgT_min); }

struct Row {
  const char* name;
  std::function<double(const AtomicRates&, int)> table;  // new table value
  std::function<double(const AtomicRates&, int)> legacy;  // legacy table value
  std::function<double(double, const Consts&)> fit;
  double rtol_paper;  // (b): paper constants vs code constants
  // (a): 1e-12, or 1e-6 where the code keeps the original's real(4)
  // exponent literals (gamma_col, alpha_rec(6): pow(x, (double)1.089f))
  double rtol_code = 1e-12;
};

}  // namespace

int main(int argc, char** argv) {
  std::string csv;
  for (int a = 1; a < argc; ++a)
    if (!std::strcmp(argv[a], "--csv") && a + 1 < argc) csv = argv[++a];

  AtomicRates L, N, L2;
  calc_atomic_rates(L);         // default: legacy
  calc_atomic_rates(N, false);  // legacy_physics = 0
  calc_atomic_rates(L2, true);

  // (d) legacy tables: default == explicit legacy; untouched slots equal.
  CHECK(L.coef_cool_ == L2.coef_cool_ && L.alpha_rec_ == L2.alpha_rec_ &&
        L.gamma_col_ == L2.gamma_col_);
  CHECK(N.gamma_col_ == L.gamma_col_);
  for (int i = 1; i <= cp::Ntemp; ++i) {
    for (int n : {1, 3, 10, 12, 14}) CHECK(N.coef_cool(n, i) == L.coef_cool(n, i));
    for (int n : {1, 3, 4}) CHECK(N.alpha_rec(n, i) == L.alpha_rec(n, i));
    CHECK(L.coef_cool(11, i) == 0 && L.coef_cool(15, i) == 0 && L.coef_cool(16, i) == 0);
  }

  const auto di = [](double T, const Consts& c) { return hg97_alpha_di(T, c); };
  // The tolerances (b) cover exp(-lambda/2)-type factors evaluated with
  // T_X differing by <= 1.4e-4 relative (largest at the lowest T).
  const std::vector<Row> rows = {
      {"P1 coef_cool(4)  HII  rec. cooling B",
       [](const AtomicRates& R, int i) { return R.coef_cool(4, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(4, i); }, hg97_recool_B_HII, 1e-3},
      {"P1 coef_cool(5)  HeII rec. cooling B",
       [](const AtomicRates& R, int i) { return R.coef_cool(5, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(5, i); }, hg97_recool_B_HeII, 1e-3},
      {"P1 coef_cool(6)  HeIII rec. cooling B",
       [](const AtomicRates& R, int i) { return R.coef_cool(6, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(6, i); }, hg97_recool_B_HeIII, 1e-3},
      {"P5 coef_cool(7)  HI  coll. ion. cooling",
       [](const AtomicRates& R, int i) { return R.coef_cool(7, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(7, i); }, hg97_cicool_HI, 3e-3, 1e-6},
      {"P5 coef_cool(8)  HeI coll. ion. cooling",
       [](const AtomicRates& R, int i) { return R.coef_cool(8, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(8, i); }, hg97_cicool_HeI, 6e-3, 1e-6},
      {"P5 coef_cool(9)  HeII coll. ion. cooling",
       [](const AtomicRates& R, int i) { return R.coef_cool(9, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(9, i); }, hg97_cicool_HeII, 1e-2, 1e-6},
      {"P2 coef_cool(11) HeI coll. exc. [cm^6]",
       [](const AtomicRates& R, int i) { return R.coef_cool(11, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(11, i); },
       [](double T, const Consts&) { return cen92_ceHeI(T); }, 1e-12},
      {"P3 coef_cool(15) HeII dielec. cooling",
       [](const AtomicRates& R, int i) { return R.coef_cool(15, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(15, i); }, hg97_dicool, 1e-2},
      {"P10 coef_cool(16) HeI(1^1S) coll. exc.",
       [](const AtomicRates& R, int i) { return R.coef_cool(16, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(16, i); },
       [](double T, const Consts&) { return gj07_ceHeI_ground(T); }, 1e-12},
      {"P4 alpha_rec(6)  HeIII rec. B",
       [](const AtomicRates& R, int i) { return R.alpha_rec(6, i); },
       [](const AtomicRates& R, int i) { return R.alpha_rec(6, i); }, hg97_alpha_B_HeIII, 1e-3, 1e-6},
      {"P6 alpha_rec(5)  HeII rec. B (+dielec.)",
       [](const AtomicRates& R, int i) { return R.alpha_rec(5, i); },
       [](const AtomicRates& R, int i) { return R.alpha_rec(5, i); },
       [di](double T, const Consts& c) { return hg97_alpha_B_HeII(T, c) + di(T, c); }, 1e-2},
      {"P7 alpha_rec(2)  HeII rec. A (+dielec.)",
       [](const AtomicRates& R, int i) { return R.alpha_rec(2, i); },
       [](const AtomicRates& R, int i) { return R.alpha_rec(2, i); },
       [di](double T, const Consts& c) { return hg97_alpha_A_HeII(T, c) + di(T, c); }, 1e-2},
      {"P7 coef_cool(2)  HeII rec. cooling A",
       [](const AtomicRates& R, int i) { return R.coef_cool(2, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(2, i); }, hg97_recool_A_HeII, 1e-3},
      {"P8 coef_cool(13) bremsstrahlung (g_ff(T))",
       [](const AtomicRates& R, int i) { return R.coef_cool(13, i); },
       [](const AtomicRates& R, int i) { return R.coef_cool(13, i); }, kwh96_brems, 1e-12},
  };

  // Table temperatures: log10 T = 4, 4.3, 4.5, 5, 5.5, 6, 7, 8. (Below
  // 1e4 K the coefficients that change are either cut off by the code's
  // lambda/2 > 100 guards or negligible, and the exp(-T_X/T) factors make
  // the comparison with the paper's T_X meaningless.)
  const int idx[] = {401, 431, 451, 501, 551, 601, 701, 801};

  std::printf("Columns: T [K]; table = legacy_physics=0 table value; fit = published fit\n"
              "with the paper's constants; table/fit; legacy = legacy_physics=1 value;\n"
              "table/legacy.\n");
  for (const Row& r : rows) {
    std::printf("\n%s\n", r.name);
    std::printf("  %10s %13s %13s %10s %13s %12s\n", "T", "table", "fit", "table/fit", "legacy",
                "table/legacy");
    for (int i : idx) {
      const double T = table_T(i);
      const double t = r.table(N, i), lg = r.legacy(L, i);
      const double f = r.fit(T, kPaper);
      // (a) transcription, with the code's constants
      CHECK_CLOSE(t, r.fit(T, kCode), r.rtol_code);
      // (b) with the paper's constants
      CHECK_CLOSE(t, f, r.rtol_paper);
      std::printf("  %10.4g %13.6e %13.6e %10.6f %13.6e %12.5g\n", T, t, f, f != 0 ? t / f : 1.0,
                  lg, lg != 0 ? t / lg : NAN);
    }
  }

  // (c) dielectronic terms vs Cen (1992)'s constants, where they matter
  std::printf("\nP3/P6 dielectronic vs Cen (1992) form (470000 K, 94000 K, 1.24e-13)\n");
  std::printf("  %10s %13s %13s %10s %13s %13s %10s\n", "T", "alpha_d", "Cen", "ratio",
              "cool_d", "Cen", "ratio");
  for (int i : {451, 501, 551, 601, 701}) {
    const double T = table_T(i);
    const double ad = N.alpha_rec(5, i) - hg97_alpha_B_HeII(T, kCode);
    const double cd = N.coef_cool(15, i);
    std::printf("  %10.4g %13.6e %13.6e %10.6f %13.6e %13.6e %10.6f\n", T, ad, cen92_alpha_di(T),
                ad / cen92_alpha_di(T), cd, cen92_dicool(T), cd / cen92_dicool(T));
    if (T >= 9.9e4) {
      CHECK_CLOSE(ad, cen92_alpha_di(T), 5e-2);
      CHECK_CLOSE(cd, cen92_dicool(T), 5e-2);
    }
  }

  // (e) atomic_rates_at reproduces the tables exactly at the nodes
  for (const AtomicRates* R : {&L, &N}) {
    const bool leg = (R == &L);
    for (int i = 1; i <= cp::Ntemp; ++i) {
      const RateSet a = atomic_rates_at(table_T(i), leg), t = rates_node(*R, i);
      CHECK(a.alpha == t.alpha && a.gcol == t.gcol && a.cool == t.cool);
    }
  }

  // (f) P9 interpolation
  // All 24 slots of a RateSet as one flat list.
  struct Slot {
    std::string name;
    std::function<double(const RateSet&)> get;
  };
  std::vector<Slot> slots;
  for (int n = 1; n <= 6; ++n)
    slots.push_back({"alpha_rec(" + std::to_string(n) + ")",
                     [n](const RateSet& r) { return r.alpha[n - 1]; }});
  for (int n = 1; n <= 3; ++n)
    slots.push_back({"gamma_col(" + std::to_string(n) + ")",
                     [n](const RateSet& r) { return r.gcol[n - 1]; }});
  for (int n = 1; n <= 16; ++n)
    slots.push_back({"coef_cool(" + std::to_string(n) + ")",
                     [n](const RateSet& r) { return r.cool[n - 1]; }});

  // Node values are returned exactly (w = 0), in both modes.
  for (const AtomicRates* R : {&L, &N})
    for (int i = 1; i < cp::Ntemp; ++i) {
      const RateSet a = rates_interp(*R, i, 0.0), t = rates_node(*R, i);
      CHECK(a.alpha == t.alpha && a.gcol == t.gcol && a.cool == t.cool);
    }
  // Lookup by temperature: T at a node, and T -> (i, w) inside an interval.
  for (int i = 2; i < cp::Ntemp; i += 7) {
    const RateSet a = rates_interp(N, table_T(i)), t = rates_node(N, i);
    for (const Slot& sl : slots) CHECK_CLOSE(sl.get(a), sl.get(t), 1e-12);
    const double w = 0.3;
    const double T = std::pow(1e1, (i - 1 + w) * cp::lgT_del + cp::lgT_min);
    const RateSet b = rates_interp(N, T), c = rates_interp(N, i, w);
    for (const Slot& sl : slots) CHECK_CLOSE(sl.get(b), sl.get(c), 1e-12);
  }
  // Table edges: the edge node, no extrapolation (also for NaN).
  {
    const RateSet lo = rates_node(N, 1), hi = rates_node(N, cp::Ntemp);
    for (double T : {1e-3, 0.5, 1.0}) {
      const RateSet r = rates_interp(N, T);
      CHECK(r.alpha == lo.alpha && r.gcol == lo.gcol && r.cool == lo.cool);
    }
    for (double T : {table_T(cp::Ntemp), 1e10, 1e12}) {
      const RateSet r = rates_interp(N, T);
      CHECK(r.alpha == hi.alpha && r.gcol == hi.gcol && r.cool == hi.cool);
    }
    const RateSet r = rates_interp(N, std::nan(""));
    CHECK(r.alpha == lo.alpha && r.gcol == lo.gcol && r.cool == lo.cool);
    // just below the last node: interpolated in the last interval
    const double Tm = std::pow(1e1, (cp::Ntemp - 1.5) * cp::lgT_del + cp::lgT_min);
    const RateSet m = rates_interp(N, Tm), mi = rates_interp(N, cp::Ntemp - 1, 0.5);
    for (const Slot& sl : slots) CHECK_CLOSE(sl.get(m), sl.get(mi), 1e-12);
  }
  // Zero nodes (below the lambda/2 > 100 cut-offs): linear in R.
  CHECK(rate_lerp(0.0, 2.0, 0.25) == 0.5 && rate_lerp(2.0, 0.0, 0.25) == 1.5);
  CHECK(rate_lerp(0.0, 0.0, 0.7) == 0.0);
  CHECK(rate_lerp(3.0, 3.0, 0.4) == 3.0);  // constants stay exact
  {
    int nzero = 0;
    for (int i = 1; i < cp::Ntemp; ++i)
      for (int n = 1; n <= 3; ++n) {
        const double a = N.gamma_col(n, i), b = N.gamma_col(n, i + 1);
        if (a == 0 && b > 0) {
          ++nzero;
          CHECK(rates_interp(N, i, 0.5).gcol[n - 1] == 0.5 * b);
        }
      }
    CHECK(nzero == 3);
  }

  // Mid-interval accuracy against the direct formula evaluation.
  std::printf("\nP9 interpolation (legacy_physics = 0 tables), relative errors vs the formula\n"
              "at mid-interval, 1e4 K <= T <= 1e8 K: max |interp - formula|, max |predicted|\n"
              "= h^2/8 |d^2 lnR/du^2| from the table, and the nearest-node lookup's max error.\n");
  std::printf("  %-14s %12s %12s %12s %8s\n", "slot", "interp", "predicted", "nearest", "checked");
  const double h = cp::lgT_del;
  for (const Slot& sl : slots) {
    double emax = 0, pmax = 0, nmax = 0;
    int nchk = 0;
    for (int i = 2; i <= cp::Ntemp - 2; ++i) {
      const double r0 = sl.get(rates_node(N, i - 1)), r1 = sl.get(rates_node(N, i)),
                   r2 = sl.get(rates_node(N, i + 1)), r3 = sl.get(rates_node(N, i + 2));
      if (!(r0 > 0 && r1 > 0 && r2 > 0 && r3 > 0)) continue;
      const double Tmid = std::pow(1e1, (i - 0.5) * h + cp::lgT_min);
      const double direct = sl.get(atomic_rates_at(Tmid, false));
      const double interp = sl.get(rates_interp(N, i, 0.5));
      const double e = std::log(interp / direct);
      // ln R second differences at nodes i and i+1, averaged to the midpoint
      const double D = (std::log(r0) - 2 * std::log(r1) + std::log(r2) + std::log(r1) -
                        2 * std::log(r2) + std::log(r3)) / 2;
      const double pred = D / 8;  // chord minus curve: +h^2/8 f''
      // 10% of the prediction for the higher-order terms (they reach 6% where a
      // steep term switches on, e.g. HeII dielectronic recombination at
      // 1.6e4 K inside alpha_rec(2)), plus round-off of the logs and pows
      const double tol = 0.1 * std::fabs(pred) + 2e-13 * (1 + std::fabs(std::log(r1)) / 50);
      CHECK(std::fabs(e - pred) <= tol);
      if (std::fabs(e - pred) > tol)
        std::fprintf(stderr, "  %s i=%d T=%g err=%.3e pred=%.3e\n", sl.name.c_str(), i, Tmid, e,
                     pred);
      ++nchk;
      if (Tmid >= 1e4 && Tmid <= 1e8) {
        emax = std::fmax(emax, std::fabs(e));
        pmax = std::fmax(pmax, std::fabs(pred));
        // nearest node: up to half an interval away, either side
        nmax = std::fmax(nmax, std::fabs(std::log(r1 / direct)));
        nmax = std::fmax(nmax, std::fabs(std::log(r2 / direct)));
      }
    }
    std::printf("  %-14s %12.3e %12.3e %12.3e %8d\n", sl.name.c_str(), emax, pmax, nmax, nchk);
    // everything that is positive above ~1e4 K is checked over the whole table
    CHECK(nchk > 400);
    if (nmax > 1e-12) CHECK(emax < 0.1 * nmax);
  }

  if (!csv.empty()) {
    FILE* fp = std::fopen(csv.c_str(), "w");
    CHECK(fp != nullptr);
    if (fp) {
      std::fprintf(fp, "logT");
      for (const char* m : {"L", "N"}) {
        for (int n = 1; n <= 16; ++n) std::fprintf(fp, ",%s_cool%d", m, n);
        for (int n = 1; n <= 6; ++n) std::fprintf(fp, ",%s_alpha%d", m, n);
        for (int n = 1; n <= 3; ++n) std::fprintf(fp, ",%s_gcol%d", m, n);
      }
      std::fprintf(fp, "\n");
      for (int i = 1; i <= cp::Ntemp; ++i) {
        std::fprintf(fp, "%.2f", (i - 1) * cp::lgT_del + cp::lgT_min);
        for (const AtomicRates* R : {&L, &N}) {
          for (int n = 1; n <= 16; ++n) std::fprintf(fp, ",%.9e", R->coef_cool(n, i));
          for (int n = 1; n <= 6; ++n) std::fprintf(fp, ",%.9e", R->alpha_rec(n, i));
          for (int n = 1; n <= 3; ++n) std::fprintf(fp, ",%.9e", R->gamma_col(n, i));
        }
        std::fprintf(fp, "\n");
      }
      std::fclose(fp);
    }
  }

  TEST_MAIN_END();
}
