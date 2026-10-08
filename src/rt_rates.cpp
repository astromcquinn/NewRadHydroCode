// Port of rt_rates.f90.
//
// Precision notes: exponents such as `x**1.503` are real(4) literals
// promoted to double, i.e. pow(x, (double)1.503f). Literals such as 0.522 in
// `(x/0.522)` are likewise (double)0.522f. `1.0+...` with 1.0 a real(4)
// literal is exact. `(8./3.)` is a real(4) division.
#include "radray/rt_rates.hpp"

#include <cmath>

#include "radray/fortran.hpp"

namespace radray {

namespace {
// single-precision literal promoted to double
constexpr double s(float v) { return static_cast<double>(v); }
}  // namespace

// rt_rates.f90: calc_atomic_rates, at one table temperature. The Fortran's
// two loops over i (rates, then cooling) are merged here; every value is
// computed by the same expression as there, so the tables are bit-identical.
RateSet atomic_rates_at(double T, bool legacy_physics) {
  using namespace pc;
  using namespace cp;
  const double gam = gamma_ideal;
  using std::exp;
  using std::pow;
  using std::sqrt;
  RateSet R{};

  {
    const double xHI = HI_K / T;
    const double xHeI = HeI_K / T;
    const double xHeII = HeII_K / T;
    double x;

    // HII recombination case A
    x = 2 * xHI;
    R.alpha[0] = 1.269e-13 * pow(x, s(1.503f)) /
                        pow(1.0 + pow(x / s(0.522f), s(0.470f)), s(1.923f));
    // HeII recombination case A
    x = 2 * xHeI;
    R.alpha[1] = 2.78e-14 * pow(x, s(0.681f));
    // HeIII recombination case A
    x = 2 * xHeII;
    R.alpha[2] = 2 * 1.269e-13 * pow(x, s(1.503f)) /
                        pow(1.0 + pow(x / s(0.522f), s(0.470f)), s(1.923f));
    // HII recombination case B
    x = 2 * xHI;
    R.alpha[3] = 2.753e-14 * pow(x, s(1.500f)) /
                        pow(1.0 + pow(x / s(2.740f), s(0.407f)), s(2.242f));
    // HeII recombination case B
    x = 2 * xHeI;
    R.alpha[4] = 1.26e-14 * pow(x, s(0.750f));
    // HeIII recombination case B
    x = 2 * xHeII;
    R.alpha[5] = 2 * 3.753e-14 * pow(x, s(1.500f)) /
                        pow(1.0 + pow(x / s(2.740f), s(0.407f)), s(2.242f));

    // legacy_physics = 0 (docs/port_notes.md 5b; harness PATCHES.md 9)
    if (!legacy_physics) {
      // P4: HeIII case B recombination. 2*3.753e-14 in the original is a
      // typo; Hui & Gnedin (1997, MNRAS 292, 27, Appendix A):
      // alpha_HeIII^B = 2*alpha_HII^B(lambda_HeII), coefficient 2.753e-14.
      x = 2 * xHeII;
      R.alpha[5] = 2 * 2.753e-14 * pow(x, s(1.500f)) /
                          pow(1.0 + pow(x / s(2.740f), s(0.407f)), s(2.242f));
      // P6: HeII dielectronic recombination, Aldrovandi & Pequignot (1973)
      // as given by Hui & Gnedin (1997, Appendix A):
      // 1.90e-3 T^-1.5 exp(-0.75 lambda_HeII/2) (1 + 0.3 exp(-0.15 lambda_HeII/2)),
      // added to the case A and case B radiative rates.
      double adi;
      if (xHeII > 100) {
        adi = 0;
      } else {
        x = 2 * xHeII;
        adi = 1.90e-3 / pow(T, 1.5) * exp(-0.375 * x) * (1 + 0.3 * exp(-0.075 * x));
      }
      // P7 (harness PATCHES.md 10): HeII case A radiative recombination,
      // Hui & Gnedin (1997, Appendix A; Burgess & Seaton fit):
      //   alpha_HeII^A = 3.0e-14 lambda_HeI^0.654 [cm^3/s],
      // replacing the original 2.78e-14 x^0.681 ("personal communication
      // with Nick"), so case A and case B (alpha_rec(5)) share one source.
      x = 2 * xHeI;
      R.alpha[1] = 3.0e-14 * pow(x, 0.654);
      R.alpha[1] = R.alpha[1] + adi;
      R.alpha[4] = R.alpha[4] + adi;
    }

    // HI collisional ionization
    if (xHI > 100) {
      R.gcol[0] = 0;
    } else {
      x = 2 * xHI;
      R.gcol[0] = 2.111e+01 / pow(T, s(1.5f)) * exp(-0.5 * x) /
                          pow(x, s(1.089f)) /
                          pow(1.0 + pow(x / s(0.354f), s(0.874f)), s(1.101f));
    }
    // HeI collisional ionization
    if (xHeI > 100) {
      R.gcol[1] = 0;
    } else {
      x = 2 * xHeI;
      R.gcol[1] = 3.238e+01 / pow(T, s(1.5f)) * exp(-0.5 * x) /
                          pow(x, s(1.146f)) /
                          pow(1.0 + pow(x / s(0.416f), s(0.987f)), s(1.056f));
    }
    // HeII collisional ionization
    if (xHeII > 100) {
      R.gcol[2] = 0;
    } else {
      x = 2 * xHeII;
      R.gcol[2] = 1.995e+01 / pow(T, s(1.5f)) * exp(-0.5 * x) /
                          pow(x, s(1.089f)) /
                          pow(1.0 + pow(x / s(0.553f), s(0.735f)), s(1.275f));
    }
  }

  {
    const double xHI = HI_K / T;
    const double xHeI = HeI_K / T;
    const double xHeII = HeII_K / T;
    double x;

    // HII recombination cooling case A
    x = 2 * xHI;
    R.cool[0] = 1.288e-13 * (k_cgs * T) * pow(x, s(1.965f)) /
                        pow(1.0 + pow(x / s(0.541f), s(0.502f)), s(2.697f));
    // HeII recombination cooling case A
    x = 2 * xHeI;
    R.cool[1] = 2.78e-14 * (k_cgs * T) * pow(x, s(0.681f));
    // HeIII recombination cooling case A
    x = 2 * xHeII;
    R.cool[2] = 8 * 1.288e-13 * (k_cgs * T) * pow(x, s(1.965f)) /
                        pow(1.0 + pow(x / s(0.541f), s(0.502f)), s(2.697f));
    // HII recombination cooling case B (same formula as case A in original)
    x = 2 * xHI;
    R.cool[3] = 1.288e-13 * (k_cgs * T) * pow(x, s(1.965f)) /
                        pow(1.0 + pow(x / s(0.541f), s(0.502f)), s(2.697f));
    // HeII recombination cooling case B
    x = 2 * xHeI;
    R.cool[4] = 2.78e-14 * (k_cgs * T) * pow(x, s(0.681f));
    // HeIII recombination cooling case B
    x = 2 * xHeII;
    R.cool[5] = 8 * 1.288e-13 * (k_cgs * T) * pow(x, s(1.965f)) /
                        pow(1.0 + pow(x / s(0.541f), s(0.502f)), s(2.697f));

    // HI collisional ionization cooling
    if (xHI > 100) {
      R.cool[6] = 0;
    } else {
      x = 2 * xHI;
      R.cool[6] = 2.111e+01 * (k_cgs * T) / pow(T, s(1.5f)) *
                          exp(-0.5 * x) / pow(x, s(1.089f)) /
                          pow(1.0 + pow(x / s(0.354f), s(0.874f)), s(1.101f));
    }
    // HeI collisional ionization cooling
    if (xHeI > 100) {
      R.cool[7] = 0;
    } else {
      x = 2 * xHeI;
      R.cool[7] = 3.238e+01 * (k_cgs * T) / pow(T, s(1.5f)) *
                          exp(-0.5 * x) / pow(x, s(1.146f)) /
                          pow(1.0 + pow(x / s(0.416f), s(0.987f)), s(1.056f));
    }
    // HeII collisional ionization cooling
    if (xHeII > 100) {
      R.cool[8] = 0;
    } else {
      x = 2 * xHeII;
      R.cool[8] = 1.995e+01 * (k_cgs * T) / pow(T, s(1.5f)) *
                          exp(-0.5 * x) / pow(x, s(1.089f)) /
                          pow(1.0 + pow(x / s(0.553f), s(0.735f)), s(1.275f));
    }
    // HI collisional excitation cooling
    if (xHI > 100) {
      R.cool[9] = 0;
    } else {
      x = 2 * xHI;
      R.cool[9] =
          7.50e-19 * exp(-0.375 * x) / (1.0 + sqrt(T / 1e5));
    }
    // HeI collisional excitation cooling (zero in original)
    R.cool[10] = 0;
    // HeII collisional excitation cooling
    if (xHeII > 100) {
      R.cool[11] = 0;
    } else {
      x = 2 * xHeII;
      R.cool[11] = 5.54e-17 * exp(-0.375 * x) / pow(T, s(0.397f)) /
                           (1.0 + sqrt(T / 1e5));
    }
    // Bremsstrahlung, Gaunt factor 1.30 (real(4) literal)
    R.cool[12] = 1.42e-27 * sqrt(T) * s(1.30f);
    // Inverse Compton: 4*sS_cgs/c_cgs*(8./3.)*3.36E-24/(gam-1)
    R.cool[13] =
        4 * sS_cgs / c_cgs * s(8.0f / 3.0f) * s(3.36e-24f) / (gam - 1);
    // coef_cool(15,:) never set in the original (stays 0, module array)

    // legacy_physics = 0 (docs/port_notes.md 5b; harness PATCHES.md 9):
    // overwrite the legacy values. radiative_physics then uses the corrected
    // cooling sum, in which coef_cool(11) and (15) enter.
    if (!legacy_physics) {
      // P1: case B recombination cooling, Hui & Gnedin (1997, Appendix A).
      // The original repeats the case A fits.
      //   HII  : 3.435e-30 T lambda_HI^1.970 / (1+(lambda_HI/2.250)^0.376)^3.720
      //   HeII : k T alpha_HeII^B = 1.26e-14 k T lambda_HeI^0.750 (radiative)
      //   HeIII: 8 x the HII fit at lambda_HeII
      x = 2 * xHI;
      R.cool[3] = 3.435e-30 * T * pow(x, 1.970) / pow(1 + pow(x / 2.250, 0.376), 3.720);
      x = 2 * xHeI;
      R.cool[4] = 1.26e-14 * (k_cgs * T) * pow(x, 0.750);

      // P7 (harness PATCHES.md 10): HeII case A recombination cooling, Hui &
      // Gnedin (1997, Appendix A): k T alpha_HeII^A = 3.0e-14 k T
      // lambda_HeI^0.654 [erg cm^3/s] (radiative part), replacing
      // 2.78e-14 k T x^0.681.
      x = 2 * xHeI;
      R.cool[1] = 3.0e-14 * (k_cgs * T) * pow(x, 0.654);

      // P8 (harness PATCHES.md 10): temperature-dependent free-free Gaunt
      // factor instead of the constant 1.30, Katz, Weinberg & Hernquist
      // (1996, ApJS 105, 19, Table 1; from Spitzer 1978, as in Cen 1992):
      //   g_ff = 1.1 + 0.34 exp(-(5.5 - log10 T)^2 / 3)   (1.1 to 1.44).
      // Rate 1.42e-27 g_ff T^0.5 [erg cm^3/s], times ne*nff as before.
      // (-(d*d))/3 is the Fortran -(d)**2/3D0 exactly.
      {
        const double d = 5.5 - std::log10(T);
        R.cool[12] = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(d * d) / 3.0));
      }

      x = 2 * xHeII;
      R.cool[5] =
          8 * 3.435e-30 * T * pow(x, 1.970) / pow(1 + pow(x / 2.250, 0.376), 3.720);

      // P5: collisional ionisation cooling = ionisation potential x rate,
      // k T_X Gamma_X (Hui & Gnedin 1997, Appendix A; Cen 1992). The
      // original uses k T, which underestimates it by T_X/T.
      R.cool[6] = (k_cgs * HI_K) * R.gcol[0];
      R.cool[7] = (k_cgs * HeI_K) * R.gcol[1];
      R.cool[8] = (k_cgs * HeII_K) * R.gcol[2];

      // P2: HeI (2^3S) collisional excitation cooling, Cen (1992, ApJS 78,
      // 341): 9.10e-27 T^-0.1687 (1+sqrt(T/1e5))^-1 exp(-13179/T)
      // [erg cm^6/s], multiplies ne^2 nHeII.
      if (13179 / T > 100) {
        R.cool[10] = 0;
      } else {
        R.cool[10] =
            9.10e-27 / pow(T, 0.1687) * exp(-13179 / T) / (1 + sqrt(T / 1e5));
      }

      // P3: coef_cool(15) = HeII dielectronic recombination cooling,
      // 0.75 k T_HeII alpha_d (Hui & Gnedin 1997, Appendix A; = Cen 1992's
      // 1.24e-13 T^-1.5 exp(-470000/T)(1+0.3 exp(-94000/T))), times ne*nHeII.
      if (xHeII > 100) {
        R.cool[14] = 0;
      } else {
        x = 2 * xHeII;
        const double adi =
            1.90e-3 / pow(T, 1.5) * exp(-0.375 * x) * (1 + 0.3 * exp(-0.075 * x));
        R.cool[14] = 0.75 * (k_cgs * HeII_K) * adi;
      }
    }
  }
  return R;
}

void calc_atomic_rates(AtomicRates& R, bool legacy_physics) {
  for (int i = 1; i <= cp::Ntemp; ++i) {
    const double T = std::pow(1e1, (i - 1) * cp::lgT_del + cp::lgT_min);
    const RateSet r = atomic_rates_at(T, legacy_physics);
    for (int n = 1; n <= 6; ++n) R.alpha_rec(n, i) = r.alpha[n - 1];
    for (int n = 1; n <= 3; ++n) R.gamma_col(n, i) = r.gcol[n - 1];
    for (int n = 1; n <= 15; ++n) R.coef_cool(n, i) = r.cool[n - 1];
  }
}

RateSet rates_node(const AtomicRates& R, int iT) {
  RateSet r;
  for (int n = 1; n <= 6; ++n) r.alpha[n - 1] = R.alpha_rec(n, iT);
  for (int n = 1; n <= 3; ++n) r.gcol[n - 1] = R.gamma_col(n, iT);
  for (int n = 1; n <= 15; ++n) r.cool[n - 1] = R.coef_cool(n, iT);
  return r;
}

// P9 (legacy_physics = 0; harness PATCHES.md 10, rt_rates.f90 rate_lerp).
// Linear interpolation of log R in log T: R = R_i exp(w log(R_{i+1}/R_i)).
// It is exact for power laws (bremsstrahlung T^0.5, the recombination fits
// away from their breaks). For a Boltzmann factor exp(-T_X/T) the relative
// error at mid-interval is (lgT_del ln10)^2/8 * T_X/T, i.e. 7.8e-4 for HI
// excitation at 1e4 K, about 12x smaller than linear interpolation of R
// (whose error is dominated by (d ln R/d log T)^2). w = 0 gives a exactly.
// If either node is 0 (below the lambda/2 > 100 cut-offs), linear in R.
double rate_lerp(double a, double b, double w) {
  if (a > 0 && b > 0) return a * std::exp(w * std::log(b / a));
  return a + w * (b - a);
}

RateSet rates_interp(const AtomicRates& R, int i, double w) {
  RateSet r;
  for (int n = 1; n <= 6; ++n) r.alpha[n - 1] = rate_lerp(R.alpha_rec(n, i), R.alpha_rec(n, i + 1), w);
  for (int n = 1; n <= 3; ++n) r.gcol[n - 1] = rate_lerp(R.gamma_col(n, i), R.gamma_col(n, i + 1), w);
  for (int n = 1; n <= 15; ++n)
    r.cool[n - 1] = rate_lerp(R.coef_cool(n, i), R.coef_cool(n, i + 1), w);
  return r;
}

// Harness rt_rates.f90 rates_interp: u = (log10 T - lgT_min)/lgT_del,
// i = 1 + int(u), w = u - (i-1). Outside the table (u <= 0 or
// u >= Ntemp-1, including NaN) the edge node is used, no extrapolation.
RateSet rates_interp(const AtomicRates& R, double tmp) {
  const double u = (std::log10(tmp) - cp::lgT_min) / cp::lgT_del;
  if (u > 0 && u < cp::Ntemp - 1) {
    const int i = 1 + static_cast<int>(u);
    const double w = u - (i - 1);
    return rates_interp(R, i, w);
  }
  return rates_node(R, u > 0 ? cp::Ntemp : 1);
}

// rt_rates.f90: photoion_sigma (Verner et al. 1996)
SigmaTable photoion_sigma(const std::array<double, kNfreq>& nu) {
  using namespace pc;
  SigmaTable sig{};
  for (int i = 0; i < kNfreq; ++i) {
    const double fq = nu[i];
    double x, y;
    // HI
    if (fq < HI_eV) {
      sig[0][i] = 0;
    } else {
      x = fq / s(0.4298f);
      y = x;
      sig[0][i] = 5.475e-14 * f::pow2(x - 1) / std::pow(y, s(4.0185f)) /
                  std::pow(1 + std::sqrt(y / s(32.88f)), s(2.963f));
    }
    // HeI: 2.136**2 and 2.039**2 are real(4) products
    if (fq < HeI_eV) {
      sig[1][i] = 0;
    } else {
      x = fq / s(13.61f) - s(0.4434f);
      y = std::sqrt(f::pow2(x) + s(2.136f * 2.136f));
      sig[1][i] = 9.492e-16 * (f::pow2(x - 1) + s(2.039f * 2.039f)) /
                  std::pow(y, s(3.9060f)) /
                  std::pow(1 + std::sqrt(y / s(1.469f)), s(3.188f));
    }
    // HeII
    if (fq < HeII_eV) {
      sig[2][i] = 0;
    } else {
      x = fq / s(1.720f);
      y = x;
      sig[2][i] = 1.369e-14 * f::pow2(x - 1) / std::pow(y, s(4.0185f)) /
                  std::pow(1 + std::sqrt(y / s(32.88f)), s(2.963f));
    }
  }
  return sig;
}

}  // namespace radray
