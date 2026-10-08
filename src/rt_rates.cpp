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

// rt_rates.f90: calc_atomic_rates
void calc_atomic_rates(AtomicRates& R, bool legacy_physics) {
  using namespace pc;
  using namespace cp;
  const double gam = gamma_ideal;
  using std::exp;
  using std::pow;
  using std::sqrt;

  for (int i = 1; i <= Ntemp; ++i) {
    const double T = pow(1e1, (i - 1) * lgT_del + lgT_min);
    const double xHI = HI_K / T;
    const double xHeI = HeI_K / T;
    const double xHeII = HeII_K / T;
    double x;

    // HII recombination case A
    x = 2 * xHI;
    R.alpha_rec(1, i) = 1.269e-13 * pow(x, s(1.503f)) /
                        pow(1.0 + pow(x / s(0.522f), s(0.470f)), s(1.923f));
    // HeII recombination case A
    x = 2 * xHeI;
    R.alpha_rec(2, i) = 2.78e-14 * pow(x, s(0.681f));
    // HeIII recombination case A
    x = 2 * xHeII;
    R.alpha_rec(3, i) = 2 * 1.269e-13 * pow(x, s(1.503f)) /
                        pow(1.0 + pow(x / s(0.522f), s(0.470f)), s(1.923f));
    // HII recombination case B
    x = 2 * xHI;
    R.alpha_rec(4, i) = 2.753e-14 * pow(x, s(1.500f)) /
                        pow(1.0 + pow(x / s(2.740f), s(0.407f)), s(2.242f));
    // HeII recombination case B
    x = 2 * xHeI;
    R.alpha_rec(5, i) = 1.26e-14 * pow(x, s(0.750f));
    // HeIII recombination case B
    x = 2 * xHeII;
    R.alpha_rec(6, i) = 2 * 3.753e-14 * pow(x, s(1.500f)) /
                        pow(1.0 + pow(x / s(2.740f), s(0.407f)), s(2.242f));

    // legacy_physics = 0 (docs/port_notes.md 5b; harness PATCHES.md 9)
    if (!legacy_physics) {
      // P4: HeIII case B recombination. 2*3.753e-14 in the original is a
      // typo; Hui & Gnedin (1997, MNRAS 292, 27, Appendix A):
      // alpha_HeIII^B = 2*alpha_HII^B(lambda_HeII), coefficient 2.753e-14.
      x = 2 * xHeII;
      R.alpha_rec(6, i) = 2 * 2.753e-14 * pow(x, s(1.500f)) /
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
      R.alpha_rec(2, i) = R.alpha_rec(2, i) + adi;
      R.alpha_rec(5, i) = R.alpha_rec(5, i) + adi;
    }

    // HI collisional ionization
    if (xHI > 100) {
      R.gamma_col(1, i) = 0;
    } else {
      x = 2 * xHI;
      R.gamma_col(1, i) = 2.111e+01 / pow(T, s(1.5f)) * exp(-0.5 * x) /
                          pow(x, s(1.089f)) /
                          pow(1.0 + pow(x / s(0.354f), s(0.874f)), s(1.101f));
    }
    // HeI collisional ionization
    if (xHeI > 100) {
      R.gamma_col(2, i) = 0;
    } else {
      x = 2 * xHeI;
      R.gamma_col(2, i) = 3.238e+01 / pow(T, s(1.5f)) * exp(-0.5 * x) /
                          pow(x, s(1.146f)) /
                          pow(1.0 + pow(x / s(0.416f), s(0.987f)), s(1.056f));
    }
    // HeII collisional ionization
    if (xHeII > 100) {
      R.gamma_col(3, i) = 0;
    } else {
      x = 2 * xHeII;
      R.gamma_col(3, i) = 1.995e+01 / pow(T, s(1.5f)) * exp(-0.5 * x) /
                          pow(x, s(1.089f)) /
                          pow(1.0 + pow(x / s(0.553f), s(0.735f)), s(1.275f));
    }
  }

  for (int i = 1; i <= Ntemp; ++i) {
    const double T = pow(1e1, (i - 1) * lgT_del + lgT_min);
    const double xHI = HI_K / T;
    const double xHeI = HeI_K / T;
    const double xHeII = HeII_K / T;
    double x;

    // HII recombination cooling case A
    x = 2 * xHI;
    R.coef_cool(1, i) = 1.288e-13 * (k_cgs * T) * pow(x, s(1.965f)) /
                        pow(1.0 + pow(x / s(0.541f), s(0.502f)), s(2.697f));
    // HeII recombination cooling case A
    x = 2 * xHeI;
    R.coef_cool(2, i) = 2.78e-14 * (k_cgs * T) * pow(x, s(0.681f));
    // HeIII recombination cooling case A
    x = 2 * xHeII;
    R.coef_cool(3, i) = 8 * 1.288e-13 * (k_cgs * T) * pow(x, s(1.965f)) /
                        pow(1.0 + pow(x / s(0.541f), s(0.502f)), s(2.697f));
    // HII recombination cooling case B (same formula as case A in original)
    x = 2 * xHI;
    R.coef_cool(4, i) = 1.288e-13 * (k_cgs * T) * pow(x, s(1.965f)) /
                        pow(1.0 + pow(x / s(0.541f), s(0.502f)), s(2.697f));
    // HeII recombination cooling case B
    x = 2 * xHeI;
    R.coef_cool(5, i) = 2.78e-14 * (k_cgs * T) * pow(x, s(0.681f));
    // HeIII recombination cooling case B
    x = 2 * xHeII;
    R.coef_cool(6, i) = 8 * 1.288e-13 * (k_cgs * T) * pow(x, s(1.965f)) /
                        pow(1.0 + pow(x / s(0.541f), s(0.502f)), s(2.697f));

    // HI collisional ionization cooling
    if (xHI > 100) {
      R.coef_cool(7, i) = 0;
    } else {
      x = 2 * xHI;
      R.coef_cool(7, i) = 2.111e+01 * (k_cgs * T) / pow(T, s(1.5f)) *
                          exp(-0.5 * x) / pow(x, s(1.089f)) /
                          pow(1.0 + pow(x / s(0.354f), s(0.874f)), s(1.101f));
    }
    // HeI collisional ionization cooling
    if (xHeI > 100) {
      R.coef_cool(8, i) = 0;
    } else {
      x = 2 * xHeI;
      R.coef_cool(8, i) = 3.238e+01 * (k_cgs * T) / pow(T, s(1.5f)) *
                          exp(-0.5 * x) / pow(x, s(1.146f)) /
                          pow(1.0 + pow(x / s(0.416f), s(0.987f)), s(1.056f));
    }
    // HeII collisional ionization cooling
    if (xHeII > 100) {
      R.coef_cool(9, i) = 0;
    } else {
      x = 2 * xHeII;
      R.coef_cool(9, i) = 1.995e+01 * (k_cgs * T) / pow(T, s(1.5f)) *
                          exp(-0.5 * x) / pow(x, s(1.089f)) /
                          pow(1.0 + pow(x / s(0.553f), s(0.735f)), s(1.275f));
    }
    // HI collisional excitation cooling
    if (xHI > 100) {
      R.coef_cool(10, i) = 0;
    } else {
      x = 2 * xHI;
      R.coef_cool(10, i) =
          7.50e-19 * exp(-0.375 * x) / (1.0 + sqrt(T / 1e5));
    }
    // HeI collisional excitation cooling (zero in original)
    R.coef_cool(11, i) = 0;
    // HeII collisional excitation cooling
    if (xHeII > 100) {
      R.coef_cool(12, i) = 0;
    } else {
      x = 2 * xHeII;
      R.coef_cool(12, i) = 5.54e-17 * exp(-0.375 * x) / pow(T, s(0.397f)) /
                           (1.0 + sqrt(T / 1e5));
    }
    // Bremsstrahlung, Gaunt factor 1.30 (real(4) literal)
    R.coef_cool(13, i) = 1.42e-27 * sqrt(T) * s(1.30f);
    // Inverse Compton: 4*sS_cgs/c_cgs*(8./3.)*3.36E-24/(gam-1)
    R.coef_cool(14, i) =
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
      R.coef_cool(4, i) = 3.435e-30 * T * pow(x, 1.970) / pow(1 + pow(x / 2.250, 0.376), 3.720);
      x = 2 * xHeI;
      R.coef_cool(5, i) = 1.26e-14 * (k_cgs * T) * pow(x, 0.750);
      x = 2 * xHeII;
      R.coef_cool(6, i) =
          8 * 3.435e-30 * T * pow(x, 1.970) / pow(1 + pow(x / 2.250, 0.376), 3.720);

      // P5: collisional ionisation cooling = ionisation potential x rate,
      // k T_X Gamma_X (Hui & Gnedin 1997, Appendix A; Cen 1992). The
      // original uses k T, which underestimates it by T_X/T.
      R.coef_cool(7, i) = (k_cgs * HI_K) * R.gamma_col(1, i);
      R.coef_cool(8, i) = (k_cgs * HeI_K) * R.gamma_col(2, i);
      R.coef_cool(9, i) = (k_cgs * HeII_K) * R.gamma_col(3, i);

      // P2: HeI (2^3S) collisional excitation cooling, Cen (1992, ApJS 78,
      // 341): 9.10e-27 T^-0.1687 (1+sqrt(T/1e5))^-1 exp(-13179/T)
      // [erg cm^6/s], multiplies ne^2 nHeII.
      if (13179 / T > 100) {
        R.coef_cool(11, i) = 0;
      } else {
        R.coef_cool(11, i) =
            9.10e-27 / pow(T, 0.1687) * exp(-13179 / T) / (1 + sqrt(T / 1e5));
      }

      // P3: coef_cool(15) = HeII dielectronic recombination cooling,
      // 0.75 k T_HeII alpha_d (Hui & Gnedin 1997, Appendix A; = Cen 1992's
      // 1.24e-13 T^-1.5 exp(-470000/T)(1+0.3 exp(-94000/T))), times ne*nHeII.
      if (xHeII > 100) {
        R.coef_cool(15, i) = 0;
      } else {
        x = 2 * xHeII;
        const double adi =
            1.90e-3 / pow(T, 1.5) * exp(-0.375 * x) * (1 + 0.3 * exp(-0.075 * x));
        R.coef_cool(15, i) = 0.75 * (k_cgs * HeII_K) * adi;
      }
    }
  }
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
