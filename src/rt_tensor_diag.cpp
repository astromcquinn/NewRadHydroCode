// The Eddington-tensor diagnostic of the rays (include/radray/rt_tensor_diag.hpp,
// docs/otvet_closure.md section 1): radhydro tensor_dump = 1, serial RT
// backend only, off by default (a null Simulation::tdiag, so byte-identical
// results). On the hydro steps that write a zdata1 snapshot, radhydro calls
// tensor_diag_begin; rt_step_1 calls tensor_diag_substep at the start of
// every RT substep; the deposits of ray_tracing.cpp (rrk::cast_1's callback
// and the thick source cell's own deposit) call tensor_diag_ray /
// tensor_diag_iso; the moment RT adds its N, F and D N into the same
// accumulator (moment_tensor_diag, src/rt_moments.cpp); and
// radhydro calls tensor_diag_write with the snapshot directory. It
// accumulates per cell the Gamma_HI-weighted angular moments S0, S1, S2 of
// the radiation (the header's formula) and writes their means over the
// substeps as tensor.bin, with tensor_odd.bin, the optically thin OTVET
// tensor of the step's sources (otvet_h.bin) and sources.txt, which
// tools/analysis/tensor_compare.py and apps/otvet_tensor.cpp read. Nothing
// here is read by the RT.
#include "radray/rt_tensor_diag.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>

#include "radray/constants.hpp"
#include "radray/kernels/rt_cell.hpp"
#include "radray/parallel.hpp"
#include "radray/rt_moments.hpp"
#include "radray/simulation.hpp"

namespace radray {

// One hydro step's accumulation: 10 doubles per RT cell in S (even substeps)
// and So (odd substeps), zeroed; the counters restart.
void tensor_diag_begin(Simulation& s) {
  if (!s.tdiag) s.tdiag = std::make_shared<TensorDiag>();
  TensorDiag& T = *s.tdiag;
  T.S.assign(10 * static_cast<std::size_t>(s.rt.ncell), 0.0);
  T.So.assign(10 * static_cast<std::size_t>(s.rt.ncell), 0.0);
  T.nsub = 0;
  T.nsub_odd = 0;
  T.dt_sum = 0;
  T.active = true;
}

// Start of a substep (rt_step_1): count it, add its length and point the
// deposits at this substep's accumulator.
void tensor_diag_substep(Simulation& s) {
  TensorDiag& T = *s.tdiag;
  T.nsub += 1;
  T.dt_sum += s.clk.dtsec_rt;
  // odd substeps accumulate into So, even ones into S; S += So at the end
  if (T.nsub % 2 == 1) {
    T.nsub_odd += 1;
    T.acc = T.So.data();
  } else {
    T.acc = T.S.data();
  }
}

namespace {
// Gamma_HI weight of a deposit: sum_g sigma_HI(nu_g) dng_g (cm^2 nH units).
// v is the cast_1 deposit: v[g] the photons dng of group g, v[5 + g] their
// dng * nu, so the group's mean frequency is v[5 + g] / v[g]; sig[0] is the
// HI cross section of rrk::photoion_sigma (the one the RT uses).
inline double gamma_weight(const double* v) {
  double nu[kNfreq], sig[3][kNfreq];
  for (int q = 0; q < kNfreq; ++q) nu[q] = v[q] > 0 ? v[5 + q] / v[q] : 0.0;
  rrk::photoion_sigma(nu, sig);
  double w = 0;
  for (int q = 0; q < kNfreq; ++q)
    if (v[q] > 0) w += sig[0][q] * v[q];
  return w * 1e18;  // O(1) numbers
}
}  // namespace

// A ray segment's deposit into cell c with the ray's direction: v[10..12]
// is dx times the summed dng, so dividing by that sum recovers the unit
// direction n. S0 += w, S1 += w n, S2 += w n n^T (xx, yy, zz, xy, xz, yz).
void tensor_diag_ray(TensorDiag& T, std::int64_t c, const double* v) {
  const double w = gamma_weight(v);
  if (!(w > 0)) return;
  double sd = 0;
  for (int q = 0; q < kNfreq; ++q) sd += v[q];
  if (!(sd > 0)) return;
  const double n[3] = {v[10] / sd, v[11] / sd, v[12] / sd};
  double* S = &T.acc[10 * c];
  S[0] += w;
  S[1] += w * n[0];
  S[2] += w * n[1];
  S[3] += w * n[2];
  S[4] += w * n[0] * n[0];
  S[5] += w * n[1] * n[1];
  S[6] += w * n[2] * n[2];
  S[7] += w * n[0] * n[1];
  S[8] += w * n[0] * n[2];
  S[9] += w * n[1] * n[2];
}

// An isotropic deposit (the optically thick source cell absorbing its own
// emission): no net flux, S2 += w I / 3.
void tensor_diag_iso(TensorDiag& T, std::int64_t c, const double* v) {
  const double w = gamma_weight(v);
  if (!(w > 0)) return;
  double* S = &T.acc[10 * c];
  S[0] += w;
  S[4] += w / 3;
  S[5] += w / 3;
  S[6] += w / 3;
}

// End of the hydro step: every rank stops accumulating, the I/O rank writes.
// tensor.bin holds the means over all substeps (S + So), tensor_odd.bin
// over the odd ones only (their difference from the even half estimates the
// rays' sampling noise); both float32, 10 per cell in the S layout.
void tensor_diag_write(Simulation& s, const std::string& dir, bool io) {
  if (!s.tdiag || !s.tdiag->active) return;
  TensorDiag& T = *s.tdiag;
  T.active = false;
  T.acc = nullptr;
  if (!io) return;
  const std::int64_t nc = s.rt.ncell;
  const double inv = T.nsub > 0 ? 1.0 / T.nsub : 0.0;
  const double invo = T.nsub_odd > 0 ? 1.0 / T.nsub_odd : 0.0;
  std::vector<float> out(10 * static_cast<std::size_t>(nc)), outo(out.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    out[i] = static_cast<float>((T.S[i] + T.So[i]) * inv);
    outo[i] = static_cast<float>(T.So[i] * invo);
  }
  auto wr = [&](const std::string& name, const std::vector<float>& v) {
    std::FILE* fp = std::fopen((dir + "/" + name).c_str(), "wb");
    if (!fp) throw std::runtime_error("tensor_diag_write: cannot open " + dir + "/" + name);
    std::fwrite(v.data(), sizeof(float), v.size(), fp);
    std::fclose(fp);
  };
  wr("tensor.bin", out);
  wr("tensor_odd.bin", outo);
  // the optically thin OTVET tensor of this step's sources, as the solver
  // computes it (rt_moments.hpp; 6 floats per cell): the closure to compare
  std::vector<float> h;
  std::int64_t nsrc = 0;
  compute_otvet_tensor(s.rt, s.units, 0, h, &nsrc);
  wr("otvet_h.bin", h);
  std::ofstream src(dir + "/sources.txt");
  src << "# i j k Luv (1-based RT cells; sfr > sfr_min); tensor.bin substeps " << T.nsub
      << " odd " << T.nsub_odd << "\n";
  {
    // the header lets the analysis put the sources' direct light in the
    // dump's units: qfac[g] = photons per substep per unit Luv in group g
    // (the emission of ray_source: fResc * pop2_init[g] * dt / vol_rt / nH_avg
    // with dt the mean substep) and sigma18[g] = 1e18 sigma_HI(freq_init[g]),
    // the weight of a fresh ray's photons in group g
    const double dt = T.nsub > 0 ? T.dt_sum / T.nsub : 0.0;
    double nu[kNfreq], sig[3][kNfreq];
    for (int g = 0; g < kNfreq; ++g) nu[g] = cp::freq_init[g];
    rrk::photoion_sigma(nu, sig);
    src.precision(10);
    src << "# qfac";
    for (int g = 0; g < kNfreq; ++g)
      src << " " << cp::fResc * cp::pop2_init[g] * dt / s.units.vol_rt / s.units.nH_avg;
    src << "\n# sigma18";
    for (int g = 0; g < kNfreq; ++g) src << " " << sig[0][g] * 1e18;
    src << "\n# redshift " << s.clk.redshift << " N " << s.rt.N << "\n";
  }
  // the source cells by the RT's own criterion (ray_source: sfr > sfr_min
  // and Luv > 0), 1-based like the Fortran, k slowest
  const RTGrid& G = s.rt;
  for (int k = 1; k <= G.N; ++k)
    for (int j = 1; j <= G.N; ++j)
      for (int i = 1; i <= G.N; ++i) {
        const std::int64_t c = G.cell(i, j, k);
        if (static_cast<double>(G.sfr[c]) > s.units.sfr_min && G.Luv[c] > 0)
          src << i << " " << j << " " << k << " " << G.Luv[c] << "\n";
      }
  // free the 160 B per cell until the next snapshot step
  T.S.clear();
  T.S.shrink_to_fit();
  T.So.clear();
  T.So.shrink_to_fit();
  T.acc = nullptr;
}

}  // namespace radray
