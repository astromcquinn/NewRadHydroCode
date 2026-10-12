// Eddington-tensor diagnostic (docs/otvet_closure.md). Not in the original
// code; a diagnostic only: it never changes the state.
//
// radhydro `tensor_dump = 1`: on the hydro steps that write a zdata1
// snapshot, every RT substep adds to each cell the HI-photoionisation-
// weighted angular moments of its radiation,
//
//   S0 = sum w,  S1 = sum w n,  S2 = sum w n n^T,   w = sum_g sigma_HI(nu_g) n_g,
//
// where for the rays every deposit of a ray segment contributes with the
// ray's direction n (an optically thick source cell's own deposit counts as
// isotropic), and for the moment methods each group's N, F and D N after the
// transport (D the closure actually used, OTVET h or M1). h_true = S2/S0 is
// the Eddington tensor and S1/S0 the reduced flux F/(c N). The sums run over
// the substeps of the hydro step and are written as their means
// (tensor.bin), with the optically thin OTVET tensor of the step's sources
// (otvet_h.bin) and the source list (sources.txt).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace radray {

struct Simulation;

struct TensorDiag {
  bool active = false;
  int nsub = 0;
  double dt_sum = 0;  // sum of dtsec_rt over the substeps
  // 10 per cell: S0, S1 (x, y, z), S2 (xx, yy, zz, xy, xz, yz)
  std::vector<double> S;
};

// Starts the accumulation for one hydro step (allocates and zeroes).
void tensor_diag_begin(Simulation& s);
// One deposit v (rrk::cast_1 layout: dng(5), dng*nu(5), dx*sum(dng)(3)) of
// a ray segment into cell c.
void tensor_diag_ray(TensorDiag& T, std::int64_t c, const double* v);
// An isotropic deposit (thick source cell), same layout (v[10..12] unused).
void tensor_diag_iso(TensorDiag& T, std::int64_t c, const double* v);
// Writes tensor.bin (float32, 10 per cell, means over the substeps),
// otvet_h.bin (6 per cell) and sources.txt (i j k Luv, 1-based) into dir,
// then stops the accumulation.
void tensor_diag_write(Simulation& s, const std::string& dir, bool io);

}  // namespace radray
