// Grid / ray-storage sizes. Mirrors the parameters of radhydro_global.f90,
// but with a runtime N (Ngas = Nrt = N per docs/validation_spec.md).
#pragma once

#include <cstdint>

namespace radray {

// Hydro solver (docs/hydro_unsplit.md): Legacy = the original's split
// moving-frame TVD scheme (src/hydro.cpp, bit-identical to the Fortran);
// Unsplit = the unsplit MUSCL-Hancock + HLLC scheme with the dual-energy
// formulation (src/hydro_unsplit.cpp, kernels/hydro_mh.hpp).
enum class HydroScheme { Legacy, Unsplit };

struct Config {
  // radhydro_global.f90: Ncpu. The reference is run with ONE thread, and the
  // port is serial; Ncpu enters ray-buffer sizes and the domain
  // decomposition, so we keep it explicit.
  int Ncpu = 1;

  int Ngas = 0;  // gas grid size per dimension
  int Nrt = 0;   // RT grid size per dimension (== Ngas here)
  int Hrt = 0;   // Nrt/2
  int Ndm = 0;   // DM particles per dimension (Npart = Ndm^3); default Ngas/2
  std::int64_t Npart = 0;  // Ndm^3

  // cosmo_parameters.f90: box (comoving Mpc/h). A compile-time parameter in
  // the original (5); a run-time value here (the Fortran harnesses set it
  // with -DBOX).
  double box = 5;
  std::int64_t Ncell = 0;  // Ngas^3
  std::int64_t Nrtc = 0;   // Nrt^3
  double ngasr = 0, nrtr = 0, hrtr = 0, nrtcr = 0, ncellr = 0;

  // N_rt_gas = Ngas/Nrt (gas cells per RT cell per dimension)
  int N_rt_gas = 1;
  int N_rt_gas_c = 1;
  double x_rt_gas = 1;

  double x_gas_box = 0;  // box/ngasr
  double x_dm_box = 0;   // box/ndmr
  double x_rt_box = 0;   // box/nrtr

  // Rays: Nray_max = nray_max_factor*Nrt^3 (64 in the original),
  // Nfree_max = Nray_max/(4*Ncpu)
  int nray_max_factor = 64;
  std::int64_t Nray_max = 0;
  std::int64_t Nfree_max = 0;

  // Runtime switch for the RT bookkeeping fixes B1-B4, B8, B9
  // (docs/port_notes.md section 5). true (params.txt legacy_rt = 1, the
  // default) reproduces the original Fortran bit-for-bit.
  bool legacy_rt = true;

  // Runtime switch for the atomic-physics fixes P1-P6 (docs/port_notes.md
  // section 5b). true (params.txt legacy_physics = 1, the default) keeps the
  // original rates and cooling sum bit-for-bit. The rate tables depend on it:
  // call calc_atomic_rates(rates, legacy_physics) after changing it.
  bool legacy_physics = true;

  // Reduced speed of light: c_hat = c_reduced * c. Enters only through
  // calc_timestep (the original's cscale): nstep_rt = nint(dl_c/dl_rt*cscale).
  // 1 (default) is the original.
  double c_reduced = 1.0;

  // Deterministic parallel ray tracing (docs/parallel.md). false (default):
  // the original sweep order (8 domains, one thread) and the Fortran ray
  // pool, bit-identical to the Fortran harness. true: ray_tracing_1/3 run
  // over 8^3-cell domains in colour-major order, in parallel within a
  // colour, with per-thread free lists; bit-identical for any thread count,
  // but not to the legacy order. Requires legacy_rt = false.
  bool det_parallel = false;

  // Maximum expansion per hydro step, da/a <= ra_max/(1-ra_max). The
  // original's real(4) 0.1 (da/a = 1/9) unless set (apps/radhydro allows a
  // different value only with legacy_rt = 0).
  double ra_max = static_cast<double>(0.1f);

  // Carry the rounding remainder of nstep_rt = nint(light cells) from one
  // hydro step to the next, so that the total number of RT substeps (light
  // travel at c_reduced) is exact over the run even when a hydro step spans
  // only a few cells of light travel. false: the original rounding.
  bool rt_substep_carry = false;

  // Rotation of the HEALPix frame of the rays a source emits, against the
  // grid-aligned photon deficit B16 (docs/port_notes.md). 0 (default): the
  // original frame, fixed to the grid. 1: one fixed rotation for the whole
  // run. 2: a new rotation every RT substep (the same for all sources).
  // 3: a new rotation every RT substep and source cell. Only emission
  // changes: a rotated ray keeps its direction in dx and the grid-frame
  // pixel containing it in ipix, so splits and merges are unchanged
  // (ray_frame_rotation in ray_tracing.hpp). Deterministic: the rotation is
  // a hash of ray_rotate_seed, the substep and the cell, independent of the
  // thread count.
  int ray_rotate = 0;
  std::uint64_t ray_rotate_seed = 0;

  // params.txt hydro_scheme = legacy (default) | unsplit, and the parameters
  // of the unsplit scheme (set_hydro_scheme_param in hydro_unsplit.hpp):
  //   unsplit_cfl: Courant number against sum_d(|v_d| + c) of the sub-cycled
  //                unsplit substeps;
  //   dual_energy: 1 = dual-energy formulation (Bryan et al. 1995), 0 = total
  //                energy only;
  //   dual_eta1, dual_eta2: its selection thresholds, the values of Bryan
  //                et al. 1995 as restated in Bryan et al. 2014 (Enzo, section
  //                4.1.1, eqs. 44-45); docs/hydro_unsplit.md section 2.
  HydroScheme hydro_scheme = HydroScheme::Legacy;
  double unsplit_cfl = 0.8;
  bool dual_energy = true;
  double dual_eta1 = 1e-3;
  double dual_eta2 = 0.1;
  //   dual_eta1_local: 1 = the eta1 test against the cell's own E (the
  //                cell-local form); 0 (default) = against the maximum of E
  //                over the cell and its 6 face neighbours, which is not
  //                fooled by stagnation cells between fast streams.
  bool dual_eta1_local = false;
  //   dual_nb_radius: the neighbourhood of that maximum, +-1 cell (6
  //                neighbours) or +-2 cells (12) along the axes. 2 (default)
  //                keeps the eta2 synchronisation from thermalising the
  //                unresolved velocity dispersion of a stagnation cell
  //                (docs/hydro_unsplit.md section 2).
  int dual_nb_radius = 2;
  //   dual_reset_rule: when the pressure comes from e, E := KE + e always
  //                (0), only if the change is at most eta1 E (1), as 1 or
  //                in thermal-dominated cells (2), only if E - KE <= eta1
  //                max_{+-1} E (3). docs/hydro_unsplit.md section 2.
  int dual_reset_rule = 0;
  //   dual_eta2_local: 1 = the eta2 test against the cell's own E (Nyx's
  //                resync); dual_p_old: 1 = the old-time cell pressure in
  //                the -p div v work term (Bryan et al. 2014 eq. 43) instead
  //                of the Hancock half-step pressure. Both for the
  //                comparison of docs/hydro_unsplit.md section 2.
  bool dual_eta2_local = false;
  bool dual_p_old = false;

  // Ngas = Nrt = N; Ndm <= 0 means the original's Ngas/2; box in Mpc/h.
  static Config make(int N, int nray_max_factor = 64, int Ndm = -1, double box = 5);
};

}  // namespace radray
