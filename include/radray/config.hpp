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

// Radiative transfer method (docs/rt_moments.md): Rays = the original's
// adaptive ray tracing (default, bit-identical to the Fortran); M1 and OTVET
// = the two-moment photon solver of src/rt_moments.cpp with the M1
// (Levermore 1984) closure or the optically thin variable Eddington tensor
// (Gnedin & Abel 2001). New methods, not in the original.
enum class RTMethod { Rays, M1, OTVET };

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

  // params.txt rt_method = rays (default) | m1 | otvet and the parameters
  // of the moment solver (set_rt_method_param in rt_moments.hpp):
  //   moment_nsub: transport subcycles per RT substep (Courant number
  //                1/moment_nsub of the explicit GLF scheme; >= 3 in 3-D,
  //                and 3 decouples odd and even cells, docs/rt_moments.md);
  //   otvet_tensor: how the OTVET tensor is summed over the sources:
  //                0 = auto, 1 = direct sum, 2 = FFT convolution.
  RTMethod rt_method = RTMethod::Rays;
  int moment_nsub = 4;
  int otvet_tensor = 0;

  // Recombinations and recombination radiation (docs/recomb_radiation.md;
  // set_recomb_param in recomb_rad.hpp). New options, not in the original:
  //   recomb_case: B (default, the original) = case B recombination rates
  //                and cooling in photo_cell (on-the-spot); A = case A.
  //   recomb_rad:  1 = the ground-state recombinations (alpha_A - alpha_B)
  //                of every ionised cell are emitted as rays (requires
  //                recomb_case = A; with the rays legacy_rt = 0; with
  //                rt_method = m1 | otvet every cell emits into the moments).
  //   recomb_rad_lines: 1 = also the helium recombination photons that
  //                ionise H (HeI 19.8/21.2 eV lines and 2-photon, HeII
  //                Lya, 2-photon and Balmer continuum; needs recomb_rad).
  //   recomb_rad_nside: HEALPix nside of the emitted rays (12 nside^2 rays
  //                per emission; default 1, i.e. 12 rays).
  //   recomb_rad_block: emission is gathered over blocks of b^3 RT cells and
  //                launched from one cell of the block, rotating every
  //                substep (default 2).
  //   recomb_rad_pool: recombination rays are only launched while the live
  //                rays stay below this fraction of Nray_max; the rest of
  //                the emission is absorbed on the spot (default 0.75).
  bool recomb_case_A = false;
  bool recomb_rad = false;
  bool recomb_rad_lines = false;
  int recomb_rad_nside = 1;
  int recomb_rad_block = 2;
  double recomb_rad_pool = 0.75;

  // Ngas = Nrt = N; Ndm <= 0 means the original's Ngas/2; box in Mpc/h.
  static Config make(int N, int nray_max_factor = 64, int Ndm = -1, double box = 5);
};

}  // namespace radray
