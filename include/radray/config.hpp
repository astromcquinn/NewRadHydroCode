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
enum class HydroBackend { Serial, AMReX };

// Radiative transfer method (docs/rt_moments.md): Rays = the original's
// adaptive ray tracing (default, bit-identical to the Fortran); M1 and OTVET
// = the two-moment photon solver of src/rt_moments.cpp with the M1
// (Levermore 1984) closure or the optically thin variable Eddington tensor
// (Gnedin & Abel 2001). New methods, not in the original.
enum class RTMethod { Rays, M1, OTVET };

// Where the ray tracing runs (docs/stage_d.md): Serial = the original's
// per-cell ray chains on the host (default, bit-identical to the Fortran);
// AMReX = the sorted-ray-array RT on an AMReX level (MPI + OpenMP),
// legacy_rt = 0 only. The ray storage of the AMReX backend: Float = the
// original's float fields (80 B of ray data, bitwise validation); RStruct =
// the compact R-struct of docs/ray_compression.md (dithered, lossy).
enum class RTBackend { Serial, AMReX };
enum class RayEncoding { Float, RStruct };

// Where the PM gravity runs (docs/stage_e_pm.md): Serial = the original's
// linked lists, float CIC and FFTW on the host (default, bit-identical to
// the Fortran); AMReX = the dark matter in an AMReX ParticleContainer, int64
// fixed-point CIC deposits and the amrex::FFT Poisson solve on an AMReX
// level of the PM mesh (MPI + OpenMP, rank- and thread-independent).
enum class PMBackend { Serial, AMReX };

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
  // PM mesh refinement (params.txt pm_refine, docs/pm_refine.md): the PM
  // mesh has pm_refine*Ngas cells per dimension; Ndm must divide it.
  // 1 (default) is the original, where the mesh is the gas grid.
  int pm_refine = 1;

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

  // Cosmological redshift of the radiation per distance travelled instead of
  // per time (docs/production_run.md, reduced speed of light). The rays
  // (ray_redshifting) and the moment RT (moment_redshift) are redshifted once
  // per hydro step by fscale = (a - da)/a. false (default): the original.
  // true with c_reduced < 1: by fscale^s, s = min(nstep_rt*dl_rt/dl_c, 1) the
  // distance the radiation actually moved in the step over the light distance
  // (redshift_fscale in cosmology.hpp); with c_reduced >= 1 it changes
  // nothing. Not tied to legacy_physics.
  bool redshift_per_distance = false;

  // RT cells a ray crosses per RT substep (docs/rt_multicell.md). 1
  // (default) is the original: nstep_rt counts the light crossings of one RT
  // cell. k > 1: nstep_rt = nint(dl_c/(k dl_rt)*cscale) substeps, each k
  // times longer, in which every ray moves k cells (up to 3k+1 segments),
  // so the light distance per hydro step is unchanged; the photo-ionisation
  // rate of a cell uses the crossing speed k*len_rt/dtsec_rt, the split
  // radius clamp is [rmin, rmin+k], and for k > 4 sources emit at the nside
  // whose split radius covers the first path (rrk::emit_nside: 768 rays for
  // k = 5..8). legacy_rt = 0 and rt_method = rays only;
  // with rt_backend = amrex the ghost layers grow to k+1. The test of
  // coarse-substep RT for AMR (Stage F): a level-l region would move rays
  // 2^l fine cells per coarse substep.
  int rt_cells_per_substep = 1;

  // B17 fix (docs/port_notes.md, docs/rt_multicell.md): photo_cell's photon
  // limit is iterated to convergence (at most rrk::kPhotonIterMax
  // iterations) instead of ending with ne after at most 5 iterations, which
  // lets photon-starved front cells ionise without the photons. false is
  // the original; part of legacy_physics = 0 since 2026-10-14 (mirrored in
  // the Fortran harnesses, PATCHES.md 14); rays and the moment RT.
  bool photon_conserve = false;
  // The params.txt key photon_conserve: -1 (default) follows legacy_physics
  // (B17 is part of legacy_physics = 0 since 2026-10-14), 0 | 1 override.
  static bool photon_conserve_effective(int key, bool legacy_physics) {
    return key < 0 ? !legacy_physics : key == 1;
  }

  // params.txt cool_suppress_overdensity (Stage F, docs/stage_f.md section
  // 2; docs/port_notes.md 5b "Other suspected physics issues"): above this
  // gas overdensity x = rho/rho_avg the original multiplies the cooling rate
  // by (200/x)^2 and the SN heating by (200/x)^1.5, an overcooling guard.
  // 200 (default) is the original; 0 turns the suppression off (dense gas
  // cools; the star-formation sink is then the guard). Not tied to
  // legacy_physics; no Fortran mirror (a new flag).
  double cool_suppress_overdensity = 200.0;

  // How a ray deposits its photons along the segments of its path
  // (rrk::cast_1/cast_3, docs/rt_multicell.md): 0 = share (the original),
  // 1 = exposure, which keeps a ray's photons for the cells further along
  // the path unless absorbed, so that an I-front reached in mid-path gets all
  // of them (with k cells per substep the share deposit smears fronts over
  // ~k cells). legacy_rt = 0 only; any rt_cells_per_substep.
  int ray_deposit = 0;

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
  //                0 = auto, 1 = direct sum, 2 = FFT convolution;
  //   moment_inject: when the source and recombination photons of an RT
  //                substep enter the moments: 0 = before the transport
  //                (default; a source cell's photons are mostly moved to its
  //                neighbours before photo_cell absorbs anything), 1 = after
  //                it, so they are absorbed in their own cell first, as the
  //                rays' photons are (docs/rt_moments.md).
  RTMethod rt_method = RTMethod::Rays;
  int moment_nsub = 4;
  int otvet_tensor = 0;
  int moment_inject = 0;
  // Alternative OTVET tensors (otvet_closure.hpp, docs/otvet_closure.md;
  // set_otvet_closure_param). New options, not in the original; the
  // defaults keep the thin minimum-image tensor:
  //   otvet_closure: 0 = thin (default), 1 = screen (e^(-r/lambda) over the
  //                periodic images), 2 = cone (HI-attenuated cones);
  //   otvet_split: 1 = the direct light carries h, the rest of the photons
  //                is isotropic (thin or screen);
  //   otvet_screen_lambda: screening length in RT cells (0 = from the
  //                photon budget);
  //   otvet_cone_nside, _coarse, _rc, _taucap, _smax, _budget, _blend: the
  //                cones (HEALPix nside, receiver grid factor, near-field
  //                radius in cells, per-cell tau cap, reach in boxes, extra
  //                cutoff from the photon budget, Gaussian sigma in cells of
  //                the blend with the split thin tensor; 0 = none).
  int otvet_closure = 0;
  int otvet_split = 0;
  double otvet_screen_lambda = 0;
  int otvet_cone_nside = 4;
  int otvet_cone_coarse = 4;
  double otvet_cone_rc = 8;
  double otvet_cone_taucap = 1;
  double otvet_cone_smax = 8;
  int otvet_cone_budget = 1;
  double otvet_cone_blend = 4;

  // params.txt rt_backend = serial (default) | amrex and the parameters of
  // the AMReX backend (set_rt_backend_param in rt_backend.hpp;
  // docs/stage_d.md):
  //   ray_encoding: float (default) | rstruct;
  //   amrex_max_grid_size: box size of the level's decomposition;
  //   amrex_ghosts: ghost layers of the optical depth / deposit FabArrays
  //                (2 is the reach of one substep, checked at run time;
  //                at least k+1 are used with rt_cells_per_substep = k);
  //   amrex_redist_local: cells a ray may move per substep for the local
  //                (neighbour-only) Redistribute; 0 = global (at least k+1
  //                are used with rt_cells_per_substep = k).
  RTBackend rt_backend = RTBackend::Serial;
  RayEncoding ray_encoding = RayEncoding::Float;
  int amrex_max_grid_size = 16;
  int amrex_ghosts = 2;
  int amrex_redist_local = 3;

  // params.txt pm_backend = serial (default) | amrex (set_pm_backend_param
  // in pm_backend.hpp; docs/stage_e_pm.md):
  //   amrex_pm_max_grid_size: box size of the PM mesh's decomposition
  //                (0 = amrex_max_grid_size; a multiple of pm_refine);
  //   amrex_pm_fft: fft of the Poisson solve: pencil (default; 1-D FFTW
  //                transforms per pencil, bitwise rank- and thread-
  //                independent) | amrex (the amrex::FFT::R2C solver).
  //   pm_gas_levels: with a refined gas hierarchy (amrex_max_level > 0,
  //                amrex_state = 1) the gas of every level enters the PM:
  //                1 (default) the mesh takes the density of the finest
  //                level covering each mesh cell and every level's cells
  //                are kicked by the mesh force over their own volume
  //                (docs/stage_f.md section 4); 0 the level-0 coupling of
  //                Stage F1 (the fine cells get the level-0 kick through
  //                the hydro's push-down). No effect without refinement.
  PMBackend pm_backend = PMBackend::Serial;
  int amrex_pm_max_grid_size = 0;
  int amrex_pm_fft = 0;
  int pm_gas_levels = 1;
  // params.txt hydro_backend = serial (default) | amrex (the unsplit step on
  // AMReX MultiFabs, MPI + OpenMP, bitwise identical to the serial unsplit
  // step; hydro_scheme = unsplit only; set_hydro_backend_param in
  // hydro_backend.hpp, docs/stage_e_hydro.md) and
  //   amrex_hydro_max_grid_size: box size of the hydro level's decomposition.
  HydroBackend hydro_backend = HydroBackend::Serial;
  int amrex_hydro_max_grid_size = 32;
  // Stage F (docs/stage_f.md): the gas hierarchy of hydro_backend = amrex.
  //   amrex_max_level: refinement levels (0 = one level, the Stage E step);
  //   amrex_refine_mass: a level-l cell is tagged when its gas mass exceeds
  //                amrex_refine_mass times the mean base-cell mass, i.e.
  //                rho/rho_avg > amrex_refine_mass * 8^l;
  //   amrex_regrid_every: regrid every n hydro calls;
  //   amrex_fine_max_grid_size, amrex_blocking_factor, amrex_n_error_buf,
  //   amrex_grid_eff: the AMReX gridding parameters of the refined levels.
  int amrex_max_level = 0;
  double amrex_refine_mass = 8.0;
  int amrex_regrid_every = 1;
  int amrex_fine_max_grid_size = 32;
  int amrex_blocking_factor = 8;
  int amrex_n_error_buf = 1;
  double amrex_grid_eff = 0.7;
  // amrex_refine_box = i0,j0,k0,i1,j1,k1 (level-0 cells, inclusive): a
  // fixed region refined to amrex_max_level in addition to the density
  // tagging (tests: a refined region around a source). Empty = none.
  int amrex_refine_box[6] = {0, 0, 0, -1, -1, -1};
  // The star-formation sink (docs/stage_f.md section 2; hydro_backend =
  // amrex): a cell above the gas overdensity sf_overdensity loses
  // min(1, sf_efficiency dt/t_ff) of its mass per hydro call, at constant
  // specific quantities, to a star particle (pm_backend = amrex with
  // amrex_state = 1; a mass loss in the hydro-only drivers). 0 = off (the
  // original).
  double sf_overdensity = 0;
  double sf_efficiency = 1;
  // params.txt amrex_state = 0 (default) | 1: the gas state lives on the
  // AMReX levels between steps on one shared decomposition
  // (amrex_max_grid_size), the host arrays only on rank 0 as a mirror for
  // the outputs (amrex_state.hpp, docs/stage_e_state.md). Needs the three
  // backends, hydro_scheme = unsplit, rt_method = rays.
  bool amrex_state = false;
  // params.txt amrex_balance = 0 (default) | 1: with amrex_state = 1, the
  // shared decomposition is rebalanced before every RT step on the cost
  // cells + amrex_balance_ray_weight * rays per box (amrex_balance_method
  // = sfc | knapsack), when the balance efficiency would improve by more
  // than amrex_balance_min_gain (negative: at every RT step;
  // docs/stage_e_state.md section 5).
  bool amrex_balance = false;
  double amrex_balance_ray_weight = 10.0;
  double amrex_balance_min_gain = 0.05;
  int amrex_balance_method = 0;  // 0 sfc, 1 knapsack

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
  //   recomb_rad_cell_onspot: tau > 0 = a cell whose own optical depth at
  //                the HI edge (group 1, where the HI recombination photons
  //                are) exceeds tau keeps its recombination photons on the
  //                spot (case B) even when its block launches; 0 (default) =
  //                the block decides for all its cells
  //                (docs/rt_methods_512.md section 6c).
  bool recomb_case_A = false;
  bool recomb_rad = false;
  bool recomb_rad_lines = false;
  int recomb_rad_nside = 1;
  int recomb_rad_block = 2;
  double recomb_rad_pool = 0.75;
  double recomb_rad_cell_onspot = 0;

  // Ngas = Nrt = N; Ndm <= 0 means the original's Ngas/2; box in Mpc/h.
  static Config make(int N, int nray_max_factor = 64, int Ndm = -1, double box = 5);
};

}  // namespace radray
