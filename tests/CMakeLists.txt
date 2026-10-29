add_executable(test_fortran_semantics test_fortran_semantics.cpp)
target_link_libraries(test_fortran_semantics PRIVATE radray)
add_test(NAME fortran_semantics
         COMMAND test_fortran_semantics ${CMAKE_CURRENT_SOURCE_DIR}/fortran_consts/expected_consts.txt)

add_executable(test_healpix test_healpix.cpp)
target_link_libraries(test_healpix PRIVATE radray)
add_test(NAME healpix COMMAND test_healpix)

add_executable(healpix_tool healpix_tool.cpp)
target_link_libraries(healpix_tool PRIVATE radray)

add_executable(test_rays test_rays.cpp)
target_link_libraries(test_rays PRIVATE radray)
add_test(NAME rays COMMAND test_rays)

# The ray encodings of the AMReX backend (no AMReX needed; docs/stage_d.md 3)
add_executable(test_ray_encoding test_ray_encoding.cpp)
target_link_libraries(test_ray_encoding PRIVATE radray)
add_test(NAME ray_encoding COMMAND test_ray_encoding)

# Moment RT (rt_method = m1 | otvet, docs/rt_moments.md)
add_executable(test_moments test_moments.cpp)
target_link_libraries(test_moments PRIVATE radray)
add_test(NAME moments COMMAND test_moments)
set_tests_properties(moments PROPERTIES LABELS moments)

# Redshift per distance travelled with a reduced speed of light
# (redshift_per_distance, docs/production_run.md)
add_executable(test_redshift_distance test_redshift_distance.cpp)
target_link_libraries(test_redshift_distance PRIVATE radray)
add_test(NAME redshift_distance COMMAND test_redshift_distance)

# The runtime cosmology and the preheating floor (Config::cosmo, T_min,
# z_floor; docs/production_run.md section 10)
add_executable(test_cosmology_runtime test_cosmology_runtime.cpp)
target_link_libraries(test_cosmology_runtime PRIVATE radray_halos)
if(TARGET radray_pm)
  target_link_libraries(test_cosmology_runtime PRIVATE radray_pm)
  target_compile_definitions(test_cosmology_runtime PRIVATE RADRAY_TEST_PM=1)
endif()
add_test(NAME cosmology_runtime COMMAND test_cosmology_runtime)
set_tests_properties(cosmology_runtime PROPERTIES LABELS "physics;cosmology")

# Halo finder and sources (docs/halos_and_sources.md)
add_executable(test_halos test_halos.cpp)
target_link_libraries(test_halos PRIVATE radray_halos)
add_test(NAME halos COMMAND test_halos)

# Abundance matching to UV LFs (docs/abundance_matching.md)
add_executable(test_abundance_matching test_abundance_matching.cpp)
target_link_libraries(test_abundance_matching PRIVATE radray_halos)
add_test(NAME abundance_matching COMMAND test_abundance_matching)
# X-ray background heating (docs/xray_heating.md)
add_executable(test_xray test_xray.cpp)
target_link_libraries(test_xray PRIVATE radray_halos)
add_test(NAME xray_background COMMAND test_xray)
set_tests_properties(xray_background PROPERTIES LABELS physics)
# thermal_init = recfast (docs/thermal_history.md)
add_executable(test_thermal test_thermal.cpp)
target_link_libraries(test_thermal PRIVATE radray)
add_test(NAME thermal_init COMMAND test_thermal)
set_tests_properties(thermal_init PROPERTIES LABELS physics)

add_executable(test_halo_catalog test_halo_catalog.cpp)
target_link_libraries(test_halo_catalog PRIVATE radray_halos)
find_program(RADRAY_GFORTRAN gfortran)
if(RADRAY_GFORTRAN AND RADRAY_HAVE_FORTRAN_REF)
  # Reference program with the original source_type layout and read_halo
  # (reference/fortran/tests, the private Fortran repository). The
  # target is named differently from its output file: Ninja rejects a target
  # that names its own output.
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/halo_catalog_ref
    COMMAND ${RADRAY_GFORTRAN} -O0 -o ${CMAKE_CURRENT_BINARY_DIR}/halo_catalog_ref
            ${RADRAY_FORTRAN_DIR}/tests/halo_catalog_ref.f90
    DEPENDS ${RADRAY_FORTRAN_DIR}/tests/halo_catalog_ref.f90)
  add_custom_target(halo_catalog_ref_build ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/halo_catalog_ref)
  set(HALO_REF ${CMAKE_CURRENT_BINARY_DIR}/halo_catalog_ref)
else()
  set(HALO_REF "")  # the test then skips the Fortran comparisons
endif()
add_test(NAME halo_catalog
         COMMAND test_halo_catalog "${HALO_REF}" ${CMAKE_CURRENT_BINARY_DIR}/halo_catalog_run)

# legacy_physics = 0 rates vs independent evaluations of the published fits
add_executable(test_rates test_rates.cpp)
target_link_libraries(test_rates PRIVATE radray)
add_test(NAME rates COMMAND test_rates)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME healpix_vs_healpy
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/healpix_vs_healpy.py
                   $<TARGET_FILE:healpix_tool> 1000000)
  set_tests_properties(healpix_vs_healpy PROPERTIES SKIP_RETURN_CODE 77)

  # The runtime cosmology's IC tool, a Planck RT run and the preheating
  # floor on a cooling test (tests/cosmology_check.py; docs/production_run.md
  # section 10)
  add_test(NAME cosmology_ics
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cosmology_check.py ics
                   ${CMAKE_CURRENT_BINARY_DIR}/cosmology_ics ${PROJECT_SOURCE_DIR}/tools/ics/make_ics.py)
  add_test(NAME cosmology_planck_rt
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cosmology_check.py rt
                   ${CMAKE_CURRENT_BINARY_DIR}/cosmology_rt $<TARGET_FILE:rt_driver>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src)
  add_test(NAME cosmology_floor_cooling
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cosmology_check.py floor
                   ${CMAKE_CURRENT_BINARY_DIR}/cosmology_floor $<TARGET_FILE:rt_driver>)
  set_tests_properties(cosmology_ics cosmology_planck_rt cosmology_floor_cooling
                       PROPERTIES TIMEOUT 900 LABELS "physics;cosmology")

  # End-to-end: run all problems (N=32) and sanity-check the outputs
  add_test(NAME end_to_end
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/e2e_check.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems
                   ${CMAKE_CURRENT_BINARY_DIR}/e2e_runs)
  set_tests_properties(end_to_end PROPERTIES TIMEOUT 3600 LABELS e2e)
  add_test(NAME end_to_end_legacy0
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/e2e_check.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems
                   ${CMAKE_CURRENT_BINARY_DIR}/e2e_runs_legacy0 --set legacy_rt=0)
  set_tests_properties(end_to_end_legacy0 PROPERTIES TIMEOUT 3600 LABELS e2e)
  # Recommended configuration: RT bookkeeping and atomic-physics fixes on
  add_test(NAME end_to_end_fixed
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/e2e_check.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems
                   ${CMAKE_CURRENT_BINARY_DIR}/e2e_runs_fixed
                   --set legacy_rt=0 --set legacy_physics=0)
  set_tests_properties(end_to_end_fixed PROPERTIES TIMEOUT 3600 LABELS e2e)

  # Strict comparison with the stored Fortran harness outputs
  # (reference/harness/outputs: stats for all steps, fields for step 10).
  foreach(p uniform_1src uniform_8src lognormal_16src uniform_8src_smallbuf
            lognormal_16src_smallbuf)
    set(refdir ${PROJECT_SOURCE_DIR}/reference/harness/outputs/${p})
    if(EXISTS ${refdir})
      add_test(NAME ref_${p}
               COMMAND ${CMAKE_COMMAND}
                 -DDRIVER=$<TARGET_FILE:rt_driver>
                 -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/${p}
                 -DOUT=${CMAKE_CURRENT_BINARY_DIR}/ref_runs/${p}
                 -DREF=${refdir}
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCOMPARE=${CMAKE_CURRENT_SOURCE_DIR}/compare_reference.py
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/run_and_compare.cmake)
      set_tests_properties(ref_${p} PROPERTIES TIMEOUT 1800 LABELS reference)
    endif()
    # Same problem with the RT bookkeeping fixes on (legacy_rt = 0), against
    # reference/harness/outputs_legacy_rt0 (same subset as above).
    set(refdir0 ${PROJECT_SOURCE_DIR}/reference/harness/outputs_legacy_rt0/${p})
    if(EXISTS ${refdir0})
      add_test(NAME ref_legacy0_${p}
               COMMAND ${CMAKE_COMMAND}
                 -DDRIVER=$<TARGET_FILE:rt_driver>
                 -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/${p}
                 -DOUT=${CMAKE_CURRENT_BINARY_DIR}/ref_runs_legacy0/${p}
                 -DREF=${refdir0}
                 -DSET=legacy_rt=0
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCOMPARE=${CMAKE_CURRENT_SOURCE_DIR}/compare_reference.py
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/run_and_compare.cmake)
      set_tests_properties(ref_legacy0_${p} PROPERTIES TIMEOUT 1800 LABELS reference)
    endif()
    # Recommended configuration, legacy_rt = 0 + legacy_physics = 0 (atomic
    # physics fixes P1-P9), against reference/harness/outputs_legacy_rt0_phys0.
    set(refdir00 ${PROJECT_SOURCE_DIR}/reference/harness/outputs_legacy_rt0_phys0/${p})
    if(EXISTS ${refdir00})
      add_test(NAME ref_fixed_${p}
               COMMAND ${CMAKE_COMMAND}
                 -DDRIVER=$<TARGET_FILE:rt_driver>
                 -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/${p}
                 -DOUT=${CMAKE_CURRENT_BINARY_DIR}/ref_runs_fixed/${p}
                 -DREF=${refdir00}
                 -DSET=legacy_rt=0,legacy_physics=0
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCOMPARE=${CMAKE_CURRENT_SOURCE_DIR}/compare_reference.py
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/run_and_compare.cmake)
      set_tests_properties(ref_fixed_${p} PROPERTIES TIMEOUT 1800 LABELS reference)
    endif()
  endforeach()
  # Ray dump (compression study) on: still bit-identical to the Fortran, and
  # the codec prototype round-trips the snapshots (docs/ray_compression.md).
  foreach(cfg "legacy;outputs;" "fixed;outputs_legacy_rt0_phys0;legacy_rt=0,legacy_physics=0")
    list(GET cfg 0 tag)
    list(GET cfg 1 refsub)
    list(GET cfg 2 sets)
    set(p lognormal_16src_smallbuf)
    set(refdir ${PROJECT_SOURCE_DIR}/reference/harness/${refsub}/${p})
    if(EXISTS ${refdir})
      add_test(NAME ray_dump_${tag}_${p}
               COMMAND ${CMAKE_COMMAND}
                 -DDRIVER=$<TARGET_FILE:rt_driver>
                 -DCODEC=$<TARGET_FILE:ray_codec>
                 -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/${p}
                 -DOUT=${CMAKE_CURRENT_BINARY_DIR}/ray_dump_runs/${tag}_${p}
                 -DREF=${refdir}
                 -DSET=${sets}
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCOMPARE=${CMAKE_CURRENT_SOURCE_DIR}/compare_reference.py
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/run_ray_dump.cmake)
      set_tests_properties(ray_dump_${tag}_${p} PROPERTIES TIMEOUT 1800 LABELS "reference;raycomp")
    endif()
  endforeach()

  # Movie tool (tools/movies, docs/movies.md) on a synthetic run directory;
  # skipped (77) without matplotlib
  add_test(NAME movie_tool
           COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/movies/test_make_movie.py
                   ${CMAKE_CURRENT_BINARY_DIR}/movie_tool)
  set_tests_properties(movie_tool PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 300)

  # Zoom ICs from MUSIC (tools/ics/music_convert.py, music_ladder.py;
  # docs/zoom_turbulence.md section 4.1): the converter and the generator on
  # a committed two-level MUSIC zoom (numpy only), and the code run from its
  # base files. With -DRADRAY_MUSIC=<MUSIC binary> (klone,
  # tools/hyak/setup_music.sh) also T0 against make_ics.py and the ladder
  # property (label music).
  add_test(NAME music_zoom
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/music_zoom_check.py
                   ${CMAKE_CURRENT_BINARY_DIR}/music_zoom ${PROJECT_SOURCE_DIR} $<TARGET_FILE:radhydro>)
  set_tests_properties(music_zoom PROPERTIES TIMEOUT 300 ENVIRONMENT OMP_NUM_THREADS=1)
  set(RADRAY_MUSIC "" CACHE FILEPATH "MUSIC binary for the music-labelled tests")
  if(RADRAY_MUSIC)
    add_test(NAME music_t0
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/music/t0_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/music_t0 ${RADRAY_MUSIC} $<TARGET_FILE:radhydro>
                     ${PROJECT_SOURCE_DIR})
    set_tests_properties(music_t0 PROPERTIES TIMEOUT 900 LABELS music ENVIRONMENT OMP_NUM_THREADS=1)
    add_test(NAME music_ladder
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/music/ladder_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/music_ladder ${RADRAY_MUSIC} ${PROJECT_SOURCE_DIR})
    set_tests_properties(music_ladder PROPERTIES TIMEOUT 900 LABELS music)
  endif()

  # I-front in a uniform medium, c_reduced = 0.1 vs 1, against the finite-c
  # (retarded-time / Shapiro et al. 2006) solution (docs/production_run.md)
  add_test(NAME ifront_rsla
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ifront_rsla_check.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                   ${CMAKE_CURRENT_BINARY_DIR}/ifront_rsla)
  set_tests_properties(ifront_rsla PROPERTIES TIMEOUT 600 LABELS physics)
  # the same with the moment RT (docs/rt_moments.md); the first-order
  # moment fronts are checked from R = 5 cells (numerical diffusion)
  foreach(m m1 otvet)
    add_test(NAME ifront_rsla_${m}
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ifront_rsla_check.py
                     $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                     ${CMAKE_CURRENT_BINARY_DIR}/ifront_rsla_${m} --set rt_method=${m} --rmin 5)
    set_tests_properties(ifront_rsla_${m} PROPERTIES TIMEOUT 600 LABELS "physics;moments")
  endforeach()

  # Plane-parallel RT-domain mode (docs/plane_parallel.md, plane_parallel = 1):
  # the Gamma_HI normalisation and the bins, slab attenuation, one front
  # across the box against the photon budget (thin and with recombinations),
  # the light-limited front, photon conservation, the domains and the thread
  # independence, on the RT-only driver (tests/plane_check.py)
  add_test(NAME plane_rt
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/plane_check.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                   ${CMAKE_CURRENT_BINARY_DIR}/plane_rt)
  set_tests_properties(plane_rt PROPERTIES TIMEOUT 1800 LABELS "physics;plane")
  # the same physics with plane_transport = sweep (the packets in flight as a
  # field per face instead of rays), and the sweep against the rays
  add_test(NAME plane_rt_sweep
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/plane_check.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                   ${CMAKE_CURRENT_BINARY_DIR}/plane_rt_sweep --only gamma slab front conserve domains threads
                   --set plane_transport=sweep)
  set_tests_properties(plane_rt_sweep PROPERTIES TIMEOUT 1800 LABELS "physics;plane")

  # rays crossing k = 2 cells per substep (docs/rt_multicell.md): the
  # uniform I-front of ifront_rsla at c_reduced = 1 against k = 1 (within
  # half a cell), the finite-c solution and the interior x_HI (the
  # crossing-speed factor of the photo-ionisation rate); and photo_cell's
  # photon budget (B17, photon_conserve)
  add_test(NAME multicell
           COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/analysis/rt_multicell_ifront.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                   ${CMAKE_CURRENT_BINARY_DIR}/multicell --k 1 2 --c 1 --check)
  set_tests_properties(multicell PROPERTIES TIMEOUT 600 LABELS "physics;multicell")
  # the exposure deposit with the converged photon limit at k = 4: the front
  # profile stays that of k = 1 (docs/rt_multicell.md)
  add_test(NAME multicell_exposure
           COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/analysis/rt_multicell_ifront.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                   ${CMAKE_CURRENT_BINARY_DIR}/multicell_exposure --k 1 4 --c 1
                   --set photon_conserve=1 --set ray_deposit=exposure --check)
  set_tests_properties(multicell_exposure PROPERTIES TIMEOUT 600 LABELS "physics;multicell")
  add_test(NAME photon_conserve
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/photon_conserve_check.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src
                   ${CMAKE_CURRENT_BINARY_DIR}/photon_conserve)
  set_tests_properties(photon_conserve PROPERTIES TIMEOUT 900 LABELS "physics;multicell")

  # Stromgren sphere: case A, case B and case A + recombination radiation
  # (docs/recomb_radiation.md section 3); ~10 min with 4 threads
  if(RADRAY_NIGHTLY_TESTS)
    add_test(NAME recomb_stromgren
             COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/analysis/recomb_stromgren.py
                     $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                     ${CMAKE_CURRENT_BINARY_DIR}/recomb_stromgren --check --threads 4)
    set_tests_properties(recomb_stromgren PROPERTIES TIMEOUT 3600
                         LABELS "nightly;physics;recomb" ENVIRONMENT OMP_NUM_THREADS=4)
    foreach(m m1 otvet)
      add_test(NAME recomb_stromgren_${m}
               COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/analysis/recomb_stromgren.py
                       $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                       ${CMAKE_CURRENT_BINARY_DIR}/recomb_stromgren --method ${m} --check --threads 4)
      set_tests_properties(recomb_stromgren_${m} PROPERTIES TIMEOUT 3600
                           LABELS "nightly;physics;recomb;moments" ENVIRONMENT OMP_NUM_THREADS=4)
    endforeach()
  endif()

  # Grid-aligned photon deficit B16 and its fix ray_rotate (docs/port_notes.md)
  add_test(NAME ray_rotate
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ray_rotate_check.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                   ${CMAKE_CURRENT_BINARY_DIR}/ray_rotate)
  set_tests_properties(ray_rotate PROPERTIES TIMEOUT 600 LABELS physics)

  # Redshift per distance travelled with a reduced speed of light
  # (redshift_per_distance, docs/production_run.md): rays and m1 smoke runs
  add_test(NAME redshift_distance_smoke
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/redshift_distance_check.py
                   $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                   ${CMAKE_CURRENT_BINARY_DIR}/redshift_distance)
  set_tests_properties(redshift_distance_smoke PROPERTIES TIMEOUT 600 LABELS physics)

  # cool_suppress_overdensity (Stage F, docs/stage_f.md section 2): 200 is
  # the original (bitwise the default); 0 switches the dense-gas cooling
  # suppression off and must run
  add_test(NAME cool_suppress_default_identity
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/backend_pair_check.py
                   ${CMAKE_CURRENT_BINARY_DIR}/cool_suppress_pair $<TARGET_FILE:radhydro>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                   "legacy_rt=0,legacy_physics=0" "legacy_rt=0,legacy_physics=0,cool_suppress_overdensity=200"
                   --quiet --set fields_every=1 --set hydro_scheme=unsplit)
  set_tests_properties(cool_suppress_default_identity PROPERTIES TIMEOUT 900 LABELS physics)
  add_test(NAME cool_suppress_off_smoke
           COMMAND radhydro ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                   ${CMAKE_CURRENT_BINARY_DIR}/cool_suppress_off --quiet --set legacy_rt=0
                   --set legacy_physics=0 --set hydro_scheme=unsplit --set cool_suppress_overdensity=0)
  set_tests_properties(cool_suppress_off_smoke PROPERTIES TIMEOUT 900 LABELS physics)

  # The runtime cosmology and the preheating floor (docs/production_run.md
  # section 10): the original's values written explicitly, T_min = 0 and a
  # floor at redshifts the run never reaches are bit-identical to the
  # default (both flag settings); a Planck run differs only where it should
  # (the ICs, the clock, the RT substeps, the units); the floor on a cooling
  # test (tests/cosmology_check.py)
  add_test(NAME cosmology_default_identity
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/backend_pair_check.py
                   ${CMAKE_CURRENT_BINARY_DIR}/cosmology_pair $<TARGET_FILE:radhydro>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                   ""
                   "omegam=0.27,omegab=0.045,hubble0=0.70,sigma8=0.80,ns=0.96,T_min=0,z_floor=1000"
                   --quiet --set fields_every=1)
  add_test(NAME cosmology_default_identity_fixes
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/backend_pair_check.py
                   ${CMAKE_CURRENT_BINARY_DIR}/cosmology_pair_fixes $<TARGET_FILE:radhydro>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                   "legacy_rt=0,legacy_physics=0"
                   "legacy_rt=0,legacy_physics=0,omegam=0.27,omegab=0.045,hubble0=0.7,T_min=1000,z_floor=5"
                   --quiet --set fields_every=1 --set hydro_scheme=unsplit)
  add_test(NAME cosmology_planck_full
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cosmology_check.py full
                   ${CMAKE_CURRENT_BINARY_DIR}/cosmology_full $<TARGET_FILE:radhydro>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke ${PROJECT_SOURCE_DIR}/tools/ics/make_ics.py)
  set_tests_properties(cosmology_default_identity cosmology_default_identity_fixes cosmology_planck_full
                       PROPERTIES TIMEOUT 900 LABELS "physics;cosmology")

  # Reduced speed of light (c_reduced = 0.1, docs/production_run.md) with the
  # recommended flags, against reference/harness/outputs_c0.1_fixed.
  foreach(p uniform_8src lognormal_16src_smallbuf)
    set(refdirc ${PROJECT_SOURCE_DIR}/reference/harness/outputs_c0.1_fixed/${p})
    if(EXISTS ${refdirc})
      add_test(NAME ref_c01_${p}
               COMMAND ${CMAKE_COMMAND}
                 -DDRIVER=$<TARGET_FILE:rt_driver>
                 -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/${p}
                 -DOUT=${CMAKE_CURRENT_BINARY_DIR}/ref_runs_c01/${p}
                 -DREF=${refdirc}
                 -DSET=c_reduced=0.1,legacy_rt=0,legacy_physics=0
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCOMPARE=${CMAKE_CURRENT_SOURCE_DIR}/compare_reference.py
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/run_and_compare.cmake)
      set_tests_properties(ref_c01_${p} PROPERTIES TIMEOUT 1800 LABELS reference)
    endif()
  endforeach()

  # Stored RT references that differ on platforms with another libm (UW Hyak
  # klone, glibc 2.28): last-place differences in the double sums of
  # rt_stats.txt and rays at the removal threshold; the gas fields are
  # byte-identical, and the gfortran harness built there shows the same
  # differences (docs/cluster.md section 3). Label platform_ref:
  # tools/hyak/build_hyak.sh leaves these out (user decision 2026-10-12);
  # CI (Ubuntu 24.04, where the references were made) runs them.
  foreach(t ref_lognormal_16src ref_lognormal_16src_smallbuf
            ref_legacy0_lognormal_16src ref_legacy0_lognormal_16src_smallbuf
            ref_fixed_lognormal_16src ref_fixed_lognormal_16src_smallbuf
            ref_c01_uniform_8src ref_c01_lognormal_16src_smallbuf
            ray_dump_legacy_lognormal_16src_smallbuf ray_dump_fixed_lognormal_16src_smallbuf)
    if(TEST ${t})
      set_property(TEST ${t} APPEND PROPERTY LABELS platform_ref)
    endif()
  endforeach()
endif()

# ---------------------------------------------------------------------------
# PM gravity (only when FFTW3f was found, see the top-level CMakeLists.txt)
if(RADRAY_HAVE_PM)
  add_executable(test_pm test_pm.cpp)
  target_link_libraries(test_pm PRIVATE radray_pm)
  add_test(NAME pm_unit COMMAND test_pm)
  set_tests_properties(pm_unit PROPERTIES LABELS pm)
  # pm_refine > 1 (docs/pm_refine.md): forces against Ewald, momentum
  add_executable(test_pm_refine test_pm_refine.cpp)
  target_link_libraries(test_pm_refine PRIVATE radray_pm)
  add_test(NAME pm_refine_unit COMMAND test_pm_refine)
  set_tests_properties(pm_refine_unit PROPERTIES LABELS pm)
  # threaded FFTW: byte-identical for any thread count
  add_executable(test_fft_threads test_fft_threads.cpp)
  target_link_libraries(test_fft_threads PRIVATE radray_fft)
  add_test(NAME fft_threads COMMAND test_fft_threads 64 1 2 3 4 8 16)
  set_tests_properties(fft_threads PROPERTIES LABELS "pm;threads")
  # the moment level's threaded OTVET FFT plans: bitwise the serial ones
  # (the guard, docs/stage_e_moments.md 4.1) at 64^3 and 128^3
  if(FFTW3F_OMP_LIBRARY AND OpenMP_CXX_FOUND)
    add_executable(test_otvet_fft_threads test_otvet_fft_threads.cpp)
    target_include_directories(test_otvet_fft_threads PRIVATE ${FFTW3F_INCLUDE_DIR})
    target_link_libraries(test_otvet_fft_threads PRIVATE radray ${FFTW3F_OMP_LIBRARY} ${FFTW3F_LIBRARY}
                          OpenMP::OpenMP_CXX)
    foreach(n 64 128)
      add_test(NAME otvet_fft_threads_${n} COMMAND test_otvet_fft_threads ${n} 2 4 8)
      set_tests_properties(otvet_fft_threads_${n} PROPERTIES LABELS "moments;threads")
    endforeach()
  endif()

  if(Python3_Interpreter_FOUND)
    set(PM_PROBLEMS zeldovich_3d growth_3d pancake_x pancake_z_nomf
                    zeldovich_3d_r2 pancake_x_r2 pancake_z_nomf_r2 pancake_x_n32)
    foreach(p ${PM_PROBLEMS})
      set(run ${CMAKE_CURRENT_BINARY_DIR}/pm_runs/${p})
      # 1. run pm_driver (fixture shared by the checks below)
      add_test(NAME pm_run_${p}
               COMMAND ${CMAKE_COMMAND} -DDRIVER=$<TARGET_FILE:pm_driver>
                       -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/${p} -DOUT=${run}
                       -P ${CMAKE_CURRENT_SOURCE_DIR}/run_pm.cmake)
      set_tests_properties(pm_run_${p} PROPERTIES FIXTURES_SETUP pm_${p} LABELS pm TIMEOUT 600)
      # 2. bit-for-bit against the stored Fortran harness outputs
      set(refdir ${PROJECT_SOURCE_DIR}/reference/harness_pm/outputs/${p})
      if(EXISTS ${refdir})
        add_test(NAME pm_ref_${p}
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_pm_reference.py
                         ${refdir} ${run})
        set_tests_properties(pm_ref_${p} PROPERTIES FIXTURES_REQUIRED pm_${p}
                             LABELS "pm;reference")
      endif()
      # 3. momentum conservation (all problems)
      add_test(NAME pm_momentum_${p}
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                       momentum ${run} ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/${p})
      set_tests_properties(pm_momentum_${p} PROPERTIES FIXTURES_REQUIRED pm_${p} LABELS pm)
    endforeach()
    # 4. analytic Zel'dovich pancakes and linear growth
    add_test(NAME pm_pancake_x
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                     pancake ${CMAKE_CURRENT_BINARY_DIR}/pm_runs/pancake_x
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/pancake_x)
    set_tests_properties(pm_pancake_x PROPERTIES FIXTURES_REQUIRED pm_pancake_x LABELS pm)
    # fixed da_sim = 0.002 from z = 300: the first step has da/a = 0.6, so the
    # time-integration error of the mode is ~1.3% (see docs/port_notes_pm.md)
    add_test(NAME pm_pancake_z_nomf
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                     pancake ${CMAKE_CURRENT_BINARY_DIR}/pm_runs/pancake_z_nomf
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/pancake_z_nomf
                     --tol-mode-x 0.025 --tol-mode-v 0.025)
    set_tests_properties(pm_pancake_z_nomf PROPERTIES FIXTURES_REQUIRED pm_pancake_z_nomf
                         LABELS pm)
    # pm_refine = 2 (docs/pm_refine.md): the same analytic pancakes with
    # Ndm = 2 Ngas on the 2 Ngas mesh. One particle per mesh cell (pancake_x
    # has one per 2): the per-particle error near the caustic is 4.0e-2,
    # exactly that of the original with Ngas = Ndm = 32 (the twin below).
    add_test(NAME pm_pancake_x_r2
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                     pancake ${CMAKE_CURRENT_BINARY_DIR}/pm_runs/pancake_x_r2
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/pancake_x_r2 --tol-x 0.045)
    set_tests_properties(pm_pancake_x_r2 PROPERTIES FIXTURES_REQUIRED pm_pancake_x_r2 LABELS pm)
    # with uniform frozen gas the refined mesh moves the DM bit for bit as
    # the original on a gas grid as fine as the mesh
    add_test(NAME pm_r2_twin_pancake_x
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                     same-dm ${CMAKE_CURRENT_BINARY_DIR}/pm_runs/pancake_x_r2
                     ${CMAKE_CURRENT_BINARY_DIR}/pm_runs/pancake_x_n32)
    set_tests_properties(pm_r2_twin_pancake_x PROPERTIES
                         FIXTURES_REQUIRED "pm_pancake_x_r2;pm_pancake_x_n32" LABELS pm)
    add_test(NAME pm_pancake_x_n32
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                     pancake ${CMAKE_CURRENT_BINARY_DIR}/pm_runs/pancake_x_n32
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/pancake_x_n32 --tol-x 0.045)
    set_tests_properties(pm_pancake_x_n32 PROPERTIES FIXTURES_REQUIRED pm_pancake_x_n32 LABELS pm)
    add_test(NAME pm_pancake_z_nomf_r2
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                     pancake ${CMAKE_CURRENT_BINARY_DIR}/pm_runs/pancake_z_nomf_r2
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/pancake_z_nomf_r2
                     --tol-x 0.045 --tol-mode-x 0.025 --tol-mode-v 0.025)
    set_tests_properties(pm_pancake_z_nomf_r2 PROPERTIES FIXTURES_REQUIRED pm_pancake_z_nomf_r2
                         LABELS pm)
    add_test(NAME pm_growth_3d
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                     growth ${CMAKE_CURRENT_BINARY_DIR}/pm_runs/growth_3d
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/growth_3d)
    set_tests_properties(pm_growth_3d PROPERTIES FIXTURES_REQUIRED pm_growth_3d LABELS pm)
  endif()
endif()

# Hydro needs FFTW3f as well (built together with PM)
if(RADRAY_HAVE_PM)
# ---------------------------------------------------------------------------
# Hydro (hydro_tools.f90 -DMF port): see docs/port_notes_hydro.md
# ---------------------------------------------------------------------------
if(Python3_Interpreter_FOUND)
  set(HYDRO_PROBLEMS ${CMAKE_CURRENT_SOURCE_DIR}/problems_hydro)
  set(HYDRO_REF ${PROJECT_SOURCE_DIR}/reference/harness_hydro/outputs)
  set(HYDRO_RUN ${CMAKE_COMMAND} -DPYTHON=${Python3_EXECUTABLE}
                -DTESTS=${CMAKE_CURRENT_SOURCE_DIR}
                -DCOMPARE=${PROJECT_SOURCE_DIR}/tools/compare/compare.py)

  # Strict comparison with the stored Fortran hydro-harness outputs
  # (stats and summaries for all steps, fields for the last step).
  foreach(p sod_x sod_y sod_z sound_x random3d random3d_freeze)
    if(EXISTS ${HYDRO_REF}/${p})
      add_test(NAME hydro_ref_${p}
               COMMAND ${HYDRO_RUN} -DMODE=ref -DDRIVER=$<TARGET_FILE:hydro_driver>
                       -DPROBLEM=${HYDRO_PROBLEMS}/${p} -DREF=${HYDRO_REF}/${p}
                       -DOUT=${CMAKE_CURRENT_BINARY_DIR}/hydro_ref_runs/${p}
                       -P ${CMAKE_CURRENT_SOURCE_DIR}/run_hydro.cmake)
      set_tests_properties(hydro_ref_${p} PROPERTIES TIMEOUT 600 LABELS "hydro;reference")
    endif()
  endforeach()

  # Physics sanity: Sod vs exact Riemann, sound speed, conservation
  add_test(NAME hydro_physics
           COMMAND ${HYDRO_RUN} -DMODE=physics -DDRIVER=$<TARGET_FILE:hydro_driver>
                   -DPROBLEMS=${HYDRO_PROBLEMS}
                   -DOUT=${CMAKE_CURRENT_BINARY_DIR}/hydro_physics_runs
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/run_hydro.cmake)
  set_tests_properties(hydro_physics PROPERTIES TIMEOUT 600 LABELS hydro)

  # Live comparison with a fresh gfortran harness build (all steps, all
  # fields), if gfortran and the Fortran reference (reference/fortran) are
  # available.
  find_program(GFORTRAN_EXE gfortran)
  if(GFORTRAN_EXE AND RADRAY_HAVE_FORTRAN_REF)
    add_test(NAME hydro_vs_fortran_live
             COMMAND ${HYDRO_RUN} -DMODE=fortran -DDRIVER=$<TARGET_FILE:hydro_driver>
                     -DPROBLEMS=${HYDRO_PROBLEMS}
                     -DHARNESS=${RADRAY_FORTRAN_DIR}/harness_hydro
                     -DOUT=${CMAKE_CURRENT_BINARY_DIR}/hydro_live_runs
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/run_hydro.cmake)
    set_tests_properties(hydro_vs_fortran_live PROPERTIES TIMEOUT 1200 LABELS "hydro;fortran")
  endif()

  # ASan/UBSan build of the hydro driver (all radray sources), run on two
  # problems for a few steps; the output must still match the reference.
  if(NOT RADRAY_SANITIZE)
    file(GLOB RADRAY_ALL_SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)
    list(APPEND RADRAY_ALL_SRC ${PROJECT_SOURCE_DIR}/src/amrex/rt_level_stub.cpp
                            ${PROJECT_SOURCE_DIR}/src/amrex/pm_level_stub.cpp)
    list(APPEND RADRAY_ALL_SRC ${PROJECT_SOURCE_DIR}/src/amrex/hydro_level_stub.cpp
                            ${PROJECT_SOURCE_DIR}/src/amrex/gas_state_stub.cpp)
    add_executable(hydro_driver_asan ${PROJECT_SOURCE_DIR}/apps/hydro_driver.cpp ${RADRAY_ALL_SRC})
    target_include_directories(hydro_driver_asan PRIVATE ${PROJECT_SOURCE_DIR}/include
                               ${FFTW3F_INCLUDE_DIR})
    target_compile_options(hydro_driver_asan PRIVATE -fsanitize=address,undefined
                           -fno-sanitize-recover=undefined -fno-omit-frame-pointer -g -O1)
    target_link_options(hydro_driver_asan PRIVATE -fsanitize=address,undefined)
    target_link_libraries(hydro_driver_asan PRIVATE ${FFTW3F_LIBRARY})
    if(RADRAY_HAVE_OPENMP)
      target_link_libraries(hydro_driver_asan PRIVATE OpenMP::OpenMP_CXX)
    endif()
    set(HYDRO_ASAN_DRIVER $<TARGET_FILE:hydro_driver_asan>)
  else()
    set(HYDRO_ASAN_DRIVER $<TARGET_FILE:hydro_driver>)
  endif()
  foreach(p sod_y random3d_freeze)
    add_test(NAME hydro_sanitize_${p}
             COMMAND ${HYDRO_RUN} -DMODE=sanitize -DDRIVER=${HYDRO_ASAN_DRIVER} -DNSTEPS=3
                     -DPROBLEM=${HYDRO_PROBLEMS}/${p} -DREF=${HYDRO_REF}/${p}
                     -DOUT=${CMAKE_CURRENT_BINARY_DIR}/hydro_asan_runs/${p}
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/run_hydro.cmake)
    set_tests_properties(hydro_sanitize_${p} PROPERTIES TIMEOUT 1200 LABELS "hydro;sanitize")
  endforeach()
endif()
endif()

# ---------------------------------------------------------------------------
# Full radhydro loop (apps/radhydro vs reference/harness_full): docs/full_run.md
# ---------------------------------------------------------------------------
if(RADRAY_HAVE_PM AND Python3_Interpreter_FOUND)
  set(FULL_PROBLEMS ${CMAKE_CURRENT_SOURCE_DIR}/problems_full)
  set(FULL_REF ${PROJECT_SOURCE_DIR}/reference/harness_full/outputs)
  set(FULL_RUN ${CMAKE_COMMAND} -DPYTHON=${Python3_EXECUTABLE}
               -DTESTS=${CMAKE_CURRENT_SOURCE_DIR}
               -DCOMPARE=${PROJECT_SOURCE_DIR}/tools/compare/compare.py)
  set(FULL_SET_rt1 "legacy_rt=1,legacy_physics=1")
  set(FULL_SET_rt0 "legacy_rt=0,legacy_physics=0")

  # Strict comparison with the stored Fortran full-harness outputs (stats of
  # every step, SHA-256 of every field file, hydro.log, P(k), particles).
  # smoke_*: seconds; legacy_* / onthefly_*: ~4 minutes each (label nightly,
  # built only with RADRAY_NIGHTLY_TESTS=ON).
  foreach(c smoke_rt1 smoke_rt0 legacy_rt1 legacy_rt0 onthefly_rt1 onthefly_rt0)
    string(REGEX REPLACE "_rt[01]$" "" p ${c})
    string(REGEX REPLACE "^.*_" "" m ${c})
    if(EXISTS ${FULL_REF}/${c})
      add_test(NAME full_ref_${c}
               COMMAND ${FULL_RUN} -DMODE=ref -DDRIVER=$<TARGET_FILE:radhydro>
                       -DPROBLEM=${FULL_PROBLEMS}/${p} -DSET=${FULL_SET_${m}}
                       -DREF=${FULL_REF}/${c} -DOUT=${CMAKE_CURRENT_BINARY_DIR}/full_ref_runs/${c}
                       -P ${CMAKE_CURRENT_SOURCE_DIR}/run_full.cmake)
      if(p STREQUAL "smoke")
        set_tests_properties(full_ref_${c} PROPERTIES TIMEOUT 600 LABELS "full;reference")
      else()
        # ~5 minutes each: label nightly, run only with RADRAY_NIGHTLY_TESTS=ON
        set_tests_properties(full_ref_${c} PROPERTIES TIMEOUT 3600
                             LABELS "full;reference;nightly")
        if(NOT RADRAY_NIGHTLY_TESTS)
          set_tests_properties(full_ref_${c} PROPERTIES DISABLED TRUE)
        endif()
      endif()
    endif()
  endforeach()

  # Live comparison with a fresh gfortran build of the full radhydro.f90
  # (smoke problem, both source modes and both flag sets), if gfortran and
  # the Fortran reference (reference/fortran) are available.
  find_program(GFORTRAN_FULL gfortran)
  if(GFORTRAN_FULL AND RADRAY_HAVE_FORTRAN_REF)
    foreach(p smoke)
      foreach(src legacy onthefly)
        foreach(m rt1 rt0)
          add_test(NAME full_vs_fortran_live_${p}_${src}_${m}
                   COMMAND ${FULL_RUN} -DMODE=fortran -DDRIVER=$<TARGET_FILE:radhydro>
                           -DPROBLEM=${FULL_PROBLEMS}/${p}
                           -DSET=source_mode=${src},${FULL_SET_${m}}
                           -DHARNESS=${RADRAY_FORTRAN_DIR}/harness_full
                           -DRADRAY_BUILD=${CMAKE_BINARY_DIR}
                           -DOUT=${CMAKE_CURRENT_BINARY_DIR}/full_live_runs/${p}_${src}_${m}
                           -P ${CMAKE_CURRENT_SOURCE_DIR}/run_full.cmake)
          set_tests_properties(full_vs_fortran_live_${p}_${src}_${m} PROPERTIES
                               TIMEOUT 1800 LABELS "full;fortran"
                               RESOURCE_LOCK harness_full_build)
        endforeach()
      endforeach()
    endforeach()
    # F3 ray pool guard (legacy_rt = 0): legacy problem with
    # nray_max_factor = 1, which overflows the pool without the guard; the
    # guard limits splitting in ~1000 of the 1669 RT substeps, identically in
    # both codes. ~3 minutes: nightly.
    add_test(NAME full_vs_fortran_live_legacy_poolguard
             COMMAND ${FULL_RUN} -DMODE=fortran -DDRIVER=$<TARGET_FILE:radhydro>
                     -DPROBLEM=${FULL_PROBLEMS}/legacy
                     -DSET=source_mode=legacy,${FULL_SET_rt0},nray_max_factor=1
                     -DHARNESS=${RADRAY_FORTRAN_DIR}/harness_full
                     -DRADRAY_BUILD=${CMAKE_BINARY_DIR}
                     -DOUT=${CMAKE_CURRENT_BINARY_DIR}/full_live_runs/legacy_poolguard
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/run_full.cmake)
    set_tests_properties(full_vs_fortran_live_legacy_poolguard PROPERTIES TIMEOUT 3600
                         LABELS "full;fortran;nightly" RESOURCE_LOCK harness_full_build)
    if(NOT RADRAY_NIGHTLY_TESTS)
      set_tests_properties(full_vs_fortran_live_legacy_poolguard PROPERTIES DISABLED TRUE)
    endif()
    # Reduced speed of light: c_reduced = 0.1 (both flag sets) and 1e-4,
    # where nstep_rt rounds to 0 and the nstep_rt >= 1 guard is used.
    foreach(cs "legacy;rt1;0.1" "onthefly;rt0;0.1" "legacy;rt0;0.0001")
      list(GET cs 0 src)
      list(GET cs 1 m)
      list(GET cs 2 cr)
      set(tname full_vs_fortran_live_smoke_${src}_${m}_c${cr})
      add_test(NAME ${tname}
               COMMAND ${FULL_RUN} -DMODE=fortran -DDRIVER=$<TARGET_FILE:radhydro>
                       -DPROBLEM=${FULL_PROBLEMS}/smoke
                       -DSET=source_mode=${src},${FULL_SET_${m}},c_reduced=${cr}
                       -DHARNESS=${RADRAY_FORTRAN_DIR}/harness_full
                       -DRADRAY_BUILD=${CMAKE_BINARY_DIR}
                       -DOUT=${CMAKE_CURRENT_BINARY_DIR}/full_live_runs/smoke_${src}_${m}_c${cr}
                       -P ${CMAKE_CURRENT_SOURCE_DIR}/run_full.cmake)
      set_tests_properties(${tname} PROPERTIES TIMEOUT 1800 LABELS "full;fortran"
                           RESOURCE_LOCK harness_full_build)
    endforeach()
  endif()

  # ASan/UBSan build of radhydro (all radray sources) on the smoke problem;
  # it must still match the reference bit for bit.
  if(NOT RADRAY_SANITIZE)
    file(GLOB RADRAY_ALL_SRC_FULL ${PROJECT_SOURCE_DIR}/src/*.cpp)
    list(APPEND RADRAY_ALL_SRC_FULL ${PROJECT_SOURCE_DIR}/src/amrex/rt_level_stub.cpp
                            ${PROJECT_SOURCE_DIR}/src/amrex/pm_level_stub.cpp)
    list(APPEND RADRAY_ALL_SRC_FULL ${PROJECT_SOURCE_DIR}/src/amrex/hydro_level_stub.cpp
                            ${PROJECT_SOURCE_DIR}/src/amrex/gas_state_stub.cpp)
    add_executable(radhydro_asan ${PROJECT_SOURCE_DIR}/apps/radhydro.cpp ${RADRAY_ALL_SRC_FULL})
    target_include_directories(radhydro_asan PRIVATE ${PROJECT_SOURCE_DIR}/include
                               ${FFTW3F_INCLUDE_DIR})
    target_compile_options(radhydro_asan PRIVATE -fsanitize=address,undefined
                           -fno-sanitize-recover=undefined -fno-omit-frame-pointer -g -O1
                           -Wno-unknown-pragmas)
    target_link_options(radhydro_asan PRIVATE -fsanitize=address,undefined)
    target_link_libraries(radhydro_asan PRIVATE ${FFTW3F_LIBRARY})
    if(RADRAY_HAVE_OPENMP)
      target_link_libraries(radhydro_asan PRIVATE OpenMP::OpenMP_CXX)
    endif()
    set(FULL_ASAN_DRIVER $<TARGET_FILE:radhydro_asan>)
  else()
    set(FULL_ASAN_DRIVER $<TARGET_FILE:radhydro>)
  endif()
  foreach(m rt1 rt0)
    if(EXISTS ${FULL_REF}/smoke_${m})
      add_test(NAME full_sanitize_smoke_${m}
               COMMAND ${FULL_RUN} -DMODE=sanitize -DDRIVER=${FULL_ASAN_DRIVER}
                       -DPROBLEM=${FULL_PROBLEMS}/smoke -DSET=${FULL_SET_${m}}
                       -DREF=${FULL_REF}/smoke_${m}
                       -DOUT=${CMAKE_CURRENT_BINARY_DIR}/full_asan_runs/smoke_${m}
                       -P ${CMAKE_CURRENT_SOURCE_DIR}/run_full.cmake)
      set_tests_properties(full_sanitize_smoke_${m} PROPERTIES TIMEOUT 1800
                           LABELS "full;sanitize")
    endif()
  endforeach()
endif()

# ---------------------------------------------------------------------------
# Thread-count independence (docs/parallel.md) and restart
if(Python3_Interpreter_FOUND)
  # deterministic_parallel ray tracing, N = 32 (4^3 domains of 8^3 cells)
  foreach(p lognormal_16src uniform_8src_smallbuf)
    add_test(NAME thread_rt_det_${p}
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/rt_det_${p} $<TARGET_FILE:rt_driver>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems/${p} --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1 --set nsteps=6)
    set_tests_properties(thread_rt_det_${p} PROPERTIES TIMEOUT 1800 LABELS threads
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
  endforeach()
  # the production RT configuration: deterministic_parallel with the per-substep
  # and per-source frame rotation ray_rotate = 3 (B16 fix, docs/port_notes.md)
  add_test(NAME thread_rt_det_rotate_lognormal_16src
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                   ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/rt_det_rotate $<TARGET_FILE:rt_driver>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src --set legacy_rt=0
                   --set legacy_physics=0 --set deterministic_parallel=1 --set ray_rotate=3
                   --set nsteps=6)
  set_tests_properties(thread_rt_det_rotate_lognormal_16src PROPERTIES TIMEOUT 1800 LABELS threads
                       ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
  # rays crossing k = 3 cells per substep (docs/rt_multicell.md): the
  # largest k with 8^3 domains (reach k+1 = d/2), vs threads
  add_test(NAME thread_rt_det_multicell_lognormal_16src
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                   ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/rt_det_multicell $<TARGET_FILE:rt_driver>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src --set legacy_rt=0
                   --set legacy_physics=0 --set deterministic_parallel=1 --set ray_rotate=3
                   --set rt_cells_per_substep=3 --set nsteps=6)
  set_tests_properties(thread_rt_det_multicell_lognormal_16src PROPERTIES TIMEOUT 1800
                       LABELS "threads;multicell" ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
  # legacy order (rt_step_2 and the RT grid loops in parallel), vs threads
  add_test(NAME thread_rt_legacy_lognormal_16src
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                   ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/rt_legacy $<TARGET_FILE:rt_driver>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src --set nsteps=4)
  set_tests_properties(thread_rt_legacy_lognormal_16src PROPERTIES TIMEOUT 1800 LABELS threads
                       ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
  # moment RT (docs/rt_moments.md): transport, OTVET tensor, physics
  foreach(m m1 otvet)
    add_test(NAME thread_rt_${m}_lognormal_16src
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/rt_${m} $<TARGET_FILE:rt_driver>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src --set nsteps=3
                     --set rt_method=${m} --set legacy_physics=0)
    set_tests_properties(thread_rt_${m}_lognormal_16src PROPERTIES TIMEOUT 1800
                         LABELS "threads;moments" ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
  endforeach()
  # the alternative OTVET tensors (otvet_closure.hpp, docs/otvet_closure.md)
  foreach(cl split screen_split cone)
    if(cl STREQUAL "split")
      set(clargs --set otvet_split=1)
    elseif(cl STREQUAL "screen_split")
      set(clargs --set otvet_closure=screen --set otvet_split=1)
    else()
      set(clargs --set otvet_closure=cone --set otvet_cone_nside=2)
    endif()
    add_test(NAME thread_rt_otvet_${cl}_lognormal_16src
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/rt_otvet_${cl} $<TARGET_FILE:rt_driver>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src --set nsteps=3
                     --set rt_method=otvet --set legacy_physics=0 ${clargs})
    set_tests_properties(thread_rt_otvet_${cl}_lognormal_16src PROPERTIES TIMEOUT 1800
                         LABELS "threads;moments" ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
  endforeach()
  # recombination radiation (docs/recomb_radiation.md): emission, block
  # rotation, pool budget and on-the-spot deposits, vs threads
  add_test(NAME thread_rt_recomb_lognormal_16src
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                   ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/rt_recomb $<TARGET_FILE:rt_driver>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src --set legacy_rt=0
                   --set legacy_physics=0 --set deterministic_parallel=1 --set nsteps=6
                   --set recomb_case=A --set recomb_rad=1)
  set_tests_properties(thread_rt_recomb_lognormal_16src PROPERTIES TIMEOUT 1800
                       LABELS "threads;recomb" ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
  # the same with the moment RT (every cell emits into the moments) and the
  # helium recombination lines
  add_test(NAME thread_rt_recomb_m1_lognormal_16src
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                   ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/rt_recomb_m1 $<TARGET_FILE:rt_driver>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src --set nsteps=3
                   --set rt_method=m1 --set legacy_physics=0
                   --set recomb_case=A --set recomb_rad=1 --set recomb_rad_lines=1)
  set_tests_properties(thread_rt_recomb_m1_lognormal_16src PROPERTIES TIMEOUT 1800
                       LABELS "threads;recomb;moments" ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
  if(RADRAY_HAVE_PM)
    # hydro, PM, summaries, FoF (N = 32, 22 steps before RT)
    add_test(NAME thread_full_legacy_hydro_pm
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/full_hydro_pm $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly --quiet
                     --set zdata1=12 --set fields_every=1)
    set_tests_properties(thread_full_legacy_hydro_pm PROPERTIES TIMEOUT 1800 LABELS threads
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
    # the unsplit hydro scheme with dual energy (docs/hydro_unsplit.md)
    add_test(NAME thread_full_unsplit_hydro_pm
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/full_unsplit $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly --quiet
                     --set zdata1=12 --set fields_every=1 --set hydro_scheme=unsplit)
    set_tests_properties(thread_full_unsplit_hydro_pm PROPERTIES TIMEOUT 1800 LABELS threads
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
    # pm_refine = 2 with Ndm = 2 N (docs/pm_refine.md): the threaded FFTW,
    # the refined mesh, FoF on the finer particles
    add_test(NAME thread_full_pm_refine
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/full_pm_refine $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly_r2 --quiet
                     --set fields_every=1 --set hydro_scheme=unsplit)
    set_tests_properties(thread_full_pm_refine PROPERTIES TIMEOUT 1800 LABELS "threads;pm"
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
    # fof_dt_myr (FoF cadence, docs/pm_refine.md): 14 source updates from
    # z = 12, of which 9 reuse the cached halos (drifted) with fof_dt_myr = 60
    add_test(NAME thread_full_fof_cadence
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/full_fof_cadence $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly_r2 --quiet
                     --set "zdata1=12 7.5" --set ra_max=0.03 --set fof_dt_myr=60
                     --set fields_every=4)
    set_tests_properties(thread_full_fof_cadence PROPERTIES TIMEOUT 1800 LABELS "threads;pm"
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
    # ASan/UBSan builds with OpenMP, thread-count independence (parallel
    # paths: N = 32)
    if(NOT RADRAY_SANITIZE AND TARGET radhydro_asan)
      file(GLOB RADRAY_ALL_SRC_RT ${PROJECT_SOURCE_DIR}/src/*.cpp)
      list(APPEND RADRAY_ALL_SRC_RT ${PROJECT_SOURCE_DIR}/src/amrex/rt_level_stub.cpp
                            ${PROJECT_SOURCE_DIR}/src/amrex/pm_level_stub.cpp)
    list(APPEND RADRAY_ALL_SRC_RT ${PROJECT_SOURCE_DIR}/src/amrex/hydro_level_stub.cpp
                            ${PROJECT_SOURCE_DIR}/src/amrex/gas_state_stub.cpp)
      add_executable(rt_driver_asan ${PROJECT_SOURCE_DIR}/apps/rt_driver.cpp ${RADRAY_ALL_SRC_RT})
      target_include_directories(rt_driver_asan PRIVATE ${PROJECT_SOURCE_DIR}/include
                                 ${FFTW3F_INCLUDE_DIR})
      target_compile_options(rt_driver_asan PRIVATE -fsanitize=address,undefined
                             -fno-sanitize-recover=undefined -fno-omit-frame-pointer -g -O1
                             -Wno-unknown-pragmas)
      target_link_options(rt_driver_asan PRIVATE -fsanitize=address,undefined)
      target_link_libraries(rt_driver_asan PRIVATE ${FFTW3F_LIBRARY})
      if(RADRAY_HAVE_OPENMP)
        target_link_libraries(rt_driver_asan PRIVATE OpenMP::OpenMP_CXX)
      endif()
      add_test(NAME sanitize_thread_rt_det
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                       ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/asan_rt_det
                       $<TARGET_FILE:rt_driver_asan> ${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src
                       --set legacy_rt=0 --set legacy_physics=0 --set deterministic_parallel=1
                       --set nsteps=3)
      add_test(NAME sanitize_thread_full_hydro_pm
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                       ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/asan_full $<TARGET_FILE:radhydro_asan>
                       ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly --quiet
                       --set zdata1=30 --set fields_every=4)
      set_tests_properties(sanitize_thread_rt_det sanitize_thread_full_hydro_pm PROPERTIES
                           TIMEOUT 1800 LABELS "threads;sanitize" RUN_SERIAL TRUE
                           ENVIRONMENT "OMP_NUM_THREADS=1;ASAN_OPTIONS=detect_leaks=1:abort_on_error=1;UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1")
    endif()
    # bit-identical restart: legacy, and deterministic_parallel + movie frames
    add_test(NAME restart_smoke_legacy
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/legacy 13)
    add_test(NAME restart_smoke_det_movie
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/det 14 --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1 --set movie_dz=0.5
                     --set movie_z0=25)
    # the unsplit hydro scheme: the internal-energy field through the
    # checkpoint and the RT coupling (docs/hydro_unsplit.md)
    add_test(NAME restart_smoke_unsplit
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/unsplit 13 --set hydro_scheme=unsplit)
    # c_reduced_schedule (docs/production_run.md section 1): restart from step
    # 13, before the switches at steps 14 and 15 (x from the checkpoint's RT grid)
    add_test(NAME restart_smoke_c_schedule
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/c_schedule 13 --set legacy_rt=0
                     --set legacy_physics=0 --set c_reduced_schedule=0.3,0.005,0.6,0.01,1)
    # xray_heating = 1 (docs/xray_heating.md): stateless, the rates of step 15
    # from the restored state (the default L_UV(M, z) of sfr_model = abundance)
    add_test(NAME restart_smoke_xray
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/xray 14 --set sfr_model=abundance
                     --set xray_heating=1 --set xray_log_lx_sfr=44 --set xray_dlog=0.01)
    # the one-speed schedule is the constant c_reduced byte for byte; the
    # switching rule, step by step
    add_test(NAME c_schedule
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/c_schedule_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/c_schedule --set legacy_rt=0 --set legacy_physics=0)
    set_tests_properties(c_schedule PROPERTIES TIMEOUT 600 LABELS physics)
    # Plane-parallel RT-domain mode in the full loop (docs/plane_parallel.md):
    # the turn-on at plane_z_re (the step before lands on it), DM20's speed
    # of light schedule, the step caps dt_max_myr and rt_substeps_max
    # (docs/zoom_turbulence.md section 4.4) on a uniform 16^3 box
    add_test(NAME plane_full
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/plane_full_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_BINARY_DIR}/plane_full)
    set_tests_properties(plane_full PROPERTIES TIMEOUT 1800 LABELS "physics;plane;full")
    # the scaled-down DM20 problem (0.256 Mpc/h, 32^3, z_re = 8, Gamma_-12 =
    # 0.3, cooling on): the post-front temperature and the clumping factor
    # (~ minutes: label nightly)
    add_test(NAME plane_dm20
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/plane_full_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_BINARY_DIR}/plane_dm20 --dm20)
    set_tests_properties(plane_dm20 PROPERTIES TIMEOUT 3600 LABELS "physics;plane;full;nightly")
    if(NOT RADRAY_NIGHTLY_TESTS)
      set_tests_properties(plane_dm20 PROPERTIES DISABLED TRUE)
    endif()
    # gammaHI_avg (docs/igm_diagnostics.md): the substep-averaged Gamma_HI in
    # fields and movie frames is the mean of the per-substep gammaHI frames,
    # thread-independent and output only; the rays and M1
    add_test(NAME gamma_avg_rays
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gamma_avg_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/gamma_avg/rays --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1)
    set_tests_properties(gamma_avg_rays PROPERTIES TIMEOUT 1200 LABELS "physics;threads")
    add_test(NAME gamma_avg_m1
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gamma_avg_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/gamma_avg/m1 --set legacy_rt=0
                     --set legacy_physics=0 --set rt_method=m1)
    set_tests_properties(gamma_avg_m1 PROPERTIES TIMEOUT 1200 LABELS "physics;threads;moments")
    # thermal_init = recfast (docs/thermal_history.md): the mean-density gas on
    # RECFAST, thread independence, the default keys inert
    add_test(NAME thermal_init_recfast
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thermal_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly
                     ${CMAKE_CURRENT_BINARY_DIR}/thermal_init)
    set_tests_properties(thermal_init_recfast PROPERTIES TIMEOUT 600 LABELS "physics;threads")
    # preheat_floor from preheat_z (docs/thermal_history.md section 5): inert
    # above preheat_z, the mean-density gas on the floor below it, threads
    add_test(NAME preheat_floor
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/preheat_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly
                     ${CMAKE_CURRENT_BINARY_DIR}/preheat_floor)
    set_tests_properties(preheat_floor PROPERTIES TIMEOUT 600 LABELS "physics;threads")
    # dual_e_self = 1 (docs/hydro_unsplit.md section 7, U7): no cell of the
    # cold unheated IGM over-cooled, the mean-density gas on its adiabat,
    # thread independence, the key inert at 0; the default shows the defect
    add_test(NAME unsplit_dual_e_self
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/dual_e_self_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly
                     ${CMAKE_CURRENT_BINARY_DIR}/dual_e_self)
    set_tests_properties(unsplit_dual_e_self PROPERTIES TIMEOUT 600 LABELS "physics;threads")
    # restart: the Compton step and the floors are stateless
    add_test(NAME restart_smoke_thermal
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/thermal 13 --set thermal_init=recfast
                     --set hydro_scheme=unsplit)
    # xray_heating = 1 (docs/xray_heating.md): thread independence, the keys
    # alone change nothing, the neutral cells heated and ionised
    add_test(NAME xray_heating
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/xray_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/xray_heating
                     ${PROJECT_SOURCE_DIR}/data/uv_lf/am_params_st_fof.txt)
    set_tests_properties(xray_heating PROPERTIES TIMEOUT 600 LABELS "physics;threads")
    # ray_rotate = 3 (B16 fix): the per-substep rotation is keyed on the clock
    add_test(NAME restart_smoke_rotate
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/rotate 14 --set legacy_rt=0
                     --set legacy_physics=0 --set ray_rotate=3)
    # Ndm > N through the checkpoint (pm_refine = 2, docs/pm_refine.md)
    add_test(NAME restart_smoke_pm_refine
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly_r2
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/pm_refine 22 --set hydro_scheme=unsplit)
    set_tests_properties(restart_smoke_pm_refine PROPERTIES TIMEOUT 900 LABELS "restart;pm")
    # the FoF cache through the checkpoint (fof_dt_myr; checkpoint version 7):
    # restart at step 77, between FoF calls with cached halos
    add_test(NAME restart_smoke_fof_cadence
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly_r2
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/fof_cadence 77
                     --set "zdata1=12 7.5" --set ra_max=0.03 --set fof_dt_myr=60)
    set_tests_properties(restart_smoke_fof_cadence PROPERTIES TIMEOUT 900 LABELS "restart;pm")
    # the moment RT state through the checkpoint (docs/rt_moments.md)
    foreach(m m1 otvet)
      add_test(NAME restart_smoke_${m}
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                       $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                       ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/${m} 13 --set rt_method=${m})
      set_tests_properties(restart_smoke_${m} PROPERTIES TIMEOUT 600 LABELS "restart;moments")
    endforeach()
    # the recombination-radiation buffers through the checkpoint
    # (docs/recomb_radiation.md)
    add_test(NAME restart_smoke_recomb
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/recomb 14 --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1
                     --set recomb_case=A --set recomb_rad=1)
    set_tests_properties(restart_smoke_recomb PROPERTIES TIMEOUT 600 LABELS "restart;recomb")
    # a checkpoint written before the RT starts (empty recombination buffers)
    add_test(NAME restart_smoke_recomb_prert
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/recomb_prert 2 --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1
                     --set recomb_case=A --set recomb_rad=1)
    set_tests_properties(restart_smoke_recomb_prert PROPERTIES TIMEOUT 600 LABELS "restart;recomb")
    add_test(NAME restart_smoke_recomb_otvet
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/recomb_otvet 13 --set rt_method=otvet
                     --set recomb_case=A --set recomb_rad=1 --set recomb_rad_lines=1)
    set_tests_properties(restart_smoke_recomb_otvet PROPERTIES TIMEOUT 600 LABELS "restart;recomb;moments")
    # moment_inject = 1: sources and recombination photons injected after the
    # transport (docs/rt_moments.md section 5a)
    add_test(NAME restart_smoke_m1_inject
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/m1_inject 13 --set rt_method=m1
                     --set recomb_case=A --set recomb_rad=1 --set moment_inject=1)
    set_tests_properties(restart_smoke_m1_inject PROPERTIES TIMEOUT 600 LABELS "restart;recomb;moments")
    # the alternative OTVET tensors (otvet_closure.hpp): the split and the
    # cone's photon budget read the checkpointed N
    add_test(NAME restart_smoke_otvet_split
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/otvet_split 13 --set rt_method=otvet
                     --set otvet_split=1)
    add_test(NAME restart_smoke_otvet_cone
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/otvet_cone 13 --set rt_method=otvet
                     --set otvet_closure=cone --set otvet_cone_nside=2 --set recomb_case=A
                     --set recomb_rad=1)
    set_tests_properties(restart_smoke_otvet_split restart_smoke_otvet_cone PROPERTIES TIMEOUT 900
                         LABELS "restart;moments")
    # tensor_dump (docs/otvet_closure.md section 1) is a pure diagnostic: the
    # other outputs byte-identical with and without it, for rays (with
    # recombination rays) and OTVET; and independent of the thread count
    add_test(NAME tensor_dump_rays
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tensor_dump_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/tensor_dump/rays $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1 --set "zdata1=30 25 20"
                     --set recomb_case=A --set recomb_rad=1)
    add_test(NAME tensor_dump_otvet
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tensor_dump_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/tensor_dump/otvet $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1 --set "zdata1=30 25 20"
                     --set rt_method=otvet)
    # ... and the diagnostic's D N with the alternative tensors
    add_test(NAME tensor_dump_otvet_cone
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tensor_dump_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/tensor_dump/otvet_cone $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1 --set "zdata1=30 25 20"
                     --set rt_method=otvet --set otvet_closure=cone --set otvet_cone_nside=2)
    add_test(NAME tensor_dump_otvet_split
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tensor_dump_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/tensor_dump/otvet_split $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1 --set "zdata1=30 25 20"
                     --set rt_method=otvet --set otvet_closure=screen --set otvet_split=1)
    set_tests_properties(tensor_dump_rays tensor_dump_otvet tensor_dump_otvet_cone
                         tensor_dump_otvet_split PROPERTIES TIMEOUT 900 LABELS moments)
    add_test(NAME thread_full_tensor_dump
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                     ${CMAKE_CURRENT_BINARY_DIR}/thread_runs/tensor_dump $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1 --set "zdata1=30 25 20"
                     --set tensor_dump=1)
    set_tests_properties(thread_full_tensor_dump PROPERTIES TIMEOUT 1800 LABELS threads
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
    set_tests_properties(restart_smoke_legacy restart_smoke_det_movie restart_smoke_unsplit
                         restart_smoke_rotate PROPERTIES TIMEOUT 600 LABELS restart)
    # slab-only movie frames (movie_slab, docs/movies.md): the slabs against
    # the full cubes, and the movie tools and reion.py on them; and through
    # a restart (slab.txt)
    add_test(NAME movie_slab
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/movie_slab_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/movie_slab --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1 --set movie_dz=0.5
                     --set movie_z0=25)
    set_tests_properties(movie_slab PROPERTIES TIMEOUT 900 LABELS movie)
    # rays crossing 2 cells per substep (docs/rt_multicell.md): the
    # checkpoint records rt_cells_per_substep and the k-cell substep carry
    add_test(NAME restart_smoke_det_multicell
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/det_multicell 14 --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1 --set rt_substep_carry=1
                     --set rt_cells_per_substep=2)
    set_tests_properties(restart_smoke_det_multicell PROPERTIES TIMEOUT 600 LABELS "restart;multicell")
    add_test(NAME restart_smoke_det_movie_slab
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${CMAKE_CURRENT_BINARY_DIR}/restart_runs/det_slab 14 --set legacy_rt=0
                     --set legacy_physics=0 --set deterministic_parallel=1 --set movie_dz=0.5
                     --set movie_z0=25 --set movie_slab=4)
    set_tests_properties(restart_smoke_det_movie_slab PROPERTIES TIMEOUT 600 LABELS "restart;movie")
  endif()
endif()

# ---------------------------------------------------------------------------
# AMReX backend (rt_backend = amrex, docs/stage_d.md; label amrex, only with
# RADRAY_AMREX=ON). The sorted-ray-array RT must be bitwise identical to
# spike/serial rt_serial for any rank count, box size and thread count, and
# agree with the serial port (legacy_rt = 0) to round-off; restarts and the
# full loop are checked as for the serial code.
if(RADRAY_HAVE_AMREX AND Python3_Interpreter_FOUND)
  # the sorted-array reference of the spike, built against the library
  add_executable(rt_serial ${PROJECT_SOURCE_DIR}/spike/serial/rt_serial.cpp)
  target_include_directories(rt_serial PRIVATE ${PROJECT_SOURCE_DIR}/spike/common)
  target_link_libraries(rt_serial PRIVATE ${RADRAY_RT_BACKEND_LIB})
  target_compile_options(rt_serial PRIVATE -ffp-contract=off -fno-fast-math)
  # The MPI launcher of the tests: mpiexec from FindMPI (set RADRAY_MPIEXEC
  # and RADRAY_MPIEXEC_FLAGS for another one, e.g. srun on a cluster)
  find_package(MPI COMPONENTS CXX QUIET)
  set(RADRAY_MPIEXEC "${MPIEXEC_EXECUTABLE}" CACHE STRING "MPI launcher of the amrex tests")
  set(RADRAY_MPIEXEC_FLAGS "--oversubscribe;--allow-run-as-root;--bind-to;none" CACHE STRING
      "launcher flags (OpenMPI defaults; empty for srun)")
  set(RADRAY_MPIEXEC_NP "-np" CACHE STRING "launcher flag for the rank count")
  set(amrex_runs ${CMAKE_CURRENT_BINARY_DIR}/amrex_runs)
  # The preheating floor on the RT backend (bitwise against the serial RT
  # on the source-free cooling test) and a Planck full loop on the three
  # backends (docs/production_run.md section 10)
  add_test(NAME cosmology_floor_amrex
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cosmology_check.py floor
                   ${CMAKE_CURRENT_BINARY_DIR}/cosmology_floor_amrex $<TARGET_FILE:rt_driver>
                   --amrex ${RADRAY_MPIEXEC} "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} 2)
  add_test(NAME cosmology_planck_full_amrex
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/cosmology_check.py full
                   ${CMAKE_CURRENT_BINARY_DIR}/cosmology_full_amrex $<TARGET_FILE:radhydro>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke ${PROJECT_SOURCE_DIR}/tools/ics/make_ics.py
                   --set legacy_rt=0 --set legacy_physics=0 --set hydro_scheme=unsplit
                   --set rt_backend=amrex --set hydro_backend=amrex --set pm_backend=amrex)
  set_tests_properties(cosmology_floor_amrex cosmology_planck_full_amrex PROPERTIES
                       LABELS "amrex;cosmology" ENVIRONMENT OMP_NUM_THREADS=1 TIMEOUT 1800 PROCESSORS 2)
  # bitwise matrix: 3 cases x (ranks, threads, max_grid_size), the matrix of
  # spike/amrex/validate.sh (N = 32: box sizes 32, 16 and 8)
  # (ln16p0: the spike's rt_serial keeps the original photon limiter, so
  # the rt_driver side runs the fixed physics with photon_conserve = 0, the
  # pre-2026-10-14 meaning of legacy_physics = 0)
  foreach(case "u8;uniform_8src;-" "ln16;lognormal_16src_smallbuf;-"
               "ln16p0;lognormal_16src_smallbuf;legacy_physics=0;photon_conserve=0")
    list(GET case 0 tag)
    list(GET case 1 prob)
    list(GET case 2 xs)  # one extra key=value for both sides, or -
    list(LENGTH case ncase)
    set(extra)
    set(xset "")
    if(NOT xs STREQUAL "-")
      set(extra --set ${xs})
      set(xset ",${xs}")
    endif()
    if(ncase GREATER 3)  # a key=value for the rt_driver side only
      list(GET case 3 xs2)
      set(xset "${xset},${xs2}")
    endif()
    add_test(NAME amrex_ref_${tag}
             COMMAND rt_serial ${CMAKE_CURRENT_SOURCE_DIR}/problems/${prob} ${amrex_runs}/ref_${tag}
                     --set legacy_rt=0 ${extra})
    set_tests_properties(amrex_ref_${tag} PROPERTIES FIXTURES_SETUP amrex_ref_${tag} LABELS amrex
                         ENVIRONMENT OMP_NUM_THREADS=2 TIMEOUT 1800)
    foreach(cfg "1;1;32" "1;2;16" "2;1;16" "2;2;8" "4;1;16" "4;2;8" "8;1;8" "8;2;16")
      list(GET cfg 0 np)
      list(GET cfg 1 nt)
      list(GET cfg 2 mgs)
      set(tname amrex_bitwise_${tag}_r${np}t${nt}_mgs${mgs})
      add_test(NAME ${tname}
               COMMAND ${CMAKE_COMMAND}
                 -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
                 -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=${np} -DNT=${nt}
                 -DDRIVER=$<TARGET_FILE:rt_driver>
                 -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/${prob}
                 -DOUT=${amrex_runs}/${tname} -DREF=${amrex_runs}/ref_${tag}
                 "-DSET=legacy_rt=0,rt_backend=amrex,amrex_max_grid_size=${mgs}${xset}"
                 -DPYTHON=${Python3_EXECUTABLE} -DCOMPARE=${PROJECT_SOURCE_DIR}/tools/compare/compare.py
                 -DSTRICT=1
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/run_amrex.cmake)
      set_tests_properties(${tname} PROPERTIES FIXTURES_REQUIRED amrex_ref_${tag} LABELS amrex
                           ENVIRONMENT OMP_NUM_THREADS=${nt} TIMEOUT 1800 PROCESSORS ${np})
    endforeach()
    # amrex_slab = 1 (docs/stage_e_moments.md section 2): z-slabs of P planes
    # in z order, with more slabs than ranks (3 ranks, 8 slabs: unequal
    # shares), thinner than the ghost and Redistribute reach (P = 2), and
    # the automatic P (0)
    foreach(cfg "1;1;0" "2;1;8" "3;1;4" "4;2;2" "8;1;0")
      list(GET cfg 0 np)
      list(GET cfg 1 nt)
      list(GET cfg 2 planes)
      set(tname amrex_bitwise_${tag}_r${np}t${nt}_slab${planes})
      add_test(NAME ${tname}
               COMMAND ${CMAKE_COMMAND}
                 -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
                 -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=${np} -DNT=${nt}
                 -DDRIVER=$<TARGET_FILE:rt_driver>
                 -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/${prob}
                 -DOUT=${amrex_runs}/${tname} -DREF=${amrex_runs}/ref_${tag}
                 "-DSET=legacy_rt=0,rt_backend=amrex,amrex_slab=1,amrex_slab_planes=${planes}${xset}"
                 -DPYTHON=${Python3_EXECUTABLE} -DCOMPARE=${PROJECT_SOURCE_DIR}/tools/compare/compare.py
                 -DSTRICT=1
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/run_amrex.cmake)
      set_tests_properties(${tname} PROPERTIES FIXTURES_REQUIRED amrex_ref_${tag} LABELS "amrex;amrex_slab"
                           ENVIRONMENT OMP_NUM_THREADS=${nt} TIMEOUT 1800 PROCESSORS ${np})
    endforeach()
  endforeach()
  # the moment RT on the level (rt_method = m1 | otvet with rt_backend =
  # amrex, docs/stage_e_moments.md): bitwise identical to the serial moment
  # RT (rt_driver, the same settings) for ranks x threads x slab planes,
  # with every sum in the serial order (amrex_moment_stats = exact); with
  # amrex_moment_stats = fixed the fields stay bitwise and the rt_stat sums
  # agree to round-off. N = 32: 3 ranks hold 3, 3, 2 slabs of 4.
  # (case; problem; settings; configurations: all, or the short list "few")
  foreach(mcase "m1;lognormal_16src;rt_method=m1,nsteps=4;all"
                "m1;uniform_8src;rt_method=m1,nsteps=4;all"
                "m1inj;lognormal_16src;rt_method=m1,nsteps=4,moment_inject=1,moment_nsub=4;all"
                "otvet;lognormal_16src;rt_method=otvet,nsteps=4;all"
                "otvetfft;lognormal_16src;rt_method=otvet,nsteps=4,otvet_tensor=fft;few"
                "split;lognormal_16src;rt_method=otvet,nsteps=4,otvet_split=1;few"
                "splitfft;lognormal_16src;rt_method=otvet,nsteps=4,otvet_split=1,otvet_tensor=fft;few"
                "screenfix;lognormal_16src;rt_method=otvet,nsteps=4,otvet_closure=screen,otvet_screen_lambda=6;few"
                "screensplit;uniform_8src;rt_method=otvet,otvet_closure=screen,otvet_split=1;few"
                "cone2;lognormal_16src;rt_method=otvet,nsteps=4,otvet_closure=cone,otvet_cone_nside=2;cone"
                "cone;uniform_8src;rt_method=otvet,otvet_closure=cone;cone"
                "conenoblend;lognormal_16src;rt_method=otvet,nsteps=4,otvet_closure=cone,otvet_cone_blend=0;cone"
                "conelamx;lognormal_16src;rt_method=otvet,otvet_closure=cone,otvet_cone_rc=1,otvet_cone_coarse=2;cone2"
                "m1recomb;lognormal_16src;rt_method=m1,nsteps=4,recomb_case=A,recomb_rad=1,recomb_rad_lines=1;few"
                "m1recombinj;lognormal_16src;rt_method=m1,nsteps=4,moment_inject=1,recomb_case=A,recomb_rad=1;few"
                "conerecomb;uniform_8src;rt_method=otvet,otvet_closure=cone,recomb_case=A,recomb_rad=1;cone")
    list(GET mcase 0 mtag)
    list(GET mcase 1 prob)
    list(GET mcase 2 xs)
    list(GET mcase 3 which)
    set(mset "legacy_rt=0,legacy_physics=0,${xs}")
    string(REPLACE "," ";" mset_list "${mset}")
    set(mref_args)
    foreach(kv IN LISTS mset_list)
      list(APPEND mref_args --set ${kv})
    endforeach()
    set(mref amrex_moment_ref_${mtag}_${prob})
    add_test(NAME ${mref}
             COMMAND rt_driver ${CMAKE_CURRENT_SOURCE_DIR}/problems/${prob} ${amrex_runs}/${mref} ${mref_args})
    set_tests_properties(${mref} PROPERTIES FIXTURES_SETUP ${mref} LABELS "amrex;moments"
                         ENVIRONMENT OMP_NUM_THREADS=2 TIMEOUT 1800)
    if(which STREQUAL "all")
      set(mcfgs "1:1:0:exact" "2:1:8:exact" "3:2:4:exact" "4:1:2:exact" "4:2:0:exact" "3:1:4:fixed")
    elseif(which STREQUAL "cone")
      # (slab planes a multiple of otvet_cone_coarse = 4)
      set(mcfgs "1:1:0:exact" "3:2:4:exact" "4:1:4:exact" "2:2:8:fixed")
    elseif(which STREQUAL "cone2")
      # (otvet_cone_coarse = 2, r_c = 1: lambda_x from the budget, the fine
      # mip level gathered)
      set(mcfgs "1:1:0:exact" "3:2:2:exact" "4:1:2:exact" "2:2:4:fixed")
    else()
      set(mcfgs "1:1:0:exact" "3:2:4:exact" "4:1:2:exact" "2:2:8:fixed")
    endif()
    foreach(cfg IN LISTS mcfgs)
      string(REPLACE ":" ";" cfg "${cfg}")
      list(GET cfg 0 np)
      list(GET cfg 1 nt)
      list(GET cfg 2 planes)
      list(GET cfg 3 stats)
      set(tname amrex_moment_${mtag}_${prob}_r${np}t${nt}_p${planes}_${stats})
      set(xtol "")
      if(stats STREQUAL "fixed")
        set(xtol "--stats-rtol;1e-12")
      endif()
      add_test(NAME ${tname}
               COMMAND ${CMAKE_COMMAND}
                 -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
                 -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=${np} -DNT=${nt}
                 -DDRIVER=$<TARGET_FILE:rt_driver>
                 -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/${prob}
                 -DOUT=${amrex_runs}/${tname} -DREF=${amrex_runs}/${mref}
                 "-DSET=${mset},rt_backend=amrex,amrex_slab=1,amrex_slab_planes=${planes},amrex_moment_stats=${stats}"
                 -DPYTHON=${Python3_EXECUTABLE} -DCOMPARE=${PROJECT_SOURCE_DIR}/tools/compare/compare.py
                 -DSTRICT=1 "-DEXTRA_TOL=${xtol}"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/run_amrex.cmake)
      set_tests_properties(${tname} PROPERTIES FIXTURES_REQUIRED ${mref} LABELS "amrex;moments;amrex_moments"
                           ENVIRONMENT OMP_NUM_THREADS=${nt} TIMEOUT 1800 PROCESSORS ${np})
    endforeach()
  endforeach()
  # round-off against the serial port (rt_driver, legacy_rt = 0) on the five
  # RT problems: tolerances of docs/framework_decision.md section 2 (dx in
  # source cells is the sum of 192 symmetric directions: float noise vs an
  # exact integer sum; fgam in a few I-front cells)
  foreach(p uniform_1src uniform_8src lognormal_16src uniform_8src_smallbuf lognormal_16src_smallbuf)
    add_test(NAME amrex_serial_ref_${p}
             COMMAND rt_driver ${CMAKE_CURRENT_SOURCE_DIR}/problems/${p} ${amrex_runs}/serial_${p}
                     --set legacy_rt=0)
    set_tests_properties(amrex_serial_ref_${p} PROPERTIES FIXTURES_SETUP amrex_serial_${p} LABELS amrex
                         TIMEOUT 1800)
    add_test(NAME amrex_roundoff_${p}
             COMMAND ${CMAKE_COMMAND}
               -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
               -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=2 -DNT=2
               -DDRIVER=$<TARGET_FILE:rt_driver>
               -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/${p}
               -DOUT=${amrex_runs}/roundoff_${p} -DREF=${amrex_runs}/serial_${p}
               "-DSET=legacy_rt=0,rt_backend=amrex"
               -DPYTHON=${Python3_EXECUTABLE} -DCOMPARE=${PROJECT_SOURCE_DIR}/tools/compare/compare.py
               -DSTRICT=0
               -P ${CMAKE_CURRENT_SOURCE_DIR}/run_amrex.cmake)
    set_tests_properties(amrex_roundoff_${p} PROPERTIES FIXTURES_REQUIRED amrex_serial_${p} LABELS amrex
                         ENVIRONMENT OMP_NUM_THREADS=2 TIMEOUT 1800 PROCESSORS 2)
  endforeach()
  # rays crossing k cells per substep (docs/rt_multicell.md): the ghost
  # layers and the local Redistribute grow to k+1; round-off against the
  # serial port with the same k
  foreach(k 2 4 4x)
    set(xs "")
    set(xsets)
    set(xtol "")
    if(k STREQUAL "4x")
      # k = 4 with the exposure deposit and the converged photon limit. The
      # exposure deposit passes photons attenuated by exp(-tau) into thick
      # gas; deposits of ~1e-12 photons per H atom are at the resolution of
      # the backend's int64 fixed point (2^-44), so the faint heating of
      # cold neutral cells near the temperature floor (3-30 K) differs by a
      # few per cent (ray statistics exact, x_HII to 1e-7)
      set(xs ",ray_deposit=exposure,photon_conserve=1")
      set(xsets --set ray_deposit=exposure --set photon_conserve=1)
      set(xtol "--field-rtol;temp=0.05;--field-rtol;dx=2")
    endif()
    string(REPLACE "x" "" kk "${k}")
    foreach(p uniform_8src lognormal_16src)
      add_test(NAME amrex_serial_ref_multicell${k}_${p}
               COMMAND rt_driver ${CMAKE_CURRENT_SOURCE_DIR}/problems/${p}
                       ${amrex_runs}/serial_multicell${k}_${p} --set legacy_rt=0
                       --set rt_cells_per_substep=${kk} ${xsets})
      set_tests_properties(amrex_serial_ref_multicell${k}_${p} PROPERTIES
                           FIXTURES_SETUP amrex_serial_multicell${k}_${p} LABELS "amrex;multicell"
                           TIMEOUT 1800)
      add_test(NAME amrex_roundoff_multicell${k}_${p}
               COMMAND ${CMAKE_COMMAND}
                 -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
                 -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=2 -DNT=2
                 -DDRIVER=$<TARGET_FILE:rt_driver>
                 -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/${p}
                 -DOUT=${amrex_runs}/roundoff_multicell${k}_${p}
                 -DREF=${amrex_runs}/serial_multicell${k}_${p}
                 "-DSET=legacy_rt=0,rt_backend=amrex,rt_cells_per_substep=${kk},amrex_max_grid_size=8${xs}"
                 -DPYTHON=${Python3_EXECUTABLE} -DCOMPARE=${PROJECT_SOURCE_DIR}/tools/compare/compare.py
                 -DSTRICT=0 "-DEXTRA_TOL=${xtol}"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/run_amrex.cmake)
      set_tests_properties(amrex_roundoff_multicell${k}_${p} PROPERTIES
                           FIXTURES_REQUIRED amrex_serial_multicell${k}_${p} LABELS "amrex;multicell"
                           ENVIRONMENT OMP_NUM_THREADS=2 TIMEOUT 1800 PROCESSORS 2)
    endforeach()
  endforeach()
  # thread-count independence of the backend (1 rank)
  add_test(NAME amrex_thread_lognormal_16src
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                   ${amrex_runs}/thread_lognormal_16src $<TARGET_FILE:rt_driver>
                   ${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src --set legacy_rt=0
                   --set legacy_physics=0 --set rt_backend=amrex --set nsteps=6)
  set_tests_properties(amrex_thread_lognormal_16src PROPERTIES TIMEOUT 1800 LABELS "amrex;threads"
                       ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
  # the R-struct encoding runs and stays close to the float layout
  add_test(NAME amrex_rstruct_lognormal_16src
           COMMAND ${CMAKE_COMMAND}
             -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
             -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=2 -DNT=1
             -DDRIVER=$<TARGET_FILE:rt_driver>
             -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems/lognormal_16src
             -DOUT=${amrex_runs}/rstruct_lognormal_16src -DREF=${amrex_runs}/serial_lognormal_16src
             "-DSET=legacy_rt=0,rt_backend=amrex,ray_encoding=rstruct"
             -DPYTHON=${Python3_EXECUTABLE} -DCOMPARE=${PROJECT_SOURCE_DIR}/tools/compare/compare.py
             -DSTRICT=lossy
             -P ${CMAKE_CURRENT_SOURCE_DIR}/run_amrex.cmake)
  set_tests_properties(amrex_rstruct_lognormal_16src PROPERTIES FIXTURES_REQUIRED amrex_serial_lognormal_16src
                       LABELS amrex ENVIRONMENT OMP_NUM_THREADS=1 TIMEOUT 1800 PROCESSORS 2)
  if(RADRAY_HAVE_PM)
    # the full loop: 1 and 2 ranks byte-identical (fields, movie frames,
    # checkpoints), and a bit-identical restart with the backend
    add_test(NAME amrex_full_ranks_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rank_check.py
                     ${amrex_runs}/full_ranks "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP}
                     "1 2" $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet
                     --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1 --set rt_backend=amrex
                     --set movie_dz=0.5 --set movie_z0=25 --set restart_every=1 --set restart_keep_latest=0)
    set_tests_properties(amrex_full_ranks_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex;full"
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
    add_test(NAME amrex_restart_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${amrex_runs}/restart 13 --set legacy_rt=0 --set legacy_physics=0
                     --set rt_backend=amrex --set movie_dz=0.5 --set movie_z0=25)
    set_tests_properties(amrex_restart_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex;restart")
    # the unsplit hydro scheme's dual-energy field through the backend
    add_test(NAME amrex_restart_smoke_unsplit
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${amrex_runs}/restart_unsplit 13 --set legacy_rt=0 --set legacy_physics=0
                     --set rt_backend=amrex --set hydro_scheme=unsplit)
    set_tests_properties(amrex_restart_smoke_unsplit PROPERTIES TIMEOUT 1800 LABELS "amrex;restart")
    # Stage E1 (docs/stage_e_hydro.md): hydro_backend = amrex, the unsplit
    # hydro step on AMReX MultiFabs, must be BITWISE identical to the serial
    # unsplit step for any rank count, box size and thread count (label
    # amrex_hydro, inside amrex)
    foreach(p sod_x random3d sound_x)
      add_test(NAME amrex_hydro_ref_${p}
               COMMAND hydro_driver ${CMAKE_CURRENT_SOURCE_DIR}/problems_hydro/${p} ${amrex_runs}/hydro_ref_${p}
                       --quiet --set hydro_scheme=unsplit)
      set_tests_properties(amrex_hydro_ref_${p} PROPERTIES FIXTURES_SETUP amrex_hydro_ref_${p}
                           LABELS "amrex;amrex_hydro" ENVIRONMENT OMP_NUM_THREADS=2 TIMEOUT 600)
      foreach(cfg "1;1;32" "1;2;16" "2;1;16" "2;2;8" "4;1;16" "4;2;8" "8;1;8")
        list(GET cfg 0 np)
        list(GET cfg 1 nt)
        list(GET cfg 2 mgs)
        set(tname amrex_hydro_bitwise_${p}_r${np}t${nt}_mgs${mgs})
        add_test(NAME ${tname}
                 COMMAND ${CMAKE_COMMAND}
                   -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
                   -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=${np} -DNT=${nt}
                   -DDRIVER=$<TARGET_FILE:hydro_driver>
                   -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems_hydro/${p}
                   -DOUT=${amrex_runs}/${tname} -DREF=${amrex_runs}/hydro_ref_${p}
                   "-DSET=hydro_scheme=unsplit,hydro_backend=amrex,amrex_hydro_max_grid_size=${mgs}"
                   -DEXTRA=--quiet
                   -DPYTHON=${Python3_EXECUTABLE} -DCOMPARE=${PROJECT_SOURCE_DIR}/tools/compare/compare.py
                   -DSTRICT=1
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/run_amrex.cmake)
        set_tests_properties(${tname} PROPERTIES FIXTURES_REQUIRED amrex_hydro_ref_${p}
                             LABELS "amrex;amrex_hydro" ENVIRONMENT OMP_NUM_THREADS=${nt} TIMEOUT 600
                             PROCESSORS ${np})
      endforeach()
    endforeach()
    # Stage F1 (docs/stage_f.md): the refined hierarchy. The result must not
    # depend on the rank count, the thread count or the box sizes (the fine
    # grids are chopped per rank), and the conserved totals of every step
    # must equal the uniform run's to round-off (periodic box, refluxing
    # and average-down). Label amr inside amrex_hydro.
    # random3d (density contrasts of 1e6) hits the dual-energy resets and the
    # temperature floor millions of times, which are not conservative in
    # energy by design: its energy is checked to 5e-3 only
    foreach(cfg "sod_x;1;0.5;" "sound_x;1;0.99;" "random3d;2;8;--energy-rtol;5e-3")
      list(GET cfg 0 p)
      list(GET cfg 1 L)
      list(GET cfg 2 M)
      list(SUBLIST cfg 3 -1 etol)  # the rest: extra check arguments, if any
      set(xtol "")
      if(etol)
        set(xtol ";${etol}")
      endif()
      set(amrset "hydro_scheme=unsplit,hydro_backend=amrex,amrex_max_level=${L},amrex_refine_mass=${M}")
      math(EXPR Lp1 "${L} + 1")
      add_test(NAME amrex_hydro_amr_ref_${p}_L${L}
               COMMAND ${CMAKE_COMMAND}
                 -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
                 -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=1 -DNT=2
                 -DDRIVER=$<TARGET_FILE:hydro_driver>
                 -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems_hydro/${p}
                 -DOUT=${amrex_runs}/amr_ref_${p}_L${L} -DREF=${amrex_runs}/hydro_ref_${p}
                 "-DSET=${amrset},amrex_hydro_max_grid_size=16"
                 -DEXTRA=--quiet
                 -DPYTHON=${Python3_EXECUTABLE} -DCOMPARE=${CMAKE_CURRENT_SOURCE_DIR}/amr_conservation_check.py
                 -DSTRICT=cons "-DEXTRA_TOL=--levels;${Lp1}${xtol}"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/run_amrex.cmake)
      set_tests_properties(amrex_hydro_amr_ref_${p}_L${L} PROPERTIES FIXTURES_REQUIRED amrex_hydro_ref_${p}
                           FIXTURES_SETUP amrex_hydro_amr_ref_${p}_L${L} LABELS "amrex;amrex_hydro;amr"
                           ENVIRONMENT OMP_NUM_THREADS=2 TIMEOUT 900)
      foreach(rc "2;1;16" "4;2;8" "1;1;32")
        list(GET rc 0 np)
        list(GET rc 1 nt)
        list(GET rc 2 mgs)
        set(tname amrex_hydro_amr_bitwise_${p}_L${L}_r${np}t${nt}_mgs${mgs})
        add_test(NAME ${tname}
                 COMMAND ${CMAKE_COMMAND}
                   -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
                   -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=${np} -DNT=${nt}
                   -DDRIVER=$<TARGET_FILE:hydro_driver>
                   -DPROBLEM=${CMAKE_CURRENT_SOURCE_DIR}/problems_hydro/${p}
                   -DOUT=${amrex_runs}/${tname} -DREF=${amrex_runs}/amr_ref_${p}_L${L}
                   "-DSET=${amrset},amrex_hydro_max_grid_size=${mgs}"
                   -DEXTRA=--quiet
                   -DPYTHON=${Python3_EXECUTABLE} -DCOMPARE=${PROJECT_SOURCE_DIR}/tools/compare/compare.py
                   -DSTRICT=1
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/run_amrex.cmake)
        set_tests_properties(${tname} PROPERTIES FIXTURES_REQUIRED amrex_hydro_amr_ref_${p}_L${L}
                             LABELS "amrex;amrex_hydro;amr" ENVIRONMENT OMP_NUM_THREADS=${nt} TIMEOUT 900
                             PROCESSORS ${np})
      endforeach()
    endforeach()
    # The star-formation sink (docs/stage_f.md section 2, label sink): the
    # smoke problem with the state on the levels and one refinement level,
    # stars above 2 x the mean density (unphysically low: every step forms
    # some). 1, 2 and 3 ranks byte-identical (the stars file is sorted by
    # formation time and position), a bit-identical restart (checkpoint
    # version 8 with the stars and the .amr levels), and the sink's
    # bookkeeping (sink_check.py).
    set(sink_sets --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1 --set rt_backend=amrex
                  --set hydro_scheme=unsplit --set hydro_backend=amrex --set amrex_hydro_max_grid_size=0
                  --set pm_backend=amrex --set amrex_state=1 --set amrex_max_level=1 --set amrex_refine_mass=4
                  --set sf_overdensity=2 --set sf_efficiency=1 --set movie_dz=0.5 --set movie_z0=25)
    add_test(NAME amrex_sink_ranks_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rank_check.py
                     ${amrex_runs}/sink_ranks "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP}
                     "1 2 3" $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet
                     ${sink_sets} --set restart_every=1 --set restart_keep_latest=0)
    set_tests_properties(amrex_sink_ranks_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex;amrex_hydro;sink;full"
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
    add_test(NAME amrex_sink_check_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/sink_check.py ${amrex_runs}/sink_ranks/np1
                     --min-stars 1)
    set_tests_properties(amrex_sink_check_smoke PROPERTIES DEPENDS amrex_sink_ranks_smoke LABELS "amrex;amrex_hydro;sink")
    add_test(NAME amrex_sink_restart_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${amrex_runs}/sink_restart 13 ${sink_sets})
    set_tests_properties(amrex_sink_restart_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex;amrex_hydro;sink;restart"
                         ENVIRONMENT OMP_NUM_THREADS=1)
    # The PM on the refined gas levels (pm_gas_levels, docs/stage_f.md
    # section 4, label amr): a 16^3 run refined everywhere to level 1 with
    # the 32^3 mesh against the uniform 32^3 run with the same particles
    # (gas, momenta and particles; 10x better than the level-0 coupling),
    # the DM + gas momentum with two levels and pm_refine = 1, 2, 4
    # (pm_gas_levels_check.py, ~30 s on 2 ranks); and with two levels and
    # pm_refine = 4, 1, 2 and 3 ranks byte-identical and a bit-identical
    # restart.
    if(Python3_Interpreter_FOUND)
      add_test(NAME amrex_pm_gas_levels_check
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_gas_levels_check.py
                       ${amrex_runs}/pm_gas_levels "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}"
                       ${RADRAY_MPIEXEC_NP} $<TARGET_FILE:radhydro> ${PROJECT_SOURCE_DIR})
      set_tests_properties(amrex_pm_gas_levels_check PROPERTIES TIMEOUT 1800 LABELS "amrex;amrex_pm;amr"
                           ENVIRONMENT OMP_NUM_THREADS=1 PROCESSORS 2)
      set(gaslev_sets --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1 --set rt_backend=amrex
                      --set hydro_scheme=unsplit --set hydro_backend=amrex --set amrex_hydro_max_grid_size=0
                      --set pm_backend=amrex --set amrex_state=1 --set amrex_max_grid_size=8
                      --set amrex_max_level=2 --set amrex_refine_mass=0.05 --set pm_refine=4
                      --set movie_dz=0.5 --set movie_z0=25)
      add_test(NAME amrex_pm_gas_levels_ranks_smoke
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rank_check.py
                       ${amrex_runs}/pm_gas_levels_ranks "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}"
                       ${RADRAY_MPIEXEC_NP} "1 2 3" $<TARGET_FILE:radhydro>
                       ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet ${gaslev_sets}
                       --set restart_every=1 --set restart_keep_latest=0)
      set_tests_properties(amrex_pm_gas_levels_ranks_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex;amrex_pm;amr;full"
                           ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
      add_test(NAME amrex_pm_gas_levels_restart_smoke
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                       $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                       ${amrex_runs}/pm_gas_levels_restart 13 ${gaslev_sets})
      set_tests_properties(amrex_pm_gas_levels_restart_smoke PROPERTIES TIMEOUT 1800
                           LABELS "amrex;amrex_pm;amr;restart" ENVIRONMENT OMP_NUM_THREADS=1)
    endif()
    # Stage G (pm_level_solve, docs/stage_g.md, label level_solve): a 16^3 run
    # refined everywhere to level 1 with the level solve against the uniform
    # 32^3 run with the same particles (and against pm_gas_levels), the DM +
    # gas momentum with two levels over the whole box (periodic level
    # solves) and over half of it (pm_level_solve_check.py, ~40 s on 2
    # ranks); two levels over half the box with pm_refine = 1 and 2: 1, 2
    # and 3 ranks byte-identical and a bit-identical restart (to z = 28: the
    # gravity steps before the smoke problem's sources appear; with the RT
    # on, the checkpoints of a rank-dependent fine-level layout order their
    # rays differently and a restart through an RT step is not bit-identical,
    # with pm_gas_levels as well: docs/stage_g.md section 3).
    if(Python3_Interpreter_FOUND)
      add_test(NAME amrex_pm_level_solve_check
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_level_solve_check.py
                       ${amrex_runs}/pm_level_solve "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}"
                       ${RADRAY_MPIEXEC_NP} $<TARGET_FILE:radhydro> ${PROJECT_SOURCE_DIR})
      set_tests_properties(amrex_pm_level_solve_check PROPERTIES TIMEOUT 1800 LABELS "amrex;amrex_pm;amr;level_solve"
                           ENVIRONMENT OMP_NUM_THREADS=1 PROCESSORS 2)
      foreach(r 1 2)
        set(levsolve_sets --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1 --set rt_backend=amrex
                          --set hydro_scheme=unsplit --set hydro_backend=amrex --set amrex_hydro_max_grid_size=0
                          --set pm_backend=amrex --set amrex_state=1 --set amrex_max_grid_size=8
                          --set amrex_max_level=2 --set amrex_refine_mass=1e9 --set amrex_refine_box=0,0,0,7,7,7
                          --set pm_level_solve=1 --set pm_refine=${r} --set movie_dz=0.5 --set movie_z0=25
                          --set "zdata1=30 28")  # the run ends at z = 27.5 (step 12), before the smoke problem's sources appear at z = 24.7 (step 13)
                                                 # (the RT costs 100 s a step here; the restart must end where the run does)
        add_test(NAME amrex_pm_level_solve_ranks_smoke_r${r}
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rank_check.py
                         ${amrex_runs}/pm_level_solve_ranks_r${r} "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}"
                         ${RADRAY_MPIEXEC_NP} "1 2 3" $<TARGET_FILE:radhydro>
                         ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet ${levsolve_sets}
                         --set restart_every=1 --set restart_keep_latest=0)
        set_tests_properties(amrex_pm_level_solve_ranks_smoke_r${r} PROPERTIES TIMEOUT 1800
                             LABELS "amrex;amrex_pm;amr;level_solve;full" ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
        add_test(NAME amrex_pm_level_solve_restart_smoke_r${r}
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                         $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                         ${amrex_runs}/pm_level_solve_restart_r${r} 9 ${levsolve_sets})
        set_tests_properties(amrex_pm_level_solve_restart_smoke_r${r} PROPERTIES TIMEOUT 1800
                             LABELS "amrex;amrex_pm;amr;level_solve;restart" ENVIRONMENT OMP_NUM_THREADS=1)
      endforeach()
      # T1 of the zoom ladder (docs/zoom_turbulence.md section 5, docs/stage_g.md
      # section 4; nightly, ~15 min on 2 ranks): 32^3 base with two levels over
      # the central eighth and the level solve, and the same hierarchy with the
      # uniform 128^3 mesh (pm_gas_levels), against the uniform 128^3 run with
      # the same 128^3 particles: the particles, the FoF halos and the density
      # PDF inside the region
      add_test(NAME amrex_pm_zoom_t1
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_zoom_t1_check.py
                       ${amrex_runs}/pm_zoom_t1 "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}"
                       ${RADRAY_MPIEXEC_NP} 2 $<TARGET_FILE:radhydro> $<TARGET_FILE:find_halos> ${PROJECT_SOURCE_DIR})
      set_tests_properties(amrex_pm_zoom_t1 PROPERTIES TIMEOUT 7200 LABELS "amrex;amrex_pm;amr;level_solve;nightly"
                           ENVIRONMENT OMP_NUM_THREADS=2 PROCESSORS 4 DISABLED $<NOT:$<BOOL:${RADRAY_NIGHTLY_TESTS}>>)
      # without refinement the level solve is the Stage E level: byte-identical
      add_test(NAME amrex_pm_level_solve_identity_smoke
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/flag_identity_check.py
                       ${amrex_runs}/pm_level_solve_identity "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}"
                       ${RADRAY_MPIEXEC_NP} 2 $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                       pm_level_solve --quiet --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1
                       --set rt_backend=amrex --set hydro_scheme=unsplit --set hydro_backend=amrex
                       --set amrex_hydro_max_grid_size=0 --set pm_backend=amrex --set amrex_state=1
                       --set amrex_max_grid_size=8 --set run_steps=6)
      set_tests_properties(amrex_pm_level_solve_identity_smoke PROPERTIES TIMEOUT 1800
                           LABELS "amrex;amrex_pm;level_solve" ENVIRONMENT OMP_NUM_THREADS=1 PROCESSORS 2)
    endif()
    # Stage F3 (docs/stage_f.md section 3): the I-front of one source in
    # uniform gas with a refined level over the central region, the front
    # crossing the level boundary, against the uniform runs at the base and
    # at twice the resolution (ifront_full.py --check). Nightly: ~2 minutes
    # on 2 ranks.
    if(Python3_Interpreter_FOUND)
      add_test(NAME amrex_rt_amr_ifront
               COMMAND ${CMAKE_COMMAND}
                 -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
                 -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=2 -DNT=2
                 -DDRIVER=$<TARGET_FILE:radhydro> -DPYTHON=${Python3_EXECUTABLE}
                 -DSRC=${PROJECT_SOURCE_DIR} -DOUT=${amrex_runs}/rt_amr_ifront -DN=32 -DHALF=4
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/rt_amr_ifront.cmake)
      set_tests_properties(amrex_rt_amr_ifront PROPERTIES TIMEOUT 3600 LABELS "amrex;amr;nightly"
                           ENVIRONMENT OMP_NUM_THREADS=2 PROCESSORS 2)
      if(NOT RADRAY_NIGHTLY_TESTS)
        set_tests_properties(amrex_rt_amr_ifront PROPERTIES DISABLED TRUE)
      endif()
      # the 8 sources of uniform_8src, one level over a box where fronts
      # cross the level boundary and overlap there (rt_amr_ifront8.slurm's
      # BOX1), against uniform 32^3 and 64^3 (ifront_multi.py --check)
      add_test(NAME amrex_rt_amr_ifront8
               COMMAND ${CMAKE_COMMAND}
                 -DLAUNCH=${RADRAY_MPIEXEC} "-DLAUNCH_FLAGS=${RADRAY_MPIEXEC_FLAGS}"
                 -DNPFLAG=${RADRAY_MPIEXEC_NP} -DNP=2 -DNT=2
                 -DDRIVER=$<TARGET_FILE:radhydro> -DPYTHON=${Python3_EXECUTABLE}
                 -DSRC=${PROJECT_SOURCE_DIR} -DOUT=${amrex_runs}/rt_amr_ifront8 -DN=32
                 -DSOURCES=8src -DBOX=0,4,4,15,19,19 -DLEVEL=1
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/rt_amr_ifront.cmake)
      set_tests_properties(amrex_rt_amr_ifront8 PROPERTIES TIMEOUT 3600 LABELS "amrex;amr;nightly"
                           ENVIRONMENT OMP_NUM_THREADS=2 PROCESSORS 2)
      if(NOT RADRAY_NIGHTLY_TESTS)
        set_tests_properties(amrex_rt_amr_ifront8 PROPERTIES DISABLED TRUE)
      endif()
      # fields_levels = 1: the plane at the finest resolution of two nested
      # refined levels on 2 ranks (the level map, the level-0 pixels equal to
      # the base field)
      add_test(NAME amrex_fields_levels
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/levels_slice_check.py
                       ${amrex_runs}/fields_levels "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}"
                       ${RADRAY_MPIEXEC_NP} 2 $<TARGET_FILE:radhydro> ${PROJECT_SOURCE_DIR})
      set_tests_properties(amrex_fields_levels PROPERTIES TIMEOUT 900 LABELS "amrex;amr"
                           ENVIRONMENT OMP_NUM_THREADS=2 PROCESSORS 4)
      # two refined levels: byte-identical for 1 rank (two fine box sizes)
      # and 2 and 3 ranks (the fine grids are chopped by the rank count)
      add_test(NAME amrex_rt_amr_layout
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/amr_layout_check.py
                       ${amrex_runs}/amr_layout "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}"
                       ${RADRAY_MPIEXEC_NP} $<TARGET_FILE:radhydro> ${PROJECT_SOURCE_DIR})
      set_tests_properties(amrex_rt_amr_layout PROPERTIES TIMEOUT 2400 LABELS "amrex;amr;nightly"
                           ENVIRONMENT OMP_NUM_THREADS=2 PROCESSORS 6)
      if(NOT RADRAY_NIGHTLY_TESTS)  # ~18 min on klone
        set_tests_properties(amrex_rt_amr_layout PROPERTIES DISABLED TRUE)
      endif()
    endif()
    # the full loop with both backends: 1, 2 and 3 ranks byte-identical to
    # each other (rank_check) and to the run with the serial hydro
    add_test(NAME amrex_hydro_full_ranks_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rank_check.py
                     ${amrex_runs}/hydro_full_ranks "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP}
                     "1 2 3" $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet
                     --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1 --set rt_backend=amrex
                     --set hydro_scheme=unsplit --set hydro_backend=amrex --set amrex_hydro_max_grid_size=8
                     --set movie_dz=0.5 --set movie_z0=25 --set restart_every=1 --set restart_keep_latest=0)
    set_tests_properties(amrex_hydro_full_ranks_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex;amrex_hydro;full"
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
    add_test(NAME amrex_hydro_full_vs_serial_hydro
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/backend_pair_check.py
                     ${amrex_runs}/hydro_full_pair $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     "hydro_backend=serial" "hydro_backend=amrex,amrex_hydro_max_grid_size=8"
                     --quiet --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1
                     --set hydro_scheme=unsplit --set movie_dz=0.5 --set movie_z0=25 --set restart_every=1
                     --set restart_keep_latest=0)
    set_tests_properties(amrex_hydro_full_vs_serial_hydro PROPERTIES TIMEOUT 1800 LABELS "amrex;amrex_hydro;full"
                         ENVIRONMENT OMP_NUM_THREADS=2)
    # dual_e_self = 1 (docs/hydro_unsplit.md section 7, U7) on the cold
    # unheated IGM, where many cells evolve e from its own value: the AMReX
    # level byte-identical to the serial step (8^3 boxes, so the stencils
    # cross box boundaries)
    add_test(NAME amrex_hydro_dual_e_self_vs_serial
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/backend_pair_check.py
                     ${amrex_runs}/hydro_dual_e_self_pair $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly
                     "hydro_backend=serial" "hydro_backend=amrex,amrex_hydro_max_grid_size=8"
                     --quiet --set zdata1=20 --set hydro_scheme=unsplit --set dual_e_self=1
                     --set thermal_floor=0.01 --set fields_every=1 --set fields=u_gas,e_gas)
    set_tests_properties(amrex_hydro_dual_e_self_vs_serial PROPERTIES TIMEOUT 1800 LABELS "amrex;amrex_hydro;full"
                         ENVIRONMENT OMP_NUM_THREADS=2)
    # amrex_state = 1 (docs/stage_e_state.md): the gas on the levels, byte-
    # identical to the three backends with the host state on every rank
    # (fields, stats, movie frames, checkpoints), for 1/2/3 ranks x 1/2
    # threads, and a bit-identical restart (label amrex_state, inside amrex)
    foreach(rt "1 1" "2 1" "3 1" "2 2")
      separate_arguments(rt)
      list(GET rt 0 np)
      list(GET rt 1 nt)
      add_test(NAME amrex_state_full_pair_r${np}t${nt}
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/backend_pair_check.py
                       --launch "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} ${np}
                       ${amrex_runs}/state_pair_r${np}t${nt} $<TARGET_FILE:radhydro>
                       ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke "amrex_state=0" "amrex_state=1"
                       --quiet --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1
                       --set rt_backend=amrex --set hydro_scheme=unsplit --set hydro_backend=amrex
                       --set amrex_hydro_max_grid_size=0 --set pm_backend=amrex --set movie_dz=0.5
                       --set movie_z0=25 --set restart_every=1 --set restart_keep_latest=0)
      set_tests_properties(amrex_state_full_pair_r${np}t${nt} PROPERTIES TIMEOUT 1800
                           LABELS "amrex;amrex_state;full" ENVIRONMENT OMP_NUM_THREADS=${nt}
                           PROCESSORS ${np})
    endforeach()
    # Plane-parallel RT-domain mode on the backend (docs/plane_parallel.md):
    # the RT-only physics checks with 2 ranks (the face rays across the
    # boxes, the far-face deletion and the fixed-point photon bookkeeping),
    # and the full loop on the three backends against the serial run
    add_test(NAME amrex_plane_rt
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/plane_check.py
                     $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                     ${amrex_runs}/plane_rt --only gamma slab front conserve domains
                     --set rt_backend=amrex --set amrex_max_grid_size=8
                     --launch "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} 2)
    set_tests_properties(amrex_plane_rt PROPERTIES TIMEOUT 1800 LABELS "amrex;plane;physics"
                         ENVIRONMENT OMP_NUM_THREADS=1 PROCESSORS 2)
    add_test(NAME amrex_plane_rt_sweep
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/plane_check.py
                     $<TARGET_FILE:rt_driver> ${CMAKE_CURRENT_SOURCE_DIR}/problems/uniform_1src
                     ${amrex_runs}/plane_rt_sweep --only gamma slab front conserve domains sweep
                     --set rt_backend=amrex --set amrex_max_grid_size=8
                     --launch "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} 2)
    set_tests_properties(amrex_plane_rt_sweep PROPERTIES TIMEOUT 1800 LABELS "amrex;plane;physics"
                         ENVIRONMENT OMP_NUM_THREADS=1 PROCESSORS 2)
    add_test(NAME amrex_plane_full
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/plane_full_check.py
                     $<TARGET_FILE:radhydro> ${amrex_runs}/plane_full --amrex $<TARGET_FILE:radhydro>
                     --launch "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} 2)
    set_tests_properties(amrex_plane_full PROPERTIES TIMEOUT 3600 LABELS "amrex;plane;physics;full"
                         ENVIRONMENT OMP_NUM_THREADS=1 PROCESSORS 2)
    # amrex_slab = 1: the shared decomposition as z-slabs (P planes; 0 = N /
    # ranks), byte-identical to the host state on the default boxes; N = 16:
    # 3 ranks hold 2, 1, 1 slabs of 4; P = 2 is thinner than the hydro's
    # ghosts and the rays' Redistribute reach
    foreach(rt "1 1 4" "2 1 2" "3 1 4" "2 2 0")
      separate_arguments(rt)
      list(GET rt 0 np)
      list(GET rt 1 nt)
      list(GET rt 2 planes)
      add_test(NAME amrex_slab_full_pair_r${np}t${nt}_p${planes}
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/backend_pair_check.py
                       --launch "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} ${np}
                       ${amrex_runs}/slab_pair_r${np}t${nt}_p${planes} $<TARGET_FILE:radhydro>
                       ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke "amrex_state=0"
                       "amrex_state=1,amrex_slab=1,amrex_slab_planes=${planes}"
                       --quiet --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1
                       --set rt_backend=amrex --set hydro_scheme=unsplit --set hydro_backend=amrex
                       --set amrex_hydro_max_grid_size=0 --set pm_backend=amrex --set movie_dz=0.5
                       --set movie_z0=25 --set restart_every=1 --set restart_keep_latest=0)
      set_tests_properties(amrex_slab_full_pair_r${np}t${nt}_p${planes} PROPERTIES TIMEOUT 1800
                           LABELS "amrex;amrex_state;amrex_slab;full" ENVIRONMENT OMP_NUM_THREADS=${nt}
                           PROCESSORS ${np})
    endforeach()
    # the full loop with the moment RT on the level and the state on the
    # slabs, byte-identical to the serial moment RT with the hydro and PM
    # backends (fields, stats, movie frames)
    foreach(rt "m1 1 1 0" "m1 2 1 4" "m1 3 2 4" "otvet 3 2 4" "split 2 2 2" "screensplit 3 1 4" "cone 2 2 4"
               "cone 4 1 4" "conerecomb 3 2 4")
      separate_arguments(rt)
      list(GET rt 0 mc)
      list(GET rt 1 np)
      list(GET rt 2 nt)
      list(GET rt 3 planes)
      if(mc STREQUAL "m1")
        set(margs --set rt_method=m1)
      elseif(mc STREQUAL "otvet")
        set(margs --set rt_method=otvet)
      elseif(mc STREQUAL "split")
        set(margs --set rt_method=otvet --set otvet_split=1)
      elseif(mc STREQUAL "cone")
        set(margs --set rt_method=otvet --set otvet_closure=cone --set otvet_cone_nside=2)
      elseif(mc STREQUAL "conerecomb")
        set(margs --set rt_method=otvet --set otvet_closure=cone --set otvet_cone_nside=2 --set recomb_case=A
                  --set recomb_rad=1)
      else()
        set(margs --set rt_method=otvet --set otvet_closure=screen --set otvet_split=1)
      endif()
      set(tname amrex_moment_full_pair_${mc}_r${np}t${nt}_p${planes})
      add_test(NAME ${tname}
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/backend_pair_check.py
                       --launch "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} ${np}
                       ${amrex_runs}/moment_pair_${mc}_r${np}t${nt}_p${planes} $<TARGET_FILE:radhydro>
                       ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke "rt_backend=serial,amrex_state=0"
                       "rt_backend=amrex,amrex_state=1,amrex_slab=1,amrex_slab_planes=${planes}"
                       --quiet ${margs} --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1
                       --set hydro_scheme=unsplit --set hydro_backend=amrex --set amrex_hydro_max_grid_size=0
                       --set pm_backend=amrex --set movie_dz=0.5 --set movie_z0=25 --set restart_every=0)
      set_tests_properties(${tname} PROPERTIES TIMEOUT 1800
                           LABELS "amrex;amrex_state;moments;amrex_moments;full" ENVIRONMENT OMP_NUM_THREADS=${nt}
                           PROCESSORS ${np})
    endforeach()
    # gammaHI_avg (the substep-averaged Gamma_HI) of the moment RT on the
    # level: the fields and a movie frame after every substep (the running
    # average) byte-identical to the serial moment RT, the state on the
    # slabs (state) or on rank 0's host arrays (host)
    foreach(rt "m1 state 2 1 4" "cone state 4 1 4" "otvet state 3 2 4" "m1 host 2 2 4")
      separate_arguments(rt)
      list(GET rt 0 mc)
      list(GET rt 1 md)
      list(GET rt 2 np)
      list(GET rt 3 nt)
      list(GET rt 4 planes)
      if(mc STREQUAL "m1")
        set(margs --set rt_method=m1)
      elseif(mc STREQUAL "otvet")
        set(margs --set rt_method=otvet)
      else()
        set(margs --set rt_method=otvet --set otvet_closure=cone --set otvet_cone_nside=2)
      endif()
      if(md STREQUAL "state")
        set(st 1)
        set(bargs --set hydro_backend=amrex --set amrex_hydro_max_grid_size=0 --set pm_backend=amrex)
      else()
        set(st 0)
        set(bargs)
      endif()
      set(tname amrex_moment_gamma_avg_${mc}_${md}_r${np}t${nt})
      add_test(NAME ${tname}
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/backend_pair_check.py
                       --launch "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} ${np}
                       ${amrex_runs}/moment_gamma_avg_${mc}_${md}_r${np}t${nt} $<TARGET_FILE:radhydro>
                       ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke "rt_backend=serial,amrex_state=0"
                       "rt_backend=amrex,amrex_state=${st},amrex_slab=1,amrex_slab_planes=${planes}"
                       --quiet ${margs} ${bargs} --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1
                       --set hydro_scheme=unsplit --set fields=fHI,gammaHI,gammaHI_avg
                       --set movie_fields=fHI,gammaHI,gammaHI_avg --set movie_dz=0.5 --set movie_z0=25
                       --set movie_every_substep=1 --set movie_slab=1
                       --set restart_every=0)
      set_tests_properties(${tname} PROPERTIES TIMEOUT 1800
                           LABELS "amrex;moments;amrex_moments;full" ENVIRONMENT OMP_NUM_THREADS=${nt}
                           PROCESSORS ${np})
    endforeach()
    # the moment RT on the level through a checkpoint (rank 0 streams the
    # moment and recombination blocks slab by slab): bit-identical restarts
    foreach(rc "m1 2 4" "split 1 0" "conerecomb 3 4")
      separate_arguments(rc)
      list(GET rc 0 mc)
      list(GET rc 1 np)
      list(GET rc 2 planes)
      if(mc STREQUAL "m1")
        set(margs --set rt_method=m1)
      elseif(mc STREQUAL "split")
        set(margs --set rt_method=otvet --set otvet_split=1)
      else()
        set(margs --set rt_method=otvet --set otvet_closure=cone --set otvet_cone_nside=2 --set recomb_case=A
                  --set recomb_rad=1)
      endif()
      add_test(NAME amrex_moment_restart_${mc}_r${np}
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                       --launch "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} ${np}
                       $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                       ${amrex_runs}/moment_restart_${mc}_r${np} 13 ${margs} --set legacy_rt=0 --set legacy_physics=0
                       --set rt_backend=amrex --set hydro_scheme=unsplit --set hydro_backend=amrex
                       --set amrex_hydro_max_grid_size=0 --set pm_backend=amrex --set amrex_state=1
                       --set amrex_slab=1 --set amrex_slab_planes=${planes} --set movie_dz=0.5 --set movie_z0=25)
      set_tests_properties(amrex_moment_restart_${mc}_r${np} PROPERTIES TIMEOUT 1800
                           LABELS "amrex;amrex_state;moments;amrex_moments;restart" ENVIRONMENT OMP_NUM_THREADS=1
                           PROCESSORS ${np})
    endforeach()
    add_test(NAME amrex_slab_restart_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${amrex_runs}/slab_restart 13 --set legacy_rt=0 --set legacy_physics=0
                     --set rt_backend=amrex --set hydro_scheme=unsplit --set hydro_backend=amrex
                     --set amrex_hydro_max_grid_size=0 --set pm_backend=amrex --set amrex_state=1
                     --set amrex_slab=1 --set amrex_slab_planes=4 --set movie_dz=0.5 --set movie_z0=25)
    set_tests_properties(amrex_slab_restart_smoke PROPERTIES TIMEOUT 1800
                         LABELS "amrex;amrex_state;amrex_slab;restart")
    # amrex_balance = 1: the mapping rebalanced every RT step, byte-identical
    # to the fixed mapping (the levels are mapping-independent), 2 and 3 ranks
    # (min_gain = -1 remaps at every RT step whatever the gain; knapsack on 2
    # ranks gives a mapping unlike the default's, sfc on 3)
    foreach(np 2 3)
      if(np EQUAL 2)
        set(bm knapsack)
      else()
        set(bm sfc)
      endif()
      add_test(NAME amrex_balance_full_pair_r${np}
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/backend_pair_check.py
                       --launch "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} ${np}
                       ${amrex_runs}/balance_pair_r${np} $<TARGET_FILE:radhydro>
                       ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke "amrex_balance=0"
                       "amrex_balance=1,amrex_balance_min_gain=-1,amrex_balance_method=${bm}"
                       --quiet --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1
                       --set rt_backend=amrex --set hydro_scheme=unsplit --set hydro_backend=amrex
                       --set amrex_hydro_max_grid_size=0 --set pm_backend=amrex --set amrex_state=1
                       --set amrex_max_grid_size=8 --set movie_dz=0.5 --set movie_z0=25
                       --set restart_every=1 --set restart_keep_latest=0)
      set_tests_properties(amrex_balance_full_pair_r${np} PROPERTIES TIMEOUT 1800
                           LABELS "amrex;amrex_state;full" ENVIRONMENT OMP_NUM_THREADS=1 PROCESSORS ${np})
    endforeach()
    add_test(NAME amrex_state_restart_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${amrex_runs}/state_restart 13 --set legacy_rt=0 --set legacy_physics=0
                     --set rt_backend=amrex --set hydro_scheme=unsplit --set hydro_backend=amrex
                     --set amrex_hydro_max_grid_size=0 --set pm_backend=amrex --set amrex_state=1
                     --set movie_dz=0.5 --set movie_z0=25)
    set_tests_properties(amrex_state_restart_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex;amrex_state;restart")
    # c_reduced_schedule on the levels: x on rank 0, broadcast; restart across
    # the switches (1 rank) and the switching rule on 2 ranks
    add_test(NAME amrex_state_restart_smoke_c_schedule
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${amrex_runs}/state_restart_c_schedule 13 --set legacy_rt=0 --set legacy_physics=0
                     --set rt_backend=amrex --set hydro_scheme=unsplit --set hydro_backend=amrex
                     --set amrex_hydro_max_grid_size=0 --set pm_backend=amrex --set amrex_state=1
                     --set movie_dz=0.5 --set movie_z0=25 --set c_reduced_schedule=0.3,0.005,0.6,0.01,1)
    set_tests_properties(amrex_state_restart_smoke_c_schedule PROPERTIES TIMEOUT 1800
                         LABELS "amrex;amrex_state;restart")
    add_test(NAME amrex_state_c_schedule_r2
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/c_schedule_check.py
                     --launch "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} 2
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${amrex_runs}/state_c_schedule_r2 --set legacy_rt=0 --set legacy_physics=0
                     --set rt_backend=amrex --set hydro_scheme=unsplit --set hydro_backend=amrex
                     --set amrex_hydro_max_grid_size=0 --set pm_backend=amrex --set amrex_state=1
                     --set amrex_max_grid_size=8)
    set_tests_properties(amrex_state_c_schedule_r2 PROPERTIES TIMEOUT 1800
                         LABELS "amrex;amrex_state" ENVIRONMENT OMP_NUM_THREADS=1 PROCESSORS 2)
    # a bit-identical restart through the hydro backend
    add_test(NAME amrex_hydro_restart_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${amrex_runs}/hydro_restart 13 --set legacy_rt=0 --set legacy_physics=0
                     --set rt_backend=amrex --set hydro_scheme=unsplit --set hydro_backend=amrex)
    set_tests_properties(amrex_hydro_restart_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex;amrex_hydro;restart")
    # nightly: the 64^3 production configuration, backend vs serial
    # (history, tau_es, per-cell fields; docs/stage_d.md section 4)
    add_test(NAME amrex_full64_vs_serial
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/amrex_full64_check.py
                     $<TARGET_FILE:radhydro> ${amrex_runs}/full64 ${PROJECT_SOURCE_DIR}/tools)
    set_tests_properties(amrex_full64_vs_serial PROPERTIES TIMEOUT 14400 LABELS "amrex;full;nightly"
                         ENVIRONMENT OMP_NUM_THREADS=4 RUN_SERIAL TRUE DISABLED $<NOT:$<BOOL:${RADRAY_NIGHTLY_TESTS}>>)
    # nightly: the same with today's production defaults (unsplit) and the
    # three backends with the state on the levels (amrex_state = 1), against
    # the serial code (docs/stage_e_pm.md section 3.1)
    add_test(NAME amrex_state_full64_vs_serial
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/amrex_full64_check.py
                     $<TARGET_FILE:radhydro> ${amrex_runs}/state_full64 ${PROJECT_SOURCE_DIR}/tools state)
    set_tests_properties(amrex_state_full64_vs_serial PROPERTIES TIMEOUT 14400
                         LABELS "amrex;amrex_state;full;nightly"
                         ENVIRONMENT OMP_NUM_THREADS=4 RUN_SERIAL TRUE DISABLED $<NOT:$<BOOL:${RADRAY_NIGHTLY_TESTS}>>)

    # ---- the PM backend (pm_backend = amrex, docs/stage_e_pm.md; label amrex_pm) ----
    # The dark matter, the CIC deposit, the Poisson solve and the kicks on
    # an AMReX level of the PM mesh. Bitwise identical for any rank count,
    # thread count and box size (int64 deposits, fixed-point statistics,
    # per-line FFTs); round-off against pm_driver's serial PM.
    set(PMC ${CMAKE_CURRENT_SOURCE_DIR}/amrex_pm_check.py)
    set(amrex_pm_runs ${CMAKE_CURRENT_BINARY_DIR}/amrex_pm_runs)
    # unit tests on the level: the single-particle force against Ewald, no
    # self-force, the gas kick, the mesh force against the serial code,
    # momentum conservation of every kick variant (1, 2 and 4 ranks, both FFTs)
    add_executable(test_amrex_pm test_amrex_pm.cpp)
    target_link_libraries(test_amrex_pm PRIVATE radray_pm ${RADRAY_RT_BACKEND_LIB})
    foreach(cfg "1;pencil;16" "2;pencil;8" "4;pencil;16" "1;amrex;16" "2;amrex;16")
      list(GET cfg 0 np)
      list(GET cfg 1 fft)
      list(GET cfg 2 mgs)
      add_test(NAME amrex_pm_unit_r${np}_${fft}
               COMMAND ${RADRAY_MPIEXEC} ${RADRAY_MPIEXEC_FLAGS} ${RADRAY_MPIEXEC_NP} ${np}
                       $<TARGET_FILE:test_amrex_pm> --fft ${fft} --mgs ${mgs})
      set_tests_properties(amrex_pm_unit_r${np}_${fft} PROPERTIES LABELS "amrex_pm;pm"
                           ENVIRONMENT OMP_NUM_THREADS=2 TIMEOUT 1800 PROCESSORS ${np})
    endforeach()
    # Stage G (pm_level_solve, docs/stage_g.md): the level Poisson solve on a
    # gas state with refined levels built by the test (no hydro): the
    # single-particle force against Ewald with level 1 everywhere, over a
    # region and with two nested levels; no self-force; level 1 everywhere
    # against the uniform pm_refine = 2 backend (the same discrete problem:
    # MLMG against the FFT); the DM + gas momentum over one step (1, 2 and 3
    # ranks; the test needs C++20 for the AMReX headers it includes)
    add_executable(test_amrex_pm_levels test_amrex_pm_levels.cpp)
    target_link_libraries(test_amrex_pm_levels PRIVATE radray_pm ${RADRAY_RT_BACKEND_LIB})
    set_target_properties(test_amrex_pm_levels PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
    foreach(np 1 2 3)
      add_test(NAME amrex_pm_levels_unit_r${np}
               COMMAND ${RADRAY_MPIEXEC} ${RADRAY_MPIEXEC_FLAGS} ${RADRAY_MPIEXEC_NP} ${np}
                       $<TARGET_FILE:test_amrex_pm_levels> --mgs 8)
      set_tests_properties(amrex_pm_levels_unit_r${np} PROPERTIES LABELS "amrex_pm;pm;amr;level_solve"
                           ENVIRONMENT OMP_NUM_THREADS=2 TIMEOUT 1800 PROCESSORS ${np})
    endforeach()
    # round-off against the serial pm_driver on the PM problems (fixtures of
    # the pm label), 2 ranks x 2 threads: every value within 1e-4 relative or
    # 1e-5 of the field's maximum (the serial float CIC sum and 3-D FFTW plan
    # against the exact deposit and the per-line transforms: 1e-7 per cell,
    # amplified by the gradient where the force is small and over 30 steps,
    # growth_3d: 5e-6 of the maximum), the stats within 1e-4 (cmax_gas of the
    # non-MF kick: 1e-5); the backend's run is the fixture of the physics
    # checks below
    foreach(p ${PM_PROBLEMS})
      add_test(NAME amrex_pm_roundoff_${p}
               COMMAND ${Python3_EXECUTABLE} ${PMC} runcmp roundoff ${CMAKE_CURRENT_BINARY_DIR}/pm_runs/${p}
                       ${amrex_pm_runs}/${p} ${RADRAY_MPIEXEC} "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} 2 2 $<TARGET_FILE:pm_driver>
                       ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/${p}
                       --rtol 1e-4 --atol-frac 1e-5 --stats-rtol 1e-4
                       -- --set pm_backend=amrex --set amrex_pm_max_grid_size=8)
      set_tests_properties(amrex_pm_roundoff_${p} PROPERTIES FIXTURES_REQUIRED pm_${p}
                           FIXTURES_SETUP amrex_pm_${p} LABELS "amrex_pm;pm" TIMEOUT 1800 PROCESSORS 2)
      add_test(NAME amrex_pm_momentum_${p}
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                       momentum ${amrex_pm_runs}/${p} ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/${p})
      set_tests_properties(amrex_pm_momentum_${p} PROPERTIES FIXTURES_REQUIRED amrex_pm_${p}
                           LABELS "amrex_pm;pm")
    endforeach()
    # the analytic pancakes, the growth test and the r = 2 twin on the backend's runs
    foreach(case "pancake_x;" "pancake_x_r2;--tol-x;0.045" "pancake_x_n32;--tol-x;0.045"
                 "pancake_z_nomf;--tol-mode-x;0.025;--tol-mode-v;0.025"
                 "pancake_z_nomf_r2;--tol-x;0.045;--tol-mode-x;0.025;--tol-mode-v;0.025")
      list(GET case 0 p)
      list(REMOVE_AT case 0)
      add_test(NAME amrex_pm_pancake_${p}
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                       pancake ${amrex_pm_runs}/${p} ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/${p} ${case})
      set_tests_properties(amrex_pm_pancake_${p} PROPERTIES FIXTURES_REQUIRED amrex_pm_${p}
                           LABELS "amrex_pm;pm")
    endforeach()
    add_test(NAME amrex_pm_growth_3d
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                     growth ${amrex_pm_runs}/growth_3d ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/growth_3d)
    set_tests_properties(amrex_pm_growth_3d PROPERTIES FIXTURES_REQUIRED amrex_pm_growth_3d LABELS "amrex_pm;pm")
    add_test(NAME amrex_pm_r2_twin_pancake_x
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pm_physics_check.py
                     same-dm ${amrex_pm_runs}/pancake_x_r2 ${amrex_pm_runs}/pancake_x_n32)
    set_tests_properties(amrex_pm_r2_twin_pancake_x PROPERTIES
                         FIXTURES_REQUIRED "amrex_pm_pancake_x_r2;amrex_pm_pancake_x_n32" LABELS "amrex_pm;pm")
    # bitwise matrix: two problems (pm_refine 1 and 2) x both FFTs x (ranks,
    # threads, box size); the reference is the 1 x 1 run with one box
    foreach(p zeldovich_3d zeldovich_3d_r2)
      foreach(fft pencil amrex)
        set(ref ${amrex_pm_runs}/bit_${p}_${fft}_ref)
        add_test(NAME amrex_pm_bitref_${p}_${fft}
                 COMMAND ${Python3_EXECUTABLE} ${PMC} run ${ref} ${RADRAY_MPIEXEC} "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} 1 1 $<TARGET_FILE:pm_driver>
                         ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/${p} --nsteps 6
                         --set pm_backend=amrex --set amrex_pm_fft=${fft} --set amrex_pm_max_grid_size=32)
        set_tests_properties(amrex_pm_bitref_${p}_${fft} PROPERTIES FIXTURES_SETUP amrex_pm_bit_${p}_${fft}
                             LABELS "amrex_pm;pm" TIMEOUT 1800)
        foreach(cfg "1;2;16" "2;1;16" "2;2;8" "4;1;16" "4;2;8" "8;1;8" "8;2;16")
          list(GET cfg 0 np)
          list(GET cfg 1 nt)
          list(GET cfg 2 mgs)
          set(tname amrex_pm_bitwise_${p}_${fft}_r${np}t${nt}_mgs${mgs})
          add_test(NAME ${tname}
                   COMMAND ${Python3_EXECUTABLE} ${PMC} runcmp bitwise ${ref} ${amrex_pm_runs}/${tname}
                           ${RADRAY_MPIEXEC} "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP} ${np} ${nt} $<TARGET_FILE:pm_driver>
                           ${CMAKE_CURRENT_SOURCE_DIR}/problems_pm/${p}
                           -- --nsteps 6 --set pm_backend=amrex --set amrex_pm_fft=${fft}
                           --set amrex_pm_max_grid_size=${mgs})
          set_tests_properties(${tname} PROPERTIES FIXTURES_REQUIRED amrex_pm_bit_${p}_${fft}
                               LABELS "amrex_pm;pm" TIMEOUT 1800 PROCESSORS ${np})
        endforeach()
      endforeach()
    endforeach()
    # the full loop with rt_backend = amrex and pm_backend = amrex: 1, 2 and
    # 4 ranks byte-identical (fields, movie frames, checkpoints), 1/2/4
    # threads byte-identical, bit-identical restarts (legacy and unsplit
    # hydro)
    add_test(NAME amrex_pm_full_ranks_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rank_check.py
                     ${amrex_pm_runs}/full_ranks "${RADRAY_MPIEXEC}" "${RADRAY_MPIEXEC_FLAGS}" ${RADRAY_MPIEXEC_NP}
                     "1 2 4" $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet
                     --set legacy_rt=0 --set legacy_physics=0 --set fields_every=1 --set rt_backend=amrex
                     --set pm_backend=amrex --set movie_dz=0.5 --set movie_z0=25 --set restart_every=1
                     --set restart_keep_latest=0)
    set_tests_properties(amrex_pm_full_ranks_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex_pm;full"
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
    add_test(NAME amrex_pm_thread_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/thread_check.py
                     ${amrex_pm_runs}/thread_smoke $<TARGET_FILE:radhydro>
                     ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke --quiet --set legacy_rt=0
                     --set legacy_physics=0 --set fields_every=1 --set rt_backend=amrex --set pm_backend=amrex
                     --set hydro_scheme=unsplit)
    set_tests_properties(amrex_pm_thread_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex_pm;threads"
                         ENVIRONMENT OMP_NUM_THREADS=1 RUN_SERIAL TRUE)
    add_test(NAME amrex_pm_restart_smoke
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${amrex_pm_runs}/restart 13 --set legacy_rt=0 --set legacy_physics=0
                     --set rt_backend=amrex --set pm_backend=amrex --set movie_dz=0.5 --set movie_z0=25)
    set_tests_properties(amrex_pm_restart_smoke PROPERTIES TIMEOUT 1800 LABELS "amrex_pm;restart")
    add_test(NAME amrex_pm_restart_smoke_unsplit
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/smoke
                     ${amrex_pm_runs}/restart_unsplit 13 --set legacy_rt=0 --set legacy_physics=0
                     --set rt_backend=amrex --set pm_backend=amrex --set hydro_scheme=unsplit)
    set_tests_properties(amrex_pm_restart_smoke_unsplit PROPERTIES TIMEOUT 1800 LABELS "amrex_pm;restart")
    # the pm_refine = 2 full loop (Ndm = 2N, on-the-fly sources) through the
    # backend and its checkpoint
    add_test(NAME amrex_pm_restart_smoke_pm_refine
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/restart_check.py
                     $<TARGET_FILE:radhydro> ${CMAKE_CURRENT_SOURCE_DIR}/problems_full/onthefly_r2
                     ${amrex_pm_runs}/restart_pm_refine 22 --set legacy_rt=0 --set legacy_physics=0
                     --set rt_backend=amrex --set pm_backend=amrex --set hydro_scheme=unsplit
                     --set deterministic_parallel=0)
    set_tests_properties(amrex_pm_restart_smoke_pm_refine PROPERTIES TIMEOUT 1800 LABELS "amrex_pm;restart")
  endif()
endif()

# ---------------------------------------------------------------------------
# OpenMP: every test runs with 2 threads unless it sets OMP_NUM_THREADS
# itself (ctest -j runs several tests at once; results do not depend on the
# thread count, see docs/parallel.md and the thread_* tests).
get_property(RADRAY_ALL_TESTS DIRECTORY PROPERTY TESTS)
foreach(t ${RADRAY_ALL_TESTS})
  get_test_property(${t} ENVIRONMENT tenv)
  if(NOT tenv MATCHES "OMP_NUM_THREADS")
    set_property(TEST ${t} APPEND PROPERTY ENVIRONMENT OMP_NUM_THREADS=2)
  endif()
endforeach()

# Hydro accuracy test suite (docs/hydro_tests.md)
if(RADRAY_HAVE_PM)
  add_subdirectory(hydro_suite)
endif()
