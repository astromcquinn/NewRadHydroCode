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

# Halo finder and sources (docs/halos_and_sources.md)
add_executable(test_halos test_halos.cpp)
target_link_libraries(test_halos PRIVATE radray_halos)
add_test(NAME halos COMMAND test_halos)

add_executable(test_halo_catalog test_halo_catalog.cpp)
target_link_libraries(test_halo_catalog PRIVATE radray_halos)
find_program(RADRAY_GFORTRAN gfortran)
if(RADRAY_GFORTRAN)
  # Reference program with the original source_type layout and read_halo
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/halo_catalog_ref
    COMMAND ${RADRAY_GFORTRAN} -O0 -o ${CMAKE_CURRENT_BINARY_DIR}/halo_catalog_ref
            ${CMAKE_CURRENT_SOURCE_DIR}/halo_catalog_ref.f90
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/halo_catalog_ref.f90)
  add_custom_target(halo_catalog_ref ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/halo_catalog_ref)
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
endif()

# ---------------------------------------------------------------------------
# PM gravity (only when FFTW3f was found, see the top-level CMakeLists.txt)
if(RADRAY_HAVE_PM)
  add_executable(test_pm test_pm.cpp)
  target_link_libraries(test_pm PRIVATE radray_pm)
  add_test(NAME pm_unit COMMAND test_pm)
  set_tests_properties(pm_unit PROPERTIES LABELS pm)

  if(Python3_Interpreter_FOUND)
    set(PM_PROBLEMS zeldovich_3d growth_3d pancake_x pancake_z_nomf)
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
  # fields), if gfortran is available.
  find_program(GFORTRAN_EXE gfortran)
  if(GFORTRAN_EXE)
    add_test(NAME hydro_vs_fortran_live
             COMMAND ${HYDRO_RUN} -DMODE=fortran -DDRIVER=$<TARGET_FILE:hydro_driver>
                     -DPROBLEMS=${HYDRO_PROBLEMS}
                     -DHARNESS=${PROJECT_SOURCE_DIR}/reference/harness_hydro
                     -DOUT=${CMAKE_CURRENT_BINARY_DIR}/hydro_live_runs
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/run_hydro.cmake)
    set_tests_properties(hydro_vs_fortran_live PROPERTIES TIMEOUT 1200 LABELS "hydro;fortran")
  endif()

  # ASan/UBSan build of the hydro driver (all radray sources), run on two
  # problems for a few steps; the output must still match the reference.
  if(NOT RADRAY_SANITIZE)
    file(GLOB RADRAY_ALL_SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)
    add_executable(hydro_driver_asan ${PROJECT_SOURCE_DIR}/apps/hydro_driver.cpp ${RADRAY_ALL_SRC})
    target_include_directories(hydro_driver_asan PRIVATE ${PROJECT_SOURCE_DIR}/include
                               ${FFTW3F_INCLUDE_DIR})
    target_compile_options(hydro_driver_asan PRIVATE -fsanitize=address,undefined
                           -fno-sanitize-recover=undefined -fno-omit-frame-pointer -g -O1)
    target_link_options(hydro_driver_asan PRIVATE -fsanitize=address,undefined)
    target_link_libraries(hydro_driver_asan PRIVATE ${FFTW3F_LIBRARY})
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
  # smoke_*: seconds; legacy_* / onthefly_*: ~4 minutes each (label long).
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
        set_tests_properties(full_ref_${c} PROPERTIES TIMEOUT 3600
                             LABELS "full;reference;long")
      endif()
    endif()
  endforeach()

  # Live comparison with a fresh gfortran build of the full radhydro.f90
  # (smoke problem, both source modes and both flag sets).
  find_program(GFORTRAN_FULL gfortran)
  if(GFORTRAN_FULL)
    foreach(p smoke)
      foreach(src legacy onthefly)
        foreach(m rt1 rt0)
          add_test(NAME full_vs_fortran_live_${p}_${src}_${m}
                   COMMAND ${FULL_RUN} -DMODE=fortran -DDRIVER=$<TARGET_FILE:radhydro>
                           -DPROBLEM=${FULL_PROBLEMS}/${p}
                           -DSET=source_mode=${src},${FULL_SET_${m}}
                           -DHARNESS=${PROJECT_SOURCE_DIR}/reference/harness_full
                           -DRADRAY_BUILD=${CMAKE_BINARY_DIR}
                           -DOUT=${CMAKE_CURRENT_BINARY_DIR}/full_live_runs/${p}_${src}_${m}
                           -P ${CMAKE_CURRENT_SOURCE_DIR}/run_full.cmake)
          set_tests_properties(full_vs_fortran_live_${p}_${src}_${m} PROPERTIES
                               TIMEOUT 1800 LABELS "full;fortran"
                               RESOURCE_LOCK harness_full_build)
        endforeach()
      endforeach()
    endforeach()
  endif()

  # ASan/UBSan build of radhydro (all radray sources) on the smoke problem;
  # it must still match the reference bit for bit.
  if(NOT RADRAY_SANITIZE)
    file(GLOB RADRAY_ALL_SRC_FULL ${PROJECT_SOURCE_DIR}/src/*.cpp)
    add_executable(radhydro_asan ${PROJECT_SOURCE_DIR}/apps/radhydro.cpp ${RADRAY_ALL_SRC_FULL})
    target_include_directories(radhydro_asan PRIVATE ${PROJECT_SOURCE_DIR}/include
                               ${FFTW3F_INCLUDE_DIR})
    target_compile_options(radhydro_asan PRIVATE -fsanitize=address,undefined
                           -fno-sanitize-recover=undefined -fno-omit-frame-pointer -g -O1
                           -Wno-unknown-pragmas)
    target_link_options(radhydro_asan PRIVATE -fsanitize=address,undefined)
    target_link_libraries(radhydro_asan PRIVATE ${FFTW3F_LIBRARY})
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

# Hydro accuracy test suite (docs/hydro_tests.md)
if(RADRAY_HAVE_PM)
  add_subdirectory(hydro_suite)
endif()
