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
