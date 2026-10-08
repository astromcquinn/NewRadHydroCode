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
  endforeach()
endif()
