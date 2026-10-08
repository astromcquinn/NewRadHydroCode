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
                               ${FFTW3_INCLUDE_DIR})
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
