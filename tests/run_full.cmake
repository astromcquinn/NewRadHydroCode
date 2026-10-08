# Full radhydro loop CTest driver (docs/full_run.md).
#   MODE=ref      run DRIVER on PROBLEM (with SET, comma-separated key=value)
#                 into OUT and compare with the stored Fortran reference REF
#                 (compare_full_reference.py: every number and byte identical)
#   MODE=fortran  run the gfortran full harness (HARNESS dir, linking the C++
#                 halo finder from RADRAY_BUILD) and DRIVER on PROBLEM with SET,
#                 compare with compare.py --rtol 0 --atol-frac 0 and
#                 compare_full_reference.py
#   MODE=sanitize same as ref with an ASan/UBSan build of DRIVER
# Other variables: PYTHON, TESTS (tests source dir), COMPARE (compare.py)
set(setargs)
if(DEFINED SET AND NOT SET STREQUAL "")
  string(REPLACE "," ";" setlist "${SET}")
  foreach(kv ${setlist})
    list(APPEND setargs --set ${kv})
  endforeach()
endif()

function(run_driver out)
  file(REMOVE_RECURSE ${out})
  file(MAKE_DIRECTORY ${out})
  execute_process(COMMAND ${DRIVER} ${PROBLEM} ${out} --quiet ${setargs}
                  OUTPUT_FILE ${out}/driver.stdout ERROR_FILE ${out}/driver.stderr
                  RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    file(READ ${out}/driver.stderr err)
    message(FATAL_ERROR "radhydro failed (${rc}) on ${PROBLEM}:\n${err}")
  endif()
endfunction()

function(compare_ref ref out)
  execute_process(COMMAND ${PYTHON} ${TESTS}/compare_full_reference.py ${ref} ${out}
                  RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "outputs of ${out} differ from ${ref}")
  endif()
endfunction()

if(MODE STREQUAL "ref")
  run_driver(${OUT})
  compare_ref(${REF} ${OUT})
elseif(MODE STREQUAL "sanitize")
  set(ENV{ASAN_OPTIONS} "detect_leaks=1:abort_on_error=1")
  set(ENV{UBSAN_OPTIONS} "halt_on_error=1:print_stacktrace=1")
  run_driver(${OUT})
  compare_ref(${REF} ${OUT})
elseif(MODE STREQUAL "fortran")
  set(fout ${OUT}/fortran)
  file(REMOVE_RECURSE ${fout})
  set(ENV{RADRAY_BUILD} ${RADRAY_BUILD})
  execute_process(COMMAND ${HARNESS}/run_problem.sh ${PROBLEM} ${fout} ${setargs}
                  RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "Fortran full harness failed on ${PROBLEM}")
  endif()
  run_driver(${OUT}/cpp)
  get_filename_component(pname ${PROBLEM} NAME)
  execute_process(COMMAND ${PYTHON} ${COMPARE} ${fout} ${OUT}/cpp --rtol 0 --atol-frac 0 -q
                  RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "C++ differs from the Fortran full harness on ${pname}")
  endif()
  compare_ref(${fout} ${OUT}/cpp)
else()
  message(FATAL_ERROR "unknown MODE ${MODE}")
endif()
