# Hydro CTest driver.
#   MODE=ref      run DRIVER on PROBLEM into OUT, compare with REF
#                 (compare_hydro_reference.py: byte/number identical)
#   MODE=physics  run DRIVER on every problem of PROBLEMS into OUT, then
#                 hydro_physics_check.py
#   MODE=fortran  run the gfortran harness (HARNESS dir) and DRIVER on every
#                 problem, compare with compare.py --rtol 0 --atol-frac 0
#   MODE=sanitize run DRIVER (ASan/UBSan build) on PROBLEM for NSTEPS steps
#                 and compare the steps it wrote with REF
# Other variables: PYTHON, TESTS (tests source dir), COMPARE (compare.py)
function(run_driver prob out)
  file(REMOVE_RECURSE ${out})
  file(MAKE_DIRECTORY ${out})
  set(extra)
  if(DEFINED NSTEPS)
    set(extra --nsteps ${NSTEPS})
  endif()
  execute_process(COMMAND ${DRIVER} ${prob} ${out} --quiet ${extra}
                  OUTPUT_FILE ${out}/driver.stdout ERROR_FILE ${out}/driver.stderr
                  RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    file(READ ${out}/driver.stderr err)
    message(FATAL_ERROR "hydro_driver failed (${rc}) on ${prob}:\n${err}")
  endif()
endfunction()

if(MODE STREQUAL "ref")
  run_driver(${PROBLEM} ${OUT})
  execute_process(COMMAND ${PYTHON} ${TESTS}/compare_hydro_reference.py ${REF} ${OUT}
                  RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "outputs differ from the Fortran reference")
  endif()
elseif(MODE STREQUAL "physics")
  file(GLOB probs LIST_DIRECTORIES true ${PROBLEMS}/*)
  foreach(p ${probs})
    if(IS_DIRECTORY ${p})
      get_filename_component(name ${p} NAME)
      run_driver(${p} ${OUT}/${name})
    endif()
  endforeach()
  execute_process(COMMAND ${PYTHON} ${TESTS}/hydro_physics_check.py ${OUT}
                  RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "hydro physics checks failed")
  endif()
elseif(MODE STREQUAL "fortran")
  file(GLOB probs LIST_DIRECTORIES true ${PROBLEMS}/*)
  foreach(p ${probs})
    if(IS_DIRECTORY ${p})
      get_filename_component(name ${p} NAME)
      file(REMOVE_RECURSE ${OUT}/fortran/${name})
      execute_process(COMMAND ${HARNESS}/run_problem.sh ${p} ${OUT}/fortran/${name}
                      RESULT_VARIABLE rc)
      if(NOT rc EQUAL 0)
        message(FATAL_ERROR "Fortran harness failed on ${name}")
      endif()
      run_driver(${p} ${OUT}/cpp/${name})
      execute_process(COMMAND ${PYTHON} ${COMPARE} ${OUT}/fortran/${name} ${OUT}/cpp/${name}
                              --rtol 0 --atol-frac 0 -q
                      RESULT_VARIABLE rc)
      if(NOT rc EQUAL 0)
        message(FATAL_ERROR "C++ differs from the Fortran harness on ${name}")
      endif()
      execute_process(COMMAND ${PYTHON} ${TESTS}/compare_hydro_reference.py
                              ${OUT}/fortran/${name} ${OUT}/cpp/${name}
                      RESULT_VARIABLE rc)
      if(NOT rc EQUAL 0)
        message(FATAL_ERROR "C++ text outputs differ from the Fortran harness on ${name}")
      endif()
    endif()
  endforeach()
elseif(MODE STREQUAL "sanitize")
  set(ENV{ASAN_OPTIONS} "detect_leaks=1:abort_on_error=1")
  set(ENV{UBSAN_OPTIONS} "halt_on_error=1:print_stacktrace=1")
  run_driver(${PROBLEM} ${OUT})
  # the sanitized build must still be bit-identical on the steps it ran
  set(sub ${OUT}.ref)
  file(REMOVE_RECURSE ${sub})
  file(GLOB steps LIST_DIRECTORIES true ${OUT}/step_*)
  file(COPY ${REF}/domain_gas.txt DESTINATION ${sub})
  foreach(s ${steps})
    get_filename_component(sn ${s} NAME)
    if(EXISTS ${REF}/${sn})
      file(COPY ${REF}/${sn}/stats.txt ${REF}/${sn}/gas_summary.txt DESTINATION ${sub}/${sn})
    endif()
  endforeach()
  execute_process(COMMAND ${PYTHON} ${TESTS}/compare_hydro_reference.py ${sub} ${OUT}
                  RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "sanitized run differs from the Fortran reference")
  endif()
else()
  message(FATAL_ERROR "unknown MODE ${MODE}")
endif()
