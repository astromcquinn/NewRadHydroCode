# Runs rt_driver (rt_backend = amrex) under the MPI launcher and compares
# with a reference run directory.
# Variables: LAUNCH LAUNCH_FLAGS NPFLAG NP NT DRIVER PROBLEM OUT REF SET
#            PYTHON COMPARE STRICT
#   SET     params.txt overrides, comma-separated (--set key=value each)
#   EXTRA   further driver arguments (list)
#   STRICT  1: bitwise (--rtol 0 --atol-frac 0); 0: round-off against the
#           serial port (rtol 1e-5 on the fields, dx and fgam excepted:
#           dx=1 (the direction field of source cells), fgam=2e-3 (I-front
#           cells), stats rtol 1e-9); lossy: the R-struct encoding
#           (rtol 0.1, integer ray stats within 5 %)
file(REMOVE_RECURSE ${OUT})
file(MAKE_DIRECTORY ${OUT})
set(extra)
if(SET)
  string(REPLACE "," ";" set_list "${SET}")
  foreach(kv IN LISTS set_list)
    list(APPEND extra --set ${kv})
  endforeach()
endif()
if(EXTRA)
  list(APPEND extra ${EXTRA})
endif()
set(launch)
if(NP GREATER 1 OR LAUNCH)
  if(NOT LAUNCH)
    message(FATAL_ERROR "no MPI launcher (RADRAY_MPIEXEC) for ${NP} ranks")
  endif()
  set(launch ${LAUNCH} ${LAUNCH_FLAGS} ${NPFLAG} ${NP})
endif()
set(ENV{OMP_NUM_THREADS} ${NT})
execute_process(COMMAND ${launch} ${DRIVER} ${PROBLEM} ${OUT} ${extra}
                OUTPUT_FILE ${OUT}.stdout ERROR_FILE ${OUT}.stderr
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "rt_driver (amrex) failed (${rc}), see ${OUT}.stdout / ${OUT}.stderr")
endif()
if(STRICT STREQUAL "1")
  set(tol --rtol 0 --atol-frac 0)
elseif(STRICT STREQUAL "lossy")
  set(tol --rtol 0.1 --atol-frac 1e-3 --stats-int-tol 0.05 --stats-rtol 0.1 --ignore dx,nu,fgam,ngam,gamma,lambda)
else()
  set(tol --rtol 1e-5 --atol-frac 1e-6 --field-rtol dx=1 --field-rtol fgam=2e-3 --stats-rtol 1e-9)
endif()
execute_process(COMMAND ${PYTHON} ${COMPARE} ${REF} ${OUT} ${tol} -q RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "outputs differ from ${REF}")
endif()
