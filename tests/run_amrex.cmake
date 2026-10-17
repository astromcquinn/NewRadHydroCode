# Runs rt_driver (rt_backend = amrex) under the MPI launcher and compares
# with a reference run directory.
# Variables: LAUNCH LAUNCH_FLAGS NPFLAG NP NT DRIVER PROBLEM OUT REF SET
#            PYTHON COMPARE STRICT
#   SET     params.txt overrides, comma-separated (--set key=value each)
#   EXTRA   further driver arguments (list)
#   STRICT  1: bitwise (--rtol 0 --atol-frac 0); cons: COMPARE is
#           tests/amr_conservation_check.py (the conserved totals of every
#           step, EXTRA_TOL e.g. --levels;2); 0: round-off against the
#           serial port: rtol 1e-3 on the fields (single cells at the
#           removal threshold or in I-fronts amplify the round-off of the
#           deposit and merge order to 1e-4), dx excepted (dx=1: the
#           direction field of source cells is a float sum of 192 symmetric
#           directions against an exact integer sum), float stats rtol 1e-6
#           (ngam_rays_total is a double sum in a different order), integer
#           stats exact; lossy: the R-struct encoding (rtol 0.1, integer ray
#           stats within 5 %)
#   EXTRA_TOL  further compare.py arguments (list), e.g. per-field rtol
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
  # srun inherits the allocation's --cpus-per-task (SLURM_CPUS_PER_TASK), which
  # multiplied by the rank count exceeds the allocation: give it the threads
  if(LAUNCH MATCHES "srun$" AND NOT LAUNCH_FLAGS MATCHES "cpus-per-task")
    list(APPEND launch --cpus-per-task=${NT})
  endif()
endif()
set(ENV{OMP_NUM_THREADS} ${NT})
execute_process(COMMAND ${launch} ${DRIVER} ${PROBLEM} ${OUT} ${extra}
                OUTPUT_FILE ${OUT}.stdout ERROR_FILE ${OUT}.stderr
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "rt_driver (amrex) failed (${rc}), see ${OUT}.stdout / ${OUT}.stderr")
endif()
if(STRICT STREQUAL "cons")
  execute_process(COMMAND ${PYTHON} ${COMPARE} ${REF} ${OUT} ${EXTRA_TOL} RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "conservation check against ${REF} failed")
  endif()
  return()
endif()
if(STRICT STREQUAL "1")
  set(tol --rtol 0 --atol-frac 0)
elseif(STRICT STREQUAL "lossy")
  set(tol --rtol 0.1 --atol-frac 1e-3 --stats-int-tol 0.05 --stats-rtol 0.1 --ignore dx,nu,fgam,ngam,gamma,lambda)
else()
  set(tol --rtol 1e-3 --atol-frac 1e-6 --field-rtol dx=1 --stats-rtol 1e-6)
endif()
if(EXTRA_TOL)
  list(APPEND tol ${EXTRA_TOL})
endif()
execute_process(COMMAND ${PYTHON} ${COMPARE} ${REF} ${OUT} ${tol} -q RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "outputs differ from ${REF}")
endif()
