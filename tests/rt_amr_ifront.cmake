# The refined I-front test (docs/stage_f.md section 3): generates the
# one-source uniform problems at N and 2N (tests/problems_full/make_ifront.py),
# runs the uniform N, the refined N (one level over the central HALF-radius
# box) and the uniform 2N cases under the MPI launcher, and checks the
# refined front against both (tools/analysis/ifront_full.py --check).
# Variables: LAUNCH LAUNCH_FLAGS NPFLAG NP NT DRIVER PYTHON SRC OUT N HALF
file(REMOVE_RECURSE ${OUT})
file(MAKE_DIRECTORY ${OUT})
math(EXPR N2 "2 * ${N}")
foreach(n ${N} ${N2})
  execute_process(COMMAND ${PYTHON} ${SRC}/tests/problems_full/make_ifront.py ${OUT}/prob${n} --N ${n}
                  RESULT_VARIABLE rc OUTPUT_QUIET)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "make_ifront.py failed for N = ${n}")
  endif()
endforeach()
set(launch)
if(NP GREATER 1 OR LAUNCH)
  if(NOT LAUNCH)
    message(FATAL_ERROR "no MPI launcher (RADRAY_MPIEXEC) for ${NP} ranks")
  endif()
  set(launch ${LAUNCH} ${LAUNCH_FLAGS} ${NPFLAG} ${NP})
  if(LAUNCH MATCHES "srun$" AND NOT LAUNCH_FLAGS MATCHES "cpus-per-task")
    list(APPEND launch --cpus-per-task=${NT})
  endif()
endif()
set(ENV{OMP_NUM_THREADS} ${NT})
math(EXPR q0 "${N} / 2 - ${HALF}")
math(EXPR q1 "${N} / 2 + ${HALF} - 1")
function(run name prob)
  execute_process(COMMAND ${launch} ${DRIVER} ${prob} ${OUT}/${name} --quiet --set run_steps=40 ${ARGN}
                  OUTPUT_FILE ${OUT}/${name}.stdout ERROR_FILE ${OUT}/${name}.stderr RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "radhydro failed (${rc}) on ${name}, see ${OUT}/${name}.stderr")
  endif()
endfunction()
run(uniform${N} ${OUT}/prob${N})
run(amr${N} ${OUT}/prob${N} --set amrex_max_level=1 --set amrex_refine_mass=1e9
    --set amrex_refine_box=${q0},${q0},${q0},${q1},${q1},${q1})
run(uniform${N2} ${OUT}/prob${N2})
execute_process(COMMAND ${PYTHON} ${SRC}/tools/analysis/ifront_full.py --N ${N}
                        --check ${OUT}/uniform${N2} ${OUT}/uniform${N} ${OUT}/amr${N} 0.25
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "the refined I-front is off (ifront_full.py --check)")
endif()
