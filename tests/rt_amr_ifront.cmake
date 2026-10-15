# The refined I-front test (docs/stage_f.md section 3): generates the
# one-source uniform problems at N and 2N (tests/problems_full/make_ifront.py),
# runs the uniform N, the refined N (one level over the central HALF-radius
# box) and the uniform 2N cases under the MPI launcher, and checks the
# refined front against both (tools/analysis/ifront_full.py --check).
# Variables: LAUNCH LAUNCH_FLAGS NPFLAG NP NT DRIVER PYTHON SRC OUT N HALF
# Optional (the 8 sources of uniform_8src): SOURCES (make_ifront.py
# --sources), BOX (amrex_refine_box, one box per level, instead of HALF;
# amrex_n_error_buf = 0), LEVEL (refined levels, default 1; the reference is
# uniform 2^LEVEL N); with SOURCES the check is tools/analysis/ifront_multi.py
# (per-source fronts, x_HII,V/M, photons in flight, no jump).
file(REMOVE_RECURSE ${OUT})
file(MAKE_DIRECTORY ${OUT})
if(NOT DEFINED LEVEL)
  set(LEVEL 1)
endif()
math(EXPR N2 "${N} << ${LEVEL}")
set(src_args)
if(SOURCES)
  set(src_args --sources ${SOURCES})
endif()
foreach(n ${N} ${N2})
  execute_process(COMMAND ${PYTHON} ${SRC}/tests/problems_full/make_ifront.py ${OUT}/prob${n} --N ${n} ${src_args}
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
if(BOX)
  set(refine --set amrex_refine_box=${BOX} --set amrex_n_error_buf=0)
else()
  math(EXPR q0 "${N} / 2 - ${HALF}")
  math(EXPR q1 "${N} / 2 + ${HALF} - 1")
  set(refine --set amrex_refine_box=${q0},${q0},${q0},${q1},${q1},${q1})
endif()
function(run name prob)
  execute_process(COMMAND ${launch} ${DRIVER} ${prob} ${OUT}/${name} --quiet --set run_steps=40 ${ARGN}
                  OUTPUT_FILE ${OUT}/${name}.stdout ERROR_FILE ${OUT}/${name}.stderr RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "radhydro failed (${rc}) on ${name}, see ${OUT}/${name}.stderr")
  endif()
endfunction()
run(uniform${N} ${OUT}/prob${N})
run(amr${N} ${OUT}/prob${N} --set amrex_max_level=${LEVEL} --set amrex_refine_mass=1e9 ${refine})
run(uniform${N2} ${OUT}/prob${N2})
if(SOURCES)
  execute_process(COMMAND ${PYTHON} ${SRC}/tools/analysis/ifront_multi.py --N ${N}
                          --sources ${OUT}/prob${N}/sources_pos.txt
                          --check ${OUT}/uniform${N2} ${OUT}/uniform${N} ${OUT}/amr${N}
                  RESULT_VARIABLE rc)
else()
  execute_process(COMMAND ${PYTHON} ${SRC}/tools/analysis/ifront_full.py --N ${N}
                          --check ${OUT}/uniform${N2} ${OUT}/uniform${N} ${OUT}/amr${N} 0.25
                  RESULT_VARIABLE rc)
endif()
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "the refined I-front is off (ifront_full.py --check)")
endif()
