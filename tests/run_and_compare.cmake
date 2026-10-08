# Runs rt_driver on one problem and compares with a Fortran reference dir.
# Variables: DRIVER PROBLEM OUT REF PYTHON COMPARE
file(REMOVE_RECURSE ${OUT})
file(MAKE_DIRECTORY ${OUT})
execute_process(COMMAND ${DRIVER} ${PROBLEM} ${OUT}
                OUTPUT_FILE ${OUT}.stdout
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "rt_driver failed (${rc}), see ${OUT}.stdout")
endif()
execute_process(COMMAND ${PYTHON} ${COMPARE} ${REF} ${OUT} RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "outputs differ from the Fortran reference")
endif()
