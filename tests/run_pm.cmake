# Runs pm_driver on one PM problem.  Variables: DRIVER PROBLEM OUT
file(REMOVE_RECURSE ${OUT})
file(MAKE_DIRECTORY ${OUT})
execute_process(COMMAND ${DRIVER} ${PROBLEM} ${OUT}
                OUTPUT_FILE ${OUT}.stdout
                ERROR_FILE ${OUT}.stderr
                RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  file(READ ${OUT}.stderr err)
  message(FATAL_ERROR "pm_driver failed (${rc}): ${err}")
endif()
