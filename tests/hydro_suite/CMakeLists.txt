# Hydro accuracy test suite (docs/hydro_tests.md). Included from
# tests/CMakeLists.txt when the hydro/PM libraries are built.
add_executable(hydro_suite_driver suite_driver.cpp)
target_link_libraries(hydro_suite_driver PRIVATE radray_pm)

# Physics accuracy tests (tests/hydro_suite/run_suite.py). The quick tier
# (label hydro_suite) runs small grids in a few minutes in total; the full
# tier (label nightly) reproduces every number of docs/hydro_tests.md
# (about 3 h per scheme on one 4-core machine). The full tier writes ~7 GB
# of runs per scheme, more than a CI runner's disk: --clean removes each
# test's runs once it is analysed (the results JSON stays).
if(Python3_Interpreter_FOUND)
  set(HS_RUN ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_suite.py)
  foreach(sec tubes galilean sedov sound pancake cosmobox kh)
    add_test(NAME hydro_suite_${sec}
             COMMAND ${HS_RUN} ${sec} --tier quick --check --jobs 2
                     --driver $<TARGET_FILE:hydro_suite_driver>
                     --rundir ${CMAKE_CURRENT_BINARY_DIR}/runs_quick)
    set_tests_properties(hydro_suite_${sec} PROPERTIES TIMEOUT 1800 LABELS hydro_suite)
    add_test(NAME hydro_suite_full_${sec}
             COMMAND ${HS_RUN} ${sec} --tier full --check --jobs 3
                     --driver $<TARGET_FILE:hydro_suite_driver>
                     --rundir ${CMAKE_CURRENT_BINARY_DIR}/runs_full --clean
                     --results ${CMAKE_CURRENT_BINARY_DIR}/results_full)
    set_tests_properties(hydro_suite_full_${sec} PROPERTIES TIMEOUT 14400 LABELS nightly)
    # the unsplit scheme (docs/hydro_unsplit.md), its own tolerances
    add_test(NAME hydro_suite_unsplit_${sec}
             COMMAND ${HS_RUN} ${sec} --tier quick --scheme unsplit --check --jobs 2
                     --driver $<TARGET_FILE:hydro_suite_driver>
                     --rundir ${CMAKE_CURRENT_BINARY_DIR}/runs_quick_unsplit)
    set_tests_properties(hydro_suite_unsplit_${sec} PROPERTIES TIMEOUT 1800 LABELS hydro_suite)
    add_test(NAME hydro_suite_full_unsplit_${sec}
             COMMAND ${HS_RUN} ${sec} --tier full --scheme unsplit --check --jobs 3
                     --driver $<TARGET_FILE:hydro_suite_driver>
                     --rundir ${CMAKE_CURRENT_BINARY_DIR}/runs_full_unsplit --clean
                     --results ${CMAKE_CURRENT_BINARY_DIR}/results_full)
    set_tests_properties(hydro_suite_full_unsplit_${sec} PROPERTIES TIMEOUT 14400 LABELS nightly)
  endforeach()
endif()
