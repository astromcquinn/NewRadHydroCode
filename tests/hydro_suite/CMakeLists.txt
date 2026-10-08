# Hydro accuracy test suite (docs/hydro_tests.md). Included from
# tests/CMakeLists.txt when the hydro/PM libraries are built.
add_executable(hydro_suite_driver suite_driver.cpp)
target_link_libraries(hydro_suite_driver PRIVATE radray_pm)
