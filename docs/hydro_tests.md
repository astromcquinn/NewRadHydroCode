# Hydro accuracy test suite

Physics tests of the C++ hydro solver (`src/hydro.cpp`, the bit-for-bit port
of `hydro_tools.f90 -DMF`, see `docs/port_notes_hydro.md`) against analytic
solutions and independent references. This is **not** a Fortran comparison:
it measures how accurate the moving-frame (MF) TVD scheme is.

RESULTS_PLACEHOLDER
