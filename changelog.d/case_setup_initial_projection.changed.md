Factored the remaining copy-pasted regression-case idioms into
`tests/case_setup.py`: `do_initial_projection` (the ~35-line incompressible
initial-projection block, duplicated verbatim across the travelling-vortex, SWE
and 3D-Coriolis cases) and `mirror_centers` (the periodic nearest-image vortex
centre, duplicated 3×). The unused `T_from_p_rho` helper (defined in 5 cases,
called in none) and the local `class obj` attribute-bag (replaced by
`types.SimpleNamespace` inside the helper) were removed, along with the imports
they orphaned. Test-only; FAST/affected cases bit-identical at tol 0.
