Fixed the incompressible + field-mode gravity ghost fill in
`cell_boundary._calculate_ghost_values`: it indexed the field-mode
`HydroState.rhoY0` (and, in the `ATMOSPHERIC_EXTENSION` branch, `p20`) with
only the vertical component `nimage[y_axs]`. That is correct for the 1D
PROFILE-mode hydrostates (vertical profiles broadcast on demand) but wrong
for the full grid-shaped FIELD-mode hydrostates a terrain / sphere run
carries — a scalar radial index then slices the wrong (longitude) axis and
mis-shapes the result, raising on the `rho = rhoY * S` broadcast (or filling
ghosts with garbage where the axes happen to match). The branch runs only
during `do_initial_projection` (which freezes the regime to incompressible),
so no prior case hit it: the SWE sphere cases that project have `grav = 0`,
so the gravity ghost fill never runs. Now indexes the whole ghost slice
tuple when the hydrostate is in field mode (via a `_hydro_at` helper), like
the sibling `metric.height[idx]` / `sol.rhoY[idx]` branches. Regression:
`test_scripts/test_field_mode_gravity_ghost.py`. This unblocks the
TC2-style initial projection on the deep compressible shell (Hughes &
Jablonowski pt 2).
