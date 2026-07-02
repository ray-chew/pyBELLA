Removed dead code from the data-assimilation layer ahead of the ModelState
repair: `localisation.py` (broken absolute import, zero callers), the legacy
`utils.ensemble` wrapper class, `set_p2_nodes`/`set_rhoY_cells`,
`params.converter`, and the `HSprojector_3t2D/2t3D` horizontal-slice
projectors (the supported 2D x-y cases are native 2D since the ModelState
refactor, so the projector guard could never fire).
