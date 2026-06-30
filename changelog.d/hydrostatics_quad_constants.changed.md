Named the terrain hydrostate quadrature grid constants in `physics/hydrostatics.py`
(`_HYDRO_QUAD_MIN_POINTS = 2048`, `_HYDRO_QUAD_POINTS_PER_CELL = 16`) with a comment,
replacing the bare `max(2048, 16 * int(elem.sc[vv]))` magic numbers. Values
unchanged; bit-identical at tol 0.
