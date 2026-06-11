Removed the unused `utils/operators/laplacian/lap2D_numba.py` (stencil-based
2D laplacian): the production 2D elliptic path uses `lap2D_manual.lap2D_gather`
exclusively (`@nb.njit`, via `implicit_euler._prepare_2d_system`), and nothing
imported the module. The 3D path's `lap3D` kernel likewise remains numba-jitted.
