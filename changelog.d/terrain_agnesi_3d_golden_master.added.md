New golden-master regression case `test_agnesi_3d`: an isolated circular 3D
Agnesi bell (h0 = 100 m, a = 10 km, N = 0.01 1/s, U = 10 m/s, Gal-Chen,
Rayleigh sponge above 12 km) on a 64x32x64 grid, 10 spin-up steps. First
case where the orography depends on BOTH horizontal coordinates, locking the
G2 != 0 terrain dynamics (second-slope contravariant fluxes, pressure-map
column, elliptic cross terms, bottom BC) as a golden master — closes
terrain-following limitation #4. Not wired into CI (full-3D elliptic solves,
~5 min/run); run locally via `pybella -ic test_agnesi_3d -N 1`.
