Fixed the full-3D (`inz > 1`) implicit/elliptic solver path, broken since the
package restructure: re-wired `lap3D` to the interior-sized `npf` arrays
(missing imports, coefficient slicing, 3D preconditioner via
`preconditioner.prepare_diag`), fixed the flat-vector memory layout
(C-order `[x, y, z]`; the old reshape silently transposed x and z on
non-cubic grids), fixed the 3D `rhs` shape mismatch and a sign error on the
x-component of the 3D nodal divergence (legacy, dating to Oct 2021), and made
the pressure diagnostic kernel dimension-agnostic. The 3D path is validated
against the 2D solver on a y-uniform quasi-2D problem to ~1e-10
(`test_scripts/test_3d_elliptic_oracle.py`), and the
`test_travelling_vortex_3d_coriolis` regression case now runs with a
committed golden-master target.
