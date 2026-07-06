"""Field-mode HydroState gravity ghost fill — the incompressible branch.

Regression for ``cell_boundary._calculate_ghost_values`` indexing the
field-mode ``HydroState.rhoY0`` / ``p20`` with only the vertical component
``nimage[y_axs]``. That is correct for the 1D PROFILE-mode hydrostates
(vertical profiles broadcast horizontally on demand), but WRONG for the
full grid-shaped FIELD-mode hydrostates a terrain / sphere run carries: a
scalar radial index then slices the wrong (longitude) axis and mis-shapes
the result, which either raises on the ``rho = rhoY * S`` broadcast or —
if the axis lengths happen to match — fills the ghosts with silent garbage.

The incompressible branch runs only during ``do_initial_projection`` (which
freezes the regime to incompressible), so no prior case hit it: the SWE
sphere cases that project have ``grav = 0``, so the gravity ghost fill never
runs. Hughes & Jablonowski's adjusted ridge balance (pt 2) may need the
projection.

Build the H&J compressible field-mode shell, freeze it incompressible, run
the gravity boundary fill, and assert:

* the fill stays finite, and
* it is CONSISTENT — the incompressible branch defines the ghost ``rhoY`` as
  the field-mode reference ``HydroState.rhoY0`` at that ghost cell, so the
  two must match exactly (the buggy single-axis index does not).

The ATMOSPHERIC_EXTENSION ``p20`` sibling (same file, same one-axis-index
pattern) is fixed by the same helper but has no field-mode case to gate: the
only ATMOSPHERIC_EXTENSION cases (Lamb waves) are profile-mode.
"""

import logging

import numpy as np

logging.disable(logging.INFO)

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.tests import test_hj_baroclinic as hj
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState


def _build_shell(nx, ny, nz):
    """H&J compressible field-mode spherical shell (no projection at build)."""
    udo = hj.UserData()
    udo.inx, udo.iny, udo.inz = nx + 1, ny + 1, nz + 1
    udo.initial_projection = False
    udo.diag = False
    udo.output_timesteps = False
    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = hj.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    return mem, ud


def test_field_mode_gravity_ghost_incompressible():
    # nx != ny on purpose: with the one-axis-index bug, rhoY0[radial_int]
    # slices the longitude axis (length nx) while S is a (nx, nz) radial
    # slice, so rho = rhoY * S raises a broadcast error outright.
    mem, ud = _build_shell(nx=16, ny=12, nz=32)
    assert mem.npf.HydroState.field_mode  # sphere shell -> field-mode hydrostates

    # freeze to incompressible exactly as do_initial_projection does, so the
    # gravity ghost fill takes the HydroState.rhoY0 branch
    ud.is_compressible = 0
    ud.compressibility = 0.0

    bdry_c.set_ghost_cells(mem, ud)

    for a in ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX"):
        assert np.all(np.isfinite(getattr(mem.sol, a))), a

    # Consistency: in the incompressible branch the gravity fill sets the
    # ghost rhoY to the field-mode reference rhoY0 at that ghost cell. Check
    # the radial (gravity-axis) ghost rows on the latitude interior (the
    # latitude-wall fill, applied after the radial fill, overwrites the phi
    # ghost corners).
    vaxis = 1  # gravity_direction
    igy = mem.elem.igs[vaxis]
    igz = mem.elem.igs[2]
    rhoY0 = mem.npf.HydroState.rhoY0
    for k in list(range(igy)) + [-(j + 1) for j in range(igy)]:  # bottom + top rows
        sl = [slice(None), slice(None), slice(igz, -igz)]
        sl[vaxis] = k
        sl = tuple(sl)
        assert np.allclose(mem.sol.rhoY[sl], rhoY0[sl]), ("row", k)
