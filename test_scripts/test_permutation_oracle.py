"""Permutation oracle: the endgame proof of axial agnosticity.

The 2D internal-long-wave reference (x-y plane, gravity on axis 1,
stratified, walls, full Coriolis) is embedded as a quasi-2D 3D state with
its axes relabeled by a CYCLIC permutation sigma (ref axis i -> twin axis
sigma[i]) and the vertical placed accordingly:

- T2: sigma = (1, 2, 0), gravity_direction = 2  (z-vertical convention)
- T0: sigma = (2, 0, 1), gravity_direction = 0  (x-vertical)

Only cyclic (even) permutations are physical relabelings: the rotation
vector is a pseudovector, so its components permute WITHOUT sign flips —
and indeed `compute_coriolis_strength` produces exactly the sigma-mapped
Omega from `omega` and `gravity_direction` alone, with no hand-tuning.

Both twins are stepped K times through the full solver (advection +
explicit + implicit/elliptic) and must reproduce the sigma-mapped 2D
reference fields plus exact uniformity along the degenerate axis. This
exercises, off the y-axis: gravity reads, buoyancy on MOMENTA[v],
hydrostatic profiles along v, wall zeroing / gravity ghost cells /
quasi-2D broadcasts on arbitrary axes, role-bound Coriolis (explicit +
H^-1), the full-tensor 3D elliptic operator, and the advection sweeps.
"""

import numpy as np
import pytest

from pybella.utils import user_data, data_structures, axes
from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import time_update as dis_time_update
from pybella.flow_solver.utils import fields, cache
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.flow_solver.physics import thermodynamics as gd_thermodynamics
from pybella.tests.test_internal_long_wave import UserData as ILW, sol_init

K_STEPS = 8
TOL = 1e-6

PROFILE_ATTRS = ("p0", "p20", "rho0", "S0", "S10", "pi0", "rhoY0", "Y0")


class _StubWriter:
    def write(self, *a, **k):
        pass

    def populate(self, *a, **k):
        pass


def _build_ref():
    ud = user_data.UserDataInit(**vars(ILW()))
    ud.diag = False
    ud.stepmax = K_STEPS

    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = gd_thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = sol_init(sol, npf, elem, node, th, ud)

    mem = data_structures.ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


def _embed(ref2d, sigma, dummy_axis, dummy_count):
    """ref (a, b) array -> twin 3D array, broadcast along the dummy axis."""
    arr = ref2d[:, :, None]  # ref axes (0, 1, dummy=2)
    arr = axes.permute_axes(arr, sigma)
    return np.repeat(arr, dummy_count, axis=dummy_axis)


def _build_twin(ref_mem, ud_ref, sigma, v):
    d = vars(ILW())
    ud = user_data.UserDataInit(**d)
    ud.diag = False
    ud.stepmax = K_STEPS
    ud.gravity_direction = v  # retriggers gravity/coriolis computation

    # relabel grid extents and counts: ref axis i -> twin axis sigma[i];
    # the former (collapsed) ref z-axis becomes the twin's degenerate axis
    ins = [None] * 3
    mins = [None] * 3
    maxs = [None] * 3
    bdry = [None] * 3
    ref_ins = (ud_ref.inx, ud_ref.iny, 2)
    ref_mins = (ud_ref.xmin, ud_ref.ymin, 0.0)
    ref_maxs = (ud_ref.xmax, ud_ref.ymax, 1.0)
    for i in range(3):
        ins[sigma[i]] = ref_ins[i]
        mins[sigma[i]] = ref_mins[i]
        maxs[sigma[i]] = ref_maxs[i]
        bdry[sigma[i]] = ud_ref.bdry_type[i]
    ud.inx, ud.iny, ud.inz = ins
    ud.xmin, ud.ymin, ud.zmin = mins
    ud.xmax, ud.ymax, ud.zmax = maxs
    for i in range(3):
        ud.bdry_type[i] = bdry[i]

    # the pseudovector check: with omega + gravity_direction alone, the
    # computed Coriolis vector must equal the sigma-mapped reference one
    expected = np.zeros(3)
    for i in range(3):
        expected[sigma[i]] = ud_ref.coriolis_strength[i]
    assert np.array_equal(ud.coriolis_strength, expected), (
        ud.coriolis_strength,
        expected,
    )

    ud.nonhydrostasy = ud_ref.nonhydrostasy
    ud.compressibility = ud_ref.compressibility

    elem, node = dis_grid.grid_init(ud)
    axes.validate(ud, elem.ndim)
    sol = fields.CellSolField(elem.sc)
    th = gd_thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)

    dummy = sigma[2]
    nc_dummy = elem.sc[dummy]
    nn_dummy = node.sc[dummy]

    # cell fields: embed sigma-mapped; momenta map componentwise by sigma
    for fld in ("rho", "rhoY", "rhoX"):
        getattr(sol, fld)[...] = _embed(
            getattr(ref_mem.sol, fld), sigma, dummy, nc_dummy
        )
    for i in range(3):
        twin_name = axes.MOMENTA[sigma[i]]
        getattr(sol, twin_name)[...] = _embed(
            getattr(ref_mem.sol, axes.MOMENTA[i]), sigma, dummy, nc_dummy
        )

    # node / cell pressure fields
    npf.p2_nodes[...] = _embed(ref_mem.npf.p2_nodes, sigma, dummy, nn_dummy)
    npf.p2_cells[...] = _embed(ref_mem.npf.p2_cells, sigma, dummy, nc_dummy)

    # hydrostatic profiles: 1D along the vertical, lengths match by design
    for attr in PROFILE_ATTRS:
        getattr(npf.HydroState, attr)[...] = getattr(ref_mem.npf.HydroState, attr)
        getattr(npf.HydroState_n, attr)[...] = getattr(ref_mem.npf.HydroState_n, attr)

    mem = data_structures.ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


def _run(mem, ud):
    return dis_time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter())


@pytest.mark.parametrize(
    "sigma,v",
    [((1, 2, 0), 2), ((2, 0, 1), 0)],
    ids=["T2-z-vertical", "T0-x-vertical"],
)
def test_permutation_twin(sigma, v):
    ref_mem, ud_ref = _build_ref()
    twin_mem, ud_twin = _build_twin(ref_mem, ud_ref, sigma, v)

    ref_mem = _run(ref_mem, ud_ref)
    twin_mem = _run(twin_mem, ud_twin)
    assert abs(ref_mem.time.t - twin_mem.time.t) == 0.0

    dummy = sigma[2]
    # interior slices; pick the first interior index on the dummy axis
    i2_ref = (slice(2, -2), slice(2, -2))

    def twin_slice(arr3, node=False):
        idx = [slice(2, -2)] * 3
        idx[dummy] = 2
        sl = arr3[tuple(idx)]
        # remaining two axes are (ref0, ref1) in sigma order; sort them back
        a, b = [sigma[i] for i in (0, 1)]
        return sl if a < b else sl.T

    failures = []

    def check(name, twin3, ref2):
        diff = np.max(np.abs(twin_slice(twin3) - ref2[i2_ref]))
        if not diff < TOL:
            failures.append(f"{name}: max|diff| = {diff:.3e}")

    check("rho", twin_mem.sol.rho, ref_mem.sol.rho)
    check("rhoY", twin_mem.sol.rhoY, ref_mem.sol.rhoY)
    check("rhoX", twin_mem.sol.rhoX, ref_mem.sol.rhoX)
    for i in range(3):
        check(
            f"momentum {axes.MOMENTA[i]} -> {axes.MOMENTA[sigma[i]]}",
            getattr(twin_mem.sol, axes.MOMENTA[sigma[i]]),
            getattr(ref_mem.sol, axes.MOMENTA[i]),
        )
    check("p2_nodes", twin_mem.npf.p2_nodes, ref_mem.npf.p2_nodes)

    # exact uniformity along the degenerate axis (interior layers)
    lo = [slice(None)] * 3
    hi = [slice(None)] * 3
    lo[dummy], hi[dummy] = 2, 3
    for fld in ("rho", "rhou", "rhov", "rhow"):
        arr = getattr(twin_mem.sol, fld)
        uni = np.max(np.abs(arr[tuple(lo)] - arr[tuple(hi)]))
        if uni != 0.0:
            failures.append(f"uniformity {fld}: {uni:.3e}")

    assert not failures, "; ".join(failures)
