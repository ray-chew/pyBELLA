"""Terrain elliptic-operator oracle: discrete div∘correction == operator.

The semi-implicit projection is consistent only if the assembled elliptic
operator (lap3D with the C_ij = wplus * (J A^T H^-1 A) tensor) is exactly
the composition of the discrete divergence with the discrete metric-
corrected momentum correction — mismatched staggering or a dropped 1/J
shows up here long before it corrupts a mountain-wave run.

Checks, on the smoke_agnesi configuration (quasi-2D x-y-vertical 3D):

1. flat baseline: the composition identity holds for the plain solver
   (guards the test's own assembly against phantom mismatches),
2. forced-flat metric == plain operator application (~1e-13),
3. terrain (witch-of-Agnesi hill): composition identity still holds.

Comparisons exclude a 3-node boundary window: ghost reconstruction inside
the lap3D kernel and production ghost filling are convention-equivalent
but not bit-equal at the box edge; interior nodes carry the full stencil.
"""

import numpy as np

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import terrain
from pybella.flow_solver.numerics import coriolis, implicit_euler
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.flow_solver.utils.boundary import node_boundary as bdry_n
from pybella.tests import smoke_agnesi
from pybella.utils import axes, user_data
from pybella.utils.operators import divergence
from pybella.utils.operators.laplacian import preconditioner


def _agnesi_hill(ud, h0_m=300.0, a_m=5000.0):
    h0 = h0_m / ud.h_ref
    a = a_m / ud.h_ref

    def h(xi1, xi2):
        return h0 * a**2 / (xi1**2 + a**2) + 0.0 * xi2

    return h


def _make_mem(orography=None):
    ud = user_data.UserDataInit(**vars(smoke_agnesi.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    # smoke_agnesi carries its own hill; the oracle controls terrain itself
    ud.orography = orography(ud) if orography is not None else None
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = smoke_agnesi.sol_init(sol, npf, elem, node, th, ud)
    from pybella.utils.data_structures import ModelState

    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


def _smooth_p_box(node):
    """Deterministic smooth pressure on the node.isc solve box."""
    x = node.x[1:-1].reshape(-1, 1, 1)
    y = node.y[1:-1].reshape(1, -1, 1)
    z = node.z[1:-1].reshape(1, 1, -1)
    Lx = node.x[-1] - node.x[0]
    return (
        np.sin(2 * np.pi * x / Lx) * np.cos(np.pi * y)
        + 0.1 * np.cos(2 * np.pi * x / Lx + 0.3) * y
        + 0.0 * z
    )


def _diag_inv_like_solver(mem, ud, dt):
    """Replicate the diag_inv construction of _prepare_3d_system."""
    hv = coriolis.compute_inverse_coefficients(mem, ud, dt)
    h_role = ((hv[0], hv[1], hv[2]), (hv[3], hv[4], hv[5]), (hv[6], hv[7], hv[8]))
    if mem.elem.metric is not None:
        h_role = terrain.elliptic_tensor(mem.elem.metric, h_role)
    rho_of = axes.role_of_axis(axes.vertical_axis(ud))
    cij = [
        [mem.npf.wplus[i] * h_role[rho_of[i]][rho_of[j]] for j in range(3)]
        for i in range(3)
    ]
    return preconditioner.prepare_diag(
        mem.npf, mem.node, cii=(cij[0][0], cij[1][1], cij[2][2])
    )


def _operator_and_composition(orography=None):
    mem, ud = _make_mem(orography)
    node = mem.node
    dt = float(ud.dtfixed)

    implicit_euler.operator_coefficients_nodes(mem, ud, dt)

    # operator side (captures coefficient copies, immune to later mutation)
    lap, _ = implicit_euler._prepare_linear_system(mem, ud, dt)
    p_box = _smooth_p_box(node)
    lhs = np.asarray(lap @ p_box.ravel()).reshape(node.isc)

    # composition side: divergence of the pure pressure correction
    diag_inv = _diag_inv_like_solver(mem, ud, dt)
    p_full = np.zeros(node.sc)
    p_full[node.i1] = p_box
    bdry_n.set_ghost_nodes(p_full, node, ud)

    mem.sol.rhou[...] = 0.0
    mem.sol.rhov[...] = 0.0
    mem.sol.rhow[...] = 0.0
    implicit_euler._correction_nodes(mem, ud, dt, p_full, 0)
    rhs = np.zeros(node.isc)
    divergence.compute_at_nodes(rhs, mem.elem, mem.sol, ud)

    comp = diag_inv * (-(1.0 / dt) * rhs + mem.npf.wcenter * p_box)
    return lhs, comp


def _window(shape):
    """3-node inset where the axis allows it (degenerate axes keep 1)."""
    return tuple(slice(3, -3) if n > 8 else slice(1, -1) for n in shape)


def _rel_err(lhs, comp):
    win = _window(lhs.shape)
    scale = np.max(np.abs(lhs[win]))
    return np.max(np.abs((lhs - comp)[win])) / scale


def test_composition_identity_flat_baseline():
    lhs, comp = _operator_and_composition()
    assert _rel_err(lhs, comp) <= 1e-12


def test_flat_metric_matches_plain_operator():
    flat = lambda ud: (lambda xi1, xi2: 0.0 * xi1 + 0.0 * xi2)
    lhs_plain, _ = _operator_and_composition()
    lhs_flat, _ = _operator_and_composition(flat)
    scale = np.max(np.abs(lhs_plain))
    assert np.max(np.abs(lhs_flat - lhs_plain)) / scale <= 1e-13


def test_composition_identity_with_terrain():
    lhs, comp = _operator_and_composition(_agnesi_hill)
    assert _rel_err(lhs, comp) <= 1e-12
