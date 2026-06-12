"""Native-2D terrain elliptic oracle: discrete div∘correction == operator.

2D sibling of ``test_terrain_elliptic_oracle.py``: the lap2D gather kernel
with the terrain-folded coefficient arrays (C = wplus ⊙ J A^T H^-1 A, the
cross terms riding the pre-existing cxy/cyx slots) must be exactly the
composition of the 2D metric divergence with the metric-corrected momentum
correction.

Checks, on the smoke_agnesi configuration collapsed to native 2D (inz = 1):

1. flat baseline: the composition identity holds for the plain 2D solver
   (the lap2D family was never proven this way — this gates the phase),
2. the same with out-of-plane Coriolis (exercises the legacy cxy/cyx path),
3. forced-flat metric == plain lap2D application (bypass contract),
4. terrain (witch-of-Agnesi hill): composition identity, with and without
   Coriolis (terrain cross terms + H^-1 off-diagonals folded together).

The 2D solve vector is the node.i2 interior box (transposed C-ravel); the
operator output lives on the rhs[node.i1] box of the convolution-shaped
node arrays. Comparisons exclude a 3-node boundary window as in 3D.
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
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState
from pybella.utils.operators import divergence
from pybella.utils.operators.laplacian import preconditioner


def _agnesi_hill(ud, h0_m=300.0, a_m=5000.0):
    h0 = h0_m / ud.h_ref
    a = a_m / ud.h_ref

    def h(xi1, xi2):
        return h0 * a**2 / (xi1**2 + a**2) + 0.0 * xi2

    return h


class _StretchedHillMap2D(terrain.CurvilinearMap):
    """Native-2D Tier-2 map: periodic x-stretch + periodic hill, Gal-Chen z."""

    def __init__(self, ud, h0_m=300.0, ax_rel=0.05):
        self.L = ud.xmax - ud.xmin
        self.h0 = h0_m / ud.h_ref
        self.ax = ax_rel * self.L
        self.eta0, self.etat = ud.ymin, ud.ymax

    def _h(self, xi1):
        return self.h0 * np.cos(np.pi * xi1 / self.L) ** 2

    def _dh(self, xi1):
        return -self.h0 * (np.pi / self.L) * np.sin(2.0 * np.pi * xi1 / self.L)

    def _decay(self, eta):
        return (self.etat - eta) / (self.etat - self.eta0)

    def coordinates(self, xi):
        xi1, eta = xi
        x = xi1 + self.ax * np.sin(2.0 * np.pi * xi1 / self.L)
        z = eta + self._h(xi1) * self._decay(eta)
        return [x + 0.0 * eta, z]

    def tangents(self, xi):
        xi1, eta = xi
        xp = 1.0 + self.ax * (2.0 * np.pi / self.L) * np.cos(2.0 * np.pi * xi1 / self.L)
        J = (1.0 - self._h(xi1) / (self.etat - self.eta0)) + 0.0 * eta
        G1 = self._dh(xi1) * self._decay(eta)
        zero = 0.0 * (J + xp)
        return [[xp + zero, G1 + zero], [zero, J + zero]]


def _make_mem(orography=None, coriolis_z=0.0, sleve=False, cmap=None):
    ud = user_data.UserDataInit(**vars(smoke_agnesi.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    ud.coriolis_strength[2] = coriolis_z
    ud.inz = 1  # collapse the degenerate axis: native 2D grid
    ud.orography = orography(ud) if orography is not None else None
    if sleve:
        hill = ud.orography
        ud.orography_smooth = lambda xi1, xi2: 0.5 * hill(xi1, xi2)
        ud.vertical_transform = terrain.SLEVETransform(
            s1=6000.0 / ud.h_ref, s2=1500.0 / ud.h_ref
        )
    elem, node = dis_grid.grid_init(ud)
    assert elem.ndim == 2
    if cmap is not None:
        # general path: override with the curvilinear-map metric before any
        # consumer (hydrostates, boundary) is built
        m = cmap(ud)
        elem.metric = terrain.build_metric_fields_from_map(elem, ud, m)
        node.metric = terrain.build_metric_fields_from_map(node, ud, m)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = smoke_agnesi.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


def _smooth_p_box(node):
    """Deterministic smooth pressure on the node.i2 solve box."""
    x = node.x[node.igx : -node.igx].reshape(-1, 1)
    y = node.y[node.igy : -node.igy].reshape(1, -1)
    Lx = node.x[-1] - node.x[0]
    return (
        np.sin(2 * np.pi * x / Lx) * np.cos(np.pi * y)
        + 0.1 * np.cos(2 * np.pi * x / Lx + 0.3) * y
    )


def _diag_inv_like_solver(mem, ud, dt):
    """Replicate the diag_inv construction of _prepare_2d_system.

    Geometric factors only — the legacy 2D preconditioner keeps H^-1 out
    of the diagonal, and the terrain branch preserves that.
    """
    del ud, dt
    if mem.elem.metric is not None:
        geo = terrain.elliptic_diag_geometric(mem.elem.metric)
        return preconditioner.prepare_diag(
            mem.npf,
            mem.node,
            cii=(mem.npf.wplus[0] * geo[0], mem.npf.wplus[1] * geo[1], None),
        )
    return preconditioner.prepare_diag(mem.npf, mem.node)


def _operator_and_composition(orography=None, coriolis_z=0.0, sleve=False, cmap=None):
    mem, ud = _make_mem(orography, coriolis_z, sleve, cmap)
    node = mem.node
    dt = float(ud.dtfixed)

    implicit_euler.operator_coefficients_nodes(mem, ud, dt)

    # operator side (the gather kernel ravels coefficient copies at build
    # time, so the linop is immune to later mutation)
    lap, _ = implicit_euler._prepare_linear_system(mem, ud, dt)
    p_box = _smooth_p_box(node)
    lhs = np.asarray(lap @ p_box.T.ravel()).reshape(node.iicy, node.iicx).T

    # composition side: divergence of the pure pressure correction
    diag_inv = _diag_inv_like_solver(mem, ud, dt)
    p_full = np.zeros((node.icx, node.icy))
    p_full[node.i2] = p_box
    bdry_n.set_ghost_nodes(p_full, node, ud)

    mem.sol.rhou[...] = 0.0
    mem.sol.rhov[...] = 0.0
    mem.sol.rhow[...] = 0.0
    implicit_euler._correction_nodes(mem, ud, dt, p_full, 0)
    rhs = np.zeros_like(mem.npf.rhs)
    divergence.compute_at_nodes(rhs, mem.elem, mem.sol, ud)

    comp = diag_inv[node.i1] * (
        -(1.0 / dt) * rhs[node.i1] + mem.npf.wcenter[node.i1] * p_box
    )
    return lhs, comp


def _rel_err(lhs, comp):
    win = (slice(3, -3), slice(3, -3))
    scale = np.max(np.abs(lhs[win]))
    return np.max(np.abs((lhs - comp)[win])) / scale


def test_composition_identity_flat_baseline():
    lhs, comp = _operator_and_composition()
    assert _rel_err(lhs, comp) <= 1e-12


def test_composition_identity_flat_with_coriolis():
    lhs, comp = _operator_and_composition(coriolis_z=0.2)
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


def test_composition_identity_terrain_and_coriolis():
    lhs, comp = _operator_and_composition(_agnesi_hill, coriolis_z=0.2)
    assert _rel_err(lhs, comp) <= 1e-12


def test_composition_identity_terrain_sleve():
    """Eta-dependent Jacobian through the native-2D elliptic assembly."""
    lhs, comp = _operator_and_composition(_agnesi_hill, sleve=True)
    assert _rel_err(lhs, comp) <= 1e-12


def test_composition_identity_general_map_2d():
    """Stretched Tier-2 map through the native-2D general N-fold, with and
    without out-of-plane Coriolis (cross terms + H^-1 off-diagonals)."""
    lhs, comp = _operator_and_composition(cmap=_StretchedHillMap2D)
    assert _rel_err(lhs, comp) <= 1e-12
    lhs_c, comp_c = _operator_and_composition(cmap=_StretchedHillMap2D, coriolis_z=0.2)
    assert _rel_err(lhs_c, comp_c) <= 1e-12
