"""Gates for the pole ghost exchange (cells + nodes).

The pole ghost fill is a PURE INDEX REMAP (lambda -> lambda + pi, phi
mirrored across the pole) with NO vector rotation, because momenta are
global Cartesian components. The theorem these gates pin: any GLOBALLY
smooth field sampled consistently at the grid's physical coordinates is a
FIXED POINT of the exchange — filling a ghost cell reproduces the far-side
interior value, which equals the field at the ghost's own physical point
to roundoff. If the exchange rotated a vector or shifted the wrong way,
the momenta would not come back.

Gates:
1. cell exchange restores an analytic Cartesian field (scalars AND a
   smooth vector field in the momenta) after the ghosts are clobbered;
2. same at a sweep-flipped metric orientation (axis bookkeeping);
3. node exchange restores an analytic scalar; the two pole-node rows are
   exactly lambda-uniform after the fill (single-valued at the pole);
4. the exchange is idempotent.
"""

import types

import numpy as np
import pytest

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import spherical
from pybella.flow_solver.utils import fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.flow_solver.utils.boundary import node_boundary as bdry_n
from pybella.utils import options as opts

_A = 2.0


class _PoleUD:
    def __init__(self, n=(32, 4, 16), a=_A):
        self.inx, self.iny, self.inz = n[0] + 1, n[1] + 1, n[2] + 1
        self.xmin, self.xmax = -np.pi, np.pi
        self.ymin, self.ymax = a * 0.95, a * 1.05
        self.zmin, self.zmax = -0.5 * np.pi, 0.5 * np.pi
        self.gravity_direction = 1
        self.bdry_type = np.array(
            [opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.POLE]
        )
        self.curvilinear_map = spherical.SphericalShellMap(a, pole=True)


def _scalar(x):
    return 2.0 + np.sin(x[0]) + 0.3 * x[1] + np.cos(2.0 * x[2])


def _vector(x):
    # a smooth Cartesian vector field (no particular symmetry) — the
    # momenta must copy component-by-component with no rotation
    return (
        np.sin(x[1]) * x[2],
        np.cos(x[0]) + 0.2 * x[2],
        x[0] * x[1],
    )


def _fill_analytic(sol, m):
    x = m.x
    sol.rho[...] = 1.2 + 0.2 * np.sin(x[0]) * np.cos(x[2])
    sol.rhoY[...] = 1.0 + 0.1 * np.cos(x[1])
    sol.rhoX[...] = _scalar(x)
    v = _vector(x)
    sol.rhou[...], sol.rhov[...], sol.rhow[...] = v[0], v[1], v[2]


def _build(n=(32, 4, 16)):
    ud = _PoleUD(n=n)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    _fill_analytic(sol, elem.metric)
    mem = types.SimpleNamespace(elem=elem, node=node, sol=sol)
    return ud, mem, elem, node, sol


# --------------------------------------------------------- cell exchange


def test_cell_pole_exchange_restores_analytic():
    ud, mem, elem, node, sol = _build()
    m = elem.metric
    ig = int(elem.igs[2])
    ncz = sol.rho.shape[2]
    ref = {
        n: getattr(sol, n).copy()
        for n in ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX")
    }

    # clobber the phi ghost slabs, then fill
    for name in ref:
        arr = getattr(sol, name)
        arr[:, :, :ig] = 1.0e9
        arr[:, :, ncz - ig :] = -1.0e9

    handler = bdry_c.CellBoundaryHandler(mem, ud)
    handler._apply_pole_boundary(sol)

    # phi ghosts, at interior lambda/r, must match the analytic reference
    ii = (slice(ig, -ig), slice(1, -1))
    scale = 10.0
    for name in ref:
        got = getattr(sol, name)
        for slab in (slice(0, ig), slice(ncz - ig, ncz)):
            a = got[ii[0], ii[1], slab]
            b = ref[name][ii[0], ii[1], slab]
            np.testing.assert_allclose(a, b, rtol=1e-11, atol=1e-11 * scale)


def test_cell_pole_exchange_is_idempotent():
    ud, mem, elem, node, sol = _build()
    handler = bdry_c.CellBoundaryHandler(mem, ud)
    handler._apply_pole_boundary(sol)
    snap = {
        n: getattr(sol, n).copy()
        for n in ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX")
    }
    handler._apply_pole_boundary(sol)
    for name, arr in snap.items():
        np.testing.assert_array_equal(getattr(sol, name), arr)


def test_cell_pole_exchange_sweep_orientation():
    """The fill must use the metric's tracked lambda/phi axes: run it in a
    flipped orientation and it must agree with the canonical fill."""
    ud, mem, elem, node, sol = _build()
    handler = bdry_c.CellBoundaryHandler(mem, ud)
    handler._apply_pole_boundary(sol)
    ref = {
        n: getattr(sol, n).copy()
        for n in ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX")
    }

    # rebuild, clobber phi ghosts, fill in a flipped orientation
    ud, mem, elem, node, sol = _build()
    ig = int(elem.igs[2])
    ncz = sol.rho.shape[2]
    for name in ref:
        arr = getattr(sol, name)
        arr[:, :, :ig] = 5.0e8
        arr[:, :, ncz - ig :] = 5.0e8
    sol.flip_forward()
    elem.metric.flip_forward()
    handler = bdry_c.CellBoundaryHandler(mem, ud)
    handler._apply_pole_boundary(sol)
    sol.flip_backward()
    elem.metric.flip_backward()

    ii = (slice(int(elem.igs[0]), -int(elem.igs[0])), slice(1, -1))
    for name in ref:
        got = getattr(sol, name)
        for slab in (slice(0, ig), slice(ncz - ig, ncz)):
            np.testing.assert_allclose(
                got[ii[0], ii[1], slab],
                ref[name][ii[0], ii[1], slab],
                rtol=1e-11,
                atol=1e-9,
            )


# --------------------------------------------------------- node exchange


def test_node_pole_exchange_restores_and_uniform_poles():
    ud, mem, elem, node, sol = _build()
    mn = node.metric
    p = np.ascontiguousarray(_scalar(mn.x))
    ref = p.copy()
    ig = int(node.igs[2])
    ncz = p.shape[2]

    # clobber phi ghost node rows
    p[:, :, :ig] = 7.0e8
    p[:, :, ncz - ig :] = 7.0e8

    bdry_n._apply_pole_nodes(p, node, dim=2)

    # ghost node rows restored (interior lambda/r)
    ii = (slice(int(node.igs[0]), -int(node.igs[0])), slice(1, -1))
    for slab in (slice(0, ig), slice(ncz - ig, ncz)):
        np.testing.assert_allclose(
            p[ii[0], ii[1], slab], ref[ii[0], ii[1], slab], rtol=1e-10, atol=1e-9
        )

    # the two pole node rows are exactly lambda-uniform (one physical point)
    for pole_row in (ig, ncz - 1 - ig):
        row = p[:, :, pole_row]  # (lambda, r)
        spread = np.max(np.abs(row - row[:1, :]), axis=0)
        assert np.max(spread) < 1e-12, (pole_row, np.max(spread))


def test_node_pole_exchange_full_set_ghost_nodes():
    """Driving the public set_ghost_nodes (lambda-periodic then pole)
    leaves the pole rows uniform and the ghost rows analytic."""
    ud, mem, elem, node, sol = _build()
    mn = node.metric
    p = np.ascontiguousarray(_scalar(mn.x))
    ref = p.copy()
    ig = int(node.igs[2])
    ncz = p.shape[2]
    p[:, :, :ig] = 0.0
    p[:, :, ncz - ig :] = 0.0

    bdry_n.set_ghost_nodes(p, node, ud)

    # interior lambda/r only: the phi x r corners are r-reflected before the
    # pole fill (set_ghost_nodes does r then phi), so they legitimately
    # differ from the analytic-at-ghost reference there
    igr = int(node.igs[1])
    ii = (slice(int(node.igs[0]), -int(node.igs[0])), slice(igr, -igr))
    for slab in (slice(0, ig), slice(ncz - ig, ncz)):
        np.testing.assert_allclose(
            p[ii[0], ii[1], slab], ref[ii[0], ii[1], slab], rtol=1e-10, atol=1e-9
        )
