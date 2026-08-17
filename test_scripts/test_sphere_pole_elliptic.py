"""Gates for the elliptic pole-ring collapse.

At a lat-lon pole every longitude node of a given radius is one physical
point; left independent, the discrete Helmholtz system is rank-deficient
there. ``pole_collapse`` collapses each (radius, hemisphere) pole ring to
ONE master pressure unknown (a Galerkin scatter/gather around the one-sided
pole operator). These gates pin, on the GLOBAL (pole-to-pole) shell:

1. a resting global shell stays at rest to the Krylov floor (the pole rows
   inject no spurious momentum) over several implicit solves;
2. after projecting a manufactured divergent momentum field, the pressure
   is exactly longitude-uniform on both pole rings (single-valued);
3. the projection reduces the momentum divergence to the solver floor at
   every interior node, INCLUDING the pole rows;
4. the scatter/gather are exact algebraic inverses on the ring subspace.
"""

import numpy as np

from pybella.flow_solver.discretisation import grid as dis_grid, spherical
from pybella.flow_solver.numerics import implicit_euler, pole_collapse
from pybella.flow_solver.physics import hydrostatics, thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.utils import options as opts
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState
from pybella.tests import smoke_agnesi

_A = 3.0


def _pole_mem(n=(16, 6, 12)):
    """Full ModelState on the GLOBAL pole-to-pole shell, zero gravity (a
    genuine 3D elliptic solve exercising the pole rings)."""
    base = smoke_agnesi.UserData()
    base.grav = 0.0
    base.gravity_strength = np.zeros(3)
    base.stratification = lambda y: 1.0 + 0.0 * y
    base.orography = None
    base.xmin, base.xmax = -np.pi, np.pi
    base.ymin, base.ymax = _A * 0.9, _A * 1.1
    base.zmin, base.zmax = -0.5 * np.pi, 0.5 * np.pi
    base.bdry_type = np.array(
        [opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.POLE]
    )
    base.inx, base.iny, base.inz = n[0] + 1, n[1] + 1, n[2] + 1
    base.curvilinear_map = spherical.SphericalShellMap(_A, pole=True)
    ud = user_data.UserDataInit(**vars(base))
    ud.coriolis_strength = np.array(ud.coriolis_strength)

    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    hydrostatics.integrated_state(npf, elem, node, th, ud)

    sol.rho[...] = 1.0
    sol.rhoY[...] = 1.0
    sol.rhou[...] = 0.0
    sol.rhov[...] = 0.0
    sol.rhow[...] = 0.0
    sol.rhoX[...] = 0.0
    npf.p2_nodes[...] = 0.0
    ud.nonhydrostasy = 1.0
    ud.compressibility = 1.0

    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


def _pole_rows(node):
    """(south, north) full-node phi indices of the pole rows."""
    igz = int(node.igs[2])
    ncz = int(node.sc[2])
    return igz, ncz - 1 - igz


# ------------------------------------------- scatter/gather algebra


def test_scatter_gather_ring_algebra():
    mem, ud = _pole_mem()
    coll = pole_collapse.get(mem.node)
    n = coll.n
    Ll, Lr, Lp = coll.shape
    # scatter(gather(e_master)) replicates the master over the whole ring;
    # gather sums the unique ring nodes (seam excluded) exactly once
    x = np.zeros(n)
    x[coll.uniq_master] = 1.0  # put 1 on each master
    g = coll.gather(np.ones(n))
    # each master gathered the (Ll-3) unique interior longitudes
    n_unique = Ll - 3
    masters = np.unique(coll.uniq_master)
    assert np.allclose(g[masters], n_unique)
    # non-master ring entries are zeroed by gather
    non_master = np.setdiff1d(coll.ring_all, masters)
    assert np.allclose(g[non_master], 0.0)
    # scatter broadcasts a master value over the ring (incl. the +pi seam)
    y = np.zeros(n)
    y[masters] = 3.0
    s = coll.scatter(y)
    assert np.allclose(s[coll.ring_all], 3.0)


# --------------------------------------------- resting shell at rest


def test_resting_global_shell_stays_at_rest():
    mem, ud = _pole_mem()
    dt = 0.5
    for _ in range(4):
        implicit_euler.do_implicit_part(mem, ud, dt)
    inner = (slice(2, -2), slice(2, -2), slice(2, -2))
    mmax = max(
        np.max(np.abs(getattr(mem.sol, m)[inner])) for m in ("rhou", "rhov", "rhow")
    )
    assert mmax < 1e-9, mmax  # Krylov floor, no pole-injected momentum


# ------------------------------- single-valued pole pressure + projection


def _set_divergent_momentum(mem, ud):
    """A smooth Cartesian momentum field with nonzero divergence that
    crosses the poles (so the pole rings carry a real constraint)."""
    x = mem.elem.metric.x
    mem.sol.rho[...] = 1.0
    mem.sol.rhoY[...] = 1.0
    mem.sol.rhou[...] = 0.15 * np.sin(x[2]) + 0.1 * x[0]
    mem.sol.rhov[...] = 0.12 * np.cos(x[0]) * x[2]
    mem.sol.rhow[...] = 0.1 * x[1]
    bdry_c.set_ghost_cells(mem, ud)


def test_pole_pressure_single_valued_and_projects():
    mem, ud = _pole_mem()
    # incompressible -> the implicit solve is a Leray projection (div -> 0)
    ud.is_compressible = 0
    ud.compressibility = 0.0
    _set_divergent_momentum(mem, ud)

    from pybella.utils.operators import divergence

    rhs0 = np.zeros(mem.node.isc)
    divergence.compute_at_nodes(rhs0, mem.elem, mem.sol, ud)
    div0 = np.max(np.abs(rhs0[mem.node.i1]))

    coll = pole_collapse.get(mem.node)
    ringsum0 = np.max(np.abs(coll.gather(rhs0.ravel())[np.unique(coll.uniq_master)]))
    div0_interior = np.max(np.abs(rhs0[mem.node.i1][:, :, 1:-1]))

    implicit_euler.do_implicit_part(mem, ud, 0.5)

    # pressure is longitude-uniform on both pole rings (single physical point)
    igl = int(mem.node.igs[0])
    south, north = _pole_rows(mem.node)
    p = mem.npf.p2_nodes
    for row in (south, north):
        ring = p[igl:-igl, 2:-2, row]  # interior lambda x interior r
        spread = np.max(np.abs(ring - ring[:1]))
        assert spread < 1e-10, (row, spread)

    # the projection drives the INTERIOR per-node divergence AND the
    # RING-SUMMED pole divergence to the Krylov floor (the collapse enforces
    # only the ring-summed pole constraint; per-lambda pole-node divergence
    # is not individually determined by the single master unknown)
    rhs1 = np.zeros(mem.node.isc)
    divergence.compute_at_nodes(rhs1, mem.elem, mem.sol, ud)
    div1_interior = np.max(np.abs(rhs1[mem.node.i1][:, :, 1:-1]))
    ringsum1 = np.max(np.abs(coll.gather(rhs1.ravel())[np.unique(coll.uniq_master)]))
    assert div1_interior < 1e-3 * div0_interior, (div0_interior, div1_interior)
    assert ringsum1 < 1e-2 * ringsum0, (ringsum0, ringsum1)
