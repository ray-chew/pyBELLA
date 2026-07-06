"""Stage-D gates: 3D compressible spherical shell with RADIAL gravity.

The resting isothermal atmosphere on the true-radius shell must stay at
rest to the solve floor. This exercises the whole radial-gravity chain
at once: `analytical_state` along radial columns (height/h_v general
branch), the general gravity ghost fill (e_up reflection), the e_up
buoyancy source in both explicit parts, the general H^-1 kernel with
nu != 0, and the elliptic projection with the spherical metric. Run at
a = 20 AND a = 5 scale heights (the latter is a strongly curved
small-planet stress: metric variation across the shell is O(1)).

Measured floors (32x8x16, 5 steps): |u| ~ 2e-12, |p2| ~ 2e-11; gates
carry ~100x headroom.
"""

import logging

import numpy as np
import pytest

logging.disable(logging.INFO)

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import spherical, time_update
from pybella.flow_solver.physics import hydrostatics, thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.tests import smoke_agnesi
from pybella.utils import options as opts
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState
from pybella.utils.io.debug import NullDebugWriter


def _spherical_hill(a_nd, h0_frac=0.06, sigma=0.5):
    """2pi-periodic hill h(lambda, phi) with analytic gradients."""
    h0 = h0_frac  # nondim (h_ref units)

    def h(lam, phi):
        return h0 * (1.0 + np.cos(lam)) / 2.0 * np.exp(-((phi / sigma) ** 2))

    def dh_dlam(lam, phi):
        return -h0 / 2.0 * np.sin(lam) * np.exp(-((phi / sigma) ** 2))

    def dh_dphi(lam, phi):
        return (
            h0
            * (1.0 + np.cos(lam))
            / 2.0
            * (-2.0 * phi / sigma**2)
            * np.exp(-((phi / sigma) ** 2))
        )

    return h, (dh_dlam, dh_dphi)


def _resting_shell(a_nd, nsteps, terrain_hill=False):
    base = smoke_agnesi.UserData()
    base.orography = None
    base.xmin, base.xmax = -np.pi, np.pi
    depth = 10000.0 / base.h_ref
    if terrain_hill:
        h, grad = _spherical_hill(a_nd)
        base.curvilinear_map = spherical.SphericalTerrainMap(a_nd, depth, h, grad)
        base.ymin, base.ymax = 0.0, depth  # radial axis carries eta
    else:
        base.curvilinear_map = spherical.SphericalShellMap(a_nd, frozen_radius=False)
        base.ymin, base.ymax = a_nd, a_nd + depth
    base.zmin, base.zmax = -1.2, 1.2
    base.bdry_type[0] = opts.BdryType.PERIODIC
    base.bdry_type[1] = opts.BdryType.WALL
    base.bdry_type[2] = opts.BdryType.WALL
    base.inx, base.iny, base.inz = 32 + 1, 8 + 1, 16 + 1
    base.stepmax = nsteps
    base.u_wind_speed = 0.0
    ud = user_data.UserDataInit(**vars(base))
    ud.coriolis_strength = np.array(ud.coriolis_strength)

    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)

    hydrostatics.analytical_state(npf, elem, node, th, ud)
    S0c = npf.HydroState.get_S0c(elem)
    sol.rhoY[...] = npf.HydroState.rhoY0
    sol.rho[...] = npf.HydroState.rhoY0 * S0c
    sol.rhou[...] = 0.0
    sol.rhov[...] = 0.0
    sol.rhow[...] = 0.0
    sol.rhoX[...] = sol.rho * (sol.rho / sol.rhoY - S0c)
    npf.p2_nodes[...] = 0.0
    ud.nonhydrostasy = 1.0
    ud.compressibility = 1.0

    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    inner = (slice(2, -2), slice(2, -2), slice(2, -2))
    rho0 = sol.rho[inner].copy()
    mem = time_update.do(
        mem, ud, tout=1.0e9, bld=None, writer=None, debug_writer=NullDebugWriter()
    )
    assert mem.time.step == nsteps
    return mem, inner, rho0


@pytest.mark.parametrize("a_nd", [20.0, 5.0], ids=["a20", "a5-small-planet"])
def test_resting_shell_with_radial_gravity(a_nd):
    mem, inner, rho0 = _resting_shell(a_nd, nsteps=10)
    moms = (mem.sol.rhou, mem.sol.rhov, mem.sol.rhow)
    umax = max(np.max(np.abs(m[inner] / mem.sol.rho[inner])) for m in moms)
    drho = np.max(np.abs(mem.sol.rho[inner] - rho0)) / np.max(np.abs(rho0))
    p2max = np.max(np.abs(mem.npf.p2_nodes))
    assert umax < 1.0e-8, umax
    assert drho < 1.0e-10, drho
    assert p2max < 1.0e-7, p2max


def test_resting_atmosphere_over_spherical_hill():
    """Stage E: the terrain-following spherical shell (Gal-Chen in r over
    a Gaussian-belt hill) preserves the resting isothermal atmosphere —
    the balanced-init contract of the planar terrain machinery, now on
    curved geometry with radial gravity."""
    mem, inner, rho0 = _resting_shell(20.0, nsteps=10, terrain_hill=True)
    moms = (mem.sol.rhou, mem.sol.rhov, mem.sol.rhow)
    umax = max(np.max(np.abs(m[inner] / mem.sol.rho[inner])) for m in moms)
    drho = np.max(np.abs(mem.sol.rho[inner] - rho0)) / np.max(np.abs(rho0))
    assert umax < 1.0e-7, umax
    assert drho < 1.0e-9, drho
