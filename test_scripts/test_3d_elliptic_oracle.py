"""Oracle: the 3D elliptic path on a y-uniform quasi-2D problem must
reproduce the 2D solver's solution.

Builds twin states carrying the same vortex — 3D in the x-z plane on a
64x1x64 grid, 2D in the x-y plane on a 64x64 grid — with zero gravity and
zero Coriolis, runs implicit_euler.do_implicit_part (a pure projection,
compressibility = 0) on both, and compares the elliptic solve output
field-by-field. Pass criterion ~ solver tolerance (ud.tol = 1e-8), NOT 1e-2.

On an nx == nz grid a transposed solve is silent in shape terms, so this
test also asserts exact y-uniformity of the 3D result and that rhov stays
zero.

Mapping: 2D (x, y, u, v) <-> 3D (x, z, u, w).
"""

import numpy as np

from pybella.utils import user_data
from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.utils import fields, cache
from pybella.flow_solver.physics import thermodynamics as gd_thermodynamics
from pybella.flow_solver.physics import hydrostatics
from pybella.flow_solver.numerics import implicit_euler
from pybella.tests import test_travelling_vortex_3d_coriolis as tv3d

# ~ solver tolerance (1e-8) with headroom for preconditioner differences
# between the 2D and 3D paths
TOL = 1e-6


class obj:
    pass


def vortex_plane(a, b):
    """rho, swirl_a, swirl_b, on a 2D (a, b) meshgrid (cell centres)."""
    rho0, del_rho, R0, fac = 0.5, 0.5, 0.4, 1024.0
    ac, bc = 0.5, 0.5
    r = np.sqrt((a - ac) ** 2 + (b - bc) ** 2)
    uth = (fac * (1.0 - r / R0) ** 6 * (r / R0) ** 6) * (r < R0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ua = np.where(r > 0, uth * (-(b - bc) / r), 0.0)
        ub = np.where(r > 0, uth * (+(a - ac) / r), 0.0)
    rho = rho0 + del_rho * (1.0 - (r / R0) ** 2) ** 6 * (r < R0)
    return rho, ua, ub


def p2_plane(a, b, th):
    """Balanced nodal pressure (centrifugal part only; f = 0)."""
    R0, fac = 0.4, 1024.0
    ac, bc = 0.5, 0.5
    coe = np.array(
        [
            1.0 / 12.0,
            -12.0 / 13.0,
            9.0 / 2.0,
            -184.0 / 15.0,
            609.0 / 32.0,
            -222.0 / 17.0,
            -38.0 / 9.0,
            54.0 / 19.0,
            783.0 / 20.0,
            -558.0 / 7.0,
            1053.0 / 22.0,
            1014.0 / 23.0,
            -1473.0 / 16.0,
            204.0 / 5.0,
            510.0 / 13.0,
            -1564.0 / 27.0,
            153.0 / 8.0,
            450.0 / 29.0,
            -269.0 / 15.0,
            174.0 / 31.0,
            57.0 / 32.0,
            -74.0 / 33.0,
            15.0 / 17.0,
            -6.0 / 35.0,
            1.0 / 72.0,
        ]
    )
    r = np.sqrt((a - ac) ** 2 + (b - bc) ** 2)
    p2n = np.zeros_like(r)
    for ip in range(25):
        p2n += fac * (coe[ip] * ((r / R0) ** (12 + ip) - 1.0))
    p2n *= r / R0 < 1.0
    return th.Gamma * fac * p2n


def build_ud(iny, inz):
    d = vars(tv3d.UserData())
    d["iny"] = iny
    d["inz"] = inz
    ud = user_data.UserDataInit(**d)
    ud.coriolis_strength = np.array([0.0, 0.0, 0.0])
    ud.gravity_strength = np.zeros(3)
    ud.nonhydrostasy = 1.0
    ud.is_compressible = 0
    ud.compressibility = 0.0
    return ud


def build_state(ud, plane):
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = gd_thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    hydrostatics.integrated_state(npf, elem, node, th, ud)

    if plane == "xz":
        a = elem.x.reshape(-1, 1, 1)
        b = elem.z.reshape(1, 1, -1)
        rho, ua, ub = vortex_plane(a + 0 * b, b + 0 * a)
        sol.rho[...] = rho
        sol.rhou[...] = rho * ua
        sol.rhov[...] = 0.0
        sol.rhow[...] = rho * ub
        an = node.x[2:-2].reshape(-1, 1, 1)
        bn = node.z[2:-2].reshape(1, 1, -1)
        p2 = p2_plane(an + 0 * bn, bn + 0 * an, th)
        npf.p2_nodes[2:-2, 2:-2, 2:-2] = p2
    else:  # xy
        a = elem.x.reshape(-1, 1)
        b = elem.y.reshape(1, -1)
        rho, ua, ub = vortex_plane(a + 0 * b, b + 0 * a)
        sol.rho[...] = rho
        sol.rhou[...] = rho * ua
        sol.rhov[...] = rho * ub
        sol.rhow[...] = 0.0
        an = node.x[2:-2].reshape(-1, 1)
        bn = node.y[2:-2].reshape(1, -1)
        p2 = p2_plane(an + 0 * bn, bn + 0 * an, th)
        npf.p2_nodes[2:-2, 2:-2] = p2

    sol.rhoY[...] = 1.0
    sol.rhoX[...] = 0.0

    mem = obj()
    mem.sol = sol
    mem.npf = npf
    mem.elem = elem
    mem.node = node
    mem.th = th
    mem.cache = cache.FlowSolverCache()
    return mem


def test_3d_elliptic_path_matches_2d():
    dt = 0.01

    ud3 = build_ud(iny=2, inz=65)
    mem3 = build_state(ud3, "xz")
    implicit_euler.do_implicit_part(mem3, ud3, dt)

    ud2 = build_ud(iny=65, inz=1)
    mem2 = build_state(ud2, "xy")
    implicit_euler.do_implicit_part(mem2, ud2, dt)

    jc, jn = 2, 2  # interior y cell / node slice of the 3D state
    checks = {
        "p2_nodes (3D xz-slice vs 2D)": (
            mem3.npf.p2_nodes[:, jn, :],
            mem2.npf.p2_nodes,
        ),
        "rhou (3D xz-slice vs 2D)": (mem3.sol.rhou[:, jc, :], mem2.sol.rhou),
        "rhow vs rhov (3D vs 2D)": (mem3.sol.rhow[:, jc, :], mem2.sol.rhov),
        "p2_nodes y-uniformity (3D)": (
            mem3.npf.p2_nodes[:, 2, :],
            mem3.npf.p2_nodes[:, 3, :],
        ),
        "rhou y-uniformity (3D)": (mem3.sol.rhou[:, 2, :], mem3.sol.rhou[:, 3, :]),
        "rhov stays zero (3D)": (mem3.sol.rhov, np.zeros_like(mem3.sol.rhov)),
    }

    for name, (lhs, rhs) in checks.items():
        diff = np.max(np.abs(lhs - rhs))
        assert diff < TOL, f"{name}: max|diff| = {diff:.3e} >= {TOL:.0e}"
