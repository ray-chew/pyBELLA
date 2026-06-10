"""Travelling vortex in the horizontal (x-z) plane with full Coriolis force.

Port of the legacy ``travelling_vortex_3D_Coriolis`` initial condition
(recoverable from the git tag ``archive/full_coriolis``,
``RKLM_Python/inputs/travelling_vortex_3D_Coriolis.py``) onto the current API.

Relative to the 2D travelling vortex (``test_travelling_vortex``), the legacy
case changes:

- the vortex lives in the x-z plane (u and w carry the swirl, v = 0) on a
  quasi-2D 64 x 1 x 64 grid over the unit cube,
- a nonzero rotation rate ``omega`` gives ``coriolis_strength = [100, 0, 100]``
  (omega * t_ref on the first and third components),
- the cyclostrophic pressure balance gains a Coriolis correction: a second
  polynomial (``ccoe``, exponents 7..25) scaled by f = coriolis_strength[0] is
  added to the centrifugal polynomial (``coe``, exponents 12..36),
- rho0 = del_rho = 0.5, no background wind, pseudo-incompressible regime,
  periodic boundaries in all three directions.

This case exercises the full-3D (``inz > 1``) implicit/elliptic solver path
(27-point Laplacian, ``utils/operators/laplacian/lap3D.py``) on a quasi-2D
grid; being y-uniform, its elliptic solve is validated against the 2D solver
to ~1e-10.
"""

import numpy as np

from ..utils import options as opts
from ..flow_solver.physics import hydrostatics
from ..flow_solver.numerics import implicit_euler
from ..flow_solver.utils import cache

from ..utils.data_structures import DiagnosticState


class UserData(object):
    grav = 0.0
    omega = 0.01 * 100.0  # => coriolis_strength[0] = [2] = omega * t_ref = 100.0

    h_ref = 10000.0
    t_ref = 100.0
    T_ref = 300.00
    p_ref = 1e5

    def __init__(self):
        self.h_ref = self.h_ref
        self.t_ref = self.t_ref
        self.T_ref = self.T_ref
        self.p_ref = self.p_ref
        self.grav = self.grav
        self.omega = self.omega

        self.xmin = 0.0
        self.xmax = 1.0
        self.ymin = 0.0
        self.ymax = 1.0
        self.zmin = 0.0
        self.zmax = 1.0

        self.u_wind_speed = 0.0
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = opts.BdryType.PERIODIC
        self.bdry_type[1] = opts.BdryType.PERIODIC
        self.bdry_type[2] = opts.BdryType.PERIODIC

        # legacy case runs in the pseudo-incompressible regime
        self.is_compressible = 0

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.95
        # legacy file used dtfixed = 2.1 * 1.200930e-2; a round 0.01 keeps the
        # regression run at exactly 100 steps (000..099) to tout = 1.0
        self.dtfixed = 0.01
        self.dtfixed0 = 0.01

        self.inx = 64 + 1
        self.iny = 1 + 1
        self.inz = 64 + 1

        self.initial_projection = True

        self.tout = [1.0]
        self.stepmax = 100

        self.stratification = self.stratification_function
        self.rhoe = self.rhoe_function
        self.output_timesteps = True

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_travelling_vortex_3d_coriolis"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.diag_state = DiagnosticState(
            test_name="test_travelling_vortex_3d_coriolis",
            file_name="target_travelling_vortex_3d_coriolis",
            Nx=self.inx - 1,
            Ny=self.iny - 1,
            steps=[self.stepmax - 1],
            # 3D fields cannot be contour-plotted by the comparison plotter
            plot_compare=False,
        )

        self.autogen_fn = False

    def stratification_function(self, y):
        if type(y) == float:
            return 1.0
        else:
            return np.ones((y.shape))

    def rhoe_function(self, rho, u, v, w, p, ud, th):
        Msq = ud.compressibility * ud.Msq
        gm1inv = th.gm1inv

        return p * gm1inv + 0.5 * Msq * rho * (u**2 + v**2 + w**2)


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed

    rotdir = 1.0

    a_rho = 1.0
    rho0 = a_rho * 0.5
    del_rho = a_rho * 0.5
    R0 = 0.4
    fac = 1.0 * 1024.0
    xc = 0.5
    zc = 0.5

    # Coriolis parameter entering the gradient-wind pressure balance
    f = ud.coriolis_strength[0]

    igs = elem.igs

    hydrostatics.integrated_state(npf, elem, node, th, ud)

    # centrifugal part of the balanced pressure (exponents 12..36)
    coe = np.zeros((25))
    coe[0] = 1.0 / 12.0
    coe[1] = -12.0 / 13.0
    coe[2] = 9.0 / 2.0
    coe[3] = -184.0 / 15.0
    coe[4] = 609.0 / 32.0
    coe[5] = -222.0 / 17.0
    coe[6] = -38.0 / 9.0
    coe[7] = 54.0 / 19.0
    coe[8] = 783.0 / 20.0
    coe[9] = -558.0 / 7.0
    coe[10] = 1053.0 / 22.0
    coe[11] = 1014.0 / 23.0
    coe[12] = -1473.0 / 16.0
    coe[13] = 204.0 / 5.0
    coe[14] = 510.0 / 13.0
    coe[15] = -1564.0 / 27.0
    coe[16] = 153.0 / 8.0
    coe[17] = 450.0 / 29.0
    coe[18] = -269.0 / 15.0
    coe[19] = 174.0 / 31.0
    coe[20] = 57.0 / 32.0
    coe[21] = -74.0 / 33.0
    coe[22] = 15.0 / 17.0
    coe[23] = -6.0 / 35.0
    coe[24] = 1.0 / 72.0

    # Coriolis part of the balanced pressure (exponents 7..25)
    ccoe = np.zeros((19))
    ccoe[0] = 1.0 / 7.0
    ccoe[1] = -3.0 / 4.0
    ccoe[2] = 4.0 / 3.0
    ccoe[3] = -1.0 / 5.0
    ccoe[4] = -45.0 / 22.0
    ccoe[5] = 3.0 / 4.0
    ccoe[6] = 9.0 / 2.0
    ccoe[7] = -36.0 / 7.0
    ccoe[8] = -11.0 / 5.0
    ccoe[9] = 55.0 / 8.0
    ccoe[10] = -33.0 / 17.0
    ccoe[11] = -4.0
    ccoe[12] = 58.0 / 19.0
    ccoe[13] = 3.0 / 5.0
    ccoe[14] = -10.0 / 7.0
    ccoe[15] = 4.0 / 11.0
    ccoe[16] = 9.0 / 46.0
    ccoe[17] = -1.0 / 8.0
    ccoe[18] = 1.0 / 50.0

    xcm = xc - (ud.xmax - ud.xmin)
    zcm = zc - (ud.zmax - ud.zmin)

    # cell-centred radius in the x-z plane, broadcast over y
    xs = elem.x.reshape(-1, 1, 1)
    zs = elem.z.reshape(1, 1, -1)
    xccs = np.zeros_like(xs)
    zccs = np.zeros_like(zs)

    xccs[...] = xc * (np.abs(xs - xc) < np.abs(xs - xcm))
    xccs[...] += xcm * (np.abs(xs - xc) > np.abs(xs - xcm))

    zccs[...] = zc * (np.abs(zs - zc) < np.abs(zs - zcm))
    zccs[...] += zcm * (np.abs(zs - zc) > np.abs(zs - zcm))

    r = np.sqrt((xs - xccs) ** 2 + (zs - zccs) ** 2)

    uth = (rotdir * fac * (1.0 - r / R0) ** 6 * (r / R0) ** 6) * (r < R0)

    u = u0 + uth * (-(zs - zccs) / r)
    v = v0 + np.zeros_like(r)
    w = w0 + uth * (+(xs - xccs) / r)

    rho = np.zeros_like(r)
    rho[...] += (rho0 + del_rho * (1.0 - (r / R0) ** 2) ** 6) * (r < R0)
    rho[...] += rho0 * (r >= R0)

    # broadcast the (icx, 1, icz) slabs over the full y extent (incl. ghosts),
    # as the legacy file did via np.repeat
    Sol.rho[...] = rho
    Sol.rhou[...] = rho * u
    Sol.rhov[...] = rho * v
    Sol.rhow[...] = rho * w

    # neutral stratification, pseudo-incompressible: rhoY = 1 everywhere
    # (the legacy file left npf.p2_cells unset)
    Sol.rhoY[...] = 1.0

    # nodal balanced pressure
    xs = node.x[igs[0] : -igs[0]].reshape(-1, 1, 1)
    zs = node.z[igs[2] : -igs[2]].reshape(1, 1, -1)
    xccs = np.zeros_like(xs)
    zccs = np.zeros_like(zs)

    xccs[np.where(np.abs(xs - xc) < np.abs(xs - xcm))] = xc
    xccs[np.where(np.abs(xs - xc) >= np.abs(xs - xcm))] = xcm

    zccs[np.where(np.abs(zs - zc) < np.abs(zs - zcm))] = zc
    zccs[np.where(np.abs(zs - zc) >= np.abs(zs - zcm))] = zcm

    r = np.sqrt((xs - xccs) ** 2 + (zs - zccs) ** 2)

    i2 = tuple(slice(igs[dim], -igs[dim]) for dim in range(elem.ndim))

    p2n = np.zeros_like(r)
    for ip in range(25):
        p2n += fac * (a_rho * coe[ip] * ((r / R0) ** (12 + ip) - 1.0) * rotdir**2)
    for ip in range(19):
        p2n += f * ccoe[ip] * ((r / R0) ** (7 + ip) - 1.0)
    p2n *= r / R0 < 1.0

    npf.p2_nodes[i2] = th.Gamma * fac * p2n

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    if ud.initial_projection == True:
        is_compressible = np.copy(ud.is_compressible)
        compressibility = np.copy(ud.compressibility)
        ud.is_compressible = 0
        ud.compressibility = 0.0

        p2aux = np.copy(npf.p2_nodes)

        Sol.rhou -= u0 * Sol.rho
        Sol.rhov -= v0 * Sol.rho
        Sol.rhow -= w0 * Sol.rho

        mem = obj()
        mem.sol = Sol
        mem.npf = npf
        mem.elem = elem
        mem.node = node
        mem.th = th
        mem.time = obj()
        mem.time.t = ud.dtfixed
        mem.time.step = 0
        mem.cache = cache.FlowSolverCache()

        implicit_euler.do_implicit_part(
            mem, ud, ud.dtfixed, writer=None, label="initial_projection"
        )

        npf.p2_nodes[...] = p2aux
        npf.dp2_nodes[...] = 0.0

        Sol.rhou += u0 * Sol.rho
        Sol.rhov += v0 * Sol.rho
        Sol.rhow += w0 * Sol.rho

        ud.is_compressible = is_compressible
        ud.compressibility = compressibility

    return Sol


def T_from_p_rho(p, rho):
    return np.divide(p, rho)


class obj(object):
    pass
