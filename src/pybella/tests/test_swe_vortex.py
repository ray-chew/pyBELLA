import numpy as np

from ..utils import options as opts
from ..flow_solver.utils.boundary import node_boundary as bdry_n
from ..flow_solver.physics import hydrostatics

from .case_setup import (
    build_bdry,
    do_initial_projection,
    make_diag_state,
    mirror_centers,
)


class UserData(object):
    """
    Balanced shallow-water vortex — demonstration-grade SWE regression case.

    Port of ``RKLM_Python/inputs/balanced_shallow_water_2D.py`` (git tag
    ``archive/full_coriolis``). The shallow-water equations are run through
    the gas-dynamics solver via the gamma = 2 equivalence::

        rho  <->  fluid depth h
        p    <->  g h^2 / 2
        rhoY <->  p^(1/gamma) = sqrt(g/2) h

    with Msq = 1 (h_ref = t_ref = T_ref = R_gas = 1) and zero vertical
    gravity; the SWE gravity g enters only through the pressure law.

    Deviations from that case (kept to demo scope):

    - grid reduced from 150x150 to 64x64, run shortened to t = 31000 s
      (31 steps of dt = 1000 s),
    - x and y boundaries set periodic (the original used walls; the mirror-image
      vortex construction assumes periodicity),
    - the nodal pressure is evaluated analytically on the node grid instead
      of cubic scattered-data interpolation of the cell field, so the golden
      master is independent of the scipy ``griddata`` implementation.
    """

    grav = 0.0
    omega = 0.0

    R_gas = 1.0
    gamm = 2.0

    h_ref = 1.0
    t_ref = 1.0
    T_ref = 1.0
    p_ref = 1.0

    def __init__(self):
        self.grav = self.grav
        self.omega = self.omega
        self.R_gas = self.R_gas
        self.gamm = self.gamm
        self.h_ref = self.h_ref
        self.t_ref = self.t_ref
        self.T_ref = self.T_ref
        self.p_ref = self.p_ref

        self.xmin = -0.5e6
        self.xmax = 0.5e6
        self.ymin = -0.5e6
        self.ymax = 0.5e6
        self.zmin = -0.5
        self.zmax = 0.5

        self.u_wind_speed = 0.0
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        # SWE gravitational acceleration; enters via the pressure law
        # p = g h^2 / 2 only (the vertical gravity ``grav`` is zero).
        self.g_swe = 9.81

        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC, opts.BdryType.PERIODIC, opts.BdryType.WALL
        )

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.9 / 2.0
        self.dtfixed = 1000.0
        self.dtfixed0 = 1000.0

        self.inx = 64 + 1
        self.iny = 64 + 1
        self.inz = 1

        self.initial_projection = True

        # the run is stopped by the stepmax cap (31 steps of dt = 1000 s,
        # i.e. t = 31000 s); tout is set beyond it as in the other tests.
        self.tout = [1.0e6]
        self.stepmax = 31

        self.stratification = self.stratification_function
        self.rhoe = self.rhoe_function
        self.output_timesteps = True

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_swe_vortex"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.diag_state = make_diag_state(
            "test_swe_vortex",
            "target_swe_vortex",
            self.inx,
            self.iny,
            self.stepmax,
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


def _depth_field(xs, ys, xc, yc, xcm, ycm, R0, fac, Frsq):
    """
    Cyclostrophically balanced depth field h(r) for the vortex with the
    tangential velocity profile uth = fac * (1 - r/R0)^6 * (r/R0)^6: the
    polynomial below is the exact radial integral of uth^2 / (g r).
    """
    coe = np.zeros((13))
    coe[0] = +1.0 / 12
    coe[1] = -12.0 / 13
    coe[2] = +33.0 / 7
    coe[3] = -44.0 / 3
    coe[4] = +495.0 / 16
    coe[5] = -792.0 / 17
    coe[6] = +154.0 / 3
    coe[7] = -792.0 / 19
    coe[8] = +99.0 / 4
    coe[9] = -220.0 / 21
    coe[10] = +3.0
    coe[11] = -12.0 / 23
    coe[12] = +1.0 / 24

    xccs = mirror_centers(xs, xc, xcm)
    yccs = mirror_centers(ys, yc, ycm)

    r = np.sqrt((xs - xccs) ** 2 + (ys - yccs) ** 2)

    rho = np.zeros_like(r)
    for i in range(12, 24 + 1):
        rho[...] += fac**2 * coe[i - 12] * (r / R0) ** i * (r < R0)

    rho *= Frsq
    rho = (rho - rho.max()) * (r < R0)
    rho += 1.0

    return rho, r, xccs, yccs


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = 0.0

    rotdir = 1.0

    R0 = 400000.0
    fac = 1.0 * 1024.0
    xc = 0.0
    yc = 0.0

    g = ud.g_swe
    Frsq = 1.0 / g

    xcm = xc - (ud.xmax - ud.xmin)
    ycm = yc - (ud.ymax - ud.ymin)

    igs = elem.igs
    igy = igs[1]

    igxn = node.igx
    igyn = node.igy

    hydrostatics.integrated_state(npf, elem, node, th, ud)

    # cell-centred depth (rho), radius and vortex centres
    xs = elem.x.reshape(-1, 1)
    ys = elem.y[igy:-igy].reshape(1, -1)
    rho, r, xccs, yccs = _depth_field(xs, ys, xc, yc, xcm, ycm, R0, fac, Frsq)

    uth = (rotdir * fac * (1.0 - r / R0) ** 6 * (r / R0) ** 6) * (r < R0)

    u = u0 + uth * (-(ys - yccs) / r)
    v = v0 + uth * (+(xs - xccs) / r)
    w = w0

    Sol.rho[:, igy:-igy] = rho
    Sol.rhou[:, igy:-igy] = rho * u
    Sol.rhov[:, igy:-igy] = rho * v
    Sol.rhow[:, igy:-igy] = rho * w

    # shallow-water pressure law p = g h^2 / 2; rhoY = p^(1/gamma)
    p = g / 2.0 * rho**2
    Sol.rhoY[:, igy:-igy] = p**th.gamminv

    # nodal pressure: same analytic depth field evaluated on the node grid
    xs_n = node.x[igxn:-igxn].reshape(-1, 1)
    ys_n = node.y[igyn:-igyn].reshape(1, -1)
    rho_n, _, _, _ = _depth_field(xs_n, ys_n, xc, yc, xcm, ycm, R0, fac, Frsq)

    p_n = g / 2.0 * rho_n**2
    npf.p2_nodes[igxn:-igxn, igyn:-igyn] = p_n**th.gamminv
    bdry_n.set_ghost_nodes(npf.p2_nodes, node, ud)

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    do_initial_projection(Sol, npf, elem, node, th, ud, u0=u0, v0=v0)

    return Sol
