"""Hughes & Jablonowski (2023) mountain baroclinic wave — pt 1: the
well-balanced sphere initial condition (NO topography) + steadiness gate.

The Ullrich et al. (2014/2016) dry balanced baroclinic base state
(:mod:`pybella.tests.ullrich_baroclinic`, gated by
``test_scripts/test_hj_background.py``) mapped onto the full-size deep
spherical shell as a pyBELLA initial condition. The crux of the Hughes &
Jablonowski case is to split the 3D balance so the discretisation carries
each part the way it expects:

* **Vertical hydrostatic balance -> the field-mode HydroState reference.**
  Built by quadrature (``hydrostatics.integrated_state`` -> the terrain
  ``_integrated_state_fields`` branch) from a z-only EQUATORIAL Ullrich
  potential-temperature column supplied through ``ud.stratification``. At
  the equator the jet vanishes (``U(0, z) == 0``, Appendix B) and the
  surface pressure is exactly ``P0``, so the equatorial column IS the
  resting reference: ``rhoY0 == 1`` at ``r = a`` and ``p2 == 0`` there.

* **Meridional (phi) pressure structure -> the Exner perturbation
  ``p2_nodes``.** The momentum pressure-gradient force uses ONLY
  ``grad(p2)`` (``explicit_euler.do_forward_step``): the reference
  contributes only buoyancy + the vertical balance. So the meridional
  pressure gradient that the Coriolis force on the jet must balance has to
  live in ``p2`` — exactly the 3D analogue of TC2 placing the geostrophic
  depth in ``p2_nodes``. ``p2 = (pi_full - pi_ref) / Msq``.

* **Full Coriolis: the constant rotation vector ``2 Omega`` in the embedded
  frame**, ``2 Omega_nd * (0, 0, +1)`` (the mirrored-embedding pseudovector
  flip; NOT the geographic pole — see
  ``SphericalShellMap.rotation_axis_cart``). The wrong sign anti-balances
  and blows up immediately.

Gradient-wind balance is a steady state of the CONTINUOUS equations;
discretely it is only approximate. In practice the truncation-level residual
is tiny — the analytic momenta are purely zonal (no meridional wind at init),
and stepping generates only an O(0.1 m/s) meridional adjustment that decays
(<< the jet, << the eventual O(10 m/s) baroclinic wave). So this case simply
ACCEPTS the small adjustment (``initial_projection = False``); the steadiness
is asserted in ``test_scripts/test_hj_baroclinic.py``.

(The TC2-style incompressible initial projection is NOT used here: it would
exercise the incompressible + field-mode-HydroState gravity ghost fill in
``cell_boundary._calculate_ghost_values`` — a path no prior case hits, since
the SWE sphere cases that project have ``grav = 0``. That path indexes the
3D field-mode ``HydroState.rhoY0`` with only the vertical component
``nimage[y_axs]`` (correct for the 1D profile mode, wrong for field mode);
projection here would need that fixed first. pt 2's adjusted ridge balance
may want it.)

Nondimensionalization mirrors the sphere gravity-wave case
(``h_ref = R_gas T_ref / grav`` -> ``gravity_strength[1] = 1``), with
``p_ref = P0 = 1e5`` so ``rhoY = (p / P0) ** (1/gamma)`` and
``Y = theta / T_ref``. The dimensional constants (``grav``, ``R_gas``,
``P0``, planet radius, rotation rate) are taken from
:mod:`pybella.tests.ullrich_baroclinic` so the reference and the full state
are built from one consistent set.
"""

import numpy as np

from ..utils import options as opts
from ..flow_solver.discretisation import spherical
from ..flow_solver.physics import hydrostatics
from ..flow_solver.utils.boundary import node_boundary as bdry_n

from . import ullrich_baroclinic as ub
from .case_setup import build_bdry, do_initial_projection, make_diag_state

# --- reference scales (shared thermodynamics with the sphere GW case) ------
_T_REF = 250.0  # [K]
_U_REF = 10.0  # [m/s]
_PHI_MAX = np.deg2rad(80.0)  # channel walls in latitude (see the global case)
_DEPTH_M = 30.0e3  # shell depth r_top - a [m] (model top ~30 km)

# Ullrich formulas use cos(phi) ** K with fractional K, so a negative cosine
# (latitude ghosts past the pole, |phi| > 90 deg) would be NaN. Clamp the
# latitude used for the ANALYTIC fill just inside the pole; those ghost rows
# are overwritten by the free-slip wall reflection before any dynamics.
_PHI_CLIP = np.deg2rad(89.5)


def _constant_coriolis_field(two_omega_nd):
    """``ud.coriolis_field`` for the FULL (non-traditional) Coriolis force:
    the constant embedded rotation vector ``2 Omega_nd * (0, 0, +1)``. The
    callable takes the three Cartesian coordinate fields and returns the
    three Cartesian rotation-vector components the H^-1 kernel consumes."""

    def field(x0, x1, x2):
        zero = 0.0 * x0
        return (zero, zero, two_omega_nd + zero)

    return field


class UserData(object):
    def __init__(self):
        self.grav = ub.G
        self.omega = 0.0  # rotation enters via coriolis_field, not omega
        self.R_gas = ub.R_D
        self.gamm = 1.4
        self.T_ref = _T_REF
        self.u_ref = _U_REF
        self.p_ref = ub.P0
        self.h_ref = self.R_gas * self.T_ref / self.grav
        self.t_ref = self.h_ref / self.u_ref
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        self.is_compressible = 1
        self.is_nonhydrostatic = 1

        # radial gravity along the y = r axis, unit magnitude (h_ref choice)
        self.gravity_direction = 1
        self.gravity_strength = np.zeros(3)
        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        self.i_gravity = np.zeros(3)
        self.i_gravity[1] = 1

        self.planet_radius = ub.A_EARTH
        a_nd = self.planet_radius / self.h_ref
        depth = _DEPTH_M / self.h_ref

        # axes: x = lambda (periodic), y = r (gravity), z = phi (walls)
        self.xmin, self.xmax = -np.pi, np.pi
        self.ymin, self.ymax = a_nd, a_nd + depth
        self.zmin, self.zmax = -_PHI_MAX, _PHI_MAX

        self.curvilinear_map = spherical.SphericalShellMap(a_nd, frozen_radius=False)

        omega_nd = ub.OMEGA * self.t_ref
        self.coriolis_field = _constant_coriolis_field(2.0 * omega_nd)

        self.u_wind_speed = 0.0
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.WALL
        )

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.45
        # advective CFL is generous on the full planet (implicit acoustics);
        # a moderate physical dt keeps the coarse smoke well inside stability
        self.dtfixed = 300.0 / self.t_ref
        self.dtfixed0 = 300.0 / self.t_ref

        self.inx = 64 + 1
        self.iny = 20 + 1
        self.inz = 32 + 1

        # accept the (truncation-small) discrete balance residual; see the
        # module docstring for why the TC2-style projection is not used here
        self.initial_projection = False

        self.tout = [1.0e9]
        self.stepmax = 20

        self.stratification = self.stratification_function
        self.rhoe = self.rhoe_function
        self.output_timesteps = True

        self.diag = False
        self.diag_updt_targets = False

        self.output_base_name = "_hj_baroclinic"
        self.output_type = "test"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.inz - 1)
        self.autogen_fn = False

        self.diag_state = make_diag_state(
            "test_hj_baroclinic",
            "target_hj_baroclinic",
            self.inx,
            self.iny,
            self.stepmax,
        )

    def stratification_function(self, z):
        """Reference (background) potential temperature theta(z) / T_ref of
        the z-only EQUATORIAL Ullrich column. ``z`` is the NONDIMENSIONAL
        height (metric.height = r - a); dimensionalise, evaluate Ullrich at
        the equator (phi = 0), convert to potential temperature with the
        standard P0 reference, and nondimensionalise by T_ref."""
        z_dim = z * self.h_ref
        kappa = (self.gamm - 1.0) / self.gamm
        T = ub.temperature(0.0, z_dim)
        p = ub.pressure(0.0, z_dim)
        theta = T * (ub.P0 / p) ** kappa
        return theta / self.T_ref

    def rhoe_function(self, rho, u, v, w, p, ud, th):
        Msq = ud.compressibility * ud.Msq
        gm1inv = th.gm1inv
        return p * gm1inv + 0.5 * Msq * rho * (u**2 + v**2 + w**2)


def _coord(grid_obj, axis):
    """1D grid coordinate along ``axis`` reshaped to broadcast over the grid."""
    shape = [1, 1, 1]
    shape[axis] = -1
    arr = (grid_obj.x, grid_obj.y, grid_obj.z)[axis]
    return arr.reshape(shape)


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    # latitude clamp for the analytic Ullrich fill. The channel case walls at
    # +-80 deg, so its ghost latitudes must be clipped just inside the pole
    # (cos^K is NaN past |phi| = 90 deg); the global (pole=True) variant folds
    # its ghosts back to interior latitudes and puts a node exactly at the
    # pole, so it overrides this to pi/2. Both leave the interior
    # untouched (interior cells / non-pole nodes are strictly inside the clip).
    phi_clip = float(getattr(ud, "phi_clip", _PHI_CLIP))

    # --- vertical hydrostatic balance: z-only equatorial reference ---------
    # integrated_state routes to the field-mode terrain branch (g != 0 and a
    # metric is present): a discretely well-balanced column from the
    # equatorial Ullrich theta profile (ud.stratification), constant in
    # (lambda, phi) since the shell height depends only on r.
    hydrostatics.integrated_state(npf, elem, node, th, ud)
    S0c = npf.HydroState.get_S0c(elem)  # reference 1/theta (cell field)

    # --- full analytic Ullrich state on the CELL grid ----------------------
    lam_c = _coord(elem, 0)
    phi_c = _coord(elem, 2)
    phi_eval = np.clip(phi_c, -phi_clip, phi_clip)
    z_c = elem.metric.height * ud.h_ref  # dimensional height [m]

    p_c = ub.pressure(phi_eval, z_c)
    rho_c = ub.density(phi_eval, z_c)
    u_c = ub.zonal_wind(phi_eval, z_c)

    rho_nd = rho_c * ud.R_gas * ud.T_ref / ud.p_ref  # rho / rho_ref
    rhoY = (p_c / ud.p_ref) ** th.gamminv  # P = (p / P0) ** (1/gamma)
    u_nd = u_c / ud.u_ref

    # zonal wind in fixed Cartesian components: u * e_lambda,
    # e_lambda = (-sin lambda, cos lambda, 0) (v = w = 0)
    sl, cl = np.sin(lam_c), np.cos(lam_c)

    shp = Sol.rho.shape
    Sol.rho[...] = np.broadcast_to(rho_nd, shp)
    Sol.rhoY[...] = np.broadcast_to(rhoY, shp)
    Sol.rhou[...] = np.broadcast_to(rho_nd * u_nd * (-sl), shp)
    Sol.rhov[...] = np.broadcast_to(rho_nd * u_nd * cl, shp)
    Sol.rhow[...] = 0.0
    Sol.rhoX[...] = Sol.rho * (Sol.rho / Sol.rhoY - S0c)

    # --- meridional pressure structure -> Exner perturbation p2 (nodes) ----
    phi_n = np.clip(_coord(node, 2), -phi_clip, phi_clip)
    z_n = node.metric.height * ud.h_ref
    pi_full_n = (ub.pressure(phi_n, z_n) / ud.p_ref) ** th.Gamma  # (p/P0)^kappa
    npf.p2_nodes[...] = pi_full_n / ud.Msq - npf.HydroState_n.p20
    bdry_n.set_ghost_nodes(npf.p2_nodes, node, ud)

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    do_initial_projection(Sol, npf, elem, node, th, ud, u0=0.0, v0=0.0)

    return Sol
