"""Williamson TC2 on the FULL pole-to-pole sphere (Stage F acid test).

The +-80 deg channel :mod:`test_sphere_swe_tc2` extended to |phi| <= pi/2:
``BdryType.POLE`` latitude walls, a pole-enabled ``SphericalShellMap`` and the
FFT-in-longitude polar filter. TC2's steady zonal flow u0 cos(phi) e_lambda
VANISHES at the poles and every field is longitude-independent, so the filter
is a no-op and the exact steady state is pole-regular — any drift is our
discretisation. This is the acid test that the pole ghost exchange (F1), the
pole-face flux closure (F2), the polar filter (F3) and the elliptic pole
collapse (F4) compose correctly in a real forecast.
"""

import numpy as np

from ..flow_solver.discretisation import spherical
from ..flow_solver.numerics import polar_filter
from ..utils import options as opts
from .case_setup import build_bdry, make_diag_state
from . import test_sphere_swe_tc2 as tc2

A_ND = tc2.A_ND
OMEGA_ND = tc2.OMEGA_ND
U0_ND = tc2.U0_ND

# reuse the analytic steady state; sol_init fills every cell analytically
# (valid pole-to-pole: u0 cos(phi) -> 0 at the poles) then projects
depth_and_velocity = tc2.depth_and_velocity
sol_init = tc2.sol_init


class UserData(tc2.UserData):
    def __init__(self):
        super().__init__()
        # extend the latitude band to the poles
        self.zmin, self.zmax = -0.5 * np.pi, 0.5 * np.pi
        self.inz = 72 + 1  # ~2.5 deg/cell, matching the channel spacing
        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.POLE
        )
        self.curvilinear_map = spherical.SphericalShellMap(
            A_ND, frozen_radius=True, pole=True
        )
        self.polar_filter = polar_filter.PolarFilter(np.deg2rad(60.0))

        self.output_base_name = "_sphere_swe_tc2_global"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.inz - 1)
        self.diag_state = make_diag_state(
            "test_sphere_swe_tc2_global",
            "target_sphere_swe_tc2_global",
            self.inx,
            self.iny,
            self.stepmax,
            plot_compare=True,
        )
