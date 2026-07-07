"""Williamson TC6 (Rossby-Haurwitz wave) on the FULL pole-to-pole sphere.

The +-80 deg channel :mod:`test_sphere_swe_tc6` extended to |phi| <= pi/2
with ``BdryType.POLE`` walls, a pole-enabled shell map and the polar filter.
Unlike TC2, TC6 has genuine longitude structure at all latitudes (the RH-4
pattern), so it exercises the pole ghost exchange and the polar filter under
real dynamics — the RH velocities still vanish at the poles (cos^k phi), so
the state is pole-regular. The 7-day phase-speed validation runs from
run_scripts; the registered regression is a short tripwire.
"""

import numpy as np

from ..flow_solver.discretisation import spherical
from ..flow_solver.numerics import polar_filter
from ..utils import options as opts
from .case_setup import build_bdry, make_diag_state
from . import test_sphere_swe_tc6 as tc6

A_ND = tc6.A_ND
OMEGA_ND = tc6.OMEGA_ND
NU_PHASE = tc6.NU_PHASE

sol_init = tc6.sol_init


class UserData(tc6.UserData):
    def __init__(self):
        super().__init__()
        self.zmin, self.zmax = -0.5 * np.pi, 0.5 * np.pi
        self.inz = 72 + 1
        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.POLE
        )
        self.curvilinear_map = spherical.SphericalShellMap(
            A_ND, frozen_radius=True, pole=True
        )
        self.polar_filter = polar_filter.PolarFilter(np.deg2rad(60.0))

        self.output_base_name = "_sphere_swe_tc6_global"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.inz - 1)
        self.diag_state = make_diag_state(
            "test_sphere_swe_tc6_global",
            "target_sphere_swe_tc6_global",
            self.inx,
            self.iny,
            self.stepmax,
            plot_compare=True,
        )
