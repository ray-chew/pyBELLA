"""Hughes & Jablonowski (2023) mountain baroclinic wave on the FULL
pole-to-pole sphere (Stage F, F8) — pt-1 flat background.

The +-80 deg channel :mod:`test_hj_baroclinic` extended to |phi| <= pi/2:
``BdryType.POLE`` latitude walls, a pole-enabled deep ``SphericalShellMap``
(``pole=True``) and the FFT-in-longitude polar filter with the onset latitude
poleward of the ridge band (phi_c = 70 deg, so the 45 deg N ridges the pt-2
variant adds are never filtered). Two channel workarounds fall away here (kept
in the channel case):

* the ``_PHI_CLIP = 89.5 deg`` Ullrich cos^K clamp -- the folded ghosts are
  interior latitudes now and a node sits EXACTLY at the pole, so the clamp is
  widened to pi/2 (via ``ud.phi_clip``); the analytic jet vanishes and the
  pressure is finite at the pole, so no NaN;
* the "phi grid >= 32" constraint -- the ghosts fold instead of overshooting
  the pole (J = r^2 cos phi >= 0 with the pole nodes snapped to J = 0), so any
  ``inz`` is valid.

Everything else is pt 1 verbatim: the well-balanced Ullrich IC written in
``metric.height`` (``hj.sol_init``), the constant embedded Coriolis vector
``2 Omega_nd (0, 0, +1)`` (pseudovector -- do NOT re-derive), option (c)
balance (no projection). Gate: ``test_scripts/test_hj_baroclinic_global.py``.
"""

import numpy as np

from ..utils import options as opts
from ..flow_solver.discretisation import spherical
from ..flow_solver.numerics import polar_filter
from .case_setup import build_bdry, make_diag_state
from . import test_hj_baroclinic as hj

# the IC is written entirely in metric.height and is valid pole-to-pole
# (the jet vanishes at the poles); reuse it verbatim
sol_init = hj.sol_init


class UserData(hj.UserData):
    def __init__(self):
        super().__init__()

        # extend the latitude band to the poles; POLE walls + pole-enabled map
        self.zmin, self.zmax = -0.5 * np.pi, 0.5 * np.pi
        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.POLE
        )
        a_nd = self.planet_radius / self.h_ref
        self.curvilinear_map = spherical.SphericalShellMap(
            a_nd, frozen_radius=False, pole=True
        )
        # filter onset poleward of the 45 deg N ridge band (pt 2), so the
        # ridges are never damped; the CFL cap buys the larger dt near the pole
        self.polar_filter = polar_filter.PolarFilter(np.deg2rad(70.0))
        # a node sits exactly at the pole now (and ghosts fold to interior
        # latitudes) -> evaluate the analytic fill up to pi/2, not 89.5 deg
        self.phi_clip = 0.5 * np.pi

        # ghosts fold instead of overshooting the pole, so the channel's
        # "phi grid pinned near 32" constraint is lifted; keep a global-scale
        # resolution
        self.inz = 32 + 1

        self.output_base_name = "_hj_baroclinic_global"
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.inz - 1)
        self.diag_state = make_diag_state(
            "test_hj_baroclinic_global",
            "target_hj_baroclinic_global",
            self.inx,
            self.iny,
            self.stepmax,
        )
