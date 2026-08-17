"""Hughes & Jablonowski (2023) mountain baroclinic wave on the FULL
pole-to-pole sphere — pt-2 with the two midlatitude ridges.

The global flat background (:mod:`test_hj_baroclinic_global`) with the two
H&J ridges (Eq. 1, h0 = 2000 m at 72E / 140E, 45N) embedded as GEOMETRY via a
pole-enabled ``SphericalTerrainMap`` (Gal-Chen radial coordinate), exactly as
the channel pt-2 (:mod:`test_hj_baroclinic_ridges`) does off the pole. The
whole reduction of pt-2 carries over: the IC is written entirely in
``metric.height``, so the pt-1 ``sol_init`` reproduces the H&J adjusted
balance (Eq. 2 = the analytic Ullrich state sampled at the terrain-following
height) unchanged once the map carries the ridges; only the map and the radial
axis (eta in [0, depth], not r) change.

The filter onset (phi_c = 70 deg, inherited) sits poleward of the 45 deg N
ridge band, so the ridges are never damped. Orography stays nondimensionalised
(/ h_ref -- the J-flips-at-the-ridge-latitudes trap, sphere analogue of
agnesi's h0/h_ref). Coriolis stays the constant ``2 Omega_nd (0, 0, +1)``.
Gate: ``test_scripts/test_hj_baroclinic_global.py`` (ridge-triggered
initiation, localised at the ridges, jet bounded).
"""

from ..flow_solver.discretisation import spherical

from . import test_hj_baroclinic_ridges as hjr
from .case_setup import make_diag_state
from .test_hj_baroclinic import _DEPTH_M
from .test_hj_baroclinic_ridges import _nondim_orography
from . import test_hj_baroclinic_global as hjg

# the IC is written entirely in metric.height; with the ridge map that height
# follows the terrain, so the flat-global IC IS the adjusted-balance state
sol_init = hjg.sol_init


class UserData(hjg.UserData):
    def __init__(self):
        super().__init__()

        a_nd = self.planet_radius / self.h_ref
        depth = _DEPTH_M / self.h_ref

        # the ridges are GEOMETRY: the H&J two-ridge orography + analytic
        # gradients (nondimensionalised / h_ref) folded into a pole-enabled
        # terrain-following radial coordinate. The radial axis now carries
        # eta in [0, depth] (Gal-Chen), so ymin/ymax change with it.
        self.ymin, self.ymax = 0.0, depth
        oro, oro_grad = _nondim_orography(self.h_ref)
        self.curvilinear_map = spherical.SphericalTerrainMap(
            a_nd, depth, oro, oro_grad, pole=True
        )

        self.output_base_name = "_hj_baroclinic_ridges_global"
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.inz - 1)
        self.diag_state = make_diag_state(
            "test_hj_baroclinic_ridges_global",
            "target_hj_baroclinic_ridges_global",
            self.inx,
            self.iny,
            self.stepmax,
        )
