import numpy as np
from ..flow_solver.physics import hydrostatics

from .case_setup import make_diag_state


class UserData(object):
    # Nsq_ref = grav * 1.3e-05

    def __init__(self):
        self.grav = 10.0  # [m/s^2]
        self.t_ref = 1000.0  # [s]

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.9
        self.dtfixed0 = 100.0
        self.dtfixed = 100.0

        self.inx = 64 + 1
        self.iny = 48 + 1
        self.inz = 1

        self.tout = [1000.0]
        self.stepmax = 31

        self.is_compressible = 1

        self.continuous_blending = False
        self.no_of_pi_initial = 1
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.initial_blending = True

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_blending_warm_bubble"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = "CFLfixed"

        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.output_timesteps = True

        # default 1e-5 tolerances on all 7 fields: the blended run reproduces
        # to ~1e-7, and looser momenta tolerances would let the psinc->comp
        # conversion be silently discarded (momenta shift ~0.1, final p2
        # increment ~1e-5)
        self.diag_state = make_diag_state(
            f"{self.output_type}_blending_warm_bubble",
            "target_blending_warm_bubble",
            self.inx,
            self.iny,
            self.stepmax,
            plot_compare=True,
            time_increment=True,
        )

        self.autogen_fn = False


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed
    delth = 2.0  # [K]

    y0 = 0.2
    r0 = 0.2

    hydrostatics.integrated_state(npf, elem, node, th, ud)

    x = elem.x
    y = elem.y

    x, y = np.meshgrid(x, y)

    # ensemble spread and truth IC of the MWR-2022 OSSE: uniformly sampled
    # bubble amplitude (git tag archive/full_coriolis,
    # RKLM_Python/inputs/rising_bubble.py); inert when seed is None and aux
    # does not contain 'truth'
    if seed is not None:
        np.random.seed(seed)
        delth += 10.0 * np.random.random()

    if "truth" in ud.aux:
        np.random.seed(1234)
        delth += 10.0 * np.random.random()

    r = np.sqrt((x) ** 2 + (y - y0) ** 2) / r0

    rhoY = npf.HydroState.rhoY0[np.newaxis, :]

    perturbation = (delth / 300.0) * (np.cos(0.5 * np.pi * r) ** 2)
    perturbation[np.where(r > 1.0)] = 0.0
    rho = rhoY / (ud.stratification(y) + perturbation.T)

    x_idx = slice(None)
    y_idx = slice(None)

    u, v, w = u0, v0, w0

    Sol.rho[x_idx, y_idx] = rho
    Sol.rhou[x_idx, y_idx] = rho * u
    Sol.rhov[x_idx, y_idx] = rho * v
    Sol.rhow[x_idx, y_idx] = rho * w
    Sol.rhoY[x_idx, y_idx] = rhoY

    # imbalanced IC: a large Exner-pressure blob (amplitude 1.0 ~ 10x the
    # slow signal), momenta untouched. The initial comp->psinc conversion
    # must absorb it in one blended step -- validated against the frozen
    # paper-era code (git tag archive/localdab). Unblended, this blob
    # detonates: p2 leaves the balanced trajectory by 2-5x the total signal,
    # which is what the p2_nodes increment tolerance below guards against.
    xn, yn = np.meshgrid(node.x, node.y, indexing="ij")
    npf.p2_nodes[...] = np.exp(-(xn**2 + (yn - 0.5) ** 2) / (2 * 0.15**2))

    return Sol
