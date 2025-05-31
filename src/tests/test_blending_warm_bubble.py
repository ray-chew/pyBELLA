import numpy as np
from ..flow_solver.physics import hydrostatics
from ..flow_solver.utils import boundary as bdry

from ..utils.data_structures import DiagnosticState


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

        # if self.is_compressible == 1:
        #     self.output_suffix = "_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.tout[-1])
        # if self.is_compressible == 0:
        #     self.output_suffix = "_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.tout[-1])
        # if self.continuous_blending == True:
        #     self.output_suffix = "_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.tout[-1])

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

        self.diag_state = DiagnosticState(
            test_name=f"{self.output_type}_blending_warm_bubble",
            file_name="target_blending_warm_bubble",
            Nx=self.inx - 1,
            Ny=self.iny - 1,
            steps=[self.stepmax - 1],
            plot_compare=True,
            time_increment=True,
            # The only thing that matters here is that
            # p2_nodes remains small
            tolerances={
                "rho": 1.0e-0,
                "rhou": 1.0e-0,
                "rhov": 1.0e-0,
                "rhow": 1.0e-0,
                "rhoY": 1.0e-0,
                "rhoX": 1.0e-0,
                "p2_nodes": 1.0e-4,
            }
        )

        self.autogen_fn = False


def sol_init(Sol, mpv, elem, node, th, ud, seed=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed
    delth = 2.0  # [K]

    y0 = 0.2
    r0 = 0.2

    hydrostatics.state(mpv, elem, node, th, ud)

    x = elem.x
    y = elem.y

    x, y = np.meshgrid(x, y)

    r = np.sqrt((x) ** 2 + (y - y0) ** 2) / r0

    p = np.repeat(mpv.HydroState.p0.reshape(1, -1), elem.icx, axis=0)
    rhoY = mpv.HydroState.rhoY0[
        np.newaxis, :
    ]  # np.repeat(mpv.HydroState.rhoY0.reshape(1,-1),elem.icx,axis=0)

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

    p = mpv.HydroState_n.p0
    rhoY = mpv.HydroState_n.rhoY0
    mpv.p2_nodes[...] = (p - mpv.HydroState_n.p0) / rhoY / ud.Msq
    # mpv.p2_nodes[...] = 1.0

    bdry.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    return Sol
