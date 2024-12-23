import numpy as np

from . import (
    user_data,
    io,
    data_structures
)

from ..flow_solver.discretisation import grid           as dis_grid
from ..flow_solver.utils import variable                as var
from ..flow_solver.utils import boundary as bdry
from ..flow_solver.physics import hydrostatics
from ..flow_solver.physics.low_mach import mpv          as lm_var
from ..flow_solver.physics.gas_dynamics import thermodynamics as gd_thermodynamics

# test module
from ..tests import diagnostics as diag

def initialise():
    ####
    # Initialise simulation state
    ####
    from . import sim_params as params

    np.set_printoptions(precision = params.print_precision)

    ##########################################################
    # Initialisation of data containers and helper classes
    ##########################################################
    # get arguments for initial condition and ensemble size
    N, UserData, sol_init, restart, ud_rewrite, dap_rewrite, r_params = io.get_args()
    if N == 1:
        params.da_debug = False

    initial_data = vars(UserData())
    ud = user_data.UserDataInit(**initial_data)
    if ud_rewrite is not None:
        ud.update_ud(ud_rewrite)
    if hasattr(ud, "rayleigh_bc"):
        ud.rayleigh_bc(ud)
    if ud.output_timesteps:
        params.output_timesteps = True
    ud.coriolis_strength = np.array(ud.coriolis_strength)

    elem, node = dis_grid.grid_init(ud)

    sol = var.Vars(elem.sc, ud)

    flux = np.empty((3), dtype=object)
    flux[0] = var.States(elem.sfx, ud)
    if elem.ndim > 1:
        flux[1] = var.States(elem.sfy, ud)
    if elem.ndim > 2:
        flux[2] = var.States(elem.sfz, ud)

    th = gd_thermodynamics.ThermodynamicalQuantities(ud)
    mpv = lm_var.MPV(elem, node, ud)
    

    io.init_logger(ud)

    # handle radiative BC
    if ud.bdry_type[1].value == "radiation":
        ud.tcy, ud.tny = bdry.get_tau_y(ud, elem, node, 0.5)

    ##########################################################
    # Initialise test module
    ##########################################################
    if ud.diag:
        diag_comparison = diag.compare_sol(ud.diag_current_run)
    else:
        diag_comparison = None


    ##########################################################
    # Populate data structures
    ##########################################################

    # member_state = data_structures.MemberState(
    #     elem=elem,
    #     node=node,
    #     Sol=Sol,
    #     flux=flux,
    #     mpv=mpv,
    #     th=th,
    # )

    ensemble_state = data_structures.EnsembleState()

    sol = sol_init(sol, mpv, elem, node, th, ud)

    ensemble_state.update_member(
                elem=elem,
                node=node,
                sol=sol,
                flux=flux,
                mpv=mpv,
                th=th,
                )

    restart_params = data_structures.RestartParameters(
        ud_rewrite=ud_rewrite,
        dap_rewrite=dap_rewrite,
        r_params=r_params,
    )

    interface_params = data_structures.InterfaceParameters()

    sim_st = data_structures.SimulationState(
        N=N,
        restart=restart,

        ud=ud,
        sol_init=sol_init,

        ensemble_state=ensemble_state,

        diag_comparison=diag_comparison,

        restart_params=restart_params,
        interface_params=interface_params
    )

    return sim_st


def overwrite_init_with_restart(sst):
    es = sst.ensemble_state
    rp = sst.restart_params

    hydrostatics.state(es.mpv, es.elem, es.node, es.th, es.ud)
    sst.ud.old_suffix = np.copy(sst.ud.output_suffix)
    sst.ud.old_suffix = "_ensemble=%i%s" % (sst.N, sst.ud.old_suffix)
    Sol0, mpv0, touts = io.sim_restart(
        rp.r_params[0], rp.r_params[1], es.elem, es.node, es.ud, es.Sol, es.mpv, rp.r_params[2]
    )
    sol_ens = [[Sol0, es.flux, mpv0, [-np.inf, sst.step]]]
    # ud.tout = touts[1:]
    sst.ud.tout = [touts[-1]]
    sst.t = touts[0]

    sst.ensemble_state = sol_ens