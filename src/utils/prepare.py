import numpy as np

from . import (
    user_data,
    io,
    data_structures
)

from ..dycore.discretisation import grid           as dis_grid
from ..dycore.utils import variable                as var
from ..dycore.utils import boundary as bdry
from ..dycore.physics.low_mach import mpv          as lm_var
from ..dycore.physics.gas_dynamics import thermodynamic as gd_thermodynamics

# test module
from ..tests import diagnostics as diag

def initialise():
    ####
    # Initialise simulation state
    ####
    from . import sim_params as params

    np.set_printoptions(precision = params.print_precision)

    step = 0
    t = 0.0

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

    Sol = var.Vars(elem.sc, ud)

    flux = np.empty((3), dtype=object)
    flux[0] = var.States(elem.sfx, ud)
    if elem.ndim > 1:
        flux[1] = var.States(elem.sfy, ud)
    if elem.ndim > 2:
        flux[2] = var.States(elem.sfz, ud)

    th = gd_thermodynamics.init(ud)
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

    model_params = data_structures.ModelParameters(
        elem=elem,
        node=node,
        Sol=Sol,
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
        step=step,
        t=t,
        N=N,
        restart=restart,

        ud=ud,
        sol_init=sol_init,

        diag_comparison=diag_comparison,

        model_params=model_params,
        restart_params=restart_params,
        interface_params=interface_params
    )

    return sim_st