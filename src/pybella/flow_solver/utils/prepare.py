import numpy as np

from ...utils import axes, user_data, io, data_structures
from ..physics import hydrostatics
from ..physics import thermodynamics as gd_thermodynamics
from ..discretisation import grid as dis_grid
from . import fields, cache
from .boundary import cell_boundary as bdry_c

# test module
from ...tests import diagnostics as diag


def initialise():
    ####
    # Initialise simulation state
    ####
    from ...utils import sim_params as params

    np.set_printoptions(precision=params.print_precision)

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
    axes.validate(ud, elem.ndim)

    sol = fields.CellSolField(elem.sc)

    th = gd_thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node)

    io.init_logger(ud)

    ##########################################################
    # Initialise test module
    ##########################################################
    if ud.diag:
        diag_comparison = diag.CompareSol(ud.diag_state)
    else:
        diag_comparison = None

    ##########################################################
    # Populate data structures
    ##########################################################

    ensemble_state = data_structures.EnsembleState()

    sol = sol_init(sol, npf, elem, node, th, ud)

    # Initialise cache and add to simulation state
    flow_cache = cache.FlowSolverCache()

    ensemble_state.update_member(
        elem=elem, node=node, sol=sol, npf=npf, th=th, cache=flow_cache
    )

    for member in ensemble_state.members:
        bdry_c.set_ghost_cells(member, ud)

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
        interface_params=interface_params,
    )

    return sim_st


def overwrite_init_with_restart(sst):
    es = sst.ensemble_state
    rp = sst.restart_params

    hydrostatics.state(es.npf, es.elem, es.node, es.th, es.ud)
    sst.ud.old_suffix = np.copy(sst.ud.output_suffix)
    sst.ud.old_suffix = "_ensemble=%i%s" % (sst.N, sst.ud.old_suffix)
    Sol0, npf0, touts = io.sim_restart(
        rp.r_params[0],
        rp.r_params[1],
        es.elem,
        es.node,
        es.ud,
        es.Sol,
        es.npf,
        rp.r_params[2],
    )
    sol_ens = [[Sol0, es.flux, npf0, [-np.inf, sst.step]]]
    # ud.tout = touts[1:]
    sst.ud.tout = [touts[-1]]
    sst.t = touts[0]

    sst.ensemble_state = sol_ens
