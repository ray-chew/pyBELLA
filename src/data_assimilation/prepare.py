
import numpy as np

from ..utils.sim_params import params

import params as da_params
import letkf as da_letkf

import logging
# to generate ensemble from one sol init instantiation
from copy import deepcopy 


def initialise(sst):
    ##########################################################
    # Initialisation of data assimilation module
    ##########################################################

    # possible da_types:
    # 1) batch_obs for the LETKF with batch observations
    # 2) rloc for LETKF with grid-point localisation
    # 3) etpf for the ETPF algorithm
    dap = da_params.init(sst.N, da_type="rloc")
    if sst.dap_rewrite is not None:
        dap.update_dap(sst.dap_rewrite)

    # if elem.ndim == 2:
    if dap.da_type == "rloc" and sst.N > 1:
        rloc = da_letkf.prepare_rloc(sst.ud, sst.elem, sst.node, dap, sst.N)

    logging.info("Generating initial ensemble...")
    sol_ens = np.zeros((sst.N), dtype=object)

    # Set random seed for reproducibility
    np.random.seed(params.random_seed)

    seeds = np.random.randint(10000, size=sst.N) if sst.N > 1 else None
    if seeds is not None and sst.restart == False:
        logging.info("Seeds used in generating initial ensemble spread = ", seeds)
        for n in range(sst.N):
            Sol0 = deepcopy(sst.Sol)
            mpv0 = deepcopy(sst.mpv)
            Sol0 = sol_init(Sol0, mpv0, elem, node, th, ud, seed=seeds[n])
            sol_ens[n] = [Sol0, deepcopy(flux), mpv0, [-np.inf, step]]
    elif restart == False:
        sol_ens = [[sol_init(Sol, mpv, elem, node, th, ud), flux, mpv, [-np.inf, step]]]
    elif restart == True:
        hydrostatic.state(mpv, elem, node, th, ud)
        ud.old_suffix = np.copy(ud.output_suffix)
        ud.old_suffix = "_ensemble=%i%s" % (N, ud.old_suffix)
        Sol0, mpv0, touts = io.sim_restart(
            r_params[0], r_params[1], elem, node, ud, Sol, mpv, r_params[2]
        )
        sol_ens = [[Sol0, flux, mpv0, [-np.inf, step]]]
        # ud.tout = touts[1:]
        ud.tout = [touts[-1]]
        t = touts[0]

    if ud.bdry_type[1].value == "radiation":
        ud.tcy, ud.tny = bdry.get_tau_y(ud, elem, node, 0.5)

    ens = da_utils.ensemble(sol_ens)

    ##########################################################
    # Load data assimilation observations
    ##########################################################

    # where are my observations?
    if N > 1:
        obs = dap.load_obs(dap.obs_path)
        # obs_mask, no calculations where entries are True
        obs_mask = da_utils.sparse_obs_selector(obs, elem, node, ud, dap)
        obs_noisy, obs_covar = da_utils.obs_noiser(obs, obs_mask, dap, rloc, elem)


    ##########################################################
    # Add ensemble info into filename
    ##########################################################
    if ud.autogen_fn:
        ud.output_suffix = io.fn_gen(ud, dap, N)
