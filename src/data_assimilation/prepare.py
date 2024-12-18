import logging
# to generate ensemble from one sol init instantiation
from copy import deepcopy 

import numpy as np

from ..dycore.physics import hydrostatics
from ..utils import sim_params as params
from ..utils import io
from ..utils import data_structures

from . import utils as da_utils
from . import params as da_params
from . import letkf as da_letkf


def initialise(sst):
    mp = sst.model_params
    rp = sst.restart_params

    ##########################################################
    # Initialisation of data assimilation module
    ##########################################################

    # possible da_types:
    # 1) batch_obs for the LETKF with batch observations
    # 2) rloc for LETKF with grid-point localisation
    # 3) etpf for the ETPF algorithm
    dap = da_params.init(sst.N, da_type="rloc")
    if rp.dap_rewrite is not None:
        dap.update_dap(rp.dap_rewrite)

    # if elem.ndim == 2:
    if dap.da_type == "rloc" and sst.N > 1:
        rloc = da_letkf.prepare_rloc(mp.ud, mp.elem, mp.node, dap, sst.N)
    else:
        rloc = None

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
            Sol0 = sst.sol_init(Sol0, mpv0, mp.elem, mp.node, mp.th, mp.ud, seed=seeds[n])
            sol_ens[n] = [Sol0, deepcopy(mp.flux), mpv0, [-np.inf, sst.step]]
    elif sst.restart == False:
        sol_ens = [[sst.sol_init(mp.Sol, mp.mpv, mp.elem, mp.node, mp.th, sst.ud), mp.flux, mp.mpv, [-np.inf, sst.step]]]
    elif sst.restart == True:
        hydrostatics.state(mp.mpv, mp.elem, mp.node, mp.th, mp.ud)
        sst.ud.old_suffix = np.copy(sst.ud.output_suffix)
        sst.ud.old_suffix = "_ensemble=%i%s" % (sst.N, sst.ud.old_suffix)
        Sol0, mpv0, touts = io.sim_restart(
            rp.r_params[0], rp.r_params[1], mp.elem, mp.node, mp.ud, mp.Sol, mp.mpv, rp.r_params[2]
        )
        sol_ens = [[Sol0, mp.flux, mpv0, [-np.inf, sst.step]]]
        # ud.tout = touts[1:]
        sst.ud.tout = [touts[-1]]
        sst.t = touts[0]

    ens = da_utils.ensemble(sol_ens)

    ##########################################################
    # Load data assimilation observations
    ##########################################################

    # where are my observations?
    if sst.N > 1:
        obs = dap.load_obs(dap.obs_path)
        # obs_mask, no calculations where entries are True
        obs_mask = da_utils.sparse_obs_selector(obs, mp.elem, mp.node, sst.ud, dap)
        obs_noisy, obs_covar = da_utils.obs_noiser(obs, obs_mask, dap, rloc, mp.elem)
    else:
        obs, obs_noisy, obs_mask, obs_covar = None, None, None, None


    ##########################################################
    # Add ensemble info into filename
    ##########################################################
    if sst.ud.autogen_fn:
        sst.ud.output_suffix = io.fn_gen(sst.ud, dap, sst.N)

    #######################
    # Populate DA params
    #######################

    sst.da_params = data_structures.DataAssimilationParameters(
        dap=dap,
        rloc=rloc,
        sol_ens=ens,
        obs=obs,
        obs_noisy=obs_noisy,
        obs_mask=obs_mask,
        obs_covar=obs_covar
    )
