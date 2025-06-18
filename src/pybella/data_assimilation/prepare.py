import logging

# to generate ensemble from one sol init instantiation
from copy import deepcopy

import numpy as np

from ..utils import sim_params as params
from ..utils import io
from ..utils import data_structures

from . import utils as da_utils
from . import params as da_params
from . import letkf as da_letkf


def initialise(sst):
    es = sst.ensemble_state
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
        rloc = da_letkf.prepare_rloc(es.ud, es.elem, es.node, dap, sst.N)
    else:
        rloc = None

    logging.info("Generating initial ensemble...")
    sol_ens = data_structures.EnsembleState()

    # Set random seed for reproducibility
    np.random.seed(params.random_seed)

    seeds = np.random.randint(10000, size=sst.N) if sst.N > 1 else None

    if sst.N > 1:
        logging.info("Seeds used in generating initial ensemble spread = ", seeds)
        for n in range(sst.N):
            Sol0 = deepcopy(sst.Sol)
            npf0 = deepcopy(sst.npf)
            Sol0 = sst.sol_init(
                Sol0, npf0, es.elem, es.node, es.th, es.ud, seed=seeds[n]
            )
            # sol_ens[n] = [Sol0, deepcopy(es.flux), npf0, [-np.inf, es.step]]
            sol_ens.update_member(
                es.elem, es.node, Sol0, npf0, deepcopy(es.flux), es.th
            )

            sst.ensembble_state = sol_ens
    # elif sst.restart == False:
    # sol_ens = [[sst.sol_init(mp.Sol, mp.npf, mp.elem, mp.node, mp.th, sst.ud), mp.flux, mp.npf, [-np.inf, sst.step]]]
    # sol_ens.update_member(mp.elem, mp.node, sst.sol_init(mp.Sol, mp.npf, mp.elem, mp.node, mp.th, sst.ud), mp.npf, deepcopy(mp.flux), mp.th)
    # for n in range(sst.N):
    #     sol_ens.get_member(n).time.t = -np.inf

    # ens = da_utils.ensemble(sol_ens)

    ##########################################################
    # Load data assimilation observations
    ##########################################################

    # where are my observations?
    if sst.N > 1:
        obs = dap.load_obs(dap.obs_path)
        # obs_mask, no calculations where entries are True
        obs_mask = da_utils.sparse_obs_selector(obs, es.elem, es.node, sst.ud, dap)
        obs_noisy, obs_covar = da_utils.obs_noiser(obs, obs_mask, dap, rloc, es.elem)
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
        # sol_ens=ens,
        obs=obs,
        obs_noisy=obs_noisy,
        obs_mask=obs_mask,
        obs_covar=obs_covar,
    )
