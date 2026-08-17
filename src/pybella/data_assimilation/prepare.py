import logging

import numpy as np

from ..utils import sim_params as params
from ..utils import io
from ..utils import data_structures
from ..flow_solver.utils import fields
from ..flow_solver.utils.boundary import cell_boundary as bdry_c

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

    elem, node = es.get_grid()
    member0 = es.get_member(0)

    # needed by the rloc analysis, and by obs_noiser for every da_type (its
    # cell/node attribute partition sizes the observation covariance arrays)
    if sst.N > 1:
        rloc = da_letkf.prepare_rloc(sst.ud, elem, node, dap, sst.N)
    else:
        rloc = None

    ##########################################################
    # Generate the initial ensemble from member 0
    ##########################################################

    # Set random seed for reproducibility
    np.random.seed(params.random_seed)

    seeds = np.random.randint(10000, size=sst.N) if sst.N > 1 else None

    if sst.N > 1:
        logging.info("Generating initial ensemble...")
        logging.info("Seeds used in generating initial ensemble spread = %s", seeds)
        sol_ens = data_structures.EnsembleState()
        for n in range(sst.N):
            # members are built from fresh containers so that sol_init runs
            # exactly once per member; member 0 was already initialised by the
            # flow-solver prepare and re-running sol_init on it would double
            # the += initialisations.
            sol0 = fields.CellSolField(elem.sc)
            npf0 = fields.NodePressureField(elem, node, sst.ud)
            sol0 = sst.sol_init(
                sol0, npf0, elem, node, member0.th, sst.ud, seed=seeds[n]
            )
            # cache=None: each member gets its own FlowSolverCache
            sol_ens.update_member(elem, node, sol0, npf0, member0.th)

        for member in sol_ens.members:
            bdry_c.set_ghost_cells(member, sst.ud)

        sst.ensemble_state = sol_ens

    ##########################################################
    # Load data assimilation observations
    ##########################################################

    # da_times may be empty (e.g. an ensemble forecast without assimilation);
    # then no observation file is needed.
    if sst.N > 1 and len(dap.da_times) > 0:
        obs = dap.load_obs(dap.obs_path)
        # obs_mask, no calculations where entries are True
        obs_mask = da_utils.sparse_obs_selector(obs, elem, node, sst.ud, dap)
        obs_noisy, obs_covar = da_utils.obs_noiser(obs, obs_mask, dap, rloc, elem)
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
        obs=obs,
        obs_noisy=obs_noisy,
        obs_mask=obs_mask,
        obs_covar=obs_covar,
    )
