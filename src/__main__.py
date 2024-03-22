import numpy as np

# dependencies of the atmospheric flow solver
from dycore.discretisation.grid import grid_init
from dycore.discretisation.time_update import time_update
from dycore.utils.boundary import (
    set_explicit_boundary_data,
    set_ghostnodes_p2,
    get_tau_y,
)
from dycore.utils.variable import States, Vars
from dycore.physics.gas_dynamics.thermodynamic import ThermodynamicInit
from dycore.physics.low_mach.mpv import MPV
from dycore.physics.hydrostatics import hydrostatic_state

# dependencies of the parallelisation by dask
from dask.distributed import Client, progress

# dependencies of the data assimilation module
from data_assimilation.params import da_params
from data_assimilation.utils import (
    ensemble,
    ensemble_inflation,
    HSprojector_2t3D,
    HSprojector_3t2D,
    sparse_obs_selector,
    obs_noiser,
)
from data_assimilation.letkf import da_interface, prepare_rloc

# from data_assimilation.letkf import analysis as letkf_analysis
from data_assimilation import etpf
from data_assimilation import blending
from data_assimilation import post_processing

# input file
from utils.user_data import UserDataInit
from utils.io import io, get_args, sim_restart, fn_gen
import utils.sim_params as gparams

# some diagnostics
from copy import deepcopy
from time import time
from termcolor import colored

# test module
from tests import diagnostics as diag

import logging


logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
    datefmt="%m-%d %H:%M",
    filename="./logs/output.log",
    filemode="w",
)
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# set a format which is simpler for console use
formatter = logging.Formatter("%(name)-12s: %(levelname)-8s %(message)s")
# tell the handler to use this format
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger().addHandler(console)

debug = gparams.debug
da_debug = gparams.da_debug
output_timesteps = False
if debug == True:
    output_timesteps = True
label_type = "TIME"
np.set_printoptions(precision=gparams.print_precision)

step = 0
t = 0.0

##########################################################
# Initialisation of data containers and helper classes
##########################################################
# get arguments for initial condition and ensemble size
N, UserData, sol_init, restart, ud_rewrite, dap_rewrite, r_params = get_args()
if N == 1:
    da_debug = False

initial_data = vars(UserData())
ud = UserDataInit(**initial_data)
if ud_rewrite is not None:
    ud.update_ud(ud_rewrite)
if hasattr(ud, "rayleigh_bc"):
    ud.rayleigh_bc(ud)
if ud.output_timesteps:
    output_timesteps = True
ud.coriolis_strength = np.array(ud.coriolis_strength)

elem, node = grid_init(ud)

Sol = Vars(elem.sc, ud)

flux = np.empty((3), dtype=object)
flux[0] = States(elem.sfx, ud)
if elem.ndim > 1:
    flux[1] = States(elem.sfy, ud)
if elem.ndim > 2:
    flux[2] = States(elem.sfz, ud)

th = ThermodynamicInit(ud)
mpv = MPV(elem, node, ud)
bld = blending.Blend(ud)

logging.info("Input file is%s" % ud.output_base_name.replace("_", " "))
##########################################################
# Initialise test module
##########################################################

if ud.diag:
    diag_comparison = diag.compare_sol(ud.diag_current_run)

##########################################################
# Initialisation of data assimilation module
##########################################################

# possible da_types:
# 1) batch_obs for the LETKF with batch observations
# 2) rloc for LETKF with grid-point localisation
# 3) etpf for the ETPF algorithm
dap = da_params(N, da_type="rloc")
if dap_rewrite is not None:
    dap.update_dap(dap_rewrite)

# if elem.ndim == 2:
if dap.da_type == "rloc" and N > 1:
    rloc = prepare_rloc(ud, elem, node, dap, N)

logging.info(colored("Generating initial ensemble...", "yellow"))
sol_ens = np.zeros((N), dtype=object)

# Set random seed for reproducibility
np.random.seed(gparams.random_seed)

seeds = np.random.randint(10000, size=N) if N > 1 else None
if seeds is not None and restart == False:
    logging.info("Seeds used in generating initial ensemble spread = ", seeds)
    for n in range(N):
        Sol0 = deepcopy(Sol)
        mpv0 = deepcopy(mpv)
        Sol0 = sol_init(Sol0, mpv0, elem, node, th, ud, seed=seeds[n])
        sol_ens[n] = [Sol0, deepcopy(flux), mpv0, [-np.inf, step]]
elif restart == False:
    sol_ens = [[sol_init(Sol, mpv, elem, node, th, ud), flux, mpv, [-np.inf, step]]]
elif restart == True:
    hydrostatic_state(mpv, elem, node, th, ud)
    ud.old_suffix = np.copy(ud.output_suffix)
    ud.old_suffix = "_ensemble=%i%s" % (N, ud.old_suffix)
    Sol0, mpv0, touts = sim_restart(
        r_params[0], r_params[1], elem, node, ud, Sol, mpv, r_params[2]
    )
    sol_ens = [[Sol0, flux, mpv0, [-np.inf, step]]]
    # ud.tout = touts[1:]
    ud.tout = [touts[-1]]
    t = touts[0]

    if ud.bdry_type[1].value == "radiation":
        ud.tcy, ud.tny = get_tau_y(ud, elem, node, 0.5)

ens = ensemble(sol_ens)

##########################################################
# Load data assimilation observations
##########################################################

# where are my observations?
if N > 1:
    obs = dap.load_obs(dap.obs_path)
    # obs_mask, no calculations where entries are True
    obs_mask = sparse_obs_selector(obs, elem, node, ud, dap)
    obs_noisy, obs_covar = obs_noiser(obs, obs_mask, dap, rloc, elem)
    # obs_noisy_interp, obs_mask = sparse_obs_selector(obs_noisy, elem, node, ud, dap)


# add ensemble info to filename
if ud.autogen_fn:
    ud.output_suffix = fn_gen(ud, dap, N)
# ud.output_suffix = '_ensemble=%i%s' %(N, ud.output_suffix)

# ud.output_suffix = '%s_%s' %(ud.output_suffix, 'nr')

##########################################################
# Start main looping
##########################################################
if __name__ == "__main__":

    ######################################################
    # Initialise writer class for I/O operations
    ######################################################
    writer = io(ud, restart)
    writer.check_jar()
    writer.jar([ud, mpv, elem, node, dap])
    # sys.exit("Let's just dill the stuff and quit!")

    writer.write_attrs()
    wrtr = None
    if N > 1:
        writer.write_da_attrs(dap)
    elif output_timesteps == True:
        wrtr = writer
    for n in range(N):  # write initial ensemble
        Sol = ens.members(ens)[n][0]
        mpv = ens.members(ens)[n][2]
        if label_type == "STEP":
            label = "ensemble_mem=%i_%.3d" % (n, step)
        else:
            label = "ensemble_mem=%i_%.3f" % (n, 0.0)
        if not restart:
            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_ic")

    if da_debug:
        # writer.jar([obs,obs_noisy,obs_noisy_interp,obs_mask,obs_covar])
        # obs = obs_noisy_interp
        writer.jar([obs, obs_noisy, obs_mask, obs_covar])

    # initialise dask parallelisation and timer
    # client = Client(threads_per_worker=1, n_workers=1)
    tic = time()

    ######################################################
    # Time looping over data assimilation windows
    ######################################################
    tout_old = -np.inf
    tout_cnt = 0
    outer_step = 0
    for tout in ud.tout:
        futures = []

        # In ensemble case, do blending for each DA window
        if N > 1:
            blend = bld if tout_old in dap.da_times else None
        else:
            blend = bld

        # initial blending?
        if ud.initial_blending == True and (outer_step == 0 or outer_step == 1):
            blend = bld

        ######################################################
        # Forecast step
        ######################################################
        logging.info("##############################################")
        logging.info(colored("Next tout = %.3f" % tout, "yellow"))
        logging.info(colored("Starting forecast...", "green"))
        mem_cnt = 0
        for mem in ens.members(ens):
            # future = client.submit(time_update, *[mem[0],mem[1],mem[2], t, tout, ud, elem, node, mem[3], th, bld, None, False])

            # handling of DA window step counter
            if N > 1:
                mem[3][0] = 0 if tout_old in dap.da_times else mem[3][0]
            if N == 1:
                mem[3][0] = mem[3][1]
            logging.info(colored("For ensemble member = %i..." % mem_cnt, "yellow"))
            future = time_update(
                mem[0],
                mem[1],
                mem[2],
                t,
                tout,
                ud,
                elem,
                node,
                mem[3],
                th,
                blend,
                wrtr,
                debug,
            )

            if ud.diag:
                diag_comparison.do(
                    future[0], future[2].p2_nodes, plot=ud.diag_plot_compare
                )

            futures.append(future)
            mem_cnt += 1

        # Dask commands, used only when parallelisation is
        # enabled
        # results = client.gather(futures)
        results = np.copy(futures)
        results = np.array(results)
        # s_res = client.scatter(results)

        ######################################################
        # Analysis step
        ######################################################
        tout = np.around(tout, 3)
        if N > 1 and tout in dap.da_times:
            futures = []

            ######################################################
            # Update ensemble with forecast
            ######################################################
            for n in range(N):
                Sol = results[n][dap.loc_c]
                set_explicit_boundary_data(Sol, elem, ud, th, mpv)
                results[n][dap.loc_c] = Sol
                p2_nodes = getattr(results[n][dap.loc_n], "p2_nodes")
                set_ghostnodes_p2(p2_nodes, node, ud)
                setattr(results[n][dap.loc_n], "p2_nodes", p2_nodes)

            ens.set_members(results, tout)

            ######################################################
            # Write output before assimilating data
            ######################################################
            logging.info(colored("Starting output...", "yellow"))
            for n in range(N):
                Sol = ens.members(ens)[n][0]
                mpv = ens.members(ens)[n][2]

                if label_type == "STEP":
                    step = outer_step
                    label = "ensemble_mem=%i_%.3d" % (n, step)
                else:
                    label = "ensemble_mem=%i_%.3f" % (n, tout)
                writer.write_all(Sol, mpv, elem, node, th, str(label) + "_before_da")

            ##################################################
            # LETKF with batch observations
            ##################################################
            if dap.da_type == "batch_obs":
                logging.info("Starting analysis... for batch observations")
                for attr in dap.obs_attributes:
                    logging.info("Assimilating %s..." % attr)
                    logging.info("Assimilating %s..." % attr)
                    # future = client.submit(da_interface, *[s_res,obs_current,dap.inflation_factor,attr,N,ud,dap.loc[attr]])
                    future = da_interface(results, dap, obs, attr, tout, N, ud)
                    futures.append(future)

                # analysis = client.gather(futures)
                analysis = futures
                # analysis = np.array(analysis)

                logging.info("Writing analysis...")
                cnt = 0
                for attr in dap.obs_attributes:
                    current = analysis[cnt]
                    for n in range(N):
                        setattr(results[:, dap.loc[attr], ...][n], attr, current[n])
                    cnt += 1

            ##################################################
            # LETKF with grid-point localisation
            ##################################################
            elif dap.da_type == "rloc":
                logging.info(
                    colored("Starting analysis... for rloc algorithm", "green")
                )
                results = HSprojector_3t2D(results, elem, dap, N)
                results = rloc.analyse(results, obs, obs_covar, obs_mask, N, tout)
                results = HSprojector_2t3D(results, elem, node, dap, N)
                # if hasattr(dap, 'converter'):
                # results = dap.converter(results, N, mpv, elem, node, th, ud)

            ##################################################
            # ETPF
            ##################################################
            elif dap.da_type == "etpf":
                ensemble_inflation(results, dap.attributes, dap.inflation_factor, N)
                results = etpf.da_interface(
                    results,
                    obs,
                    dap.obs_attributes,
                    dap.rejuvenation_factor,
                    dap.da_times,
                    tout,
                    N,
                )

            ##################################################
            # Post-processing
            ##################################################
            elif dap.da_type == "pprocess":
                results = post_processing.interface()

            else:
                assert 0, "DA type not implemented: use 'rloc', 'batch_obs' or 'etpf'."

        ######################################################
        # Update ensemble with analysis
        ######################################################
        for n in range(N):
            Sol = results[n][dap.loc_c]
            set_explicit_boundary_data(Sol, elem, ud, th, mpv)
            results[n][dap.loc_c] = Sol
            p2_nodes = getattr(results[n][dap.loc_n], "p2_nodes")
            set_ghostnodes_p2(p2_nodes, node, ud)
            setattr(results[n][dap.loc_n], "p2_nodes", p2_nodes)

        ens.set_members(results, tout)

        ######################################################
        # Write output at tout
        ######################################################
        logging.info(colored("Starting output...", "yellow"))
        for n in range(N):
            Sol = ens.members(ens)[n][0]
            mpv = ens.members(ens)[n][2]

            if label_type == "STEP":
                step = outer_step
                label = "ensemble_mem=%i_%.3d" % (n, step)
            else:
                label = "ensemble_mem=%i_%.3f" % (n, tout)
            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_after_full_step")

        # synchronise_variables(mpv, Sol, elem, node, ud, th)
        t = tout
        tout_old = np.copy(tout)
        logging.info(colored("tout = %.3f" % tout, "yellow"))

        tout_cnt += 1
        outer_step += 1
        if outer_step > ud.stepmax:
            break

    toc = time()
    logging.info(colored("Time taken = %.6f" % (toc - tic), "yellow"))

    writer.close_everything()