"""Command-line argument parsing and logger setup."""

import os
import errno
import logging
import argparse
from datetime import datetime

import numpy as np
import yaml

from ...interfaces.ic_config import IC_MODULES


def get_args():
    """
    Argument parser for initial conditions and ensemble size.

    """

    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        "-N",
        action="store",
        dest="N",
        help="<Optional> Set ensemble size, if none is given N=1 is used.",
        required=False,
        type=int,
    )

    parser.add_argument(
        "-ic",
        "--initial_conditions",
        action="store",
        dest="ic",
        help="<Required> Set initial conditions",
        required=True,
        choices=set(IC_MODULES.keys()),  # Use the keys from IC_MODULES
    )

    subparsers = parser.add_subparsers(dest="subcommand")

    restart = subparsers.add_parser("restart")
    restart.add_argument(
        "-p",
        "--path",
        action="store",
        dest="path",
        help="path to data for simulation restart.",
        required=True,
        type=str,
    )
    restart.add_argument(
        "-n",
        "--name",
        action="store",
        dest="name",
        help="name of datasets for simulation restart.",
        required=True,
        type=str,
    )
    restart.add_argument(
        "-t",
        "--time",
        nargs="*",
        help="time outputs for simulation restart in format [start,stop,interval). Use None for ud.tout settings.",
        type=float,
        required=False,
        default=None,
    )

    queue = subparsers.add_parser("queue")
    queue.add_argument(
        "-w", "--rewrite", nargs="*", help="", required=True, type=yaml.safe_load
    )

    args = parser.parse_args()  # collect cmd line args
    ic = args.ic

    # Import the appropriate module
    if ic in IC_MODULES:
        module_name = IC_MODULES[ic]
        try:
            module = __import__(module_name, fromlist=["UserData", "sol_init"])
            UserData = getattr(module, "UserData")
            sol_init = getattr(module, "sol_init")
        except ImportError as e:
            raise ImportError(f"Failed to import {module_name}: {e}")
    else:
        raise ValueError(f"Unknown initial condition: {ic}")

    if UserData is None or sol_init is None:
        assert 0, "Initial condition file is not well defined."

    if args.N is None:
        N = 1
    else:
        N = args.N

    if args.subcommand == "restart":
        rstrt = True
        if args.time is not None:
            t_vals = args.time
            time = np.arange(t_vals[0], t_vals[1], t_vals[2])
        params = [args.path, args.name, time]
    else:
        rstrt = False
        params = None

    if args.subcommand == "queue":
        ud = args.rewrite[0]
        dap = args.rewrite[1]
    else:
        ud = None
        dap = None

    return N, UserData, sol_init, rstrt, ud, dap, params


def mkdir_p(path):
    """http://stackoverflow.com/a/600612/190597 (tzot)"""
    try:
        os.makedirs(path, exist_ok=True)  # Python>3.2
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def init_logger(ud):
    now = datetime.now()
    date = now.strftime("%d%m%y")
    time = now.strftime("%H%M%S")

    input_filename = "%s%s" % (ud.output_type, ud.output_base_name)
    logger_filename = "./logs/%s_%s_%s.log" % (input_filename, date, time)

    mkdir_p(os.path.dirname(logger_filename))

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%m-%d %H:%M",
        filename=logger_filename,
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

    # Suppress library specific debug outputs
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
    logging.getLogger("numba").setLevel(logging.WARNING)
    logging.getLogger("numba.core").setLevel(logging.WARNING)

    logging.getLogger().setLevel(logging.INFO)

    logging.info("Input file is %s" % input_filename)
