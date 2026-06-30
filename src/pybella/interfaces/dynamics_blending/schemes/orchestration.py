import logging

import numpy as np

from ....flow_solver.physics import eos as gd_eos
from .comp_psinc import do_comp_to_psinc_conv, do_psinc_to_comp_conv
from .swe_lake import do_swe_to_lake_conv, do_lake_to_swe_conv
from .hydro_nonhydro import do_nonhydro_to_hydro_conv, do_hydro_to_nonhydro_conv


######################################################
# Blending calls from data.py
######################################################
def blending_before_timestep(
    mem,
    ud,
    bld,
    label,
    writer,
    step,
    window_step,
    t,
    dt,
    swe_to_lake,
    debug,
):
    ######################################################
    # Blending : Do full regime to limit regime conversion
    ######################################################
    # do unpacking
    elem, node, sol, npf, th, _, _ = mem

    # these make sure that we are the correct window step
    if bld is not None and window_step == 0:
        # these make sure that blending switches are on
        if (bld.bb or bld.cb) and ud.blending_conv is not None:
            # these distinguish between SWE and Euler blending
            if ud.blending_conv == "swe":
                do_swe_to_lake_conv(sol, npf, elem, node, ud, th, writer, label, debug)
                swe_to_lake = True
            else:
                mem = do_comp_to_psinc_conv(mem, bld, ud, label, writer)

    ######################################################
    # Blending : Do full steps or transition steps?
    ######################################################
    if bld is not None:
        c_init = bld.criterion_init(window_step)
    else:
        c_init = False

    ######################################################
    # Blending : If full blending steps...
    ######################################################
    # check that blending switches are on
    if c_init and bld.cb and ud.blending_conv is not None:
        # distinguish between Euler and SWE blending
        if ud.blending_conv != "swe":
            do_psinc_to_comp_conv(
                mem,
                ud,
                bld,
                ud,
                label,
                writer,
                step,
                t + dt,
            )

    ######################################################
    # Initial Blending
    ######################################################
    # Is initial blending switch on, and if yes, are we in the 0th time-step?
    if ud.initial_blending == True and step < 1 and bld is not None:
        # Distinguish between SWE and Euler blendings
        if ud.blending_conv != "swe":
            if bld.psinc_init > 0:
                ud.is_compressible = 0
                ud.compressibility = 0.0
                mem = do_comp_to_psinc_conv(mem, bld, ud, label, writer)
            elif bld.hydro_init > 0:
                sol, npf, t = do_nonhydro_to_hydro_conv(
                    sol,
                    flux,
                    npf,
                    bld,
                    elem,
                    node,
                    th,
                    ud,
                    label,
                    writer,
                    step,
                    window_step,
                    t,
                    dt,
                )
                ud.is_nonhydrostatic = 0
                ud.nonhydrostasy = 0.0
        else:
            do_swe_to_lake_conv(sol, npf, elem, node, ud, th, writer, label, debug)
            swe_to_lake = True
            ud.is_compressible = 0
            ud.compressibility = 0.0

    # Elif, is initial blending switch on and are we on the 1st time-step?
    # If we are on the first time-step, do we do comp-psinc blending?
    elif (
        ud.initial_blending == True and step == ud.no_of_pi_initial and bld is not None
    ):
        # Distinguish between SWE and Euler blendings
        if ud.blending_conv != "swe":
            do_psinc_to_comp_conv(
                mem,
                ud,
                bld,
                label,
                writer,
                step,
                t + dt,
            )
            ud.is_compressible = 1
            ud.compressibility = 1.0
    # Else, do we do nonhydrostatic-hydrostatic blending?
    elif (
        ud.initial_blending == True and step == ud.no_of_hy_initial and bld is not None
    ):
        if ud.blending_conv != "swe":
            sol, npf = do_hydro_to_nonhydro_conv(
                sol,
                flux,
                npf,
                bld,
                elem,
                node,
                th,
                ud,
                label,
                writer,
                step,
                window_step,
                t,
                dt,
            )
            ud.is_nonhydrostatic = 1
            ud.nonhydrostasy = 1.0
    else:
        ud.is_compressible = gd_eos.is_compressible(ud, window_step)
        ud.compressibility = gd_eos.compressibility(ud, t, window_step)
        ud.is_nonhydrostatic = gd_eos.is_nonhydrostatic(ud, window_step)
        ud.nonhydrostasy = gd_eos.nonhydrostasy(ud, t, window_step)

    return swe_to_lake, sol, npf, t


def blending_after_timestep(
    sol,
    npf,
    bld,
    elem,
    node,
    th,
    ud,
    label,
    writer,
    step,
    window_step,
    t,
    dt,
    swe_to_lake,
    debug,
):
    ######################################################
    # Blending : Are we in the lake regime? And is this
    #            the window step where we go back to SWE?
    ######################################################
    if bld is not None and swe_to_lake and step > 0:
        initialise_lake_to_swe_conv = bld.criterion_init(window_step + 1)

    ######################################################
    # Blending : If we are in the lake regime, is blending
    #            on? If yes, do lake-to-swe conversion.
    ######################################################
    if (
        ud.blending_conv == "swe"
        and swe_to_lake
        and initialise_lake_to_swe_conv
        and bld is not None
    ):
        tmp_CFL = np.copy(ud.CFL)
        ud.CFL = 0.8
        sol, npf = do_lake_to_swe_conv(
            sol,
            flux,
            npf,
            elem,
            node,
            ud,
            th,
            writer,
            label,
            debug,
            step,
            window_step,
            t,
            dt,
        )
        ud.CFL = tmp_CFL
        ud.is_compressible = 1
        ud.compressibility = 1.0

    return sol, npf


def prepare_blending(
    mem,
    ud,
    bld,
    label,
    writer,
    step,
    window_step,
    t,
    dt,
    swe_to_lake,
    debug,
):

    if check_and_apply_initial_hydrostatic_conversion(step, ud, bld):
        ud.is_nonhydrostatic = 0

    swe_to_lake, sol, npf, t = blending_before_timestep(
        mem,
        ud,
        bld,
        label,
        writer,
        step,
        window_step,
        t,
        dt,
        swe_to_lake,
        debug,
    )

    return swe_to_lake, sol, npf, t


def check_and_apply_initial_hydrostatic_conversion(step, ud, bld):
    """
    Check and apply initial blending conversion if needed.

    Returns:
        bool: True if conversion was applied, False otherwise
    """
    if step != 0 or bld is None or "imbal" not in ud.aux:
        return False

    hydrostatic_case = ud.is_nonhydrostatic == 0
    nonhydrostatic_case = ud.is_nonhydrostatic == 1 and ud.initial_blending == True

    if hydrostatic_case or nonhydrostatic_case:
        logging.info("nonhydrostatic to hydrostatic conversion...")
        ud.is_nonhydrostatic = 0
        return True

    return False
