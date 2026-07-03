import logging

import numpy as np

from ....flow_solver.physics import eos as gd_eos
from .comp_psinc import do_comp_to_psinc_conv, do_psinc_to_comp_conv
from .swe_lake import do_swe_to_lake_conv, do_lake_to_swe_conv


######################################################
# Blending calls from time_update.do
######################################################
def _window_start_conversion_due(bld, ud, window_step):
    """First step of a blending window with a to-limit conversion configured."""
    return (
        bld is not None
        and window_step == 0
        and (bld.bb or bld.cb)
        and ud.blending_conv is not None
    )


def _full_blend_due(bld, ud, window_step):
    """Window step at which the scheduled limit -> full conversion is due."""
    return (
        bld is not None
        and bld.criterion_init(window_step)
        and bld.cb
        and ud.blending_conv is not None
    )


def _initial_blend_phase(ud, bld, step):
    """Which leg of the initial full -> limit -> full blend this step is on.

    Returns "to_limit" on step 0, "to_full" on step ``no_of_pi_initial``,
    else None (the eos regime schedule applies instead).
    """
    if not ud.initial_blending or bld is None:
        return None
    if step < 1:
        return "to_limit"
    if step == ud.no_of_pi_initial:
        return "to_full"
    return None


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
    lake_to_swe_pending,
):
    ######################################################
    # Blending : Do full regime to limit regime conversion
    ######################################################
    if _window_start_conversion_due(bld, ud, window_step):
        if ud.blending_conv == "swe":
            do_swe_to_lake_conv(mem, ud, writer, label)
            swe_to_lake = True
        else:
            mem = do_comp_to_psinc_conv(mem, bld, ud, label, writer)

    ######################################################
    # Blending : scheduled limit regime to full regime conversion
    ######################################################
    if _full_blend_due(bld, ud, window_step) and ud.blending_conv != "swe":
        do_psinc_to_comp_conv(mem, ud, bld, label, writer, step, t + dt)

    ######################################################
    # Initial Blending
    ######################################################
    phase = _initial_blend_phase(ud, bld, step)
    if phase == "to_limit":
        if ud.blending_conv == "swe":
            do_swe_to_lake_conv(mem, ud, writer, label)
            swe_to_lake = True
            # release the lid again at the end of this same step (the
            # initial blend is one lake step); the window-driven path
            # instead re-derives this from bld.criterion_init per step
            lake_to_swe_pending = True
            ud.is_compressible = 0
            ud.compressibility = 0.0
        elif bld.psinc_init > 0:
            ud.is_compressible = 0
            ud.compressibility = 0.0
            mem = do_comp_to_psinc_conv(mem, bld, ud, label, writer)
    elif phase == "to_full":
        if ud.blending_conv != "swe":
            do_psinc_to_comp_conv(mem, ud, bld, label, writer, step, t + dt)
            ud.is_compressible = 1
            ud.compressibility = 1.0
    else:
        ud.is_compressible = gd_eos.is_compressible(ud, window_step)
        ud.compressibility = gd_eos.compressibility(ud, t, window_step)
        ud.is_nonhydrostatic = gd_eos.is_nonhydrostatic(ud, window_step)
        ud.nonhydrostasy = gd_eos.nonhydrostasy(ud, t, window_step)

    return swe_to_lake, lake_to_swe_pending


def blending_after_timestep(
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
    lake_to_swe_pending,
):
    ######################################################
    # Blending : Are we in the lake regime? And is this
    #            the window step where we go back to SWE?
    ######################################################
    if bld is not None and swe_to_lake and step > 0:
        lake_to_swe_pending = bld.criterion_init(window_step + 1)

    ######################################################
    # Blending : If we are in the lake regime, is blending
    #            on? If yes, do lake-to-swe conversion.
    ######################################################
    if (
        ud.blending_conv == "swe"
        and swe_to_lake
        and lake_to_swe_pending
        and bld is not None
    ):
        tmp_CFL = np.copy(ud.CFL)
        ud.CFL = 0.8
        do_lake_to_swe_conv(mem, ud, label, writer, step, t + dt)
        ud.CFL = tmp_CFL
        ud.is_compressible = 1
        ud.compressibility = 1.0
        lake_to_swe_pending = False

    return lake_to_swe_pending


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
    lake_to_swe_pending,
):

    if check_and_apply_initial_hydrostatic_conversion(step, ud, bld):
        ud.is_nonhydrostatic = 0

    return blending_before_timestep(
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
        lake_to_swe_pending,
    )


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
