import copy
import logging

import numpy as np

from ...utils import options as opts

# dependencies of the flow solver subpackage
from ..utils import boundary as bdry
from ..physics.gas_dynamics import (
    numerical_flux as gd_flux,
    eos as gd_eos,
    cfl as gd_cfl,
)
from ..physics.gas_dynamics import explicit as gd_explicit
from ..physics.low_mach import second_projection as lm_sp

# for blending module
from ...interfaces.dynamics_blending import schemes
from ...interfaces.time_stepper import prestep


def do(
    mem,
    ud,
    tout,
    bld=None,
    writer=None,
    debug_writer=None,
):
    """
    For more details, refer to the write-up :ref:`time-stepping`.

    Does a time-step for the flow solver.

    """
    swe_to_lake = False

    while (mem.time.t < tout) and (mem.time.step < ud.stepmax):
        bdry.set_explicit_boundary_data(mem.sol, mem.elem, ud, mem.th, mem.mpv)

        label = "%.3d" % mem.time.step

        if mem.time.step == 0 and writer != None:
            writer.write_all(mem, str(label) + "_ic")

        dt, cfl, cfl_ac = gd_cfl.dynamic_timestep(mem.sol, mem.time.t, tout, mem.elem, ud, mem.th, mem.time.step)

        dt = prestep.apply_modifcations(dt, ud, mem.time.step)

        ######################################################
        # Blending : Do blending before timestep
        ######################################################
        swe_to_lake, mem.sol, mem.mpv, mem.time.t = schemes.prepare_blending(
            mem,
            ud,
            bld,
            label,
            writer,
            mem.time.step,
            mem.time.window_step,
            mem.time.t,
            dt,
            swe_to_lake,
            debug_writer,
        )

        ud.is_nonhydrostatic = gd_eos.is_nonhydrostatic(ud, mem.time.window_step)
        ud.nonhydrostasy = gd_eos.nonhydrostasy(ud, mem.time.t, mem.time.window_step)

        if ud.continuous_blending or ud.initial_blending:
            logging.info(f"step = {mem.time.step}, window_step = {mem.time.window_step}")

        logging.info(
            f"""
                    -------
                    is_compressible = {ud.is_compressible}, is_nonhydrostatic = {ud.is_nonhydrostatic}
                    compressibility = {ud.compressibility:.3f}, nonhydrostasy = {ud.nonhydrostasy:.3f}
                    -------
                    """
        )

        Sol0 = copy.deepcopy(mem.sol)

        debug_writer.write(f"{label}_before_flux")

        gd_flux.recompute_advective_fluxes(mem)

        debug_writer.populate_flux_components(f"{label}_before_advect", mem.flux, mem.elem)
        debug_writer.write(f"{label}_before_advect")

        if ud.do_advection:
            gd_explicit.advect_rk(
                mem,
                ud,
                0.5 * dt,
            )

        debug_writer.write(f"{label}_after_advect")
        debug_writer.populate(f"{label}_after_full_step", "p2_nodes", mem.mpv.p2_nodes)

        mem.mpv.p2_nodes0[...] = mem.mpv.p2_nodes

        lm_sp.euler_backward_non_advective_expl_part(mem, ud, 0.5 * dt)

        debug_writer.write(f"{label}_after_ebnaexp")

        Sol0_increment = Sol0 if ud.is_compressible == 0 else None

        lm_sp.euler_backward_non_advective_impl_part(
            mem.sol,
            mem.mpv,
            mem.elem,
            mem.node,
            ud,
            mem.th,
            mem.time.t,
            0.5 * dt,
            mem,
            Sol0=Sol0_increment,
            label=f"{label}_after_ebnaimp",
            writer=writer,
        )

        if ud.bdry_type[1] == opts.BdryType.RAYLEIGH:
            # top rayleight damping
            bdry.rayleigh_damping(mem.sol, mem.mpv, ud, mem.elem, mem.node)

        bdry.apply_rayleigh_forcing(
            mem.sol,
            mem.mpv,
            ud,
            mem.elem,
            mem.node,
            mem.time.t,
            mem.time.step,
            dt,
            mem.th,
            bdry,
        )

        debug_writer.write(f"{label}_after_ebnaimp")

        gd_flux.recompute_advective_fluxes(mem)

        debug_writer.populate_flux_components(f"{label}_after_half_step", mem.flux, mem.elem)
        debug_writer.write(f"{label}_after_half_step")

        Sol_half_new = copy.deepcopy(mem.sol)
        mpv_half_new = copy.deepcopy(mem.mpv)
        mem.mpv.p2_nodes_half = np.copy(mem.mpv.p2_nodes)

        if ud.is_nonhydrostatic == 0 or (
            ud.is_compressible == 1 and ud.is_nonhydrostatic == 1
        ):
            mem.mpv.p2_nodes[...] = mem.mpv.p2_nodes0

        mem.sol = copy.deepcopy(Sol0)

        lm_sp.euler_forward_non_advective(
            mem,
            ud,
            0.5 * dt,
            writer=writer,
            label=str(label) + "_after_efna",
        )

        debug_writer.write(f"{label}_after_efna")

        if ud.do_advection:
            gd_explicit.advect(
                mem,
                ud,
                dt,
                mem.time.step % 2,
                str(label) + "_full",
                writer,
            )

        debug_writer.write(f"{label}_after_full_advect")

        lm_sp.euler_backward_non_advective_expl_part(mem, ud, 0.5 * dt)

        debug_writer.write(f"{label}_after_full_ebnaexp")

        lm_sp.euler_backward_non_advective_impl_part(
            mem.sol,
            mem.mpv,
            mem.elem,
            mem.node,
            ud,
            mem.th,
            mem.time.t,
            0.5 * dt,
            mem,
            writer=writer,
            label=str(label) + "_after_full_step",
        )

        if ud.bdry_type[1] == opts.BdryType.RAYLEIGH:
            # top rayleight damping
            bdry.rayleigh_damping(mem.sol, mem.mpv, ud, mem.elem, mem.node)

        # bottom rayleigh forcing
        bdry.apply_rayleigh_forcing(
            mem.sol,
            mem.mpv,
            ud,
            mem.elem,
            mem.node,
            mem.time.t,
            mem.time.step,
            dt,
            mem.th,
            bdry,
            half=False,
            Sol_half_new=Sol_half_new,
            mpv_half_new=mpv_half_new,
        )

        ######################################################
        # Blending : Do blending after timestep
        ######################################################
        mem.sol, mem.mpv = schemes.blending_after_timestep(
            mem.sol,
            mem.flux,
            mem.mpv,
            bld,
            mem.elem,
            mem.node,
            mem.th,
            ud,
            label,
            writer,
            mem.time.step,
            mem.time.window_step,
            mem.time.t,
            dt,
            swe_to_lake,
            debug_writer,
        )

        if writer != None:
            writer.time = mem.time.t
            writer.write_all(mem, str(label) + "_after_full_step")

        logging.info(
            "###############################################################################################"
        )
        logging.info(
            "step %i done, t = %.12f, dt = %.12f, CFL = %.8f, CFL_ac = %.8f"
            % (mem.time.step, mem.time.t, dt, cfl, cfl_ac)
        )
        logging.info(
            "###############################################################################################"
        )

        mem.time.t += dt
        mem.time.step += 1
        mem.time.window_step += 1

    return mem
