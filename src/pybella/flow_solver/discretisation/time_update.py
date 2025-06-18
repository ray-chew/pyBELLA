import copy
import logging
import numpy as np

# dependencies from pybella common
from ...utils import options as opts

# dependencies of the flow solver subpackage
from ..utils.boundary import rayleigh_boundary as bdry_r
from ..physics import cfl, eos
from ..numerics.explicit_advection import advective_flux, compute_advection
from ..numerics import explicit_euler, implicit_euler

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
        label = "%.3d" % mem.time.step

        if mem.time.step == 0 and writer != None:
            writer.write_all(mem, str(label) + "_ic")

        dt, cfl_adv, cfl_acs = cfl.dynamic_timestep(
            mem.sol, mem.time.t, tout, mem.elem, ud, mem.th, mem.time.step
        )

        dt = prestep.apply_modifcations(dt, ud, mem.time.step)

        ######################################################
        # Blending : Do blending before timestep
        ######################################################
        swe_to_lake, mem.sol, mem.npf, mem.time.t = schemes.prepare_blending(
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

        ud.is_nonhydrostatic = eos.is_nonhydrostatic(ud, mem.time.window_step)
        ud.nonhydrostasy = eos.nonhydrostasy(ud, mem.time.t, mem.time.window_step)

        if ud.continuous_blending or ud.initial_blending:
            logging.info(
                f"step = {mem.time.step}, window_step = {mem.time.window_step}"
            )

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

        advective_flux.recompute(mem)

        debug_writer.write(f"{label}_before_advect")

        if ud.do_advection:
            compute_advection.first_order_runge_kutta(
                mem,
                ud,
                0.5 * dt,
            )

        debug_writer.write(f"{label}_after_advect")
        debug_writer.populate(f"{label}_after_full_step", "p2_nodes", mem.npf.p2_nodes)

        mem.npf.p2_nodes0[...] = mem.npf.p2_nodes

        implicit_euler.do_explicit_part(mem, ud, 0.5 * dt)

        debug_writer.write(f"{label}_after_ebnaexp")

        Sol0_increment = Sol0 if ud.is_compressible == 0 else None

        implicit_euler.do_implicit_part(
            mem,
            ud,
            0.5 * dt,
            Sol0=Sol0_increment,
            label=f"{label}_after_ebnaimp",
            writer=writer,
        )

        if ud.bdry_type[1] == opts.BdryType.RAYLEIGH:
            # top rayleight damping
            bdry_r.rayleigh_damping(mem.sol, mem.npf, ud, mem.elem, mem.node)

        bdry_r.apply_rayleigh_forcing(
            mem, ud, dt
        )

        debug_writer.write(f"{label}_after_ebnaimp")

        advective_flux.recompute(mem)

        debug_writer.write(f"{label}_after_half_step")

        Sol_half_new = copy.deepcopy(mem.sol)
        npf_half_new = copy.deepcopy(mem.npf)
        mem.npf.p2_nodes_half = np.copy(mem.npf.p2_nodes)

        if ud.is_nonhydrostatic == 0 or (
            ud.is_compressible == 1 and ud.is_nonhydrostatic == 1
        ):
            mem.npf.p2_nodes[...] = mem.npf.p2_nodes0

        mem.sol = copy.deepcopy(Sol0)

        explicit_euler.do_forward_step(
            mem,
            ud,
            0.5 * dt,
            writer=writer,
            label=str(label) + "_after_efna",
        )

        debug_writer.write(f"{label}_after_efna")

        if ud.do_advection:
            compute_advection.strange_splitting(
                mem,
                ud,
                dt,
                mem.time.step % 2,
                str(label) + "_full",
                writer,
            )

        debug_writer.write(f"{label}_after_full_advect")

        implicit_euler.do_explicit_part(mem, ud, 0.5 * dt)

        debug_writer.write(f"{label}_after_full_ebnaexp")

        implicit_euler.do_implicit_part(
            mem,
            ud,
            0.5 * dt,
            writer=writer,
            label=str(label) + "_after_full_step",
        )

        if ud.bdry_type[1] == opts.BdryType.RAYLEIGH:
            # top rayleight damping
            bdry_r.rayleigh_damping(mem.sol, mem.npf, ud, mem.elem, mem.node)

        # bottom rayleigh forcing
        bdry_r.apply_rayleigh_forcing(
            mem, ud, dt,
            half=False,
            Sol_half_new=Sol_half_new,
            npf_half_new=npf_half_new,
        )

        ######################################################
        # Blending : Do blending after timestep
        ######################################################
        mem.sol, mem.npf = schemes.blending_after_timestep(
            mem.sol,
            mem.npf,
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
            % (mem.time.step, mem.time.t, dt, cfl_adv, cfl_acs)
        )
        logging.info(
            "###############################################################################################"
        )

        mem.time.t += dt
        mem.time.step += 1
        mem.time.window_step += 1

    return mem
