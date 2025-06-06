import copy
import logging

import numpy as np

# dependencies of the pyBELLA package
from ...utils import io

# dependencies of the flow solver subpackage
from . import grid as dis_grid
from ..utils import boundary as bdry, options as opts
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


def data_init(ud):
    """
    Helper function to initialise the `elem` and `node` grids, corresponding to the cell and node grids, from a given user iniital data file.

    Parameters
    ----------
    ud : :class:`inputs.user_data.UserDataInit`
        Data container for the initial conditions.

    Returns
    -------
    elem : :class:`discretization.kgrid.ElemSpaceDiscr`
        Cells grid.
    node : :class:`discretization.kgrid.NodeSpaceDiscr`
        Nodes grid.

    """
    inx = ud.inx
    iny = ud.iny
    inz = ud.inz
    x0 = ud.xmin
    x1 = ud.xmax
    y0 = ud.ymin
    y1 = ud.ymax
    z0 = ud.zmin
    z1 = ud.zmax

    grid = dis_grid.Grid(inx, iny, inz, x0, x1, y0, y1, z0, z1)

    elem = dis_grid.ElemSpaceDiscr(grid, ud)
    node = dis_grid.NodeSpaceDiscr(grid, ud)

    return elem, node


def do(
    sst,
    mem,
    tout,
    bld=None,
    writer=None,
    debug_writer=None,
):
    """
    For more details, refer to the write-up :ref:`time-stepping`.

    Does a time-step for the atmospheric solver.

    Parameters
    ----------
    Sol : :class:`management.variable.Vars`
        Solution data container.
    flux : :class:`management.variable.States`
        Data container for the fluxes.
    mpv : :class:`physics.low_mach.mpv.MPV`
        Variables relating to the elliptic solver.
    t : float
        Current time
    tout : float
        Next output time
    ud : :class:`inputs.user_data.UserDataInit`
        Data container for the initial conditions
    elem : :class:`discretization.kgrid.ElemSpaceDiscr`
        Cells grid.
    node : :class:`discretization.kgrid.NodeSpaceDiscr`
        Nodes grid.
    step : int
        Current step.
    th : :class:`physics.gas_dynamics.thermodynamic.init`
        Thermodynamic variables of the system
    bld : :class:`data_assimilation.blending.Blend()`
        Blending class used to initalise interface blending methods.
    writer : :class:`management.io.io`, optional
        `default == None`. If given, output after each time-step will be written in the hdf5 format.
    debug : boolean, optional
        `default == False`. If `True`, then writer will output `Sol`:
            1. before flux calculation
            2. before advection routine
            3. after advection routine
            4. after explicit solver
            5. after implicit solver

        during both the half-step for the prediction of advective flux and the full-step.

    Returns
    -------
    list
        A list of `[Sol,flux,mpv,[window_step,step]]` data containers at time `tout`.
    """
    ud = sst.ud
    elem, node, Sol, flux, mpv, th, time = mem

    window_step = time.window_step
    step = time.step
    t = time.t
    swe_to_lake = False

    while (t < tout) and (step < ud.stepmax):
        bdry.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

        label = "%.3d" % step

        if step == 0 and writer != None:
            writer.write_all(mem, str(label) + "_ic")

        dt, cfl, cfl_ac = gd_cfl.dynamic_timestep(Sol, t, tout, elem, ud, th, step)

        dt = prestep.apply_modifcations(dt, ud, step)

        ######################################################
        # Blending : Do blending before timestep
        ######################################################
        swe_to_lake, Sol, mpv, t = schemes.prepare_blending(
            sst,
            mem,
            bld,
            label,
            writer,
            step,
            window_step,
            t,
            dt,
            swe_to_lake,
            debug_writer,
        )

        ud.is_nonhydrostatic = gd_eos.is_nonhydrostatic(ud, window_step)
        ud.nonhydrostasy = gd_eos.nonhydrostasy(ud, t, window_step)

        if ud.continuous_blending or ud.initial_blending:
            logging.info(f"step = {step}, window_step = {window_step}")

        logging.info(
            f"""
                    -------
                    is_compressible = {ud.is_compressible}, is_nonhydrostatic = {ud.is_nonhydrostatic}
                    compressibility = {ud.compressibility:.3f}, nonhydrostasy = {ud.nonhydrostasy:.3f}
                    -------
                    """
            )

        Sol0 = copy.deepcopy(Sol)

        debug_writer.write(f"{label}_before_flux")

        gd_flux.recompute_advective_fluxes(flux, Sol)

        debug_writer.populate_flux_components(f"{label}_before_advect", flux, elem)
        debug_writer.write(f"{label}_before_advect")

        if ud.do_advection:
            gd_explicit.advect_rk(
                Sol,
                flux,
                0.5 * dt,
                elem,
                step % 2,
                ud,
                th,
                mpv,
                node,
                str(label) + "_half",
                writer,
            )

        debug_writer.write(f"{label}_after_advect")
        debug_writer.populate(f"{label}_after_full_step", "p2_nodes", mpv.p2_nodes)

        mpv.p2_nodes0[...] = mpv.p2_nodes

        lm_sp.euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5 * dt, ud, th)

        debug_writer.write(f"{label}_after_ebnaexp")

        Sol0_increment = Sol0 if ud.is_compressible == 0 else None

        lm_sp.euler_backward_non_advective_impl_part(
            Sol,
            mpv,
            elem,
            node,
            ud,
            th,
            t,
            0.5 * dt,
            1.0,
            Sol0=Sol0_increment,
            label=f"{label}_after_ebnaimp",
            writer=writer,
        )

        if ud.bdry_type[1] == opts.BdryType.RAYLEIGH:
            # top rayleight damping
            bdry.rayleigh_damping(Sol, mpv, ud, elem, node)

        bdry.apply_rayleigh_forcing(
            Sol,
            mpv,
            ud,
            elem,
            node,
            t,
            step,
            dt,
            th,
            bdry,
        )

        debug_writer.write(f"{label}_after_ebnaimp")

        gd_flux.recompute_advective_fluxes(flux, Sol)

        debug_writer.populate_flux_components(f"{label}_after_half_step", flux, elem)
        debug_writer.write(f"{label}_after_half_step")

        Sol_half_new = copy.deepcopy(Sol)
        mpv_half_new = copy.deepcopy(mpv)
        mpv.p2_nodes_half = np.copy(mpv.p2_nodes)

        if ud.is_nonhydrostatic == 0 or (
            ud.is_compressible == 1 and ud.is_nonhydrostatic == 1
        ):
            mpv.p2_nodes[...] = mpv.p2_nodes0

        Sol = copy.deepcopy(Sol0)

        lm_sp.euler_forward_non_advective(
            Sol,
            mpv,
            elem,
            node,
            0.5 * dt,
            ud,
            th,
            writer=writer,
            label=str(label) + "_after_efna",
        )

        debug_writer.write(f"{label}_after_efna")

        if ud.do_advection:
            gd_explicit.advect(
                Sol,
                flux,
                dt,
                elem,
                step % 2,
                ud,
                th,
                mpv,
                node,
                str(label) + "_full",
                writer,
            )

        debug_writer.write(f"{label}_after_full_advect")

        lm_sp.euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5 * dt, ud, th)

        debug_writer.write(f"{label}_after_full_ebnaexp")

        lm_sp.euler_backward_non_advective_impl_part(
            Sol,
            mpv,
            elem,
            node,
            ud,
            th,
            t,
            0.5 * dt,
            2.0,
            writer=writer,
            label=str(label) + "_after_full_step",
        )

        if ud.bdry_type[1] == opts.BdryType.RAYLEIGH:
            # top rayleight damping
            bdry.rayleigh_damping(Sol, mpv, ud, elem, node)

        # bottom rayleigh forcing
        bdry.apply_rayleigh_forcing(
            Sol,
            mpv,
            ud,
            elem,
            node,
            t,
            step,
            dt,
            th,
            bdry,
            half=False,
            Sol_half_new=Sol_half_new,
            mpv_half_new=mpv_half_new,
        )

        ######################################################
        # Blending : Do blending after timestep
        ######################################################
        Sol, mpv = schemes.blending_after_timestep(
            Sol,
            flux,
            mpv,
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
            debug_writer,
        )

        mem.sol = Sol
        mem.flux = flux
        mem.mpv = mpv

        if writer != None:
            writer.time = t
            writer.write_all(mem, str(label) + "_after_full_step")

        logging.info(
            "###############################################################################################"
        )
        logging.info(
            "step %i done, t = %.12f, dt = %.12f, CFL = %.8f, CFL_ac = %.8f"
            % (step, t, dt, cfl, cfl_ac)
        )
        logging.info(
            "###############################################################################################"
        )

        t += dt
        step += 1
        window_step += 1

    mem.time.t = t
    mem.time.step = step
    mem.time.window_step = window_step

    return mem
    # return [Sol, flux, mpv, [window_step, step]]
