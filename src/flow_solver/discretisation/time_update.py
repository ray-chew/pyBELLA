import copy
import logging

import numpy as np
import termcolor

# dependencies of the pyBELLA package
from ...utils import io

# dependencies of the flow solver subpackage
from . import grid as dis_grid
from ..utils import (
    boundary as bdry, 
    options as opts
)
from ..physics.gas_dynamics import (
    numerical_flux as gd_flux,
    eos as gd_eos,
    cfl as gd_cfl
)
from ..physics.gas_dynamics import explicit as gd_explicit
from ..physics.low_mach import second_projection as lm_sp

# for blending module
from ...interfaces.dynamics_blending import schemes


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
    debug=False,
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
            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_ic")

        dt, cfl, cfl_ac = gd_cfl.dynamic_timestep(Sol, t, tout, elem, ud, th, step)

        if "CFLfixed" in ud.aux:
            if step < 2:
                dt = 21.69 / ud.t_ref

        c1 = (
            step == 0
            and ud.is_nonhydrostatic == 0
            and bld is not None
            and "imbal" in ud.aux
        )
        c2 = (
            step == 0
            and ud.is_nonhydrostatic == 1
            and ud.initial_blending == True
            and bld is not None
            and "imbal" in ud.aux
        )

        if c1 or c2:
            logging.info("nonhydrostatic to hydrostatic conversion...")
            ud.is_nonhydrostatic = 0
            if test_hydrob == False:
                dt *= 0.5


        ######################################################
        # Blending : Do blending before timestep
        ######################################################
        swe_to_lake, Sol, mpv, t = schemes.blending_before_timestep(
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
            debug,
        )

        ud.is_nonhydrostatic = gd_eos.is_nonhydrostatic(ud, window_step)
        ud.nonhydrostasy = gd_eos.nonhydrostasy(ud, t, window_step)
        
        if ud.continuous_blending or ud.initial_blending:
            logging.info(f"step = {step}, window_step = {window_step}")

        logging.info(f"""
                    -------
                    is_compressible = {ud.is_compressible}, is_nonhydrostatic = {ud.is_nonhydrostatic}
                    compressibility = {ud.compressibility:.3f}, nonhydrostasy = {ud.nonhydrostasy:.3f}
                    -------
                    """)

        Sol0 = copy.deepcopy(Sol)

        if debug == True:
            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_before_flux")

        gd_flux.recompute_advective_fluxes(flux, Sol)

        if debug:
            writer.populate(f"{label}_before_advect", "rhoYu", flux[0].rhoY)
            writer.populate(f"{label}_before_advect", "rhoYv", flux[1].rhoY)
            if elem.ndim == 3:
                writer.populate(f"{label}_before_advect", "rhoYw", flux[2].rhoY)
            writer.write_all(Sol, mpv, elem, node, th, f"{label}_before_advect")


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

        if debug:
            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_after_advect")
            writer.populate(str(label) + "_after_full_step", "p2_nodes0", mpv.p2_nodes)

        mpv.p2_nodes0[...] = mpv.p2_nodes

        lm_sp.euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5 * dt, ud, th)

        if debug == True:
            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_after_ebnaexp")

        Sol0 = Sol0 if ud.is_compressible == 0 else None

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
            Sol0=Sol0,
            label=f"{label}_after_ebnaimp",
            writer=writer,
        )

        if ud.bdry_type[1] == opts.BdryType.RAYLEIGH:
            # top rayleight damping
            bdry.rayleigh_damping(Sol, mpv, ud, elem, node)

        # bottom rayleigh forcing
        if hasattr(ud, "rayleigh_forcing"):
            if ud.rayleigh_forcing:

                if ud.rayleigh_forcing_type == "file":
                    reader = io.read_input(
                        ud.rayleigh_forcing_fn, ud.rayleigh_forcing_path
                    )

                    Sol_half_new = copy.deepcopy(Sol)
                    mpv_half_new = copy.deepcopy(mpv)

                    time_tag = "%.3d_after_full_step" % step
                    reader.get_data(Sol_half_new, mpv_half_new, time_tag, half=True)

                    # assuming constant background state
                    up = Sol_half_new.rhou / Sol_half_new.rho
                    vp = Sol_half_new.rhov / Sol_half_new.rho
                    Yp = (
                        Sol_half_new.rhoY / Sol_half_new.rho
                        - mpv.HydroState.Y0.reshape(1, -1)
                    )

                    pi = mpv_half_new.p2_nodes

                    bdry.rayleigh_damping(
                        Sol, mpv, ud, elem, node, [up, vp, Yp, pi, t + 0.5 * dt]
                    )

                elif ud.rayleigh_forcing_type == "func":
                    # boundary.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

                    s = 5.0e-3 + 1e-4 + 0e-5
                    ud.rf_bot.eigenfunction((t + 0.5 * dt), s)
                    up, vp, Yp, pi = ud.rf_bot.dehatter(th)

                    ud.rf_bot.eigenfunction((t + 0.5 * dt), s, grid="n")
                    _, _, _, pi_n = ud.rf_bot.dehatter(th, grid="n")

                    bdry.rayleigh_damping(
                        Sol, mpv, ud, elem, node, [up, vp, Yp, pi_n, t + 0.5 * dt]
                    )
                bdry.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

        if debug == True:
            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_after_ebnaimp")


        flux_half_new = copy.deepcopy(flux)

        gd_flux.recompute_advective_fluxes(flux, Sol)

        if debug:
            writer.populate(f"{label}_after_half_step", "rhoYu", flux[0].rhoY)
            writer.populate(f"{label}_after_half_step", "rhoYv", flux[1].rhoY)
            if elem.ndim == 3:
                writer.populate(f"{label}_after_half_step", "rhoYw", flux[2].rhoY)
            writer.write_all(Sol, mpv, elem, node, th, f"{label}_after_half_step")

        Sol_half_new = copy.deepcopy(Sol)
        mpv_half_new = copy.deepcopy(mpv)
        rho_half = np.copy(Sol.rho)
        rhou_half = np.copy(Sol.rhou)
        rhov_half = np.copy(Sol.rhov)
        rhow_half = np.copy(Sol.rhow)
        rhoX_half = np.copy(Sol.rhoX)
        rhoY_half = np.copy(Sol.rhoY)
        # pwchi = np.copy(Sol.pwchi)
        p2_nodes_half = np.copy(mpv.p2_nodes)


        if ud.is_nonhydrostatic == 0 or (ud.is_compressible == 1 and ud.is_nonhydrostatic == 1):
            mpv.p2_nodes[...] = mpv.p2_nodes0

        # Sol = copy.deepcopy(Sol0)

        # Sol.rhov0 = np.copy(Sol.rhov)
        Sol.rho_half = rho_half
        Sol.rhou_half = rhou_half
        Sol.rhov_half = rhov_half
        Sol.rhow_half = rhow_half
        Sol.rhoX_half = rhoX_half
        Sol.rhoY_half = rhoY_half
        mpv.p2_nodes_half = p2_nodes_half
        # Sol.pwchi = pwchi

        # Sol.rho_half = Sol.rho
        # Sol.rhou_half = Sol.rhou
        # Sol.rhov_half = Sol.rhov
        # Sol.rhow_half = Sol.rhow
        # Sol.rhoX_half = Sol.rhoX
        # Sol.rhoY_half = Sol.rhoY
        # mpv.p2_nodes_half = mpv.p2_nodes
        Sol_half_old = copy.deepcopy(Sol_half_new)
        flux_half_old = copy.deepcopy(flux_half_new)
        mpv_half_old = copy.deepcopy(mpv_half_new)

        # mpv.p2_nodes[...] = ud.compressibility * mpv.p2_nodes0 + (1.0-ud.compressibility) * mpv.p2_nodes

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

        if debug == True:
            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_after_efna")

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

        if debug == True:
            writer.write_all(
                Sol, mpv, elem, node, th, str(label) + "_after_full_advect"
            )

        lm_sp.euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5 * dt, ud, th)

        if debug == True:
            writer.write_all(
                Sol, mpv, elem, node, th, str(label) + "_after_full_ebnaexp"
            )

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
        if hasattr(ud, "rayleigh_forcing"):
            if ud.rayleigh_forcing:

                if ud.rayleigh_forcing_type == "file":

                    reader = io.read_input(
                        ud.rayleigh_forcing_fn, ud.rayleigh_forcing_path
                    )

                    # misusing hydrostatic blending data containers
                    time_tag = "%.3d_after_full_step" % step
                    reader.get_data(Sol_half_new, mpv_half_new, time_tag)

                    # assuming constant background state
                    up = Sol_half_new.rhou / Sol_half_new.rho
                    vp = Sol_half_new.rhov / Sol_half_new.rho
                    Yp = (
                        Sol_half_new.rhoY / Sol_half_new.rho
                        - mpv.HydroState.Y0.reshape(1, -1)
                    )
                    # vp = 0.0
                    # Yp = 0.0
                    pi = mpv_half_new.p2_nodes

                    bdry.rayleigh_damping(
                        Sol, mpv, ud, elem, node, [up, vp, Yp, pi, t + dt]
                    )

                elif ud.rayleigh_forcing_type == "func":

                    s = 5.0e-3 + 1e-4 + 0e-5
                    ud.rf_bot.eigenfunction((t + dt), s)
                    up, vp, Yp, pi = ud.rf_bot.dehatter(th)

                    ud.rf_bot.eigenfunction((t + dt), s, grid="n")
                    _, _, _, pi_n = ud.rf_bot.dehatter(th, grid="n")

                    bdry.rayleigh_damping(
                        Sol, mpv, ud, elem, node, [up, vp, Yp, pi_n, t + dt]
                    )
                bdry.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

        # if writer is not None: writer.populate(str(label)+'_after_full_step','p2_half',mpv.p2_nodes_half)

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
            debug,
        )

        if c1 or c2:
            logging.info(termcolor.colored("hydrostatic to nonhydrostatic conversion...", "blue"))

            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_half_full")
            writer.populate(str(label) + "_ic", "pwchi", Sol.pwchi)

            if test_hydrob == False:
                Sol = copy.deepcopy(Sol_half_old)
                # mpv = copy.deepcopy(mpv_half_old)

                logging.info(termcolor.colored("test_hydrob == False", "red"))
                writer.write_all(Sol, mpv, elem, node, th, str(label) + "_quarter")

                writer.populate(str(label) + "_quarter", "pwchi", Sol.pwchi)

                logging.info("quarter dt = %.8f" % (dt * 0.5))

                ret = do(
                    Sol_half_old,
                    flux_half_old,
                    mpv_half_old,
                    dt - 0.5 * dt,
                    dt + 0.5 * dt,
                    ud,
                    elem,
                    node,
                    [0, 0],
                    th,
                    bld=None,
                    writer=None,
                    debug=False,
                )

                Sol_tu = copy.deepcopy(ret[0])
                # mpv_tu = copy.deepcopy(ret[2])
                Sol.rho[...] = Sol_tu.rho_half
                Sol.rhou[...] = Sol_tu.rhou_half
                Sol.rhov[...] = Sol_tu.rhov_half
                Sol.rhow[...] = Sol_tu.rhow_half
                Sol.rhoX[...] = Sol_tu.rhoX_half
                Sol.rhoY[...] = Sol_tu.rhoY_half
                Sol.pwchi[...] = Sol_tu.pwchi

                # mpv.p2_nodes[...] = mpv_tu.p2_nodes_half

                writer.write_all(Sol, mpv, elem, node, th, str(label) + "_half")

                writer.populate(str(label) + "_half", "pwchi", Sol.pwchi)

                ret = do(
                    Sol,
                    flux,
                    mpv,
                    dt,
                    2.0 * dt,
                    ud,
                    elem,
                    node,
                    [0, 0],
                    th,
                    bld=None,
                    writer=None,
                    debug=False,
                )

                Sol = copy.deepcopy(ret[0])
                flux = copy.deepcopy(ret[1])
                mpv = copy.deepcopy(ret[2])

            if test_hydrob == True:
                Sol = copy.deepcopy(Sol_half_old)
                # mpv = copy.deepcopy(mpv_half_old)

                logging.info(termcolor.colored("test_hydrob == False", "red"))
                writer.write_all(Sol, mpv, elem, node, th, str(label) + "_quarter")

                # writer.populate(str(label)+'_quarter', 'pwchi', Sol.pwchi)

                logging.info("quarter dt = %.8f" % (dt * 0.5))

                ret = do(
                    Sol_half_old,
                    flux_half_old,
                    mpv_half_old,
                    dt - 0.5 * dt,
                    dt + 0.5 * dt,
                    ud,
                    elem,
                    node,
                    [0, 0],
                    th,
                    bld=None,
                    writer=None,
                    debug=False,
                )

                Sol_tu = copy.deepcopy(ret[0])
                # mpv_tu = copy.deepcopy(ret[2])
                Sol.rho[...] = Sol_tu.rho_half
                Sol.rhou[...] = Sol_tu.rhou_half
                Sol.rhov[...] = Sol_tu.rhov_half
                Sol.rhow[...] = Sol_tu.rhow_half
                Sol.rhoX[...] = Sol_tu.rhoX_half
                Sol.rhoY[...] = Sol_tu.rhoY_half
                Sol.pwchi[...] = Sol_tu.pwchi

                # mpv.p2_nodes[...] = mpv_tu.p2_nodes_half

                # writer.write_all(Sol,mpv,elem,node,th,str(label)+'_half')

                # writer.populate(str(label)+'_half', 'pwchi', Sol.pwchi)

                ret = do(
                    Sol,
                    flux,
                    mpv,
                    dt,
                    2.0 * dt,
                    ud,
                    elem,
                    node,
                    [0, 0],
                    th,
                    bld=None,
                    writer=None,
                    debug=False,
                )

                Sol = copy.deepcopy(ret[0])
                flux = copy.deepcopy(ret[1])
                mpv = copy.deepcopy(ret[2])
                # writer.write_all(Sol,mpv,elem,node,th,str(label)+'_half')
                # writer.populate(str(label)+'_half', 'pwchi', Sol.pwchi)

                logging.info(termcolor.colored("test_hydrob == True", "red"))

            if test_hydrob == False:
                dt *= 2.0
            if c2:
                ud.is_nonhydrostatic = 1

        t += dt

        if writer != None:
            writer.time = t
            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_after_full_step")
            # writer.populate(str(label)+'_after_full_step', 'pwchi', Sol.pwchi)
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

        step += 1
        window_step += 1

        time.step = step
        time.window_step = window_step

    return mem
    # return [Sol, flux, mpv, [window_step, step]]
