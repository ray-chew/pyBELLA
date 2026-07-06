import copy
import logging
import numpy as np

# dependencies from pybella common
from ...utils import axes
from ...utils import options as opts

# dependencies of the flow solver subpackage
from ..utils.boundary import rayleigh_boundary as bdry_r
from ..physics import cfl, eos, surface_constraint
from ..numerics.explicit_advection import advective_flux, compute_advection
from ..numerics import diffusion, explicit_euler, implicit_euler

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
    from ...backends import is_device_backend

    if is_device_backend(ud):
        # device-resident inner loop: state pushed once, one dt scalar sync
        # per step, full state pulled at tout (backends/jax_ops/device_step);
        # blend windows segment at the conversion steps (Phase D4)
        from ...backends.jax_ops import device_step

        return device_step.run_window(mem, ud, tout, bld=bld, writer=writer)

    swe_to_lake = False
    lake_to_swe_pending = False

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
        # the conversion routines mutate/rebind mem.sol and mem.npf in place;
        # assigning pre-conversion aliases back here would discard the
        # psinc -> comp conversion (the pre-refactor code threaded the
        # CONVERTED Sol/mpv through, cf. 7f0b676~1 schemes.py)
        swe_to_lake, lake_to_swe_pending = schemes.prepare_blending(
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
            lake_to_swe_pending,
        )

        ud.is_nonhydrostatic = eos.is_nonhydrostatic(ud, mem.time.window_step)
        ud.nonhydrostasy = eos.nonhydrostasy(ud, mem.time.t, mem.time.window_step)

        if ud.continuous_blending or ud.initial_blending:
            logging.info(
                f"step = {mem.time.step}, window_step = {mem.time.window_step}"
            )

        logging.info(f"""
                    -------
                    is_compressible = {ud.is_compressible}, is_nonhydrostatic = {ud.is_nonhydrostatic}
                    compressibility = {ud.compressibility:.3f}, nonhydrostasy = {ud.nonhydrostasy:.3f}
                    -------
                    """)

        if ud.is_nonhydrostatic == 0:
            # Hydrostatic regime (alpha_w = 0): the full thesis sec. 4.2.3-4.2.4
            # one-step blend, realised on the SI-midpoint predictor/corrector.
            # Fig. 4.2: over n -> n+1/2 run, in parallel, the two first-order
            # updates (solid arrows, for a balanced Pw^{n+1/2} independent of the
            # input imbalance) AND a second-order pi update (dotted arrow, from the
            # n+1/4 midpoint); hand the balanced Pw and the second-order pi to the
            # second-order corrector n+1/2 -> n+1.
            sol0_n, ps1_half, ps1_npf = predictor_half_step(
                mem, ud, 0.5 * dt, writer, debug_writer, label
            )  # n -> n+1/4 (first-order); fluxes now at n+1/4
            n14_sol = copy.deepcopy(mem.sol)
            n14_npf = copy.deepcopy(mem.npf)
            # innovation 2 (sec. 4.2.4): second-order pi update over [n, n+1/2]
            # using the n+1/4 midpoint -> pi^{n+1/2}_2nd
            corrector_full_step(
                mem,
                ud,
                0.5 * dt,
                (sol0_n, ps1_half, ps1_npf),
                writer,
                debug_writer,
                label,
            )
            pi_2nd = np.copy(mem.npf.p2_nodes)
            # innovation 1 (sec. 4.2.3): restore n+1/4, second first-order update
            # -> balanced Pw^{n+1/2}
            mem.sol = n14_sol
            mem.npf = n14_npf
            predictor_half_step(mem, ud, 0.5 * dt, writer, debug_writer, label)
            # combine: balanced Pw^{n+1/2} (mem) + second-order pi^{n+1/2}
            mem.npf.p2_nodes[...] = pi_2nd
            corrector_full_step(
                mem,
                ud,
                dt,
                (sol0_n, copy.deepcopy(mem.sol), copy.deepcopy(mem.npf)),
                writer,
                debug_writer,
                label,
            )
        else:
            # Nonhydrostatic (alpha_w = 1): predictor (first-order, n -> n+1/2)
            # then corrector (second-order, n+1/2 -> n+1). Bit-identical to the
            # thesis-era stepper (RKLM_Python data.py).
            predictor_state = predictor_half_step(
                mem, ud, dt, writer, debug_writer, label
            )
            corrector_full_step(
                mem, ud, dt, predictor_state, writer, debug_writer, label
            )

        ######################################################
        # Blending : Do blending after timestep
        ######################################################
        lake_to_swe_pending = schemes.blending_after_timestep(
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
            lake_to_swe_pending,
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


def predictor_half_step(mem, ud, dt, writer=None, debug_writer=None, label=""):
    """First-order predictor: advance ``mem`` in place from t^n to the half-time
    state t^{n+1/2}.

    This is also the hydrostatic-balance recovery step (thesis eq. 4.35): in the
    hydrostatic regime (``alpha_w = 0``) the implicit solve here reconstructs a
    balanced Exner pressure and the diagnostic vertical momentum from the other
    quantities, provided the hydrostatic background is stratification-consistent
    (see ``dev_notes/hydrostatic_blending.md`` "ROOT CAUSE"). The hydro regime in
    :func:`do` invokes this twice per step (thesis sec. 4.2.3 two first-order
    updates) plus a second-order ``pi`` update (sec. 4.2.4), not a single
    predictor+corrector; the eos schedule flips ``alpha_w`` and no separate
    conversion routine is needed.

    On return ``mem`` is at t^{n+1/2} and carries the half-time advective fluxes
    the corrector's Strang sweep consumes. Returns the auxiliary snapshots the
    corrector needs: ``(sol0, sol_half_new, npf_half_new)`` — the t^n cell state
    and the half-time cell/node-pressure copies (for the second-order pass and
    the bottom Rayleigh forcing). ``debug_writer`` must be a real writer (use
    ``io.NullDebugWriter()`` when driving this outside ``do``).
    """
    sol0 = copy.deepcopy(mem.sol)

    debug_writer.write(f"{label}_before_flux")

    advective_flux.recompute(mem, ud)

    debug_writer.write(f"{label}_before_advect")

    if ud.do_advection:
        compute_advection.first_order_runge_kutta(
            mem,
            ud,
            0.5 * dt,
        )
    surface_constraint.apply(mem, ud)

    debug_writer.write(f"{label}_after_advect")
    debug_writer.populate(f"{label}_after_full_step", "p2_nodes", mem.npf.p2_nodes)

    mem.npf.p2_nodes0[...] = mem.npf.p2_nodes

    implicit_euler.do_explicit_part(mem, ud, 0.5 * dt)
    surface_constraint.apply(mem, ud)

    debug_writer.write(f"{label}_after_ebnaexp")

    sol0_increment = sol0 if ud.is_compressible == 0 else None

    implicit_euler.do_implicit_part(
        mem,
        ud,
        0.5 * dt,
        sol0=sol0_increment,
        label=f"{label}_after_ebnaimp",
        writer=writer,
    )

    if ud.bdry_type[axes.vertical_axis(ud)] == opts.BdryType.RAYLEIGH:
        # top rayleight damping
        bdry_r.rayleigh_damping(mem.sol, mem.npf, ud)

    bdry_r.apply_rayleigh_forcing(mem, ud, dt)
    surface_constraint.apply(mem, ud)

    debug_writer.write(f"{label}_after_ebnaimp")

    advective_flux.recompute(mem, ud)

    debug_writer.write(f"{label}_after_half_step")

    sol_half_new = copy.deepcopy(mem.sol)
    npf_half_new = copy.deepcopy(mem.npf)
    mem.npf.p2_nodes_half = np.copy(mem.npf.p2_nodes)

    return sol0, sol_half_new, npf_half_new


def corrector_full_step(
    mem, ud, dt, predictor_state, writer=None, debug_writer=None, label=""
):
    """Second-order corrector: from the half-time state (``mem`` carries the
    half-time fluxes) and the ``predictor_state`` snapshots from
    :func:`predictor_half_step`, advance ``mem`` in place to t^{n+1}.
    """
    sol0, sol_half_new, npf_half_new = predictor_state

    # Seed the corrector pressure from p2_nodes0 (the start-of-step Exner
    # pressure) for the compressible-nonhydrostatic path. The hydrostatic
    # (alpha_w = 0) case is intentionally EXCLUDED: with the Phase H1a
    # operator fix the hydrostatic predictor reconstructs a balanced Exner
    # pressure, and resetting it here would discard that reconstruction
    # (re-seeding from the imbalanced p2_nodes0). Keeping it lets the
    # hydrostatic regime carry the recovered pressure into the corrector.
    # Bit-identical for alpha_w = 1. See dev_notes/hydrostatic_blending.md.
    if ud.is_compressible == 1 and ud.is_nonhydrostatic == 1:
        mem.npf.p2_nodes[...] = mem.npf.p2_nodes0

    mem.sol = copy.deepcopy(sol0)

    explicit_euler.do_forward_step(
        mem,
        ud,
        0.5 * dt,
        writer=writer,
        label=str(label) + "_after_efna",
    )

    surface_constraint.apply(mem, ud)

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

    surface_constraint.apply(mem, ud)

    debug_writer.write(f"{label}_after_full_advect")

    implicit_euler.do_explicit_part(mem, ud, 0.5 * dt)
    surface_constraint.apply(mem, ud)

    debug_writer.write(f"{label}_after_full_ebnaexp")

    implicit_euler.do_implicit_part(
        mem,
        ud,
        0.5 * dt,
        writer=writer,
        label=str(label) + "_after_full_step",
    )

    if ud.bdry_type[axes.vertical_axis(ud)] == opts.BdryType.RAYLEIGH:
        # top rayleight damping
        bdry_r.rayleigh_damping(mem.sol, mem.npf, ud)

    # bottom rayleigh forcing
    bdry_r.apply_rayleigh_forcing(
        mem,
        ud,
        dt,
        half=False,
        sol_half_new=sol_half_new,
        npf_half_new=npf_half_new,
    )
    surface_constraint.apply(mem, ud)

    if ud.diffusion:
        diffusion.apply(mem, ud, dt)
