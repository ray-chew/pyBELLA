"""Device-resident host loop (`ud.backend = "jax-device"`).

The window driver split out of `device_kernels`: host-side CFL/dt control,
the support guard, func-mode forcing evaluation, the per-window compile
cache, and :func:`run_window` itself. The traced substeps and the compiled
``step`` live in `device_kernels`; see `device_step` for the overall design
notes."""

import logging

import numpy as np
import jax
import jax.numpy as jnp

from pybella.flow_solver.physics import cfl as cfl_np

from .device_config import build_device_config
from .device_kernels import make_step
from .device_state import to_device, write_back


@jax.jit
def _cfl_maxima_plain(rho, rhou, rhov, rhow, rhoY, gamm, Msq):
    p = rhoY**gamm
    c = jnp.sqrt(gamm * p / rho) / jnp.sqrt(Msq)
    u = jnp.abs(rhou / rho)
    v = jnp.abs(rhov / rho)
    w = jnp.abs(rhow / rho)
    return jnp.stack(
        [
            u.max(),
            v.max(),
            w.max(),
            (u + c).max(),
            (v + c).max(),
            (w + c).max(),
        ]
    )


def _cfl_maxima(s, cfg):
    if not cfg.terrain:
        return _cfl_maxima_plain(
            s["rho"], s["rhou"], s["rhov"], s["rhow"], s["rhoY"], cfg.gamm, cfg.Msq
        )
    m = cfg.metric
    rho, rhoY = s["rho"], s["rhoY"]
    p = rhoY**cfg.gamm
    c = jnp.sqrt(cfg.gamm * p / rho) / jnp.sqrt(cfg.Msq)
    u = jnp.abs(s["rhou"] / rho)
    v = jnp.abs(s["rhov"] / rho)
    w = jnp.abs(s["rhow"] / rho)
    moms = (s["rhou"], s["rhov"], s["rhow"])
    # contravariant speeds |N_a . m|/(rho J) and signal bounds c |N_a|/J on
    # every sweep axis, mirroring the numpy cfl.dynamic_timestep
    vels = [u, v, w]
    cs = [c, c, c]
    cv, (ch1, ch2) = m.cart_v, m.cart_haxes
    for a in range(cfg.ndim):
        Na = m.N[a]
        contra = Na[cv] * moms[cv] + Na[ch1] * moms[ch1]
        norm_sq = Na[cv] ** 2 + Na[ch1] ** 2
        if ch2 is not None:
            contra = contra + Na[ch2] * moms[ch2]
            norm_sq = norm_sq + Na[ch2] ** 2
        vels[a] = jnp.abs(contra / rho) * m.ooJ
        cs[a] = c * jnp.sqrt(norm_sq) * m.ooJ
    u, v, w = vels
    return jnp.stack(
        [
            u.max(),
            v.max(),
            w.max(),
            (u + cs[0]).max(),
            (v + cs[1]).max(),
            (w + cs[2]).max(),
        ]
    )


def _host_dt(maxima, mem, ud, tout):
    eps = np.finfo(float).eps
    u_max, v_max, w_max, upc, vpc, wpc = (
        max(float(x), eps) for x in np.asarray(maxima)
    )
    return cfl_np._calculate_advective_timestep(
        ud.CFL,
        mem.elem,
        u_max,
        v_max,
        w_max,
        upc,
        vpc,
        wpc,
        mem.time.t,
        tout,
        ud,
        mem.time.step,
        eps,
    )


def _check_supported(mem, ud, writer):
    problems = []
    blending = getattr(ud, "continuous_blending", False) or getattr(
        ud, "initial_blending", False
    )
    # comp<->psinc blending is supported via window segmentation (Phase D4:
    # host conversions at the known blend steps, one compiled step variant
    # per regime); the SWE/lake conversions and the 'imbal' initial
    # HYDROSTATIC conversion (alpha_w = 0 step path) are not implemented
    if blending and getattr(ud, "blending_conv", None) == "swe":
        problems.append("SWE/lake blending conversions")
    if blending and "imbal" in getattr(ud, "aux", ""):
        problems.append("initial hydrostatic conversion ('imbal' aux)")
    if getattr(ud, "is_ArakawaKonor", 0):
        problems.append("Arakawa-Konor")
    if getattr(ud, "acoustic_timestep", 0) == 1:
        problems.append("acoustic timestep")
    if (
        getattr(ud, "rayleigh_forcing", False)
        and getattr(ud, "rayleigh_forcing_type", "func") == "file"
    ):
        problems.append("file-based rayleigh forcing")
    from pybella.utils import sim_params

    if getattr(sim_params, "debug", False):
        problems.append("debug writers (sim_params.debug)")
    if problems:
        raise NotImplementedError(
            "backend='jax-device' does not support: "
            + ", ".join(problems)
            + " — run with backend='jax' (hybrid) instead"
        )


def _eval_forcing(mem, ud, t_offset):
    """Host-side func-mode forcing arrays (eigenfunction at host-known t)."""
    s_par = 5.0e-3 + 1e-4 + 0e-5
    ud.rf_bot.eigenfunction(t_offset, s_par)
    up, vp, Yp, _ = ud.rf_bot.dehatter(mem.th)
    ud.rf_bot.eigenfunction(t_offset, s_par, grid="n")
    _, _, _, pi_n = ud.rf_bot.dehatter(mem.th, grid="n")
    return (
        jnp.asarray(up),
        jnp.asarray(vp),
        jnp.asarray(Yp),
        jnp.asarray(pi_n),
    )


_DUMMY_FORCING = (0.0, 0.0, 0.0, 0.0)

_WINDOW_CACHE = {}


def _get_window_cache(mem, ud):
    """Config + compiled step functions, persistent across output windows
    (each time_update.do call) for the same grid/ud."""
    key = (id(mem.elem), id(ud))
    entry = _WINDOW_CACHE.get(key)
    if entry is None or entry[0] is not mem.elem or entry[1] is not ud:
        cfg = build_device_config(mem, ud)
        _WINDOW_CACHE[key] = (mem.elem, ud, cfg, {})
        entry = _WINDOW_CACHE[key]
    return entry[2], entry[3]


def _blend_conversion_due(bld, ud, step, window_step):
    """True iff prepare_blending will touch the STATE this step (a comp<->
    psinc conversion) — the device must round-trip through host memory."""
    from pybella.interfaces.dynamics_blending.schemes import orchestration

    return (
        orchestration._window_start_conversion_due(bld, ud, window_step)
        or orchestration._full_blend_due(bld, ud, window_step)
        or orchestration._initial_blend_phase(ud, bld, step) is not None
    )


def run_window(mem, ud, tout, bld=None, writer=None):
    """Device-resident replacement for time_update.do's step loop.

    Blending (Phase D4): the window is segmented at the blend steps — the
    per-step regime control is the SAME `schemes.prepare_blending` the numpy
    loop runs (it is pure host bookkeeping except on conversion steps, where
    the state round-trips through host memory and the comp<->psinc routines
    run unchanged); each regime compiles its own step variant (the regime
    ints are static trace structure)."""
    from pybella.flow_solver.physics import eos
    from pybella.interfaces.dynamics_blending import schemes
    from pybella.interfaces.time_stepper import prestep

    _check_supported(mem, ud, writer)

    # regime fields must exist before the config build (its use_cross probe
    # evaluates the coriolis coefficients); the loop refreshes them per step
    ud.is_compressible = eos.is_compressible(ud, mem.time.window_step)
    ud.compressibility = eos.compressibility(ud, mem.time.t, mem.time.window_step)
    ud.is_nonhydrostatic = eos.is_nonhydrostatic(ud, mem.time.window_step)
    ud.nonhydrostasy = eos.nonhydrostasy(ud, mem.time.t, mem.time.window_step)

    cfg, step_fns = _get_window_cache(mem, ud)

    s = to_device(mem)
    compile_count = 0
    if writer is not None:
        logging.info(
            "jax-device: per-step output writer active — the state is "
            "pulled to host every step (disable output_timesteps for "
            "device-resident performance)"
        )

    while (mem.time.t < tout) and (mem.time.step < ud.stepmax):
        label = "%.3d" % mem.time.step
        if mem.time.step == 0 and writer is not None:
            writer.write_all(mem, str(label) + "_ic")

        maxima = _cfl_maxima(s, cfg)
        dt, cfl_adv, cfl_acs = _host_dt(maxima, mem, ud, tout)
        dt = prestep.apply_modifcations(dt, ud, mem.time.step)

        # the SAME per-step blending/regime control as the numpy loop: on
        # conversion steps the state round-trips through host memory and the
        # comp<->psinc routines run unchanged; otherwise prepare_blending is
        # the plain eos regime refresh
        conversion = _blend_conversion_due(bld, ud, mem.time.step, mem.time.window_step)
        if conversion:
            write_back(s, mem)
        schemes.prepare_blending(
            mem, ud, bld, label, writer, mem.time.step,
            mem.time.window_step, mem.time.t, dt, False, False,
        )
        if conversion:
            s = to_device(mem)
        ud.is_nonhydrostatic = eos.is_nonhydrostatic(ud, mem.time.window_step)
        ud.nonhydrostasy = eos.nonhydrostasy(ud, mem.time.t, mem.time.window_step)
        if int(ud.is_nonhydrostatic) == 0:
            raise NotImplementedError(
                "hydrostatic regime (alpha_w = 0) step path is not "
                "implemented on jax-device — run with backend='jax' (hybrid)"
            )

        parity = mem.time.step % 2
        key = (parity, int(ud.is_nonhydrostatic), int(ud.is_compressible))
        if key not in step_fns:
            step_fns[key] = make_step(
                cfg, parity, int(ud.is_nonhydrostatic), int(ud.is_compressible)
            )
            compile_count += 1

        if cfg.has_forcing:
            forcing_half = _eval_forcing(mem, ud, mem.time.t + 0.5 * dt)
            forcing_full = _eval_forcing(mem, ud, mem.time.t + dt)
        else:
            forcing_half = forcing_full = _DUMMY_FORCING

        out = step_fns[key](
            s,
            dt,
            float(ud.nonhydrostasy),
            float(ud.compressibility),
            forcing_half,
            forcing_full,
        )
        if key[2] == 0:
            # psinc step: keep the predictor half-time pressure current on
            # the host container (numpy sets npf.p2_nodes_half every step;
            # the psinc->comp blend conversion reads it)
            s, p2_half = out
            mem.npf.p2_nodes_half = np.asarray(p2_half)
        else:
            s = out

        if writer is not None:
            write_back(s, mem)
            writer.time = mem.time.t
            writer.write_all(mem, str(label) + "_after_full_step")

        logging.info(
            "device step %i done, t = %.12f, dt = %.12f, CFL = %.8f, CFL_ac = %.8f",
            mem.time.step,
            mem.time.t,
            dt,
            cfl_adv,
            cfl_acs,
        )
        mem.time.t += dt
        mem.time.step += 1
        mem.time.window_step += 1

    write_back(s, mem)
    mem._device_compile_count = compile_count
    return mem
