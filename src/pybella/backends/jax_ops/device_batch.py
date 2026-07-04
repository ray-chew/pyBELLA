"""Ensemble-batched device window driver (`ens_run_strategy: batch` + JAX).

One assimilation-window integration advances ALL ensemble members in
lockstep on the device: the batch state is the same 7-leaf dict pytree as
`device_state` with a leading member axis, and the full compiled step —
bicgstab included — is vmapped over that axis (`lax.while_loop` under vmap
select-freezes converged lanes, so the solve keeps per-member convergence
with a batch-max iteration count). Design rationale and the dt-policy
decision live in dev_notes/nedas_interface.md (Phase D2).

dt policy: batch-min, host-controlled. Per step one vmapped CFL-maxima
kernel yields (K, 6) maxima in a single D2H pull; the host computes every
member's dt through the same `_host_dt` path as the single-member loop and
the batch advances with the min. This keeps every member inside its CFL
bound and the whole batch on one step/parity/tout-clip schedule; it is NOT
bitwise-reproducible against the per-member-dt scheduler path (statistical
gate only — the cross-member dt spread is ~3e-4 relative on the TV OSSE).

``mode="loop"`` runs the identical driver — same batch-min dt sequence —
with the unbatched compiled step per member: the vmap gate comparator
(differences are then pure vmap/jit reordering).

Blending stays host-side and unsupported here exactly as in `device_loop`
(`_check_supported`); the regime scalars are batch-uniform by construction
(shared ``ud``, lockstep ``window_step``).
"""

import logging

import numpy as np
import jax
import jax.numpy as jnp

from pybella.flow_solver.physics import eos

from .device_config import _SOL_FIELDS
from .device_kernels import build_step, make_step
from .device_loop import (
    _DUMMY_FORCING,
    _blend_conversion_due,
    _cfl_maxima_plain,
    _check_supported,
    _get_window_cache,
    _host_dt,
)
from .device_state import to_device, write_back

_STATE_FIELDS = _SOL_FIELDS + ("p2_nodes",)

_batch_cfl_maxima = jax.jit(
    jax.vmap(_cfl_maxima_plain, in_axes=(0, 0, 0, 0, 0, None, None))
)


def make_batch_step(cfg, parity, is_nonhydrostatic, is_compressible):
    """Vmapped-over-members twin of `device_kernels.make_step`.

    The member axis is axis 0 of every state leaf; dt, the regime scalars
    and the forcing tuples are batch-uniform (in_axes None).
    """
    step = build_step(cfg, parity, is_nonhydrostatic, is_compressible)
    return jax.jit(
        jax.vmap(step, in_axes=(0, None, None, None, None, None)),
        donate_argnums=(0,),
    )


def _assert_lockstep(mems):
    t0, s0, w0 = mems[0].time.t, mems[0].time.step, mems[0].time.window_step
    for k, mem in enumerate(mems):
        if (mem.time.t, mem.time.step, mem.time.window_step) != (t0, s0, w0):
            raise ValueError(
                f"batch members not in lockstep at window entry: member {k} "
                f"at (t={mem.time.t}, step={mem.time.step}, "
                f"window_step={mem.time.window_step}) vs member 0 at "
                f"(t={t0}, step={s0}, window_step={w0})"
            )


def _stack(states):
    return {name: jnp.stack([s[name] for s in states]) for name in _STATE_FIELDS}


def _shard_members(sb, n_devices):
    """Lay the member axis out over the first `n_devices` GPUs.

    One process, one NamedSharding over a 1-D device mesh; the jitted
    vmapped step runs on the sharded batch directly and its outputs stay
    sharded (jit specialises per input sharding — no cache-key change).
    """
    devices = jax.devices()[:n_devices]
    if len(devices) < n_devices:
        raise ValueError(
            f"ens_devices={n_devices} but only {len(devices)} JAX devices visible"
        )
    mesh = jax.sharding.Mesh(np.array(devices), ("mem",))
    sharding = jax.sharding.NamedSharding(mesh, jax.sharding.PartitionSpec("mem"))
    return {name: jax.device_put(v, sharding) for name, v in sb.items()}


def _unstack(sb, k):
    return {name: sb[name][k] for name in _STATE_FIELDS}


def run_window_batch(mems, ud, tout, mode="vmap", n_devices=1, bld=None):
    """Advance all members of `mems` (list of ModelState) to `tout` in
    lockstep on the device; the batch twin of `device_loop.run_window`.

    `n_devices > 1` (vmap mode only) shards the member axis across GPUs
    (`len(mems)` must divide evenly). Blending windows segment at the
    conversion steps exactly as in `run_window` — the schedule is
    batch-uniform (lockstep step/window_step, shared ud), so all members
    convert on host together and re-enter the device in the new regime.
    No per-step writer support: the batch path is for DA-member forecasts
    (NEDAS windows), which never write per-step output.
    """
    from pybella.interfaces.dynamics_blending import schemes
    from pybella.interfaces.time_stepper import prestep

    if mode not in ("vmap", "loop"):
        raise ValueError(f"unknown batch mode {mode!r} (expected vmap|loop)")
    if n_devices > 1:
        if mode != "vmap":
            raise ValueError("member-axis sharding requires ens_batch_mode 'vmap'")
        if len(mems) % n_devices:
            raise ValueError(
                f"nens={len(mems)} does not divide over ens_devices={n_devices}"
            )
    mem0 = mems[0]
    _check_supported(mem0, ud, None)
    if getattr(ud, "rayleigh_forcing", False):
        raise NotImplementedError(
            "func-mode rayleigh forcing is host-evaluated per window and "
            "not wired into the batch driver yet"
        )
    _assert_lockstep(mems)

    # regime fields must exist before the config build; batch-uniform by
    # construction (shared ud, lockstep window_step) — see run_window
    ud.is_compressible = eos.is_compressible(ud, mem0.time.window_step)
    ud.compressibility = eos.compressibility(ud, mem0.time.t, mem0.time.window_step)
    ud.is_nonhydrostatic = eos.is_nonhydrostatic(ud, mem0.time.window_step)
    ud.nonhydrostasy = eos.nonhydrostasy(ud, mem0.time.t, mem0.time.window_step)

    cfg, step_fns = _get_window_cache(mem0, ud)
    if cfg.terrain:
        raise NotImplementedError(
            "terrain CFL maxima are not batched yet (2D x-y DA cases only)"
        )

    K = len(mems)
    states = [to_device(mem) for mem in mems]
    sb = _stack(states) if mode == "vmap" else None
    if n_devices > 1:
        sb = _shard_members(sb, n_devices)

    while (mem0.time.t < tout) and (mem0.time.step < ud.stepmax):
        # ONE CFL kernel for both modes (loop stacks a transient view): the
        # batch-min dt is then bitwise-identical vmap<->loop by construction,
        # so the vmap gate isolates the vmapped STEP's reordering alone —
        # per-member `_cfl_maxima_plain` reduces in a different order than the
        # vmapped kernel on GPU and would drift dt at the last ULP
        sb_cfl = sb if mode == "vmap" else _stack(states)
        maxima = np.asarray(
            _batch_cfl_maxima(
                sb_cfl["rho"], sb_cfl["rhou"], sb_cfl["rhov"], sb_cfl["rhow"],
                sb_cfl["rhoY"], cfg.gamm, cfg.Msq,
            )
        )
        # every member's dt through the unchanged host path, then lockstep min
        dts = [_host_dt(maxima[k], mem, ud, tout)[0] for k, mem in enumerate(mems)]
        dt = prestep.apply_modifcations(min(dts), ud, mem0.time.step)

        # blending/regime control, batch-uniform (lockstep step/window_step,
        # shared ud): on conversion steps ALL members round-trip through host
        # memory and convert with the unchanged numpy routines
        label = "%.3d" % mem0.time.step
        conversion = _blend_conversion_due(
            bld, ud, mem0.time.step, mem0.time.window_step
        )
        if conversion:
            if mode == "vmap":
                for k, mem in enumerate(mems):
                    write_back(_unstack(sb, k), mem)
            else:
                for s, mem in zip(states, mems):
                    write_back(s, mem)
            for mem in mems:
                schemes.prepare_blending(
                    mem, ud, bld, label, None, mem0.time.step,
                    mem0.time.window_step, mem0.time.t, dt, False, False,
                )
            states = [to_device(mem) for mem in mems]
            sb = _stack(states) if mode == "vmap" else None
            if n_devices > 1:
                sb = _shard_members(sb, n_devices)
        else:
            # regime bookkeeping only (no state access) — once, shared ud
            schemes.prepare_blending(
                mem0, ud, bld, label, None, mem0.time.step,
                mem0.time.window_step, mem0.time.t, dt, False, False,
            )
        ud.is_nonhydrostatic = eos.is_nonhydrostatic(ud, mem0.time.window_step)
        ud.nonhydrostasy = eos.nonhydrostasy(ud, mem0.time.t, mem0.time.window_step)
        if int(ud.is_nonhydrostatic) == 0:
            raise NotImplementedError(
                "hydrostatic regime (alpha_w = 0) step path is not "
                "implemented on jax-device — run with backend='jax' (hybrid)"
            )

        parity = mem0.time.step % 2
        args = (
            dt,
            float(ud.nonhydrostasy),
            float(ud.compressibility),
            _DUMMY_FORCING,
            _DUMMY_FORCING,
        )
        regime = (parity, int(ud.is_nonhydrostatic), int(ud.is_compressible))
        psinc = regime[2] == 0
        if mode == "vmap":
            key = ("batch", K) + regime
            if key not in step_fns:
                step_fns[key] = make_batch_step(cfg, *regime)
            out = step_fns[key](sb, *args)
            if psinc:
                # keep npf.p2_nodes_half current per member (see run_window)
                sb, p2_half = out
                for k, mem in enumerate(mems):
                    mem.npf.p2_nodes_half = np.asarray(p2_half[k])
            else:
                sb = out
        else:
            if regime not in step_fns:
                step_fns[regime] = make_step(cfg, *regime)
            outs = [step_fns[regime](s, *args) for s in states]
            if psinc:
                states = [o[0] for o in outs]
                for mem, o in zip(mems, outs):
                    mem.npf.p2_nodes_half = np.asarray(o[1])
            else:
                states = outs

        logging.info(
            "device batch step %i done (K=%i, mode=%s), t = %.12f, dt = %.12f "
            "(member dt spread %.3e)",
            mem0.time.step, K, mode, mem0.time.t, dt, max(dts) - min(dts),
        )
        t_new = mem0.time.t + dt
        for mem in mems:
            mem.time.t = t_new
            mem.time.step += 1
            mem.time.window_step += 1

    for k, mem in enumerate(mems):
        write_back(_unstack(sb, k) if mode == "vmap" else states[k], mem)
    return mems
