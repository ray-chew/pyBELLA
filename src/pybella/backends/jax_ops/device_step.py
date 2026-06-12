"""Device-resident jitted time step (`ud.backend = "jax-device"`).

One full semi-implicit step — RK half-advection, the two explicit/implicit
pairs, Strang advection, Rayleigh damping/forcing, optional diffusion, and
every ghost fill — traced into a single ``step(state, dt, ...)`` jit.
The host loop syncs one dt scalar per step and pulls the full state only at
output times (see :func:`run_window`).

State is a plain dict pytree with exactly 7 leaves
(rho, rhou, rhov, rhow, rhoY, rhoX, p2_nodes). Everything else —
geometry, thermodynamic constants, hydrostatic profiles, metric arrays in
every sweep orientation, boundary fill plans, laplacian gather plans,
sponge profiles — is frozen into a :class:`DeviceConfig` at window start
and closure-captured by ``make_step`` (XLA constants). Static jit
structure: Strang parity (2 compiled variants) and the regime ints, which
are constant per blending window (blending is guarded out).

Faithfulness notes:
- the numpy step's ``deepcopy`` saves (sol0, sol_half) become free
  reference holds on immutable arrays;
- ``npf`` scratch (p2_nodes0/rhs/wcenter/wplus/u,v,w) are trace locals;
- ``sol.pwchi`` (RK advection) has no live consumer and is dropped;
- the laplacian operators are rebuilt per solve from traced coefficients
  using static gather/mask plans — same arithmetic as the hybrid twins;
- bicgstab runs inside the jit (jax.scipy, scipy-matching semantics).
"""

import logging

import numpy as np
import jax
import jax.numpy as jnp

from pybella.utils import axes
from pybella.utils import options as opts
from pybella.utils.slices import (
    get_inner_slice,
    get_interface_indices,
    get_neighbor_indices,
)
from pybella.utils.operators.convolution import (
    get_averaging_kernel,
    get_flux_kernels,
)
from pybella.flow_solver.discretisation import terrain as terrain_mod
from pybella.flow_solver.physics import cfl as cfl_np

from . import advection as jax_advection
from . import boundary as jax_boundary
from . import convolution as jax_convolution
from . import coriolis as jax_coriolis
from . import diffusion as jax_diffusion
from . import divergence as jax_divergence
from . import gradient as jax_gradient
from . import rayleigh as jax_rayleigh
from .laplacian import lap2D as jax_lap2D
from .laplacian import lap3D as jax_lap3D

_SOL_FIELDS = ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX")
_MOMENTA = ("rhou", "rhov", "rhow")


class _Namespace:
    """Attribute bag (metric-like views over jnp arrays, plans, ...)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


# =========================================================================
# DeviceConfig
# =========================================================================


def _orient_perm(ndim, split):
    """Cyclic permutation placing axis `split` last (sweep orientation)."""
    return tuple((split + 1 + i) % ndim for i in range(ndim))


def build_device_config(mem, ud):
    """Freeze all static data for the jitted step (host-side, once)."""
    elem, node, th, npf = mem.elem, mem.node, mem.th, mem.npf
    ndim = elem.ndim
    v = axes.vertical_axis(ud)
    cfg = _Namespace()

    cfg.ndim = ndim
    cfg.v_phys = v
    cfg.role_perm = axes.role_perm(v)
    cfg.elem_dxyz = tuple(float(d) for d in elem.dxyz)
    cfg.node_dxyz = tuple(float(d) for d in node.dxyz)
    cfg.iisc = tuple(int(i) for i in elem.iisc[:ndim])
    cfg.node_i1 = node.i1
    cfg.node_i2 = node.i2
    cfg.node_sc = tuple(int(x) for x in npf.p2_nodes.shape)
    cfg.node_isc_shape = npf.rhs.shape

    # thermodynamics + scalars
    cfg.gamm, cfg.gm1, cfg.gm1inv = th.gamm, th.gm1, th.gm1inv
    cfg.Gamma, cfg.Gammainv = th.Gamma, th.Gammainv
    cfg.Msq = float(ud.Msq)
    cfg.g = float(ud.gravity_strength[v])
    cfg.coriolis = tuple(float(c) for c in ud.coriolis_strength)
    cfg.winds = (
        float(ud.u_wind_speed),
        float(ud.v_wind_speed),
        float(ud.w_wind_speed),
    )
    cfg.is_compressible = int(ud.is_compressible)
    cfg.atmosphere = hasattr(ud, "ATMOSPHERIC_EXTENSION")
    cfg.do_advection = bool(ud.do_advection)
    cfg.diffusion = bool(ud.diffusion)
    cfg.diffusion_coeff = float(ud.diffusion_coeff)
    cfg.rayleigh_bdry = ud.bdry_type[v] == opts.BdryType.RAYLEIGH
    cfg.has_forcing = bool(getattr(ud, "rayleigh_forcing", False))
    cfg.tol = float(ud.tol)
    cfg.max_iterations = int(ud.max_iterations)

    # hydrostatic profiles (force lazy inits; never traced)
    cfg.S0c = jnp.asarray(npf.HydroState.get_S0c(elem))
    cfg.dSdy = jnp.asarray(npf.HydroState_n.get_dSdy(elem, node))
    if npf.HydroState.field_mode:
        Ybar = npf.HydroState.Y0
    else:
        Ybar = jax_rayleigh._vertical_profile(npf.HydroState.Y0, ndim, v)
    cfg.Ybar = jnp.asarray(Ybar)

    # sponge profiles
    if cfg.rayleigh_bdry:
        cfg.tcy = jnp.asarray(jax_rayleigh._vertical_profile(ud.tcy, ndim, v))
    else:
        cfg.tcy = jnp.asarray(0.0)
    if cfg.has_forcing:
        cfg.tcy_f = jnp.asarray(jax_rayleigh._vertical_profile(ud.forcing_tcy, ndim, v))
        cfg.tny_f = jnp.asarray(ud.forcing_tny)

    # terrain metric: canonical + per-split sweep orientations
    cfg.terrain = elem.metric is not None
    if cfg.terrain:
        m = elem.metric
        assert m.vaxis == v, "config must be built from the canonical metric"

        def orient(arr, perm):
            return None if arr is None else jnp.asarray(np.transpose(arr, perm))

        ident = tuple(range(ndim))
        cfg.metric = _Namespace(
            J=orient(m.J, ident),
            ooJ=orient(m.ooJ, ident),
            G1=orient(m.G1, ident),
            G2=orient(m.G2, ident),
            # general (Klein) area normals: outer index = array axis,
            # inner = Cartesian component (see terrain.MetricFields)
            N=[[orient(c, ident) for c in Na] for Na in m.N],
            vaxis=m.vaxis,
            haxes=m.haxes,
            cart_v=m.cart_v,
            cart_haxes=m.cart_haxes,
        )
        cfg.node_metric_J_i1 = jnp.asarray(np.asarray(node.metric.J)[node.i1])
        cfg.node_metric_ooJ_i1 = jnp.asarray(np.asarray(node.metric.ooJ)[node.i1])
        cfg.ooJ_split = tuple(orient(m.ooJ, _orient_perm(ndim, s)) for s in range(ndim))
    else:
        cfg.metric = None

    # boundary fill plans (canonical + sweep orientations)
    cfg.boundary = jax_boundary.get_boundary_config(mem, ud)
    bdry_ints = tuple(
        (
            jax_boundary._PERIODIC
            if ud.bdry_type[d] == opts.BdryType.PERIODIC
            else jax_boundary._WALL
        )
        for d in range(ndim)
    )
    degen = tuple((int(d), int(node.sc[d])) for d in axes.degenerate_axes(node))
    cfg.node_fill = jax_boundary._node_fill_fn(
        ndim, tuple(int(i) for i in node.igs), bdry_ints, degen
    )

    # wall-node scaling masks (scale_wall_node_values as multiplicative plan)
    cfg.wall_scale_idx = []
    igs = node.igs
    for dim in range(ndim):
        if ud.bdry_type[dim] in (opts.BdryType.WALL, opts.BdryType.RAYLEIGH):
            for boundary_idx in [igs[dim], -igs[dim] - 1]:
                idx = [slice(igs[d], -igs[d]) for d in range(ndim)]
                idx[dim] = boundary_idx
                cfg.wall_scale_idx.append(tuple(idx))

    # flux container shapes (from the canonical cache containers)
    flux_np = mem.cache.get_flux_containers(elem)
    cfg.flux_shapes = tuple(flux_np[d].rhoY.shape for d in range(ndim))
    cfg.flux_kernels = get_flux_kernels(ndim)
    cfg.avg_kernel2 = get_averaging_kernel(ndim, width=2)

    # laplacian plans
    if ndim == 2:
        cfg.lap2d = _build_lap2d_plan(node, ud)
    else:
        cfg.lap3d = _build_lap3d_plan(mem, ud)

    return cfg


# =========================================================================
# laplacian plans (static structure; coefficients gathered traced per solve)
# =========================================================================


def _build_lap2d_plan(node, ud):
    """Static gather/mask structure of jax_ops.laplacian.lap2D.get_linop."""
    p = _Namespace()
    p.y_atmosphere = bool(
        hasattr(ud, "ATMOSPHERIC_EXTENSION") and ud.ATMOSPHERIC_EXTENSION
    )
    _wall = (opts.BdryType.WALL, opts.BdryType.RAYLEIGH)
    p.x_wall = ud.bdry_type[0] in _wall
    p.y_wall = ud.bdry_type[1] in _wall

    iicxn, iicyn = node.iicx, node.iicy
    N = iicxn * iicyn
    idx = np.arange(N)
    cnt_y, cnt_x = np.divmod(idx, iicxn)

    stencil = {
        "topleft": idx - iicxn - 1,
        "midleft": idx - 1,
        "botleft": idx + iicxn - 1,
        "topmid": idx - iicxn,
        "midmid": idx.copy(),
        "botmid": idx + iicxn,
        "topright": idx - iicxn + 1,
        "midright": idx + 1,
        "botright": idx + iicxn + 1,
    }
    left = cnt_x == 0
    right = cnt_x == iicxn - 1
    top = cnt_y == 0
    bot = cnt_y == iicyn - 1
    for name in ("topleft", "midleft", "botleft"):
        stencil[name][left] += iicxn - 1
    for name in ("topright", "midright", "botright"):
        stencil[name][right] -= iicxn - 1
    y_wrap = 2 * iicxn if p.y_atmosphere else iicxn * (iicyn - 1)
    for name in ("topleft", "topmid", "topright"):
        stencil[name][top] += y_wrap
    for name in ("botleft", "botmid", "botright"):
        stencil[name][bot] -= y_wrap

    ne_idx = cnt_y * (iicxn + 1) + cnt_x
    p.ne = {
        "tl": ne_idx,
        "tr": ne_idx + 1,
        "bl": ne_idx + (iicxn + 1),
        "br": ne_idx + (iicxn + 1) + 1,
    }
    # wall zero masks per corner (1.0 keep / 0.0 zero), matching the numpy
    # in-place zeroing
    masks = {c: np.ones(N) for c in ("tl", "tr", "bl", "br")}
    if p.x_wall:
        masks["tl"][left] = 0.0
        masks["bl"][left] = 0.0
        masks["tr"][right] = 0.0
        masks["br"][right] = 0.0
    if p.y_wall and not p.y_atmosphere:
        masks["tl"][top] = 0.0
        masks["tr"][top] = 0.0
        masks["bl"][bot] = 0.0
        masks["br"][bot] = 0.0
    p.stencil = {k: jnp.asarray(v) for k, v in stencil.items()}
    p.hp_masks = {c: jnp.asarray(m) for c, m in masks.items()}
    p.iicx, p.iicy = iicxn, iicyn
    p.oodx, p.oody = 1.0 / node.dx, 1.0 / node.dy
    return p


def _lap2d_matvec_factory(plan, hplusx, hplusy, cor4, hcenter, dinv):
    """Traced twin of jax_ops.laplacian.lap2D.get_linop's coefficient
    gathering; returns the matvec on the flat interior-node vector."""
    hpx = {c: hplusx[plan.ne[c]] * plan.hp_masks[c] for c in plan.ne}
    hpy = {c: hplusy[plan.ne[c]] * plan.hp_masks[c] for c in plan.ne}
    cxx, cyy, cxy, cyx = cor4
    cor = {
        "cxx": {c: cxx[plan.ne[c]] for c in plan.ne},
        "cyy": {c: cyy[plan.ne[c]] for c in plan.ne},
        "cxy": {c: cxy[plan.ne[c]] for c in plan.ne},
        "cyx": {c: cyx[plan.ne[c]] for c in plan.ne},
    }

    def matvec(pvec):
        return jax_lap2D._lap2D_gather(
            pvec,
            plan.stencil,
            hpx,
            hpy,
            cor,
            hcenter,
            dinv,
            plan.oodx,
            plan.oody,
        )

    return matvec


def _build_lap3d_plan(mem, ud):
    """Static structure of jax_ops.laplacian.lap3D.get_linop: periodicity
    permutations, wall slab masks, scales, and the (static) use_cross."""
    elem, node = mem.elem, mem.node
    ndim = elem.ndim
    p = _Namespace()
    oodxyz = 1.0 / (node.dxyz**2)
    p.scales = (
        float(oodxyz[0]),
        float(oodxyz[1]),
        float(oodxyz[2]),
        1.0 / node.dx,
        1.0 / node.dy,
        1.0 / node.dz,
    )
    p.i1 = (slice(1, -1), slice(1, -1), slice(1, -1))
    periodic = [ud.bdry_type[d] == opts.BdryType.PERIODIC for d in range(ndim)]

    hshape = np.asarray(mem.npf.wcenter)[p.i1].shape
    p.padded = tuple(sz + 2 for sz in hshape)

    # wall slab multiplicative mask on the CELL coefficient box (the C_ij
    # boxes are cell-shaped sc sliced by i1, one larger than the node box)
    cell_box = tuple(int(sz) - 2 for sz in np.asarray(mem.sol.rho).shape)
    mask = np.ones(cell_box)
    for dim in range(ndim):
        if not periodic[dim]:
            lo = [slice(None)] * 3
            hi = [slice(None)] * 3
            lo[dim] = 0
            hi[dim] = -1
            mask[tuple(lo)] = 0.0
            mask[tuple(hi)] = 0.0
    p.wall_mask = jnp.asarray(mask)

    perms = []
    for dim in range(3):
        n = p.padded[dim]
        perm = np.arange(n)
        if dim < ndim and periodic[dim]:
            perm[[0, 1, -2, -1]] = [n - 3, n - 2, 1, 2]
        perms.append(jnp.asarray(perm))
    p.perms = tuple(perms)

    # use_cross is structural (terrain / off-diagonal H^-1): decide once on
    # the initial state with the numpy criterion
    from pybella.flow_solver.numerics import coriolis as coriolis_np

    dt0 = float(getattr(ud, "dtfixed", 0.0)) or 1.0
    hv = coriolis_np.compute_inverse_coefficients(mem, ud, dt0)
    h_role = ((hv[0], hv[1], hv[2]), (hv[3], hv[4], hv[5]), (hv[6], hv[7], hv[8]))
    if elem.metric is not None:
        h_role = terrain_mod.elliptic_tensor(elem.metric, h_role)
    rho_of = axes.role_of_axis(axes.vertical_axis(ud))
    # npf.wplus is only populated once operator_coefficients_nodes has run;
    # on a fresh state it is zeros, so compute the coefficient directly from
    # sol (wplus[i] == Gammainv * rhoY * Y for every i)
    Y0 = np.asarray(mem.sol.rhoY) / np.asarray(mem.sol.rho)
    coeff0 = mem.th.Gammainv * np.asarray(mem.sol.rhoY) * Y0
    offmax = 0.0
    for i in range(3):
        for j in range(3):
            if i != j:
                cij = coeff0 * np.asarray(h_role[rho_of[i]][rho_of[j]])
                offmax = max(offmax, float(np.max(np.abs(cij[p.i1]))))
    p.use_cross = bool(offmax > 0.0)
    p.rho_of = rho_of
    return p


# =========================================================================
# state transfer
# =========================================================================


def to_device(mem):
    s = {name: jnp.asarray(getattr(mem.sol, name)) for name in _SOL_FIELDS}
    s["p2_nodes"] = jnp.asarray(mem.npf.p2_nodes)
    return s


def write_back(s, mem):
    for name in _SOL_FIELDS:
        getattr(mem.sol, name)[...] = np.asarray(s[name])
    mem.npf.p2_nodes[...] = np.asarray(s["p2_nodes"])


# =========================================================================
# functional substeps (all traced; s = state dict)
# =========================================================================


def _sol_tuple(s):
    return tuple(s[n] for n in _SOL_FIELDS)


def _ghost_fill(s, cfg, split=None):
    """set_ghost_cells twin on the state dict (canonical or sweep-oriented).

    With split=None the fields are canonical and all dims are processed;
    with split=k the fields are sweep-oriented (axis k last) and only the
    last array axis is filled (current_step semantics = k).
    """
    bcfg = cfg.boundary
    ndim = cfg.ndim
    dims = range(ndim) if split is None else [ndim - 1]
    for dim in dims:
        current_step = split if split is not None else dim
        arrays = _sol_tuple(s)
        if bcfg.gravity_on[current_step]:
            orientation = "sweep" if split is not None else "phys"
            out = bcfg.gravity_fill[orientation](*arrays)
        else:
            out = bcfg.no_gravity_fill[(dim, current_step)](*arrays)
        for name, val in zip(_SOL_FIELDS, out):
            s[name] = val
    return s


def _advective_flux(s, cfg):
    """advective_flux.recompute twin: per-dim rhoY interface fluxes
    (sweep-oriented containers, ghost faces zero as in the fresh cache)."""
    ndim = cfg.ndim
    inner_idx = get_inner_slice(ndim)
    components = ["u", "v"] if ndim == 2 else ["u", "v", "w"]
    rho_components = ["rhou", "rhov"] if ndim == 2 else ["rhou", "rhov", "rhow"]
    m = cfg.metric

    flux_rhoY = []
    for i, (comp, rho_comp) in enumerate(zip(components, rho_components)):
        momentum = s[rho_comp]
        if m is not None and i == m.vaxis:
            a_h1, a_h2 = m.haxes
            momentum = momentum - m.G1 * s[rho_components[a_h1]]
            if m.G2 is not None:
                momentum = momentum - m.G2 * s[rho_components[a_h2]]
            rhoY_vel = s["rhoY"] * momentum / s["rho"]
        elif m is not None:
            rhoY_vel = m.J * s["rhoY"] * momentum / s["rho"]
        else:
            rhoY_vel = s["rhoY"] * momentum / s["rho"]
        conv = jax_convolution.apply_directional_convolution(
            rhoY_vel, cfg.flux_kernels[comp], comp, ndim
        )
        flux_rhoY.append(
            jnp.zeros(cfg.flux_shapes[i], dtype=jnp.float64).at[inner_idx].set(conv)
        )
    return flux_rhoY


def _flip(s, k=1):
    """k forward flips (moveaxis 0 -> last) of all sol fields."""
    for _ in range(k):
        for name in _SOL_FIELDS:
            s[name] = jnp.moveaxis(s[name], 0, -1)
    return s


def _sweep_flux(s, cfg, flux_rhoY_d, lmbda_rec, split):
    """recovery + HLL for one sweep (fields sweep-oriented)."""
    ndim = cfg.ndim
    rho = s["rho"]
    u = s["rhou"] / rho
    v = s["rhov"] / rho
    w = s["rhow"] / rho
    Y = s["rhoY"] / rho
    X = s["rhoX"] / rho
    vel = (u, v, w)[split]
    ooJ = cfg.ooJ_split[split] if cfg.terrain else jnp.asarray(1.0)

    zeros = jnp.zeros_like(flux_rhoY_d)
    return jax_advection._recovery_hll(
        u,
        v,
        w,
        X,
        Y,
        s["rhoY"],
        vel,
        flux_rhoY_d,
        zeros,
        zeros,
        zeros,
        zeros,
        zeros,
        ooJ,
        lmbda_rec,
        cfg.gamm,
        ndim,
        cfg.terrain,
    )


def _flux_update(s, cfg, flux, flux_rhoY_d, lmbda, split):
    """_update_solution_variables twin (fields sweep-oriented)."""
    lefts_idx, rights_idx = get_neighbor_indices(cfg.ndim)
    flux_rho, flux_rhou, flux_rhov, flux_rhow, flux_rhoX = flux
    by_name = {
        "rho": flux_rho,
        "rhou": flux_rhou,
        "rhov": flux_rhov,
        "rhow": flux_rhow,
        "rhoX": flux_rhoX,
        "rhoY": flux_rhoY_d,
    }
    ooJ = cfg.ooJ_split[split] if cfg.terrain else None
    for name in _SOL_FIELDS:
        diff = by_name[name][lefts_idx] - by_name[name][rights_idx]
        if ooJ is not None:
            diff = ooJ * diff
        s[name] = s[name] + lmbda * diff
    return s


def _advect_rk(s, cfg, dt, flux_rhoY):
    """first_order_runge_kutta twin (dt here is the half step, as called)."""
    ndim = cfg.ndim
    fluxes = [None] * ndim
    for split in range(ndim):
        s = _flip(s)
        if cfg.iisc[split] > 1:
            s = _ghost_fill(s, cfg, split=split)
            fluxes[split] = _sweep_flux(s, cfg, flux_rhoY[split], 0.0, split)
    for split in range(ndim):
        lmbda = dt / cfg.elem_dxyz[split]
        s = _flip(s)
        flux = fluxes[split]
        if flux is None:
            # numpy still applies the update with the zero-initialised flux
            # container components and the recomputed rhoY flux
            z = jnp.zeros_like(flux_rhoY[split])
            flux = (z, z, z, z, z)
        s = _flux_update(s, cfg, flux, flux_rhoY[split], lmbda, split)
    s = _ghost_fill(s, cfg)
    return s


def _advect_strang(s, cfg, dt, parity, flux_rhoY):
    """strange_splitting twin; parity is static (two jit variants)."""
    ndim = cfg.ndim
    time_step = 0.5 * dt
    sweeps = [not bool(parity), bool(parity)]
    for reverse in sweeps:
        dim_range = range(ndim - 1, -1, -1) if reverse else range(ndim)
        for split in dim_range:
            lmbda = time_step / cfg.elem_dxyz[split]
            if reverse:
                if cfg.iisc[split] > 1:
                    s = _ghost_fill(s, cfg, split=split)
                    flux = _sweep_flux(s, cfg, flux_rhoY[split], lmbda, split)
                    s = _flux_update(s, cfg, flux, flux_rhoY[split], lmbda, split)
                s = _flip(s, ndim - 1)  # flip_backward == ndim-1 forward flips
            else:
                s = _flip(s)
                if cfg.iisc[split] > 1:
                    s = _ghost_fill(s, cfg, split=split)
                    flux = _sweep_flux(s, cfg, flux_rhoY[split], lmbda, split)
                    s = _flux_update(s, cfg, flux, flux_rhoY[split], lmbda, split)
    s = _ghost_fill(s, cfg)
    return s


def _coriolis_inputs(s, cfg, dt, nonhydro):
    wdt_h1 = dt * cfg.coriolis[cfg.role_perm[0]]
    wdt_v = dt * cfg.coriolis[cfg.role_perm[1]]
    wdt_h2 = dt * cfg.coriolis[cfg.role_perm[2]]
    Y = s["rhoY"] / s["rho"]
    nu = -(dt**2) * (cfg.g / cfg.Msq) * cfg.dSdy * Y
    return wdt_h1, wdt_v, wdt_h2, nu, nonhydro


def _apply_hinv(fields3, s, cfg, dt, nonhydro):
    """H^-1 @ (axis-indexed component triple), role plumbing as numpy."""
    ax_h1, ax_v, ax_h2 = cfg.role_perm
    wh1, wv, wh2, nu, _ = _coriolis_inputs(s, cfg, dt, nonhydro)
    U, V, W = fields3[ax_h1], fields3[ax_v], fields3[ax_h2]
    u, v, w = jax_coriolis.apply_inverse(U, V, W, wh1, wh2, wv, nu, nonhydro)
    out = [None, None, None]
    out[ax_h1], out[ax_v], out[ax_h2] = u, v, w
    return tuple(out)


def _explicit_part(s, cfg, dt, nonhydro):
    """implicit_euler.do_explicit_part twin."""
    vmom = _MOMENTA[cfg.v_phys]
    dbuoy = s["rhoY"] * (s["rhoX"] / s["rho"])
    s[vmom] = (nonhydro * s[vmom]) - dt * (cfg.g / cfg.Msq) * dbuoy

    u0, v0, w0 = cfg.winds
    for name, w_ in zip(_MOMENTA, (u0, v0, w0)):
        s[name] = s[name] + (-1.0) * w_ * s["rho"]

    moms = _apply_hinv((s["rhou"], s["rhov"], s["rhow"]), s, cfg, dt, nonhydro)
    s["rhou"], s["rhov"], s["rhow"] = moms

    for name, w_ in zip(_MOMENTA, (u0, v0, w0)):
        s[name] = s[name] + (+1.0) * w_ * s["rho"]
    return s


def _divergence_rhs(s, cfg):
    """divergence.compute_at_nodes twin; adopts the wall-zeroed momenta."""
    momenta = (
        s["rhou"],
        s["rhov"],
        s["rhow"] if cfg.ndim == 3 else None,
    )
    if cfg.terrain:
        m = cfg.metric
        metric = (
            (m.J, m.G1, None, None, None)
            if cfg.ndim == 2
            else (m.J, m.G1, m.G2, m.vaxis, m.haxes)
        )
    else:
        metric = None
    wall_dims = (
        ()
        if cfg.atmosphere
        else tuple(
            d for d in range(cfg.ndim) if cfg.boundary.bdry_int[d] == jax_boundary._WALL
        )
    )
    rhs, momenta_out = jax_divergence.compute_at_nodes(
        momenta,
        s["rho"],
        s["rhoY"],
        cfg.ndim,
        cfg.elem_dxyz,
        wall_dims=wall_dims,
        metric=metric,
    )
    s["rhou"] = momenta_out[0]
    s["rhov"] = momenta_out[1]
    if cfg.ndim == 3:
        s["rhow"] = momenta_out[2]
    return s, rhs


def _scale_wall_nodes(rhs, cfg, factor):
    for idx in cfg.wall_scale_idx:
        rhs = rhs.at[idx].multiply(factor)
    return rhs


def _forward_step(s, cfg, dt, nonhydro, compressibility):
    """explicit_euler.do_forward_step twin."""
    ndim = cfg.ndim
    ax_h1, ax_v, ax_h2 = cfg.role_perm
    corr = cfg.coriolis
    corr_h1, corr_v, corr_h2 = corr[ax_h1], corr[ax_v], corr[ax_h2]
    u0, v0, w0 = cfg.winds
    Ginv = cfg.Gammainv

    s, rhs = _divergence_rhs(s, cfg)
    if not cfg.atmosphere:
        rhs = _scale_wall_nodes(rhs, cfg, 2.0)

    dpidP = (cfg.gm1 / cfg.Msq) * jax_convolution.apply_convolution_kernel(
        s["rhoY"] ** (cfg.gamm - 2.0), cfg.avg_kernel2
    )

    rhoYovG = Ginv * s["rhoY"]
    dbuoy = s["rhoY"] * (s["rhoX"] / s["rho"])

    dpd = list(jax_gradient.compute_at_nodes(s["p2_nodes"], ndim, cfg.node_dxyz))
    if cfg.terrain:
        dpd = terrain_mod.apply_gradient_map(cfg.metric, dpd)

    rho = s["rho"]
    drho = (s["rhou"] - u0 * rho, s["rhov"] - v0 * rho, s["rhow"] - w0 * rho)
    dm_h1, dm_v, dm_h2 = drho[ax_h1], drho[ax_v], drho[ax_h2]
    dp_h1, dp_v, dp_h2 = dpd[ax_h1], dpd[ax_v], dpd[ax_h2]

    mom = {n: s[n] for n in _MOMENTA}
    names = _MOMENTA
    vel_v = mom[names[ax_v]] / rho

    mom[names[ax_h1]] = mom[names[ax_h1]] - dt * (
        rhoYovG * dp_h1 - corr_h2 * dm_v + corr_v * dm_h2
    )
    mom[names[ax_v]] = (
        mom[names[ax_v]]
        - dt
        * (
            rhoYovG * dp_v
            + (cfg.g / cfg.Msq) * dbuoy * nonhydro
            - corr_h1 * dm_h2
            + corr_h2 * dm_h1
        )
        * 1.0
    )  # (1 - is_ArakawaKonor), guarded to 0
    mom[names[ax_h2]] = mom[names[ax_h2]] - dt * (
        rhoYovG * dp_h2 - corr_v * dm_h1 + corr_h1 * dm_v
    )
    for n in _MOMENTA:
        s[n] = mom[n]

    s["rhoX"] = (rho * (rho / s["rhoY"] - cfg.S0c)) - dt * (vel_v * cfg.dSdy) * rho

    dp2n = jnp.zeros_like(s["p2_nodes"])
    if cfg.terrain:
        interior = -dt * dpidP * (rhs * cfg.node_metric_ooJ_i1)
    else:
        interior = -dt * dpidP * rhs
    dp2n = dp2n.at[cfg.node_i1].add(interior)
    s["p2_nodes"] = s["p2_nodes"] + compressibility * dp2n
    s["p2_nodes"] = cfg.node_fill(s["p2_nodes"])
    return s


def _operator_coefficients(s, cfg, dt):
    """operator_coefficients_nodes twin: returns (wplus_coeff, wcenter)."""
    ccenter = -cfg.Msq * cfg.gm1inv / (dt**2)
    cexp = 2.0 - cfg.gamm
    Y = s["rhoY"] / s["rho"]
    coeff = cfg.Gammainv * s["rhoY"] * Y

    wcenter = ccenter * jax_convolution.apply_convolution_kernel(
        s["rhoY"] ** cexp, cfg.avg_kernel2
    )
    if cfg.terrain:
        wcenter = wcenter * cfg.node_metric_J_i1
    if not cfg.atmosphere:
        wcenter = _scale_wall_nodes(wcenter, cfg, 0.5)
    return coeff, wcenter


def _correction_nodes(s, cfg, dt, p, updt_chi, nonhydro):
    """implicit_euler._correction_nodes twin."""
    Gammainv = cfg.Gammainv
    dpd = list(jax_gradient.compute_at_nodes(p, cfg.ndim, cfg.node_dxyz))
    if cfg.terrain:
        dpd = terrain_mod.apply_gradient_map(cfg.metric, dpd)
    Dpx, Dpy, Dpz = dpd

    thinv = s["rho"] / s["rhoY"]
    Y = s["rhoY"] / s["rho"]
    coeff = Gammainv * s["rhoY"] * Y

    pu = -dt * coeff * Dpx
    pv = -dt * coeff * Dpy
    pw = -dt * coeff * Dpz
    pu, pv, pw = _apply_hinv((pu, pv, pw), s, cfg, dt, nonhydro)

    s["rhou"] = s["rhou"] + thinv * pu
    s["rhov"] = s["rhov"] + thinv * pv
    s["rhow"] = s["rhow"] + thinv * pw
    vmom = _MOMENTA[cfg.v_phys]
    s["rhoX"] = s["rhoX"] + (-updt_chi) * dt * cfg.dSdy * s[vmom]
    return s


def _prepare_diag(cfg, wplus_list, wcenter, cii=None):
    """preconditioner.prepare_diag twin (traced)."""
    ndim = cfg.ndim
    w = cii if cii is not None else (wplus_list + [None])[:3]
    coeff = 0.75 if ndim == 2 else 0.0625
    dx, dy, dz = cfg.node_dxyz
    inv_dx2, inv_dy2 = 1.0 / (dx**2), 1.0 / (dy**2)
    diag = wcenter
    diag = diag - coeff * inv_dx2 * jax_convolution.apply_convolution_kernel(
        w[0], cfg.avg_kernel2
    )
    diag = diag - coeff * inv_dy2 * jax_convolution.apply_convolution_kernel(
        w[1], cfg.avg_kernel2
    )
    if ndim == 2:
        inv_dxdy = 1.0 / (dx * dy)
        diag = diag - coeff * inv_dxdy * jax_convolution.apply_convolution_kernel(
            w[0], cfg.avg_kernel2
        )
        diag = diag - coeff * inv_dxdy * jax_convolution.apply_convolution_kernel(
            w[1], cfg.avg_kernel2
        )
    else:
        inv_dz2 = 1.0 / (dz**2)
        diag = diag - coeff * inv_dz2 * jax_convolution.apply_convolution_kernel(
            w[2], cfg.avg_kernel2
        )
    return 1.0 / diag


def _fravel(a):
    """np.ravel(a, order='F') for 2D arrays."""
    return a.T.ravel()


def _coriolis_h_fields(s, cfg, dt, nonhydro):
    wh1, wv, wh2, nu, _ = _coriolis_inputs(s, cfg, dt, nonhydro)
    return jax_coriolis.compute_coefficients(wh1, wh2, wv, nu, nonhydro)


def _implicit_part(s, cfg, dt, nonhydro, compressibility, sol0=None):
    """implicit_euler.do_implicit_part twin."""
    ndim = cfg.ndim

    # numpy fills sol_for_boundary = sol0 (incompressible) or sol; the
    # caller passes an already-filled sol0 in the incompressible branch
    if sol0 is None:
        s = _ghost_fill(s, cfg)

    wplus_coeff, wcenter = _operator_coefficients(s, cfg, dt)

    s = _correction_nodes(s, cfg, dt, s["p2_nodes"], 0, nonhydro)
    s = _ghost_fill(s, cfg)

    s, rhs = _divergence_rhs(s, cfg)
    rhs = rhs / dt

    # compressibility correction (ArakawaKonor guarded out)
    if cfg.is_compressible == 0:
        wcenter = wcenter * compressibility
    else:
        wcenter = wcenter * compressibility

    h = _coriolis_h_fields(s, cfg, dt, nonhydro)
    h11, h12, h13, h21, h22, h23, h31, h32, h33, _ = h

    if ndim == 2:
        coriolis_params = (h11.T, h22.T, h12.T, h21.T)
        if cfg.terrain:
            h11_t, h22_t, h12_t, h21_t = coriolis_params
            h2x2 = ((h11_t.T, h12_t.T), (h21_t.T, h22_t.T))
            M = terrain_mod.elliptic_tensor_2d(cfg.metric, h2x2)
            coriolis_params = (M[0][0].T, M[1][1].T, M[0][1].T, M[1][0].T)
            met = cfg.metric
            diag_inv = _prepare_diag(
                cfg,
                None,
                wcenter,
                cii=[
                    wplus_coeff * met.J,
                    wplus_coeff * (1.0 + met.G1 * met.G1) * met.ooJ,
                    None,
                ],
            )
        else:
            diag_inv = _prepare_diag(cfg, [wplus_coeff, wplus_coeff], wcenter)
        rhs = rhs * diag_inv

        plan = cfg.lap2d
        coeff_slc = (slice(1, -1), slice(1, -1))
        hplusx = _fravel(wplus_coeff[coeff_slc])
        hplusy = _fravel(wplus_coeff[coeff_slc])
        hcenter = _fravel(wcenter[cfg.node_i1])
        cor4 = tuple(c[coeff_slc].ravel() for c in coriolis_params)
        dinv = _fravel(diag_inv[cfg.node_i1])
        matvec = _lap2d_matvec_factory(plan, hplusx, hplusy, cor4, hcenter, dinv)
        rhs_inner = rhs[cfg.node_i1].T.ravel()
    else:
        h_role = ((h11, h12, h13), (h21, h22, h23), (h31, h32, h33))
        if cfg.terrain:
            h_role = terrain_mod.elliptic_tensor(cfg.metric, h_role)
        rho_of = cfg.lap3d.rho_of
        cij = [
            [wplus_coeff * h_role[rho_of[i]][rho_of[j]] for j in range(3)]
            for i in range(3)
        ]
        diag_inv = _prepare_diag(
            cfg, None, wcenter, cii=[cij[0][0], cij[1][1], cij[2][2]]
        )
        rhs = rhs * diag_inv

        plan = cfg.lap3d
        i1 = plan.i1
        C = tuple(
            tuple(cij[i][j][i1] * plan.wall_mask for j in range(3)) for i in range(3)
        )
        hcenter = wcenter[i1]

        def matvec(pvec):
            return jax_lap3D._lap3D(
                pvec,
                C,
                hcenter,
                plan.scales,
                plan.padded,
                plan.perms,
                diag_inv,
                plan.use_cross,
            )

        rhs_inner = jnp.zeros_like(rhs).at[cfg.node_i1].set(rhs[cfg.node_i1]).ravel()

    A = lambda x: jnp.reshape(matvec(x), x.shape)
    p2, _ = jax.scipy.sparse.linalg.bicgstab(
        A, rhs_inner, tol=1e-5, atol=cfg.tol, maxiter=cfg.max_iterations
    )

    # reshape into the full node box
    p2_full = jnp.zeros(cfg.node_sc, dtype=jnp.float64)
    if ndim == 2:
        p2_full = p2_full.at[cfg.node_i2].set(p2.reshape(rhs[cfg.node_i1].T.shape).T)
    else:
        p2_full = p2_full.at[cfg.node_i1].set(p2.reshape(rhs.shape))

    s = _correction_nodes(s, cfg, dt, p2_full, 1, nonhydro)
    s["p2_nodes"] = s["p2_nodes"] + p2_full
    s["p2_nodes"] = cfg.node_fill(s["p2_nodes"])
    s = _ghost_fill(s, cfg)
    return s


def _rayleigh_damp(s, cfg, forcing=None):
    """rayleigh_damping twin on the state dict."""
    if forcing is not None:
        tcy = jnp.asarray(0.0)
        tcy_f, tny_f = cfg.tcy_f, cfg.tny_f
        u_f, v_f, Y_f, pi_f = forcing
        mfac, c_f = 1.0, 1.0
        has_forcing = True
    else:
        tcy = cfg.tcy
        tcy_f = tny_f = jnp.asarray(0.0)
        u_f = v_f = Y_f = pi_f = jnp.asarray(0.0)
        mfac, c_f = 0.0, 0.0
        has_forcing = False

    out = jax_rayleigh._damp(
        s["rho"],
        s["rhou"],
        s["rhov"],
        s["rhow"],
        s["rhoY"],
        s["p2_nodes"],
        tcy,
        tcy_f,
        tny_f,
        cfg.Ybar,
        cfg.winds[0],
        cfg.winds[1],
        cfg.winds[2],
        u_f,
        v_f,
        Y_f,
        pi_f,
        mfac,
        c_f,
        cfg.ndim,
        has_forcing,
    )
    s["rhou"], s["rhov"], rhow, s["rhoY"], p2 = out
    if cfg.ndim == 3:
        s["rhow"] = rhow
    if forcing is not None:
        s["p2_nodes"] = p2
        s = _ghost_fill(s, cfg)  # apply_rayleigh_forcing refills ghosts
    return s


def _diffuse(s, cfg, dt):
    rho, rhou, rhov, rhow = jax_diffusion._diffuse(
        s["rho"],
        s["rhou"],
        s["rhov"],
        s["rhow"],
        s["rhoY"],
        1.0 / cfg.S0c,
        dt * cfg.diffusion_coeff,
        cfg.elem_dxyz,
        cfg.ndim,
    )
    s["rho"], s["rhou"], s["rhov"] = rho, rhou, rhov
    if cfg.ndim == 3:
        s["rhow"] = rhow
    s = _ghost_fill(s, cfg)
    return s


# =========================================================================
# the step + CFL + window driver
# =========================================================================


def make_step(cfg, parity, is_nonhydrostatic):
    """Build the jitted full step for one Strang parity + regime structure."""

    def step(s, dt, nonhydro, compressibility, forcing_half, forcing_full):
        sol0 = dict(s)  # free reference hold (incl. p2 for the fill helper)

        flux_rhoY = _advective_flux(s, cfg)
        if cfg.do_advection:
            s = _advect_rk(s, cfg, 0.5 * dt, flux_rhoY)

        p2_nodes0 = s["p2_nodes"]

        s = _explicit_part(s, cfg, 0.5 * dt, nonhydro)
        if cfg.is_compressible == 0:
            sol0 = _ghost_fill(sol0, cfg)  # numpy fills the held copy
        s = _implicit_part(
            s,
            cfg,
            0.5 * dt,
            nonhydro,
            compressibility,
            sol0=sol0 if cfg.is_compressible == 0 else None,
        )

        if cfg.rayleigh_bdry:
            s = _rayleigh_damp(s, cfg)
        if cfg.has_forcing:
            s = _rayleigh_damp(s, cfg, forcing=forcing_half)

        # the Strang advective flux is computed from the POST-implicit
        # half-time state (time_update line 123), before the sol restore
        flux_rhoY_half = _advective_flux(s, cfg)

        # p2 reset branch (static regime structure)
        if is_nonhydrostatic == 0 or (
            cfg.is_compressible == 1 and is_nonhydrostatic == 1
        ):
            s["p2_nodes"] = p2_nodes0

        # restore sol to t^n for the full forward pass
        for n in _SOL_FIELDS:
            s[n] = sol0[n]

        s = _forward_step(s, cfg, 0.5 * dt, nonhydro, compressibility)

        if cfg.do_advection:
            s = _advect_strang(s, cfg, dt, parity, flux_rhoY_half)

        s = _explicit_part_post(s, cfg, dt, nonhydro, compressibility)

        if cfg.rayleigh_bdry:
            s = _rayleigh_damp(s, cfg)
        if cfg.has_forcing:
            s = _rayleigh_damp(s, cfg, forcing=forcing_full)

        if cfg.diffusion:
            s = _diffuse(s, cfg, dt)
        return s

    return jax.jit(step, donate_argnums=(0,))


def _explicit_part_post(s, cfg, dt, nonhydro, compressibility):
    s = _explicit_part(s, cfg, 0.5 * dt, nonhydro)
    s = _implicit_part(s, cfg, 0.5 * dt, nonhydro, compressibility, sol0=None)
    return s


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
    contra = moms[m.vaxis] - m.G1 * moms[m.haxes[0]]
    slope_sq = m.G1**2
    if m.G2 is not None:
        contra = contra - m.G2 * moms[m.haxes[1]]
        slope_sq = slope_sq + m.G2**2
    vels = [u, v, w]
    vels[m.vaxis] = jnp.abs(contra / rho) * m.ooJ
    u, v, w = vels
    c_vert = c * jnp.sqrt(1.0 + slope_sq) * m.ooJ
    cs = [c, c, c]
    cs[m.vaxis] = c_vert
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
    if getattr(ud, "continuous_blending", False) or getattr(
        ud, "initial_blending", False
    ):
        problems.append("dynamics blending")
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
    if "CFLfixed" in getattr(ud, "aux", ""):
        problems.append("CFLfixed prestep override")
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


def run_window(mem, ud, tout, writer=None):
    """Device-resident replacement for time_update.do's step loop."""
    from pybella.flow_solver.physics import eos

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

        # the non-blending prepare_blending path sets all four regime
        # fields per step via eos; replicate (blending itself is guarded)
        ud.is_compressible = eos.is_compressible(ud, mem.time.window_step)
        ud.compressibility = eos.compressibility(ud, mem.time.t, mem.time.window_step)
        ud.is_nonhydrostatic = eos.is_nonhydrostatic(ud, mem.time.window_step)
        ud.nonhydrostasy = eos.nonhydrostasy(ud, mem.time.t, mem.time.window_step)
        assert (
            int(ud.is_compressible) == cfg.is_compressible
        ), "is_compressible changed mid-window — unsupported on jax-device"

        parity = mem.time.step % 2
        key = (parity, int(ud.is_nonhydrostatic))
        if key not in step_fns:
            step_fns[key] = make_step(cfg, parity, int(ud.is_nonhydrostatic))
            compile_count += 1

        if cfg.has_forcing:
            forcing_half = _eval_forcing(mem, ud, mem.time.t + 0.5 * dt)
            forcing_full = _eval_forcing(mem, ud, mem.time.t + dt)
        else:
            forcing_half = forcing_full = _DUMMY_FORCING

        s = step_fns[key](
            s,
            dt,
            float(ud.nonhydrostasy),
            float(ud.compressibility),
            forcing_half,
            forcing_full,
        )

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
