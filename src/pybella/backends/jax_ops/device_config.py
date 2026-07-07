"""Device-resident step (ud.backend = "jax-device"): host-side configuration.

Builds the per-window DeviceConfig (a `_Namespace` of frozen geometry,
thermodynamic constants, hydrostatic profiles, per-sweep metric arrays,
boundary-fill plans and laplacian gather/mask plans) that the jitted step in
`device_kernels` closure-captures. Host-side; runs once per window."""

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
    # role-ordered (h1, v, h2) Coriolis components: scalars on the legacy
    # path, per-cell fields when ud.coriolis_field is set (f(phi) e_r on the
    # sphere). Reuse the numpy role_components so the field build + cache is
    # bit-identical to the hybrid path.
    from pybella.flow_solver.numerics import coriolis as coriolis_np

    w_role = coriolis_np.role_components(mem, ud)
    cfg.coriolis_role = tuple(
        jnp.asarray(c) if not np.isscalar(c) else float(c) for c in w_role
    )
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

    # general (non-vertical-line, spherical) metric: e_up buoyancy, general
    # H^-1, general free-slip walls. The tangent-plane surface constraint is
    # active only for thin-shell SWE-on-sphere cases.
    cfg.general = elem.metric is not None and not elem.metric.vertical_line
    cfg.constrain_to_surface = bool(getattr(ud, "constrain_to_surface", False))

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
            # Tier-3 (sphere) data; inert for vertical-line maps
            height=orient(m.height, ident),
            h_v=orient(m.h_v, ident),
            e_up=(None if m.e_up is None else [orient(c, ident) for c in m.e_up]),
            vertical_line=m.vertical_line,
        )
        cfg.node_metric_J_i1 = jnp.asarray(np.asarray(node.metric.J)[node.i1])
        cfg.node_metric_ooJ_i1 = jnp.asarray(np.asarray(node.metric.ooJ)[node.i1])
        cfg.ooJ_split = tuple(orient(m.ooJ, _orient_perm(ndim, s)) for s in range(ndim))
    else:
        cfg.metric = None

    # boundary fill plans (canonical + sweep orientations)
    cfg.boundary = jax_boundary.get_boundary_config(mem, ud)
    bdry_ints = tuple(jax_boundary._bdry_int(ud.bdry_type[d]) for d in range(ndim))
    degen = tuple((int(d), int(node.sc[d])) for d in axes.degenerate_axes(node))
    cfg.node_fill = jax_boundary._node_fill_fn(
        ndim, tuple(int(i) for i in node.igs), bdry_ints, degen
    )

    # pole-ring collapse index maps (Stage F): the Galerkin scatter/gather on
    # the flat 3D solve vector, mirroring numerics.pole_collapse.PoleCollapse.
    # ring_complement zeroes the ring by a mask multiply (bicgstab's
    # custom_linear_solve double-transposes the operator; an integer
    # scatter-SET is not double-transposable, a mask multiply + scatter-ADD is)
    from pybella.flow_solver.numerics import pole_collapse

    cfg.pole_collapse = None
    if pole_collapse.pole_axis_present(ud):
        pc = pole_collapse.get(node)
        ring_complement = np.ones(pc.n)
        ring_complement[pc.ring_all] = 0.0
        cfg.pole_collapse = _Namespace(
            scatter_src=jnp.asarray(pc.scatter_src),
            ring_complement=jnp.asarray(ring_complement),
            uniq_mem=jnp.asarray(pc.uniq_mem),
            uniq_master=jnp.asarray(pc.uniq_master),
        )

    # polar filter (Stage F): host-side static transfer factors + J weights +
    # the longitude-CFL cap, mirroring numerics.polar_filter.apply / cfl_cap
    from pybella.flow_solver.numerics import polar_filter as polar_filter_np

    cfg.polar_filter = getattr(ud, "polar_filter", None)
    cfg.filter_plan = None
    cfg.filter_cap = None
    if cfg.polar_filter is not None:
        m = elem.metric
        assert m is not None and not m.vertical_line, "polar filter needs a sphere"
        lam_axis = int(m.cart_haxes[0])
        phi_axis = int(m.cart_haxes[1])
        v_axis = int(m.cart_v)
        sc = [int(s) for s in elem.sc]
        igl, igp, igr = (int(elem.igs[a]) for a in (lam_axis, phi_axis, v_axis))
        Nlam = sc[lam_axis] - 2 * igl
        sl = [slice(None)] * ndim
        sl[lam_axis] = slice(igl, sc[lam_axis] - igl)
        sl[phi_axis] = slice(igp, sc[phi_axis] - igp)
        sl[v_axis] = slice(igr, sc[v_axis] - igr)
        sl = tuple(sl)
        phi_coords = axes.coords_along(elem, phi_axis)[igp : sc[phi_axis] - igp]
        r = polar_filter_np.transfer(
            Nlam, np.cos(phi_coords), cfg.polar_filter.phi_c, cfg.polar_filter.p
        )
        bidx = [None] * ndim
        bidx[lam_axis] = slice(None)
        bidx[phi_axis] = slice(None)
        cfg.filter_plan = _Namespace(
            sl=sl,
            lam_axis=lam_axis,
            N=Nlam,
            r=jnp.asarray(r[tuple(bidx)]),
            J=jnp.asarray(np.asarray(m.J)[sl]),
        )
        cap = polar_filter_np.cfl_cap(elem, ud)
        cfg.filter_cap = None if cap is None else jnp.asarray(cap)

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
