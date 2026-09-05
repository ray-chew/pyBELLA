"""Device-resident step: the traced substeps and the compiled step.

Split out of `device_step` (host loop in `device_loop`); see that module
for the overall design notes."""

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

from .device_config import (
    _SOL_FIELDS,
    _MOMENTA,
    _lap2d_matvec_factory,
)


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
        elif bcfg.bdry_int[current_step] == jax_boundary._POLE:
            # lat-lon pole fold: pure index remap of the phi ghost slabs,
            # reached canonically or during the phi sweep
            orientation = "sweep" if split is not None else "phys"
            out = bcfg.pole_cell_fill[orientation](*arrays)
        elif dim == current_step and current_step in bcfg.general_wall_fill:
            # general (spherical) free-slip wall — only reached at canonical
            # orientation (see jax_boundary.set_ghost_cells for the argument)
            out = bcfg.general_wall_fill[current_step](*arrays)
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
        if m is not None:
            # general curvilinear mass flux rhoY * (N_i . m) / rho
            # (vertical-first contraction, mirroring the numpy
            # advective_flux._normal_momentum assembly)
            Ni = m.N[i]
            cv, (ch1, ch2) = m.cart_v, m.cart_haxes
            momentum = Ni[cv] * s[rho_components[cv]] + Ni[ch1] * s[rho_components[ch1]]
            if ch2 is not None:
                momentum = momentum + Ni[ch2] * s[rho_components[ch2]]
            rhoY_vel = s["rhoY"] * momentum / s["rho"]
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
    # role-ordered (h1, v, h2) rotation components: scalars, or coriolis_field
    # per-cell arrays (cfg.coriolis_role)
    wdt_h1 = dt * cfg.coriolis_role[0]
    wdt_v = dt * cfg.coriolis_role[1]
    wdt_h2 = dt * cfg.coriolis_role[2]
    Y = s["rhoY"] / s["rho"]
    nu = -(dt**2) * (cfg.g / cfg.Msq) * cfg.dSdy * Y
    return wdt_h1, wdt_v, wdt_h2, nu, nonhydro


def _e_role(cfg):
    """Role-ordered (e_h1, e_v, e_h2) up-direction components, or None on
    vertical-line/no-metric runs (the scalar, axis-aligned H^-1 applies)."""
    if not cfg.general:
        return None
    e = cfg.metric.e_up
    ax_h1, ax_v, ax_h2 = cfg.role_perm
    return e[ax_h1], e[ax_v], e[ax_h2]


def _apply_hinv(fields3, s, cfg, dt, nonhydro):
    """H^-1 @ (axis-indexed component triple), role plumbing as numpy."""
    ax_h1, ax_v, ax_h2 = cfg.role_perm
    wh1, wv, wh2, nu, _ = _coriolis_inputs(s, cfg, dt, nonhydro)
    U, V, W = fields3[ax_h1], fields3[ax_v], fields3[ax_h2]
    e_role = _e_role(cfg)
    if e_role is None:
        u, v, w = jax_coriolis.apply_inverse(U, V, W, wh1, wh2, wv, nu, nonhydro)
    else:
        e1, e2, e3 = e_role
        u, v, w = jax_coriolis.apply_inverse_general(
            U, V, W, wh1, wh2, wv, e1, e2, e3, nu, nonhydro
        )
    out = [None, None, None]
    out[ax_h1], out[ax_v], out[ax_h2] = u, v, w
    return tuple(out)


def _explicit_part(s, cfg, dt, nonhydro):
    """implicit_euler.do_explicit_part twin."""
    dbuoy = s["rhoY"] * (s["rhoX"] / s["rho"])
    if cfg.general:
        # general map (sphere): the alpha_w discard and the buoyancy kick act
        # on the e_up-PARALLEL momentum, m <- m - (1-alpha)(m.e)e - dt(g/Msq)*dbuoy*e
        e = cfg.metric.e_up
        m_dot_e = s["rhou"] * e[0] + s["rhov"] * e[1] + s["rhow"] * e[2]
        kick = dt * (cfg.g / cfg.Msq) * dbuoy
        for name, ek in zip(_MOMENTA, e):
            s[name] = s[name] + (nonhydro - 1.0) * m_dot_e * ek - kick * ek
    else:
        vmom = _MOMENTA[cfg.v_phys]
        s[vmom] = (nonhydro * s[vmom]) - dt * (cfg.g / cfg.Msq) * dbuoy

    u0, v0, w0 = cfg.winds
    for name, w_ in zip(_MOMENTA, (u0, v0, w0)):
        s[name] = s[name] + (-1.0) * w_ * s["rho"]

    moms = _apply_hinv((s["rhou"], s["rhov"], s["rhow"]), s, cfg, dt, nonhydro)
    s["rhou"], s["rhov"], s["rhow"] = moms

    for name, w_ in zip(_MOMENTA, (u0, v0, w0)):
        s[name] = s[name] + (+1.0) * w_ * s["rho"]
    return s


def _surface_constraint(s, cfg):
    """surface_constraint.apply twin: m <- m - (m.e_up) e_up.

    Active only for thin-shell SWE-on-sphere cases (cfg.constrain_to_surface);
    a no-op otherwise, so every non-SWE device step stays bit-identical.
    Canonical orientation (all call sites are outside the advection sweeps)."""
    if not cfg.constrain_to_surface:
        return s
    e = cfg.metric.e_up
    m_dot_e = s["rhou"] * e[0] + s["rhov"] * e[1] + s["rhow"] * e[2]
    for name, ek in zip(_MOMENTA, e):
        s[name] = s[name] - m_dot_e * ek
    return s


def _polar_filter(s, cfg):
    """polar_filter.apply twin: damp the CFL-violating zonal modes near the
    poles (FFT-in-longitude, J-weighted, k=0 kept -> ring-conserving).

    A no-op unless ``ud.polar_filter`` is set (cfg.filter_plan is None), so
    every non-sphere / channel device step stays bit-identical. Canonical
    orientation (called at the end of the step, mirroring time_update.do)."""
    plan = cfg.filter_plan
    if plan is None:
        return s
    sl, lam, N, r, J = plan.sl, plan.lam_axis, plan.N, plan.r, plan.J
    for name in _SOL_FIELDS:
        f = s[name]
        Jf = J * f[sl]
        F = jnp.fft.rfft(Jf, axis=lam) * r
        s[name] = f.at[sl].set(jnp.fft.irfft(F, n=N, axis=lam) / J)
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
        # general curvilinear contract: area normals + fixed Cartesian axes
        metric = (
            tuple(tuple(Na) for Na in m.N),
            m.cart_v,
            m.cart_haxes,
        )
    else:
        metric = None
    wall_dims = (
        ()
        if cfg.atmosphere
        else tuple(
            d
            for d in range(cfg.ndim)
            if cfg.boundary.bdry_int[d] in (jax_boundary._WALL, jax_boundary._POLE)
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
    # role-ordered rotation components (scalars, or coriolis_field arrays)
    corr_h1, corr_v, corr_h2 = cfg.coriolis_role
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
    mom_h1, mom_v, mom_h2 = mom[names[ax_h1]], mom[names[ax_v]], mom[names[ax_h2]]

    if cfg.general:
        # general map (sphere): buoyancy acts along the LOCAL up e_up and the
        # nonhydro (alpha_w) factor applies to the e-PARALLEL part of the
        # WHOLE tendency (fac_par; the e-perpendicular part is never
        # alpha_w-suppressed). is_ArakawaKonor is guarded to 0 on the device.
        e = cfg.metric.e_up
        e_h1, e_v, e_h2 = e[ax_h1], e[ax_v], e[ax_h2]
        vel_up = (mom_h1 * e_h1 + mom_v * e_v + mom_h2 * e_h2) / rho
        buoy = (cfg.g / cfg.Msq) * dbuoy
        T_h1 = rhoYovG * dp_h1 + buoy * e_h1 - corr_h2 * dm_v + corr_v * dm_h2
        T_v = rhoYovG * dp_v + buoy * e_v - corr_h1 * dm_h2 + corr_h2 * dm_h1
        T_h2 = rhoYovG * dp_h2 + buoy * e_h2 - corr_v * dm_h1 + corr_h1 * dm_v
        T_dot_e = T_h1 * e_h1 + T_v * e_v + T_h2 * e_h2
        fac_par = nonhydro - 1.0
        mom[names[ax_h1]] = mom_h1 - dt * (T_h1 + fac_par * T_dot_e * e_h1)
        mom[names[ax_v]] = mom_v - dt * (T_v + fac_par * T_dot_e * e_v)
        mom[names[ax_h2]] = mom_h2 - dt * (T_h2 + fac_par * T_dot_e * e_h2)
        for n in _MOMENTA:
            s[n] = mom[n]
        s["rhoX"] = (rho * (rho / s["rhoY"] - cfg.S0c)) - dt * (vel_up * cfg.dSdy) * rho
    else:
        vel_v = mom_v / rho
        mom[names[ax_h1]] = mom_h1 - dt * (
            rhoYovG * dp_h1 - corr_h2 * dm_v + corr_v * dm_h2
        )
        mom[names[ax_v]] = (
            mom_v
            - dt
            * (
                rhoYovG * dp_v
                + (cfg.g / cfg.Msq) * dbuoy * nonhydro
                - corr_h1 * dm_h2
                + corr_h2 * dm_h1
            )
            * 1.0
        )  # (1 - is_ArakawaKonor), guarded to 0
        mom[names[ax_h2]] = mom_h2 - dt * (
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
    if cfg.general:
        # general map: stratification couples to the e_up-parallel momentum
        e = cfg.metric.e_up
        m_up = s["rhou"] * e[0] + s["rhov"] * e[1] + s["rhow"] * e[2]
        s["rhoX"] = s["rhoX"] + (-updt_chi) * dt * cfg.dSdy * m_up
    else:
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
    e_role = _e_role(cfg)
    if e_role is None:
        return jax_coriolis.compute_coefficients(wh1, wh2, wv, nu, nonhydro)
    e1, e2, e3 = e_role
    return jax_coriolis.compute_coefficients_general(
        wh1, wh2, wv, e1, e2, e3, nu, nonhydro
    )


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
            geo = terrain_mod.elliptic_diag_geometric(cfg.metric)
            diag_inv = _prepare_diag(
                cfg,
                None,
                wcenter,
                cii=[wplus_coeff * geo[0], wplus_coeff * geo[1], None],
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

        def base_matvec(pvec):
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

        # Pole-ring collapse: the lap3D pole rows are already one-sided
        # (the wall_mask slab-zeroes the non-periodic phi axis);
        # wrap the matvec + gather the rhs with the Galerkin scatter/gather so
        # each pole ring solves as one master. Same recipe as the hybrid
        # jax_lap3D.wrap_pole_collapse (mask-multiply, not integer scatter-SET).
        pc = cfg.pole_collapse
        if pc is not None:

            def matvec(pvec):
                y = jnp.reshape(base_matvec(pvec[pc.scatter_src]), (-1,))
                contrib = y[pc.uniq_mem]
                y = y * pc.ring_complement
                return y.at[pc.uniq_master].add(contrib)

            rhs_c = rhs_inner[pc.uniq_mem]
            rhs_inner = (rhs_inner * pc.ring_complement).at[pc.uniq_master].add(rhs_c)
        else:
            matvec = base_matvec

    A = lambda x: jnp.reshape(matvec(x), x.shape)
    p2, _ = jax.scipy.sparse.linalg.bicgstab(
        A, rhs_inner, tol=cfg.rtol, atol=cfg.tol, maxiter=cfg.max_iterations
    )

    # pole collapse: scatter the master value back over its ring so the
    # pressure is single-valued at the pole (3D only; 2D never sees poles)
    if ndim == 3 and cfg.pole_collapse is not None:
        p2 = p2[cfg.pole_collapse.scatter_src]

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


def build_step(cfg, parity, is_nonhydrostatic, is_compressible):
    """The raw (untraced) full-step closure for one Strang parity + regime
    structure — jitted by :func:`make_step`, vmapped over the member axis by
    ``device_batch.make_batch_step``.

    Both regime ints are STATIC trace structure: the psinc leg of a blending
    window compiles its own variant rather than branching at runtime — the
    validated compressible graph is untouched.

    The psinc variant (is_compressible == 0) returns ``(s, p2_half)`` — the
    predictor half-time nodal pressure that numpy stores as
    ``npf.p2_nodes_half`` (predictor_half_step) and the psinc->comp blend
    conversion reads; the compressible variant returns ``s`` alone,
    unchanged."""

    def step(s, dt, nonhydro, compressibility, forcing_half, forcing_full):
        sol0 = dict(s)  # free reference hold (incl. p2 for the fill helper)

        flux_rhoY = _advective_flux(s, cfg)
        if cfg.do_advection:
            s = _advect_rk(s, cfg, 0.5 * dt, flux_rhoY)
        s = _surface_constraint(s, cfg)

        p2_nodes0 = s["p2_nodes"]

        s = _explicit_part(s, cfg, 0.5 * dt, nonhydro)
        s = _surface_constraint(s, cfg)
        if is_compressible == 0:
            sol0 = _ghost_fill(sol0, cfg)  # numpy fills the held copy
        s = _implicit_part(
            s,
            cfg,
            0.5 * dt,
            nonhydro,
            compressibility,
            sol0=sol0 if is_compressible == 0 else None,
        )

        if cfg.rayleigh_bdry:
            s = _rayleigh_damp(s, cfg)
        if cfg.has_forcing:
            s = _rayleigh_damp(s, cfg, forcing=forcing_half)
        s = _surface_constraint(s, cfg)

        # the Strang advective flux is computed from the POST-implicit
        # half-time state (time_update line 123), before the sol restore
        flux_rhoY_half = _advective_flux(s, cfg)

        # numpy: predictor_half_step stores p2_nodes_half = copy(p2_nodes)
        # at exactly this point (after the half-time flux recompute)
        p2_half = s["p2_nodes"] if is_compressible == 0 else None

        # p2 reset branch (static regime structure); hydrostatic (alpha_w = 0)
        # excluded so the hydrostatic predictor's balanced-Exner
        # reconstruction is not discarded — mirrors time_update.py. Device
        # rejects the hydrostatic regime, so alpha_w = 0 never reaches here;
        # kept consistent for the numpy<->jax-device contract.
        if is_compressible == 1 and is_nonhydrostatic == 1:
            s["p2_nodes"] = p2_nodes0

        # restore sol to t^n for the full forward pass
        for n in _SOL_FIELDS:
            s[n] = sol0[n]

        s = _forward_step(s, cfg, 0.5 * dt, nonhydro, compressibility)
        s = _surface_constraint(s, cfg)

        if cfg.do_advection:
            s = _advect_strang(s, cfg, dt, parity, flux_rhoY_half)
        s = _surface_constraint(s, cfg)

        s = _explicit_part_post(s, cfg, dt, nonhydro, compressibility)

        if cfg.rayleigh_bdry:
            s = _rayleigh_damp(s, cfg)
        if cfg.has_forcing:
            s = _rayleigh_damp(s, cfg, forcing=forcing_full)
        s = _surface_constraint(s, cfg)

        if cfg.diffusion:
            s = _diffuse(s, cfg, dt)

        # Polar filter: damp the CFL-violating zonal modes near the poles
        # once per step, then re-apply the tangent-plane surface
        # constraint (the filter mixes Cartesian momentum components per ring,
        # nudging them off the local tangent plane) and refill the ghosts so
        # the pole exchange sees the filtered interior. Mirrors the numpy
        # attach point in time_update.do (filter -> surface -> ghost refill),
        # which device.run_window returns BEFORE. No-op unless ud.polar_filter.
        if cfg.polar_filter is not None:
            s = _polar_filter(s, cfg)
            s = _surface_constraint(s, cfg)
            s = _ghost_fill(s, cfg)

        if is_compressible == 0:
            return s, p2_half
        return s

    return step


def make_step(cfg, parity, is_nonhydrostatic, is_compressible):
    """Build the jitted full step for one Strang parity + regime structure."""
    return jax.jit(
        build_step(cfg, parity, is_nonhydrostatic, is_compressible),
        donate_argnums=(0,),
    )


def _explicit_part_post(s, cfg, dt, nonhydro, compressibility):
    s = _explicit_part(s, cfg, 0.5 * dt, nonhydro)
    s = _surface_constraint(s, cfg)
    s = _implicit_part(s, cfg, 0.5 * dt, nonhydro, compressibility, sol0=None)
    return s
