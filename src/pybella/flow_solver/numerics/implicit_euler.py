import numpy as np
import scipy as sp

from ...utils import axes
from ...backends import is_jax_backend
from ...utils.operators import convolution, divergence, gradient
from ...utils.operators.laplacian import preconditioner, lap2D_manual, lap3D
from ..discretisation import terrain
from ..utils.boundary import cell_boundary as bdry_c
from ..utils.boundary import node_boundary as bdry_n
from ..utils.boundary import common as bdry
from . import coriolis


def _jax_backend(ud):
    """True when the elliptic solve runs on the JAX backend (ud.backend).

    Coefficient assembly, preconditioning and the rhs stay on the numpy
    path either way — only the operator application and the Krylov
    iteration are swapped (see backends/jax_ops/elliptic_solve.py).
    """
    return is_jax_backend(ud)


class solver_counter(object):
    """
    taken from https://stackoverflow.com/questions/33512081/getting-the-number-of-iterations-of-scipys-gmres-iterative-method

    """

    def __init__(self):
        self.niter = 0

    def __call__(self, rk=None):
        self.niter += 1
        self.rk = rk


def do_explicit_part(mem, ud, dt):
    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[axes.vertical_axis(ud)]
    Msq = ud.Msq

    dbuoy = mem.sol.rhoY * (mem.sol.rhoX / mem.sol.rho)
    vmom = axes.vertical_momentum(ud)
    setattr(mem.sol, vmom, (nonhydro * getattr(mem.sol, vmom)) - dt * (g / Msq) * dbuoy)

    mem.sol.mod_bg_wind(ud, -1.0)

    coriolis.multiply_inverse_terms(mem.sol, mem, ud, dt)

    mem.sol.mod_bg_wind(ud, +1.0)


def do_implicit_part(
    mem,
    ud,
    dt,
    sol0=None,
    writer=None,
    label=None,
    debug=False,
):
    """
    Optimized version with reduced redundancy and improved structure.
    """
    # Early return optimization - disable writer if not debugging
    if not debug:
        writer = None

    nc = mem.node.sc

    # Helper function to reduce writer boilerplate
    def write_debug_data(key, data):
        if writer is not None:
            writer.populate(str(label), key, data)

    # Initial debug output
    write_debug_data("p2_initial", mem.npf.p2_nodes)

    # Set boundary data and compute operator coefficients (consolidated)
    sol_for_boundary = sol0 if sol0 is not None else mem.sol
    bdry_c.set_ghost_cells(mem, ud, sol=sol_for_boundary)
    operator_coefficients_nodes(mem, ud, dt)

    # Debug output for w components
    if writer is not None:
        write_debug_data("hcenter", mem.npf.wcenter)
        write_debug_data("wplusx", mem.npf.wplus[0])
        write_debug_data("wplusy", mem.npf.wplus[1])

        # Handle 3D case more cleanly
        wplusz_data = (
            mem.npf.wplus[2] if mem.elem.ndim == 3 else np.zeros_like(mem.npf.wplus[0])
        )
        write_debug_data("wplusz", wplusz_data)

    # Boundary and correction operations
    # bdry.set_ghost_nodes(mem.npf.p2_nodes, mem.node, ud)
    _correction_nodes(mem, ud, dt, mem.npf.p2_nodes, 0)
    bdry_c.set_ghost_cells(mem, ud)

    # Compute RHS
    mem.npf.rhs[...] = divergence.compute_at_nodes(mem.npf.rhs, mem.elem, mem.sol, ud)
    write_debug_data("rhs", mem.npf.rhs)

    mem.npf.rhs /= dt

    # Handle compressibility - simplified logic
    _apply_compressibility_correction(mem, ud)

    write_debug_data("rhs_nodes", mem.npf.rhs)

    # Prepare and solve linear system
    lap, rhs_inner = _prepare_linear_system(mem, ud, dt)

    # Solve using BiCGSTAB
    if _jax_backend(ud):
        from ...backends.jax_ops import elliptic_solve

        p2 = elliptic_solve.bicgstab(
            lap, rhs_inner, atol=ud.tol, maxiter=ud.max_iterations
        )
    else:
        counter = solver_counter()
        p2, _ = sp.sparse.linalg.bicgstab(
            lap, rhs_inner, atol=ud.tol, maxiter=ud.max_iterations, callback=counter
        )

    # Reshape solution and apply
    p2_full = _reshape_solution(p2, mem, ud, nc)
    write_debug_data("p2_full", p2_full)

    # # Final boundary and correction operations
    # bdry.set_ghost_nodes(p2_full, mem.node, ud)
    _correction_nodes(mem, ud, dt, p2_full, 1)

    mem.npf.p2_nodes[...] += p2_full
    bdry_n.set_ghost_nodes(mem.npf.p2_nodes, mem.node, ud)
    bdry_c.set_ghost_cells(mem, ud)


def _correction_nodes(mem, ud, dt, p, updt_chi):
    Gammainv = mem.th.Gammainv

    dSdy = mem.npf.HydroState_n.get_dSdy(mem.elem, mem.node)

    Dpx, Dpy, Dpz = gradient.compute_at_nodes(p, mem.elem.ndim, mem.node.dxyz)
    if mem.elem.metric is not None:
        # physical gradients via the terrain map A — the same correction the
        # elliptic operator's C_ij coefficients encode, so the projection
        # annihilates exactly the divergence it measures
        Dpx, Dpy, Dpz = terrain.apply_gradient_map(mem.elem.metric, [Dpx, Dpy, Dpz])

    thinv = mem.sol.rho / mem.sol.rhoY

    Y = mem.sol.rhoY / mem.sol.rho
    coeff = Gammainv * mem.sol.rhoY * Y

    mem.npf.u[...] = -dt * coeff * Dpx
    mem.npf.v[...] = -dt * coeff * Dpy
    mem.npf.w[...] = -dt * coeff * Dpz

    coriolis.multiply_inverse_terms(mem.npf, mem, ud, dt, attrs=["u", "v", "w"])

    mem.sol.rhou += thinv * mem.npf.u
    mem.sol.rhov += thinv * mem.npf.v
    # the w-row applies in 2D too: H^-1 rotates the pressure correction into
    # the out-of-plane momentum whenever Coriolis is active. Restricting it
    # to ndim == 3 dropped that component in 2D runs (implicit-side sibling
    # of the explicit-step defect fixed 2026-06-09; quantified at 1.4e-4 by
    # the 3D-vs-2D full-Coriolis oracle).
    mem.sol.rhow += thinv * mem.npf.w
    mem.sol.rhoX += -updt_chi * dt * dSdy * getattr(mem.sol, axes.vertical_momentum(ud))


def operator_coefficients_nodes(mem, ud, dt):
    Gammainv = mem.th.Gammainv
    ndim = mem.node.ndim

    ccenter = -ud.Msq * mem.th.gm1inv / (dt**2)
    cexp = 2.0 - mem.th.gamm

    Y = mem.sol.rhoY / mem.sol.rho
    coeff = Gammainv * mem.sol.rhoY * Y

    for dim in range(ndim):
        mem.npf.wplus[dim][...] = coeff

    kernel = convolution.get_averaging_kernel(ndim, width=2)

    mem.npf.wcenter = ccenter * convolution.apply_convolution_kernel(
        mem.sol.rhoY**cexp, kernel
    )

    if mem.node.metric is not None:
        # with terrain the solved equation is J * (Helmholtz): the rhs
        # carries J*div F, the C_ij tensor carries J — the pointwise center
        # term needs its J too (node-exact, not averaged)
        mem.npf.wcenter *= mem.node.metric.J[mem.node.i1]

    if not hasattr(ud, "ATMOSPHERIC_EXTENSION"):
        bdry.scale_wall_node_values(mem.npf.wcenter, mem.node, ud)


def _apply_compressibility_correction(mem, ud):
    """
    Extracted compressibility logic for better readability.
    """
    if ud.is_compressible == 0:
        if ud.is_ArakawaKonor:
            mem.npf.rhs -= mem.npf.wcenter * mem.npf.dp2_nodes
            mem.npf.wcenter[...] = 0.0
        else:
            # Simplified - this was redundant multiplication
            mem.npf.wcenter[...] *= ud.compressibility
    else:
        mem.npf.wcenter *= ud.compressibility


def _prepare_linear_system(mem, ud, dt):
    """
    Prepare the linear system components based on dimensionality.
    """
    if mem.elem.ndim == 2:
        return _prepare_2d_system(mem, ud, dt)
    else:  # 3D case
        return _prepare_3d_system(mem, ud, dt)


def _prepare_2d_system(mem, ud, dt):
    """Prepare 2D linear system.

    The lap2D kernel's effective coefficients are wplus ⊙ c (c = the four
    coefficient arrays passed as ``coriolis_params``), so with terrain the
    2x2 tensor M = J A^T H^-1 A is folded into c — mirroring the 3D
    ``cij = wplus[i] * h_role[i][j]`` — while wplus stays separate to keep
    the kernel's wall coefficient-zeroing intact. ``multiply_inverse_terms``
    returns the (h1, v) block TRANSPOSED (C-ravel of the transpose == the
    F-ravel the wplus arrays get), so the fold un-transposes around
    ``elliptic_tensor_2d``.
    """
    Vec = mem.npf
    coriolis_params = coriolis.multiply_inverse_terms(
        Vec, mem, ud, dt, attrs=("u", "v", "w"), get_coeffs=True
    )

    if mem.elem.metric is not None:
        h11_t, h22_t, h12_t, h21_t = coriolis_params
        h2x2 = ((h11_t.T, h12_t.T), (h21_t.T, h22_t.T))
        M = terrain.elliptic_tensor_2d(mem.elem.metric, h2x2)
        coriolis_params = (M[0][0].T, M[1][1].T, M[0][1].T, M[1][0].T)
        # diag: fold only the geometric factors (the legacy 2D preconditioner
        # keeps H^-1 out of the diagonal — preserved here so a forced-flat
        # metric preconditions bit-identically to the plain path)
        geo = terrain.elliptic_diag_geometric(mem.elem.metric)
        diag_inv = preconditioner.prepare_diag(
            mem.npf,
            mem.node,
            cii=(mem.npf.wplus[0] * geo[0], mem.npf.wplus[1] * geo[1], None),
        )
    else:
        diag_inv = preconditioner.prepare_diag(mem.npf, mem.node)
    mem.npf.rhs *= diag_inv

    rhs_inner = mem.npf.rhs[mem.node.i1].T.ravel()

    if _jax_backend(ud):
        from ...backends.jax_ops.laplacian import lap2D as jax_lap2D

        return (
            jax_lap2D.get_linop(mem.npf, mem.node, coriolis_params, diag_inv, ud),
            rhs_inner,
        )

    p2 = mem.npf.p2_nodes[mem.node.i2].T
    lap = lap2D_manual.get_linop(mem.npf, mem.node, coriolis_params, diag_inv, ud)
    sh = p2.shape[0] * p2.shape[1]

    lap = sp.sparse.linalg.LinearOperator((sh, sh), lap)

    return lap, rhs_inner


def _prepare_3d_system(mem, ud, dt):
    """Prepare 3D linear system.

    The solve vector is the full node.isc box (interior nodes plus one
    ghost layer per side) in C order. The ghost ring carries zero operator
    rows and zero rhs entries, so it stays exactly zero through BiCGSTAB.

    The operator carries the full H^-1 tensor coefficients C_ij =
    (Gamma^-1 P Theta) * h[role(i), role(j)] — the same H^-1 applied by
    _correction_nodes — so the elliptic solve is consistent with the
    momentum correction (the legacy operator had only ad-hoc x-z cross
    terms). With no rotation and no buoyancy H^-1 is the identity and the
    operator reduces bit-exactly to the plain Laplacian.
    """
    hv = coriolis.compute_inverse_coefficients(mem, ud, dt)
    h_role = ((hv[0], hv[1], hv[2]), (hv[3], hv[4], hv[5]), (hv[6], hv[7], hv[8]))
    if mem.elem.metric is not None:
        # terrain: the operator tensor is J A^T H^-1 A — the metric map A of
        # the momentum correction composed with the J-weighted divergence
        h_role = terrain.elliptic_tensor(mem.elem.metric, h_role)
    rho_of = axes.role_of_axis(axes.vertical_axis(ud))
    cij = [
        [mem.npf.wplus[i] * h_role[rho_of[i]][rho_of[j]] for j in range(3)]
        for i in range(3)
    ]

    diag_inv = preconditioner.prepare_diag(
        mem.npf, mem.node, cii=(cij[0][0], cij[1][1], cij[2][2])
    )
    mem.npf.rhs *= diag_inv

    rhs_inner = np.zeros_like(mem.npf.rhs)
    rhs_inner[mem.node.i1] = mem.npf.rhs[mem.node.i1]

    if _jax_backend(ud):
        from ...backends.jax_ops.laplacian import lap3D as jax_lap3D

        return (
            jax_lap3D.get_linop(mem.elem, mem.node, mem.npf, ud, diag_inv, dt, cij),
            rhs_inner.ravel(),
        )

    lap = lap3D.get_linop(mem.elem, mem.node, mem.npf, ud, diag_inv, dt, cij)
    sh = mem.npf.rhs.size

    lap = sp.sparse.linalg.LinearOperator((sh, sh), lap, dtype=np.float64)

    return lap, rhs_inner.ravel()


def _reshape_solution(p2, mem, ud, nc):
    """
    Reshape the solution vector back to the appropriate format.
    """
    p2_full = np.zeros(nc).squeeze()

    if mem.elem.ndim == 2:
        p2_full[mem.node.i2] = p2.reshape(mem.npf.rhs[mem.node.i1].T.shape).T
    else:  # 3D case: solution vector is the C-ordered node.isc box
        p2_full[mem.node.i1] = p2.reshape(mem.npf.rhs.shape)

    return p2_full
