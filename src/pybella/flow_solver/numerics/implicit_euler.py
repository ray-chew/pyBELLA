import numpy as np
import scipy as sp

from ...utils.operators import convolution, divergence, gradient
from ...utils.operators.laplacian import preconditioner, lap2D_manual, lap3D
from ..utils.boundary import cell_boundary as bdry_c
from ..utils.boundary import node_boundary as bdry_n
from ..utils.boundary import common as bdry
from . import coriolis


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
    g = ud.gravity_strength[1]
    Msq = ud.Msq

    dbuoy = mem.sol.rhoY * (mem.sol.rhoX / mem.sol.rho)
    mem.sol.rhov = (nonhydro * mem.sol.rhov) - dt * (g / Msq) * dbuoy

    mem.sol.mod_bg_wind(ud, -1.0)

    coriolis.multiply_inverse_terms(mem.sol, mem, ud, dt)

    mem.sol.mod_bg_wind(ud, +1.0)

    bdry_c.set_ghost_cells(mem, ud)


def do_implicit_part(
    mem,
    ud,
    dt,
    Sol0=None,
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
    sol_for_boundary = Sol0 if Sol0 is not None else mem.sol
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
    # bdry.set_ghostnodes_p2(mem.npf.p2_nodes, mem.node, ud)
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
    counter = solver_counter()
    p2, _ = sp.sparse.linalg.bicgstab(
        lap, rhs_inner, atol=ud.tol, maxiter=ud.max_iterations, callback=counter
    )

    # Reshape solution and apply
    p2_full = _reshape_solution(p2, mem, ud, nc)
    write_debug_data("p2_full", p2_full)

    # # Final boundary and correction operations
    # bdry.set_ghostnodes_p2(p2_full, mem.node, ud)
    _correction_nodes(mem, ud, dt, p2_full, 1)

    mem.npf.p2_nodes[...] += p2_full
    bdry_n.set_ghostnodes_p2(mem.npf.p2_nodes, mem.node, ud)
    bdry_c.set_ghost_cells(mem, ud)


def _correction_nodes(mem, ud, dt, p, updt_chi):
    ndim = mem.node.ndim
    Gammainv = mem.th.Gammainv

    dSdy = mem.npf.HydroState_n.get_dSdy(mem.elem, mem.node)

    Dpx, Dpy, Dpz = gradient.compute_at_nodes(p, mem.elem.ndim, mem.node.dxyz)

    thinv = mem.sol.rho / mem.sol.rhoY

    Y = mem.sol.rhoY / mem.sol.rho
    coeff = Gammainv * mem.sol.rhoY * Y

    mem.npf.u[...] = -dt * coeff * Dpx
    mem.npf.v[...] = -dt * coeff * Dpy
    mem.npf.w[...] = -dt * coeff * Dpz

    coriolis.multiply_inverse_terms(mem.npf, mem, ud, dt, attrs=["u", "v", "w"])

    mem.sol.rhou += thinv * mem.npf.u
    mem.sol.rhov += thinv * mem.npf.v
    mem.sol.rhow += thinv * mem.npf.w if ndim == 3 else 0.0
    mem.sol.rhoX += -updt_chi * dt * dSdy * mem.sol.rhov


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
    """Prepare 2D linear system."""
    Vec = mem.npf
    coriolis_params = coriolis.multiply_inverse_terms(
        Vec, mem, ud, dt, attrs=("u", "v", "w"), get_coeffs=True
    )

    diag_inv = preconditioner.prepare_diag(mem.npf, mem.node)
    mem.npf.rhs *= diag_inv

    p2 = mem.npf.p2_nodes[mem.node.i2].T
    lap = lap2D_manual.get_linop(mem.npf, mem.node, coriolis_params, diag_inv, ud)
    sh = p2.shape[0] * p2.shape[1]

    lap = sp.sparse.linalg.LinearOperator((sh, sh), lap)
    rhs_inner = mem.npf.rhs[mem.node.i1].T.ravel()

    return lap, rhs_inner


def _prepare_3d_system(mem, ud, dt):
    """Prepare 3D linear system."""
    # Note: diag_inv appears to be used but not defined in 3D case
    # This might be a bug in the original code
    diag_inv = None  # TODO: Verify if this should be computed for 3D

    lap = lap3D.get_linop(mem.elem, mem.node, mem.npf, ud, diag_inv, dt)
    p2 = mem.npf.p2_nodes  # Define p2 for 3D case
    sh = p2.reshape(-1).shape[0]

    lap = sp.sparse.linalg.LinearOperator((sh, sh), lap)
    rhs_inner = mem.npf.rhs[mem.node.i1].ravel()

    return lap, rhs_inner, sh


def _reshape_solution(p2, mem, ud, nc):
    """
    Reshape the solution vector back to the appropriate format.
    """
    p2_full = np.zeros(nc).squeeze()

    if mem.elem.ndim == 2:
        p2_full[mem.node.i2] = p2.reshape(mem.npf.rhs[mem.node.i1].T.shape).T
    else:  # 3D case
        p2_full[mem.node.i1] = p2.reshape(ud.inx + 2, ud.iny + 2, ud.inz + 2)

    return p2_full
