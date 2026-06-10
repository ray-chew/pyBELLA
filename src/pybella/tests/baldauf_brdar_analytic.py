"""Linear reference solution for the Baldauf & Brdar (2013) IGW channel test.

Builds a numerically-exact solution of the *linearised* compressible Euler
equations about the isothermal hydrostatic background of
``test_igw_baldauf_brdar`` (f-plane Coriolis included) and evolves the
simulation's own saved initial condition, providing an independent physics
oracle for the nonlinear solver: golden masters catch *change*, this catches
*wrongness* (sign, dispersion, rotation errors).

Method (instead of transcribing B&B's closed form, which hinges on delicate
rigid-lid basis choices):

1.  Bretherton transform: with rho0(z) = rho_s exp(-z/H), the substitution
    (u, vo, w, b) = exp(z/2H) (U, Vo, W, B),  p' = rho_s exp(-z/2H) P
    turns the linear system into one with CONSTANT coefficients:

        dU/dt  = -ik P - f Vo
        dVo/dt = +f U
        dW/dt  = -dP/dz + alpha P + B
        dB/dt  = -N^2 W
        dP/dt  = -c^2 ik U - c^2 dW/dz - c^2 alpha W

    with alpha = 1/(2H) - g/c^2, N^2 = (gamma-1) g^2 / (gamma R T0), and
    b = -g (rho' - p'/c^2) / rho0 the buoyancy. (Coriolis signs follow
    pyBELLA's convention: du/dt = -f vo, dvo/dt = +f u, rotation about the
    vertical y-axis; vo is pyBELLA's out-of-plane w.)

2.  x is periodic: exact FFT decomposition; each k evolves independently.

3.  z: staggered collocation on a grid ``refine``x finer than the simulation
    (W, B on interior interfaces -- the rigid lids W=0 are imposed by
    construction; U, Vo, P at cell centres; 2nd-order staggered derivatives,
    energy-neutral discretisation).

4.  Time: exact, via the matrix exponential exp(A t) per mode.

The comparison is restricted to the x-mean-free (k != 0) wave content: the
k = 0 column is a hydrostatically balanced state that is frozen in the
linear dynamics, and subtracting the instantaneous x-mean on both sides
makes the comparison immune to the (discrete-vs-continuum) background
residuals of the initial condition.
"""

import numpy as np
import scipy as sp


class IGWParams:
    """SI parameters of the test_igw_baldauf_brdar setup, derived from ud."""

    def __init__(self, ud):
        self.gamma = ud.gamm
        self.R = ud.R_gas
        self.T0 = ud.T_ref
        self.g = ud.grav
        self.p_s = ud.p_ref
        self.f = ud.omega
        self.h_ref = ud.h_ref
        self.t_ref = ud.t_ref
        self.u_ref = ud.u_ref
        self.rho_ref = ud.p_ref / (ud.R_gas * ud.T_ref)

        self.c2 = self.gamma * self.R * self.T0
        self.H = self.R * self.T0 / self.g  # density/pressure scale height
        self.N2 = (self.gamma - 1.0) * self.g**2 / (self.gamma * self.R * self.T0)
        self.alpha = 1.0 / (2.0 * self.H) - self.g / self.c2

    def rho0(self, z):
        return self.rho_ref * np.exp(-z / self.H)

    def p0(self, z):
        return self.p_s * np.exp(-z / self.H)


def _build_operator(k, zc, zi, par):
    """A(k) for the transformed state [U, Vo, P (centres); W, B (interfaces)].

    zc: cell-centre heights (nc,), zi: interior interface heights (ni = nc-1,).
    Returns the complex (3nc + 2ni) square matrix.
    """
    nc, ni = len(zc), len(zi)
    dz = zc[1] - zc[0]
    n = 3 * nc + 2 * ni
    A = np.zeros((n, n), dtype=complex)

    iU = slice(0, nc)
    iV = slice(nc, 2 * nc)
    iP = slice(2 * nc, 3 * nc)
    iW = slice(3 * nc, 3 * nc + ni)
    iB = slice(3 * nc + ni, n)

    I_c = np.eye(nc)
    # interface <- centre difference / average (ni x nc)
    D_ic = np.zeros((ni, nc))
    M_ic = np.zeros((ni, nc))
    for i in range(ni):
        D_ic[i, i] = -1.0 / dz
        D_ic[i, i + 1] = 1.0 / dz
        M_ic[i, i] = 0.5
        M_ic[i, i + 1] = 0.5
    # centre <- interface difference / average (nc x ni), W = 0 at the lids
    D_ci = np.zeros((nc, ni))
    M_ci = np.zeros((nc, ni))
    for j in range(nc):
        if j - 1 >= 0:
            D_ci[j, j - 1] = -1.0 / dz
            M_ci[j, j - 1] = 0.5
        if j < ni:
            D_ci[j, j] = 1.0 / dz
            M_ci[j, j] = 0.5

    ik = 1j * k

    # dU/dt = -ik P - f Vo
    A[iU, iP] = -ik * I_c
    A[iU, iV] = -par.f * I_c
    # dVo/dt = +f U
    A[iV, iU] = +par.f * I_c
    # dP/dt = -c^2 ik U - c^2 (D_ci W) - c^2 alpha (M_ci W)
    A[iP, iU] = -par.c2 * ik * I_c
    A[iP, iW] = -par.c2 * (D_ci + par.alpha * M_ci)
    # dW/dt = -(D_ic P) + alpha (M_ic P) + B
    A[iW, iP] = -D_ic + par.alpha * M_ic
    A[iW, iB] = np.eye(ni)
    # dB/dt = -N^2 W
    A[iB, iW] = -par.N2 * np.eye(ni)

    return A, (iU, iV, iP, iW, iB)


def evolve_linear(fields_ic, x_len, zc_sim, t_end, par, refine=2):
    """Evolve the linear system from the sim's IC perturbations.

    fields_ic: dict with keys 'u', 'vo', 'w', 'p', 'rho' — SI perturbation
    fields on the sim's interior cell grid (nx, nz), x-mean NOT yet removed.
    x_len: domain length [m]; zc_sim: sim cell-centre heights (nz,) [m].
    Returns dict of the same keys evaluated on the sim grid at t_end
    (k != 0 content only).
    """
    nx, nz = fields_ic["u"].shape
    dz_sim = zc_sim[1] - zc_sim[0]

    # fine staggered z grids
    ncf = refine * nz
    dzf = dz_sim / refine
    zcf = (np.arange(ncf) + 0.5) * dzf
    zif = (np.arange(1, ncf)) * dzf

    # buoyancy from rho', p'
    b_ic = (
        -par.g
        * (fields_ic["rho"] - fields_ic["p"] / par.c2)
        / par.rho0(zc_sim)[None, :]
    )

    # remove x-mean (k=0 column is frozen in the linear dynamics)
    def demean(q):
        return q - q.mean(axis=0, keepdims=True)

    u = demean(fields_ic["u"])
    vo = demean(fields_ic["vo"])
    w = demean(fields_ic["w"])
    p = demean(fields_ic["p"])
    b = demean(b_ic)

    # Bretherton transform on the sim grid, then x-FFT
    wgt = np.exp(-zc_sim / (2.0 * par.H))[None, :]
    U0 = np.fft.rfft(u * wgt, axis=0)
    V0 = np.fft.rfft(vo * wgt, axis=0)
    W0 = np.fft.rfft(w * wgt, axis=0)
    B0 = np.fft.rfft(b * wgt, axis=0)
    P0 = np.fft.rfft(p / (par.rho_ref * wgt), axis=0)

    ks = 2.0 * np.pi * np.fft.rfftfreq(nx, d=x_len / nx)
    nk = len(ks)

    # spline-interpolate the (smooth) transformed profiles to the fine grids
    def to_fine(Q, ztarget):
        out = np.empty((nk, len(ztarget)), dtype=complex)
        for m in range(nk):
            re = sp.interpolate.CubicSpline(zc_sim, Q[m].real, bc_type="natural")
            im = sp.interpolate.CubicSpline(zc_sim, Q[m].imag, bc_type="natural")
            out[m] = re(ztarget) + 1j * im(ztarget)
        return out

    Uf, Vf, Pf = to_fine(U0, zcf), to_fine(V0, zcf), to_fine(P0, zcf)
    Wf, Bf = to_fine(W0, zif), to_fine(B0, zif)

    ncf_, nif = len(zcf), len(zif)
    UT = np.zeros_like(Uf)
    VT = np.zeros_like(Vf)
    PT = np.zeros_like(Pf)
    WT = np.zeros_like(Wf)
    BT = np.zeros_like(Bf)

    # skip spectrally-empty modes (the Gaussian envelope kills high k)
    mode_amp = np.zeros(nk)
    for m in range(nk):
        mode_amp[m] = max(
            np.abs(Uf[m]).max(),
            np.abs(Vf[m]).max(),
            np.abs(Wf[m]).max(),
            np.abs(Bf[m]).max() / par.N2**0.5,
            np.abs(Pf[m]).max() / par.c2**0.5,
        )
    # 1e-6 of the peak: the Gaussian envelope's spectrum is ~2e-8 of peak by
    # mode 80; everything above the cut carries no metric-relevant energy
    amp_cut = 1e-6 * mode_amp.max()

    energy_drift = 0.0
    for m in range(1, nk):  # skip k = 0
        if mode_amp[m] < amp_cut:
            continue
        A, (iU, iV, iP, iW, iB) = _build_operator(ks[m], zcf, zif, par)
        q0 = np.concatenate([Uf[m], Vf[m], Pf[m], Wf[m], Bf[m]])
        # exact-in-time evolution via eigendecomposition (A is neutrally
        # stable: similar to skew-Hermitian under the energy weights)
        lam, S = np.linalg.eig(A)
        qt = S @ (np.exp(lam * t_end) * np.linalg.solve(S, q0))

        # energy conservation check (discretisation is energy-neutral)
        def energy(q):
            return (
                np.sum(np.abs(q[iU]) ** 2)
                + np.sum(np.abs(q[iV]) ** 2)
                + np.sum(np.abs(q[iP]) ** 2) / par.c2
                + np.sum(np.abs(q[iW]) ** 2)
                + np.sum(np.abs(q[iB]) ** 2) / par.N2
            )

        e0 = energy(q0)
        if e0 > 0.0:
            energy_drift = max(energy_drift, abs(energy(qt) / e0 - 1.0))

        UT[m], VT[m], PT[m] = qt[iU], qt[iV], qt[iP]
        WT[m], BT[m] = qt[iW], qt[iB]

    # back to physical space on the fine grid
    def from_modes(Q):
        return np.fft.irfft(Q, n=nx, axis=0)

    u_f = from_modes(UT) * np.exp(zcf / (2.0 * par.H))[None, :]
    vo_f = from_modes(VT) * np.exp(zcf / (2.0 * par.H))[None, :]
    p_f = from_modes(PT) * (par.rho_ref * np.exp(-zcf / (2.0 * par.H)))[None, :]
    w_f = from_modes(WT) * np.exp(zif / (2.0 * par.H))[None, :]
    b_f = from_modes(BT) * np.exp(zif / (2.0 * par.H))[None, :]

    # evaluate on the sim cell centres
    def to_sim(qf, zsrc):
        out = np.empty((nx, len(zc_sim)))
        for i in range(nx):
            out[i] = sp.interpolate.CubicSpline(zsrc, qf[i], bc_type="natural")(zc_sim)
        return out

    u_s = to_sim(u_f, zcf)
    vo_s = to_sim(vo_f, zcf)
    p_s = to_sim(p_f, zcf)
    w_s = to_sim(w_f, zif)
    b_s = to_sim(b_f, zif)

    # rho' back from (b, p'):  b = -g (rho' - p'/c^2)/rho0
    rho_s = p_s / par.c2 - par.rho0(zc_sim)[None, :] * b_s / par.g

    return (
        {"u": u_s, "vo": vo_s, "w": w_s, "p": p_s, "rho": rho_s},
        {"energy_drift": energy_drift},
    )


def sim_fields_SI(h5file, tag, par, interior=(slice(2, -2), slice(2, -2))):
    """Extract SI perturbation fields from a pyBELLA igw output file.

    tag: e.g. '000_ic' or '030_after_full_step'. Background = isothermal
    analytic profiles; the comparison layer removes x-means anyway.
    Returns dict u, vo, w, p, rho on the interior cell grid plus zc [m].
    """
    import h5py

    with h5py.File(h5file, "r") as h:
        rho = h["rho"][f"rho_{tag}"][...][interior]
        rhou = h["rhou"][f"rhou_{tag}"][...][interior]
        rhov = h["rhov"][f"rhov_{tag}"][...][interior]
        rhow = h["rhow"][f"rhow_{tag}"][...][interior]
        rhoY = h["rhoY"][f"rhoY_{tag}"][...][interior]

    nz = rho.shape[1]
    zc = (np.arange(nz) + 0.5) * (10000.0 / nz)

    rho_SI = rho * par.rho_ref
    p_SI = rhoY**par.gamma * par.p_s

    return {
        "u": (rhou / rho) * par.u_ref,
        "vo": (rhow / rho) * par.u_ref,  # pyBELLA w = out-of-plane velocity
        "w": (rhov / rho) * par.u_ref,  # pyBELLA v = vertical velocity
        "p": p_SI - par.p0(zc)[None, :],
        "rho": rho_SI - (par.p0(zc) / (par.R * par.T0))[None, :],
    }, zc


class _StubWriter:
    def write(self, *a, **k):
        pass

    def populate(self, *a, **k):
        pass


def run_sim(dt_factor=1, omega=None, nx_factor=1):
    """Run test_igw_baldauf_brdar in-process at dt = 500 s / dt_factor.

    Returns (ud, ic_fields, end_fields, zc, t_end_SI) with fields in SI
    perturbation form on the interior cell grid (see sim_fields_SI).
    The run length is held fixed at 31 * 500 s = 15500 s.
    omega overrides the case's Coriolis parameter (e.g. 0.0 to switch it off
    consistently in sim and reference).
    """
    from .test_igw_baldauf_brdar import UserData, sol_init
    from ..utils import user_data, data_structures
    from ..flow_solver.discretisation import grid as dis_grid
    from ..flow_solver.discretisation import time_update as dis_time_update
    from ..flow_solver.utils import fields as fs_fields
    from ..flow_solver.utils import cache as fs_cache
    from ..flow_solver.utils.boundary import cell_boundary as bdry_c
    from ..flow_solver.physics import thermodynamics as gd_thermodynamics

    ud = user_data.UserDataInit(**vars(UserData()))
    if omega is not None:
        ud.omega = omega
    ud.coriolis_strength = np.array([0.0, ud.omega * ud.t_ref, 0.0])
    if nx_factor != 1:
        # keep inx even (initial_pressure parity assert): 602 -> 601 cells
        ud.inx = 301 * nx_factor
    ud.dtfixed /= dt_factor
    ud.dtfixed0 /= dt_factor
    ud.stepmax = 31 * dt_factor
    ud.diag = False

    elem, node = dis_grid.grid_init(ud)
    sol = fs_fields.CellSolField(elem.sc)
    th = gd_thermodynamics.ThermodynamicalQuantities(ud)
    npf = fs_fields.NodePressureField(elem, node)
    sol = sol_init(sol, npf, elem, node, th, ud)

    mem = data_structures.ModelState(
        elem, node, sol, npf, th, fs_cache.FlowSolverCache()
    )
    bdry_c.set_ghost_cells(mem, ud)

    par = IGWParams(ud)

    def extract(s):
        i2 = (slice(2, -2), slice(2, -2))
        rho = s.rho[i2]
        rhoY = s.rhoY[i2]
        nz = rho.shape[1]
        zc = (np.arange(nz) + 0.5) * (10000.0 / nz)
        return {
            "u": (s.rhou[i2] / rho) * par.u_ref,
            "vo": (s.rhow[i2] / rho) * par.u_ref,
            "w": (s.rhov[i2] / rho) * par.u_ref,
            "p": rhoY**par.gamma * par.p_s - par.p0(zc)[None, :],
            "rho": rho * par.rho_ref - (par.p0(zc) / (par.R * par.T0))[None, :],
        }, zc

    ic, zc = extract(mem.sol)
    mem = dis_time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter())
    end, _ = extract(mem.sol)

    return ud, ic, end, zc, mem.time.t * ud.t_ref


def compare(
    h5file, ud, t_end, tag_ic="000_ic", tag_end="030_after_full_step", refine=2
):
    """Run the full comparison; returns (metrics, sim_end, ref_end, zc)."""
    par = IGWParams(ud)
    L = (ud.xmax - ud.xmin) * ud.h_ref

    ic, zc = sim_fields_SI(h5file, tag_ic, par)
    end, _ = sim_fields_SI(h5file, tag_end, par)

    ref, diag = evolve_linear(ic, L, zc, t_end, par, refine=refine)

    def demean(q):
        return q - q.mean(axis=0, keepdims=True)

    metrics = {"energy_drift": diag["energy_drift"]}
    sim_w = {}
    for key in ("u", "vo", "w", "p", "rho"):
        s = demean(end[key])
        r = ref[key]
        denom = np.linalg.norm(r)
        metrics[key] = np.linalg.norm(s - r) / denom if denom > 0 else np.inf
        sim_w[key] = s

    return metrics, sim_w, ref, zc
