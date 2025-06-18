import copy
import numpy as np
from ....utils import io
from ....utils import options as opts
from . import cell_boundary as bdry_c

def get_tau_y(ud, elem, node, alpha):
    tauc_y = np.zeros_like(elem.y)
    taun_y = np.zeros_like(node.y)

    ud.bny = node.y[-ud.inbcy - 3]

    c1n = node.y <= ud.bny
    ccn = (node.y[:-2] - ud.bny) / (node.y[:-2][-1] - ud.bny)
    c2n = np.logical_and(ccn >= 0.0, ccn <= 0.5)
    c3n = np.logical_and(ccn > 0.5, ccn <= 1.0)

    taun_y[np.where(c1n)] = 0.0
    taun_y[np.where(c2n)] = (
        -alpha
        / 2.0
        * (
            1.0
            - np.cos((node.y[np.where(c2n)] - ud.bny) / (node.y[-1] - ud.bny) * np.pi)
        )
    )
    taun_y[np.where(c3n)] = (
        -alpha
        / 2.0
        * (
            1.0
            + ((node.y[np.where(c3n)] - ud.bny) / (node.y[-1] - ud.bny) - 0.5) * np.pi
        )
    )

    taun_y[-2:] = -np.abs(taun_y).max()
    tauc_y[...] = np.interp(elem.y, node.y, taun_y)

    return tauc_y, taun_y


def get_bottom_tau_y(ud, elem, node, alpha, cutoff=0.5):
    tauc_y = np.zeros_like(elem.y)
    taun_y = np.zeros_like(node.y)

    assert ud.ymax > cutoff, "rayleigh forcing boundary below minimum domain extent"
    idx = (np.abs(elem.y - (ud.ymax - cutoff))).argmin()

    ud.forcing_bny = node.y[idx]

    c1n = node.y <= ud.forcing_bny
    ccn = (node.y[:-3] - ud.forcing_bny) / (node.y[:-3][-1] - ud.forcing_bny)
    c2n = np.logical_and(ccn >= 0.0, ccn <= 0.5)
    c3n = np.logical_and(ccn > 0.5, ccn <= 1.0)

    taun_y[np.where(c1n)] = 0.0
    taun_y[np.where(c2n)] = (
        -alpha
        / 2.0
        * (
            1.0
            - np.cos(
                (node.y[np.where(c2n)] - ud.forcing_bny)
                / (node.y[-1] - ud.forcing_bny)
                * np.pi
            )
        )
    )
    taun_y[np.where(c3n)] = (
        -alpha
        / 2.0
        * (
            1.0
            + (
                (node.y[np.where(c3n)] - ud.forcing_bny) / (node.y[-1] - ud.forcing_bny)
                - 0.5
            )
            * np.pi
        )
    )

    taun_y[-3:] = -np.abs(taun_y).max()
    taun_y[...] = taun_y[::-1]
    tauc_y = np.interp(elem.y, node.y, taun_y)

    dd = 1.0
    tauc_y = dd * tauc_y / np.abs(tauc_y).max()
    taun_y = dd * taun_y / np.abs(taun_y).max()

    # tauc_y = tauc_y[2:-2]
    # taun_y = taun_y[2:-2]

    return tauc_y, taun_y


def apply_rayleigh_forcing(
    Sol,
    npf,
    ud,
    elem,
    node,
    t,
    step,
    dt,
    th,
    half=True,
    Sol_half_new=None,
    npf_half_new=None,
):
    """Apply Rayleigh forcing boundary condition (file or function based)."""
    if not (hasattr(ud, "rayleigh_forcing") and ud.rayleigh_forcing):
        return

    t_offset = 0.5 * dt if half else dt

    if ud.rayleigh_forcing_type == "file":
        reader = io.read_input(ud.rayleigh_forcing_fn, ud.rayleigh_forcing_path)

        if Sol_half_new is None or npf_half_new is None:
            Sol_half_new = copy.deepcopy(Sol)
            npf_half_new = copy.deepcopy(npf)

        time_tag = "%.3d_after_full_step" % step
        reader.get_data(Sol_half_new, npf_half_new, time_tag, half=half)

        up = Sol_half_new.rhou / Sol_half_new.rho
        vp = Sol_half_new.rhov / Sol_half_new.rho
        Yp = Sol_half_new.rhoY / Sol_half_new.rho - npf.HydroState.Y0.reshape(1, -1)
        pi = npf_half_new.p2_nodes

        rayleigh_damping(Sol, npf, ud, elem, node, [up, vp, Yp, pi, t + t_offset])

    elif ud.rayleigh_forcing_type == "func":
        s = 5.0e-3 + 1e-4 + 0e-5
        ud.rf_bot.eigenfunction(t + t_offset, s)
        up, vp, Yp, pi = ud.rf_bot.dehatter(th)

        ud.rf_bot.eigenfunction(t + t_offset, s, grid="n")
        _, _, _, pi_n = ud.rf_bot.dehatter(th, grid="n")

        rayleigh_damping(
            Sol, npf, ud, elem, node, [up, vp, Yp, pi_n, t + t_offset]
        )

    bdry_c.set_explicit_boundary_data(Sol, elem, ud, th, npf)

def rayleigh_damping(Sol, npf, ud, elem, node, forcing=None):
    u = Sol.rhou / Sol.rho  # [elem.i2]
    v = Sol.rhov / Sol.rho  # [elem.i2]
    Y = Sol.rhoY / Sol.rho  # [elem.i2]
    rho = Sol.rho  # [elem.i2]

    if ud.bdry_type[1] == opts.BdryType.RAYLEIGH:
        tcy, tny = ud.tcy, ud.tny
    else:
        tcy, tny = 0.0, 0.0

    if forcing is not None:
        tcy_f, tny_f = ud.forcing_tcy, ud.forcing_tny
        tcy, tny = 0.0, 0.0

        u_f, v_f, Y_f, pi_f, t = forcing

        if ud.rayleigh_forcing_type == "file":
            G = np.sqrt(9.0 / 40.0)
            N = np.sqrt(ud.Nsq_ref)
            C = ud.Cs * ud.u_ref
            Gam = N * G / C
            Om = ud.coriolis_strength[2] / 2.0 / ud.t_ref
            # Om = 7.292 * 1e-5 / ud.t_ref
            growth_rate = np.sqrt(Om * C * Gam)
            mfac = np.exp(growth_rate * t * ud.t_ref)

            # if np.all(ud.coriolis_strength) == 0.0:
            # if ud.trad_forcing:
            # mfac = 1.0
        else:
            mfac = 1.0

        npf.p2_nodes[...] += tny_f * (npf.p2_nodes) + np.abs(tny_f) * mfac * pi_f
        c_f = 1.0

    else:
        u_f, v_f, Y_f = 0.0, 0.0, 0.0
        c_f = 0.0
        tcy_f, tny_f = 0.0, 0.0
        mfac = 0.0

    # assuming 2D vertical slice - not dimension agnostic
    u += tcy * (u - ud.u_wind_speed) + c_f * (
        tcy_f * (u - ud.u_wind_speed) + np.abs(tcy_f) * mfac * u_f
    )
    v += tcy * (v - ud.v_wind_speed) + c_f * (
        tcy_f * (v - ud.v_wind_speed) + np.abs(tcy_f) * mfac * v_f
    )

    Ybar = npf.HydroState.Y0.reshape(1, -1)
    Y += tcy * (Y - Ybar) + c_f * (tcy_f * (Y - Ybar) + np.abs(tcy_f) * mfac * Y_f)

    Sol.rhou[...] = rho * u
    Sol.rhov[...] = rho * v
    Sol.rhoY[...] = rho * Y

