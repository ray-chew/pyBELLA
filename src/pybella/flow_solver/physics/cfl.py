import numpy as np


def dynamic_timestep(Sol, time, time_output, elem, ud, th, step):
    """
    Calculate dynamic timestep for CFD simulation.

    Documentation to be homogenised.

    Args:
        Sol: Solution object containing flow variables
        time: Current simulation time
        time_output: Target output time
        elem: Element object with grid spacing
        ud: User data object with CFL and timestep parameters
        th: Thermodynamic properties
        step: Current step number

    Returns:
        If acoustic_timestep == 1: dt (float)
        Else: tuple of (dt, cfl, cfl_ac)
    """
    machine_epsilon = np.finfo(float).eps

    # Calculate thermodynamic properties
    p = Sol.rhoY**th.gamm
    c = np.sqrt(th.gamm * p / Sol.rho) / np.sqrt(ud.Msq)

    # Calculate velocity components
    u = np.abs(Sol.rhou / Sol.rho)
    v = np.abs(Sol.rhov / Sol.rho)
    w = np.abs(Sol.rhow / Sol.rho)

    # terrain: the coordinate velocity along every sweep axis is
    # xi_a-dot = (N_a . m) / (rho J) and the signal speed gains the
    # face-area/Jacobian factor c |N_a| / J (vertical-line reduction:
    # eta_dot = (w - G.u_h)/J with c sqrt(1 + G^2)/J on the vertical,
    # u with c on the horizontals); the metric is in the unflipped
    # orientation here (called between steps)
    vels = [u, v, w]
    cs = [c, c, c]
    if elem.metric is not None:
        m = elem.metric
        moms = (Sol.rhou, Sol.rhov, Sol.rhow)
        cv = m.cart_v
        ch1, ch2 = m.cart_haxes
        for a in range(elem.ndim):
            Na = m.N[a]
            contra = Na[cv] * moms[cv] + Na[ch1] * moms[ch1]
            norm_sq = Na[cv] ** 2 + Na[ch1] ** 2
            if ch2 is not None:
                contra = contra + Na[ch2] * moms[ch2]
                norm_sq = norm_sq + Na[ch2] ** 2
            vels[a] = np.abs(contra / Sol.rho) * m.ooJ
            cs[a] = c * np.sqrt(norm_sq) * m.ooJ
        # polar filter (Stage F): the longitude modes that survive the
        # filter have an effective speed scaled by ~cos(phi)/cos(phi_c), so
        # the longitude CFL is relieved by that factor near the poles
        from ..numerics import polar_filter as _pf  # local: physics->numerics

        cap = _pf.cfl_cap(elem, ud)
        if cap is not None:
            lam_axis = m.cart_haxes[0]
            vels[lam_axis] = vels[lam_axis] * cap
            cs[lam_axis] = cs[lam_axis] * cap
    u, v, w = vels

    # Find maximum velocities (with minimum threshold)
    u_max = max(u.max(), machine_epsilon)
    v_max = max(v.max(), machine_epsilon)
    w_max = max(w.max(), machine_epsilon)

    # Calculate acoustic velocities
    upc_max = max((u + cs[0]).max(), machine_epsilon)
    vpc_max = max((v + cs[1]).max(), machine_epsilon)
    wpc_max = max((w + cs[2]).max(), machine_epsilon)

    if ud.acoustic_timestep == 1:
        return _calculate_acoustic_timestep(
            ud.CFL,
            elem,
            upc_max,
            vpc_max,
            wpc_max,
            time,
            time_output,
            ud,
            step,
            machine_epsilon,
        )
    else:
        return _calculate_advective_timestep(
            ud.CFL,
            elem,
            u_max,
            v_max,
            w_max,
            upc_max,
            vpc_max,
            wpc_max,
            time,
            time_output,
            ud,
            step,
            machine_epsilon,
        )


def _calculate_acoustic_timestep(
    CFL, elem, upc_max, vpc_max, wpc_max, time, time_output, ud, step, machine_epsilon
):
    """Calculate timestep based on acoustic CFL condition."""
    # Calculate directional timesteps
    dt_directions = [
        CFL * elem.dx / upc_max,
        CFL * elem.dy / vpc_max,
        CFL * elem.dz / wpc_max,
    ]
    dt_cfl = min(dt_directions)

    # Apply ramping factor
    ramp_factor = ud.dtfixed0 + min(step, 1.0) * (ud.dtfixed - ud.dtfixed0)
    dt = min(dt_cfl, ramp_factor)

    # Ensure we don't overshoot the output time
    remaining_time = time_output - time
    if dt > remaining_time:
        dt = remaining_time + machine_epsilon

    return dt


def _calculate_advective_timestep(
    CFL,
    elem,
    u_max,
    v_max,
    w_max,
    upc_max,
    vpc_max,
    wpc_max,
    time,
    time_output,
    ud,
    step,
    machine_epsilon,
):
    """Calculate timestep based on convective CFL condition."""
    # Calculate directional timesteps
    dt_directions = [
        CFL * elem.dx / u_max,
        CFL * elem.dy / v_max,
        CFL * elem.dz / w_max,
    ]
    dt_cfl = min(dt_directions)

    # Apply ramping and step scaling for positive steps
    if step >= 0:
        ramp_factor = ud.dtfixed0 + min(step, 1.0) * (ud.dtfixed - ud.dtfixed0)
        dt = min(dt_cfl, ramp_factor)
        dt *= min(float(step + 1), 1.0)
    else:
        dt = dt_cfl

    # Ensure we don't overshoot the output time
    remaining_time = time_output - time
    if dt > remaining_time:
        dt = remaining_time

    # Calculate CFL numbers for output
    cfl = CFL * dt / dt_cfl
    cfl_ac = max(dt * upc_max / elem.dx, dt * vpc_max / elem.dy, dt * wpc_max / elem.dz)

    return dt, cfl, cfl_ac
