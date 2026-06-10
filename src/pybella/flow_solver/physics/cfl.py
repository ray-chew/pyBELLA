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

    # terrain: the vertical coordinate velocity is eta_dot = (w - G.u_h)/J
    # and the vertical signal speed gains the slope/Jacobian factor; the
    # metric is in the unflipped orientation here (called between steps)
    c_vert = c
    if elem.metric is not None:
        m = elem.metric
        moms = (Sol.rhou, Sol.rhov, Sol.rhow)
        contra = moms[m.vaxis] - m.G1 * moms[m.haxes[0]]
        slope_sq = m.G1**2
        if m.G2 is not None:
            contra = contra - m.G2 * moms[m.haxes[1]]
            slope_sq = slope_sq + m.G2**2
        vels = [u, v, w]
        vels[m.vaxis] = np.abs(contra / Sol.rho) * m.ooJ
        u, v, w = vels
        c_vert = c * np.sqrt(1.0 + slope_sq) * m.ooJ

    # Find maximum velocities (with minimum threshold)
    u_max = max(u.max(), machine_epsilon)
    v_max = max(v.max(), machine_epsilon)
    w_max = max(w.max(), machine_epsilon)

    # Calculate acoustic velocities
    cs = [c, c, c]
    if elem.metric is not None:
        cs[elem.metric.vaxis] = c_vert
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
