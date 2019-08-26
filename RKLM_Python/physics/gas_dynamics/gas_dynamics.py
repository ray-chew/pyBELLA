import numpy as np

machine_epsilon = np.finfo(float).eps

def dynamic_timestep(Sol, elem, ud, th, step):
    global machine_epsilon

    gamm = th.gamm

    Minv = 1.0 / np.sqrt(ud.Msq)
    CFL = ud.CFL

    u_max = machine_epsilon
    v_max = machine_epsilon
    w_max = machine_epsilon

    upc_max = machine_epsilon
    vpc_max = machine_epsilon
    wpc_max = machine_epsilon

    p = Sol.rhoY**gamm
    c = np.sqrt(gamm * p / Sol.rho) * Minv
    u = np.abs(Sol.rhou / Sol.rho)
    v = np.abs(Sol.rhov / Sol.rho)
    w = np.abs(Sol.rhow / Sol.rho)

    u_max = max(u.max(), u_max)
    v_max = max(v.max(), v_max)
    w_max = max(w.max(), w_max)

    upc_max = max((u+c).max(), upc_max)
    vpc_max = max((v+c).max(), vpc_max)
    wpc_max = max((w+c).max(), wpc_max)

    if (ud.acoustic_timestep == 1):
        dtx = CFL * elem.dx / upc_max
        dty = CFL * elem.dy / vpc_max
        dtz = CFL * elem.dz / wpc_max

        dt_cfl = min(min(dtx, dty), dtz)
        dt = min(dt_cfl, ud.dtfixed0 + min(step, 1) * (ud.dtfixed - ud.dtfixed0))

        return dt
    else:
        dtx = CFL * elem.dx / u_max
        dty = CFL * elem.dy / v_max
        dtz = CFL * elem.dz / w_max

        dt_cfl = min(min(dtx, dty), dtz)
        dt = min(dt_cfl, ud.dtfixed0 + min(step, 1) * (ud.dtfixed - ud.dtfixed0))
        dt *= min(float(step+1), 1.0)

        return dt