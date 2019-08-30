import numpy as np

machine_epsilon = np.finfo(float).eps

def dynamic_timestep(Sol, time, time_output, elem, ud, th, step):
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
        dt = min(dt_cfl, ud.dtfixed0 + min(step, 1.) * (ud.dtfixed - ud.dtfixed0))

        if (2.0*dt >= time_output - time):
            dt = 0.5 * (time_output - time) + machine_epsilon


        return dt
    else:
        dtx = CFL * elem.dx / u_max
        dty = CFL * elem.dy / v_max
        dtz = CFL * elem.dz / w_max

        dt_cfl = min(min(dtx, dty), dtz)
        dt = min(dt_cfl, ud.dtfixed0 + min(step, 1.) * (ud.dtfixed - ud.dtfixed0))
        dt *= min(float(step+1), 1.0)

        if ((2.0*dt) >= (time_output - time)):
            dt = 0.5 * (time_output - time) + machine_epsilon

        # if step == 38:
        #     dt = 0.0000668354011612
        # if step == 39:
        #     dt = 0.0000668354011612

        # if step == 1:
        #     dt = 71.5337121335972483
        # if step == 2:
        #     dt = 71.3741788666122119
            
        return dt