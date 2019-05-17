from inputs.enum_bdry import BdryType

from scipy import signal
import numpy as np

def euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, dt, alpha_diff):
    nc = node.nc
    p2 = mpv.dp2_nodes

    x_periodic = ud.bdry_type[0] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[1] == BdryType.PERIODIC
    z_periodic = ud.bdry_type[2] == BdryType.PERIODIC

    print(x_periodic)

    

def operator_coefficients_nodes(elem, node, Sol, mpv, ud, th, dt):
    g = ud.gravity_strength[1]
    Msq = ud.Msq
    Gammainv = th.Gammainv

    ndim = node.ndim
    nonhydr = ud.nonhydrostasy

    time_offset = 3.0 - ud.acoustic_order

    coriolis = ud.coriolis_strength[0]

    print(ud.bdry_type[0] == BdryType.PERIODIC)