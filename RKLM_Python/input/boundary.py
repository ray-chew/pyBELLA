from management.enumerator import HillShapes
from input.enum_bdry import BdryType

import numpy as np 

class InitializeBdry(object):
    def __init__(self, elem, ud):
        icx = elem.icx
        icz = elem.icz
        igx = elem.igx
        # igz = elem.igz

        if elem.ndim == 1:
            print("this boundary condition setting makes no sense in 1D")

        elif elem.ndim == 2:
            self.wall_slope = slanted_wall_slope(elem.x[igx:int(icx-igx)], ud)
            slope_sum = np.sum(np.abs(self.wall_slope))

            if slope_sum <= np.sqrt(np.finfo(np.float).eps):
                self.wall_relative_slope = 0.0
            else:
                slope_sum_inv = 1.0 / slope_sum
                self.wall_relative_slope = slope_sum_inv * np.abs(self.wall_slope)

                slope_sum = np.sum(self.wall_relative_slope)

                print("relative slope sum = %e" %slope_sum)

        elif elem.ndim == 3:
            self.wall_slope = slanted_wall_slope_3d(elem.x[igx:int(icx-igx)], icx * icz, ud)
            slope_sum = np.sum(np.abs(self.wall_slope))

            if slope_sum <= np.sqrt(np.finfo(np.float).eps):
                self.wall_relative_slope = 0.0
            else:
                slope_sum_inv = 1.0 / slope_sum
                self.wall_relative_slope = slope_sum_inv * np.abs(self.wall_slope)

                slope_sum = np.sum(self.wall_relative_slope)
                print("relative slope sum = %e" %slope_sum)

            
def slanted_wall_slope_3d(x, size, ud):
    hill_height = ud.hill_height
    length_scale_inv = 1.0 / ud.hill_length_scale
    x_sc = x * length_scale_inv

    return np.repeat((-hill_height * 2.0 * x_sc / ((1.0 + x_sc * x_sc)**2) * length_scale_inv), size)
    
def slanted_wall_slope(x, ud):
    if (ud.hill_shape == HillShapes.SCHULTOW):
        kx = 2.0 * 2.0 * np.pi / (ud.xmax - ud.xmin)
        kz = np.sqrt(ud.Nsq / ud.wind_speed / ud.wind_speed + kx * kx)
        q = 0.25

        xi = kx * x
        y = 0.0

        for __ in range(10):
            y = q * np.cos(xi + y)

        yp = -q * np.sin(xi + y) / (1. + q * np.sin(xi + y))
        yp = kx * yp / kz
        return yp

    else:
        hill_height = ud.hill_height
        length_scale_inv = 1.0 / ud.hill_length_scale
        x_sc = x * length_scale_inv
        return (-hill_height * 2.0 * x_sc / ((1.0 + x_sc * x_sc)**2 * length_scale_inv)) * np.ones((x.shape))

def set_wall_rhoYflux(bdry, Sol0, mpv, elem, ud):
    icx = elem.icx
    icy = elem.icy
    icz = elem.icz
    igx = elem.igx
    igy = elem.igy
    igz = elem.igz

    if elem.ndim == 1:
        print("wall flux in 1D makes no sense")
    elif elem.ndim == 2:
        if (ud.bdrytype_min[0] == BdryType.PERIODIC):
            is_x_periodic = 1
        else:
            is_x_periodic = 0

        nstart = elem.igy * elem.icx
        idx = np.arange(igx,icx-igx)
        n = nstart + idx
        bdry.wall_rhoYflux = wall_rhoYflux(elem.x[igx:icx-igy],elem.y[igx:icx-igy],np.divide(Sol0.rhou[n],Sol0.rho[n]),np.divide(Sol0.rhov[n],Sol0.rho[n]),mpv.HydroState_n.rhoY0[igy],ud)
        wall_flux_balance = np.sum(bdry.wall_rhoYflux[igx:icx-igx])

        # correction for zero net flux
        bdry.wall_rhoYflux[igx:icx-igx] -= wall_flux_balance * bdry.wall_relative_slope

        bdry.wall_rhoYflux[:igx] = bdry.wall_rhoYflux[icx-1-igx:icx-1][::-1] = 0.0

        if (is_x_periodic):
            bdry.wall_rhoYflux[:igx] = bdry.wall_rhoYflux[icx-2*igx:icx-1*igx]
            bdry.wall_rhoYflux[icx-igx:icx] = bdry.wall_rhoYflux[igx:2*igx]

        # flux_sum = np.sum(bdry.wall_rhoYflux[igx:icx-igx])
    elif elem.ndim == 3:
        if (ud.bdrytype_min[0] == BdryType.PERIODIC):
            is_x_periodic = 1
        else:
            is_x_periodic = 0

        if (ud.bdrytype_min[2] == BdryType.PERIODIC):
            is_z_periodic = 1
        else:
            is_z_periodic = 0

        nstart = igy*icx

        k_idx = np.arange(igz,icz-igy)
        i_idx = np.arange(igx,icx-igx)

        #########################
        # shit-code-alert: rewrite code on a clearer-head day.
        #
        nijk = []
        nik = []

        for k in k_idx:
            njk = nstart + k * icx * icy
            nk = k * icx
            for i in i_idx:
                nijk.append(njk + i)
                nik.append(nk + i)
        #
        #########################

        bdry.wall_rhoYflux[nik] = wall_rhoYflux(elem.x[i_idx],elem.z[k_idx], np.divide(Sol0.rhou[nijk],Sol0.rho[nijk]), np.divide(Sol0.rhow[nijk],Sol0.rho[nijk]), mpv.HydroState_n.rhoY0[igy],ud)
        wall_flux_balance = np.sum(bdry.wall_rhoYflux)

        #########################
        #
        nik_left = []
        nik_right = []
        for k in range(icz):
            nk = k * icx
            for i in range(igx):
                nik_left.append(nk+i)
                nik.right = nk + icx-1-i
        bdry.wall_rhoYflux[nik_left] = 0.0
        bdry.wall_rhoYflux[nik_right] = 0.0
        #
        #########################

        #########################
        # correction for zero net flux
        nik = []
        for k in range(igz,icz-igx):
            nk = k*icx
            for i in range(igx,icx-igx):
                nik.append(nk+i)

        bdry.wall_rhoYflux[nik] -= wall_flux_balance * bdry.wall_relative_slope[nik]
        #
        #########################

        #########################
        #
        if (is_x_periodic):
            niklt = []
            nikrs = []
            nikrt = []
            nikls = []
            for k in range(icz):
                nk = k*icx
                for i in range(igx):
                    niklt.append(nk + i)
                    nikrs.append(nk + icx - 2*igx + i)
                    nikrt.append(nk + icx - igx + i)
                    nikls.append(nk + igx + i)
            bdry.wall_rhoYflux[niklt] = bdry.wall_rhoYflux[nikrs]
            bdry.wall_rhoYflux[nikrt] = bdry.wall_rhoYflux[nikls]

        if (is_z_periodic):
            niklt = []
            nikrs = []
            nikrt = []
            nikls = []
            for i in range(icx):
                ni = i
                for k in range(igz):
                    niklt.append(ni + k*icz)
                    nikrs.append(ni + (icz - 2*igz + k)*icx)
                    nikrt.append(ni + (icz - igz + k) * icx)
                    nikls.append(ni + (igz + k) * icx)
            bdry.wall_rhoYflux[niklt] = bdry.wall_rhoYflux[nikrs]
            bdry.wall_rhoYflux[nikrt] = bdry.wall_rhoYflux[nikls]
        #
        #########################

        #########################
        #
        nik = []
        for k in range(igz, icz-igx):
            nk = k * icx
            for i in range(igx,icx-igx):
                nik.append(nk + i)
        flux_sum = np.sum(bdry.wall_rhoYflux[nik])

        print("wall flux sum = %e" %flux_sum)


def wall_rhoYflux(x,z,wind_speed_x,wind_speed_z,rhoY0,ud):
    return slanted_wall_slope(x,ud) * wind_speed_x * rhoY0

def set_explicit_boundary_data(Sol, elem):
    for split_step in range(elem.ndim):
        None