import numpy as np
from scipy import signal
import scipy.sparse.linalg as la

from dycore.physics.gas_dynamics.eos import nonhydrostasy, compressibility, is_compressible, is_nonhydrostatic

from termcolor import colored
from copy import deepcopy

class Blend(object):
    """
    Class that takes care of the blending interface.
    """

    def __init__(self,ud):
        self.bb = False
        self.cb = ud.continuous_blending
        self.psinc_init = ud.no_of_pi_initial
        # self.psinc_trans = ud.no_of_pi_transition
        self.hydro_init = ud.no_of_hy_initial

        if self.psinc_init > 0 and self.cb:
            self.bb = True

        self.c_init = self.criterion_init
        # self.c_trans = self.criterion_trans

        self.fac = ud.Msq
    
    def criterion_init(self, step):
        return step == (self.psinc_init) and self.cb and self.bb

    # def criterion_trans(self, step):
    #     return step <= self.psinc_trans and self.cb and not self.bb
        
    def convert_p2n(self, p2n):
        ndim = p2n.ndim
        dp2n = p2n - p2n.mean()

        self.kernel = np.ones([2]*ndim)
        dp2c = signal.fftconvolve(dp2n,self.kernel, mode='valid') / self.kernel.sum()

        # self.dp2n = dp2n - dp2n.mean()
        # self.dp2c = dp2c - dp2c.mean()

        self.dp2n = dp2n
        self.dp2c = dp2c

        return dp2c

    def update_Sol(self, Sol, elem, node, th, ud, mpv, sgn, label=None,writer=None):
        if writer != None: writer.populate(str(label)+'_before_blending', 'dp2n', self.dp2n)
        if writer != None: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_before_blending')

        if sgn == 'bef':
            sign = -1.0
        elif sgn == 'aft':
            sign = +1.0
        else:
            assert 0, "sgn == bef or sgn == aft"

        rho = np.copy(Sol.rho)
        rhoY = np.copy(Sol.rhoY)

        Y = rhoY / rho

        if ud.blending_mean == 'rhoY':
            rhoYc = (rhoY**th.gm1 + sign * self.fac * self.dp2c)**(th.gm1inv)
        elif ud.blending_mean == '1.0':
            rhoYc = (1.0 + sign * self.fac * self.dp2c)**(th.gm1inv)

        alpha = rhoYc / Sol.rhoY

        if ud.blending_conv == 'rho':
        ### keep theta, convert rho
            Sol.rho[...] = rho*alpha
            Sol.rhoY[...] = (Sol.rho * Y)

            rho_fac = Sol.rho / rho
            Sol.rhou[...] *= rho_fac
            Sol.rhov[...] *= rho_fac
            Sol.rhow[...] *= rho_fac
            Sol.rhoX[...] *= rho_fac

        elif ud.blending_conv == 'theta':
        ### keep rho, convert theta
            Yc = Y * alpha
            Sol.rhoY[...] = rho * Yc
            Sol.rhoX[...] = rho * (1.0/ Yc - mpv.HydroState.S0.reshape(1,-1))
        else:
            assert(0, "ud.blending_conv undefined.")

        if writer != None: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_blending')


    def update_p2n(self,Sol, mpv, node, th, ud):
        mpv.p2_nodes = self.dp2n



######################################################
# COMP - PSINC blending
######################################################

def do_comp_to_psinc_conv(Sol, mpv, bld, elem, node, th, ud, label, writer):
    print(colored("Converting COMP to PSINC",'blue'))
    dp2n = mpv.p2_nodes
    bld.convert_p2n(dp2n)
    bld.update_Sol(Sol,elem,node,th,ud,mpv,'bef',label=label,writer=writer)
    bld.update_p2n(Sol,mpv,node,th,ud)

    return Sol, mpv


def do_psinc_to_comp_conv(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt):
    from management import data
    print(colored("Blending... step = %i" %step,'blue'))
    Sol_freeze = deepcopy(Sol)
    mpv_freeze = deepcopy(mpv)

    ret = data.time_update(Sol,flux,mpv, t, t+dt, ud, elem, node, [0,step-1], th, bld=None, writer=None, debug=False)

    fac_old = ud.blending_weight
    fac_new = 1.0 - fac_old
    dp2n_0 = (fac_new * ret[2].p2_nodes_half + fac_old * mpv_freeze.p2_nodes_half)
    dp2n_1 = (fac_new * ret[2].p2_nodes + fac_old * mpv_freeze.p2_nodes)

    if ud.blending_type == 'half':
        dp2n = dp2n_0
    elif ud.blending_type == 'full':
        dp2n = dp2n_1
    else:
        assert 0, "incorrect ud.blending_type"

    if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_start', mpv_freeze.p2_nodes)
    if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_end', ret[2].p2_nodes)
    Sol = Sol_freeze
    mpv = mpv_freeze

    if writer != None: writer.populate(str(label)+'_after_full_step', 'dp2n', dp2n)
    print(colored("Converting PSINC to COMP",'blue'))
    bld.convert_p2n(dp2n)
    bld.update_Sol(Sol,elem,node,th,ud,mpv,'aft',label=label,writer=writer)
    bld.update_p2n(Sol,mpv,node,th,ud)

    return Sol, mpv 


######################################################
# SWE - Lake blending
######################################################

def do_swe_to_lake_conv(Sol, mpv, elem, node, ud, th, writer, label, debug):
    print(colored("swe to lake conversion...",'blue'))

    H1 = Sol.rho[:,2:-2:,][:,0,:]
    # setattr(ud,'mean_val',H1.mean())

    H10 = mpv.p2_nodes[:,2:-2,:].mean(axis=1)
    H10 -= H10.mean()

    # define 2D kernel
    kernel = np.ones((2,2))
    kernel /= kernel.sum()

    # do node-to-cell averaging
    H10 = signal.convolve(H10, kernel, mode='valid')

    # H1 = (H1 - ud.mean_val)
    H1 = H1 - ud.Msq * H10
    H1 = np.expand_dims(H1, axis=1)
    H1 = np.repeat(H1, elem.icy, axis=1)
    setattr(ud,'mean_val',H1)

    Sol.rhou[...] = Sol.rhou / Sol.rho * ud.mean_val
    Sol.rhov[...] = Sol.rhov / Sol.rho * ud.mean_val
    Sol.rhow[...] = Sol.rhow / Sol.rho * ud.mean_val
    Sol.rhoY[...] = Sol.rhoY / Sol.rho * ud.mean_val
    Sol.rho[...] = ud.mean_val

    # boundary.set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_swe_to_lake')


def do_lake_to_swe_conv(Sol, flux, mpv, elem, node, ud, th, writer, label, debug, step, window_step, t, dt):
    from management import data
    if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_lake_time_step')

    Sol_freeze = deepcopy(Sol)
    mpv_freeze = deepcopy(mpv)

    print(colored("doing lake-to-swe time-update...",'blue'))
    ret = data.time_update(Sol,flux,mpv, t, t+dt, ud, elem, node, [0,step], th, bld=None, writer=None, debug=False)

    fac_old = ud.blending_weight
    fac_new = 1.0 - fac_old

    dp2n_0 = (fac_new * ret[2].p2_nodes_half + fac_old * mpv_freeze.p2_nodes_half)
    dp2n_1 = (fac_new * ret[2].p2_nodes + fac_old * mpv_freeze.p2_nodes)

    if ud.blending_type == 'half':
        dp2n = dp2n_0
    elif ud.blending_type == 'full':
        dp2n = dp2n_1
    else:
        assert 0, "incorrect ud.blending_type"

    Sol = deepcopy(Sol_freeze)
    mpv = deepcopy(mpv_freeze)
    
    mpv.p2_nodes[...] = dp2n

    H10 = mpv.p2_nodes[:,2:-2,:].mean(axis=1)
    print(colored("lake to swe conversion...",'blue'))
    H10 -= H10.mean()

    # define 2D kernel
    kernel = np.ones((2,2))
    kernel /= kernel.sum()

    # do node-to-cell averaging
    H1 = signal.convolve(H10, kernel, mode='valid')
    # H1 = ud.mean_val + ud.Msq * H1
    # print(colored(H1.max(), 'red'))

    # project H1 back to horizontal slice with ghost cells
    H1 = np.expand_dims(H1, axis=1)
    H1 = np.repeat(H1, elem.icy, axis=1)
    H1 = ud.mean_val + ud.Msq * H1

    Sol.rho[...] = H1
    Sol.rhou[...] = Sol.rhou / ud.mean_val * Sol.rho
    Sol.rhov[...] = Sol.rhov / ud.mean_val * Sol.rho
    Sol.rhow[...] = Sol.rhow / ud.mean_val * Sol.rho
    Sol.rhoY[...] = Sol.rhoY / ud.mean_val * Sol.rho

    if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_lake_to_swe')
    return Sol, mpv


######################################################
# Nonhydrostatic - Hydrostatic blending
######################################################
def do_nonhydro_to_hydro_conv(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt):
    from management import data
    print(colored("nonhydrostatic to hydrostatic conversion...", 'blue'))
    # bld.convert_p2n(mpv.p2_nodes)
    # bld.update_Sol(Sol,elem,node,th,ud,mpv,'bef',label=label,writer=writer)
    # Sol.rhov = Sol.rhov_half

    # Sol_tmp = deepcopy(Sol)
    # flux_tmp = deepcopy(flux)
    # mpv_tmp = deepcopy(mpv)

    # nonhydro to hydro blending incomplete.
    # ret = data.time_update(Sol,flux,mpv, t, t+1*dt, ud, elem, node, [0,0], th, bld=None, writer=None, debug=False)
    
    # Sol = Sol_tmp
    # flux = flux_tmp
    # mpv = mpv_tmp
    # Sol = ret[0]
    # flux = ret[1]
    # mpv = ret[2]
    # Sol = deepcopy(ret[0])
    # mpv = deepcopy(ret[2])
    # Sol.rhov[...] = Sol.rhov_half
    # t += 0.5*dt
    # t += 1*dt
    return Sol, mpv, t
    

def do_hydro_to_nonhydro_conv(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt):
    from management import data
    print(colored("hydrostatic to nonhydrostatic conversion...", 'blue'))
    print(colored("Blending... step = %i" %step,'blue'))

    # Sol_tmp = deepcopy(Sol)
    # flux_tmp = deepcopy(flux)
    # mpv_tmp = deepcopy(mpv)

    # ret = data.time_update(Sol,flux,mpv, t, t+dt, ud, elem, node, [0,step-1], th, bld=None, writer=None, debug=False)

    # Sol = Sol_tmp
    # flux = flux_tmp
    # mpv = mpv_tmp

    # retv_half = ret[0].rhov_half / ret[0].rho_half
    # retv_full = ret[0].rhov / ret[0].rho

    # solv_half = Sol.rhov_half / Sol.rho_half
    # solv_full = Sol.rhov / Sol.rho

    # fac_full = 0.5
    # fac_half = 1.0 - fac_full

    # # print(np.sum(solv_full))
    # # print(np.sum(retv_half))
    # # print(np.sum((fac_full * solv_full + fac_half * retv_half)))
    # # print(np.sum(fac_half * retv_half))

    # fac_full = 0.5
    # fac_half = 0.5

    # # print(np.sum(solv_full))
    # # print(np.sum(retv_half))
    # # print(np.sum((fac_full * solv_full + fac_half * retv_half)))
    # # print(np.sum(fac_half * retv_half))

    # if writer != None: writer.populate(str(label)+'_after_full_step', 'ret_half', ret[0].rhov_half)
    # if writer != None: writer.populate(str(label)+'_after_full_step', 'ret_full', ret[0].rhov)


    # if writer != None: writer.populate(str(label)+'_after_full_step', 'solv_half', Sol.rhov_half)
    # if writer != None: writer.populate(str(label)+'_after_full_step', 'solv_full', Sol.rhov)

    # Sol.rhov = Sol.rho * (fac_full * solv_full + fac_half * retv_half)
    # if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_end', ret[2].p2_nodes)

    # fac_mpv_half = 0.5
    # fac_mpv_full = 1.0 - fac_mpv_half
    # mpv.p2_nodes = fac_mpv_half * mpv.p2_nodes + fac_mpv_full * ret[2].p2_nodes
    # dp2n = ret[2].p2_nodes_half
    # bld.convert_p2n(dp2n)
    # bld.update_Sol(Sol,elem,node,th,ud,mpv,'aft',label=label,writer=writer)
    # bld.update_p2n(Sol,mpv,node,th,ud)
# 
    return Sol, mpv


######################################################
# Blending calls from data.py
######################################################
def blending_before_timestep(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt, swe_to_lake, debug):
    ######################################################
    # Blending : Do full regime to limit regime conversion
    ######################################################
    # these make sure that we are the correct window step
    if bld is not None and window_step == 0: 
        # these make sure that blending switches are on
        if (bld.bb or bld.cb) and ud.blending_conv is not None:
            # these distinguish between SWE and Euler blending
            if ud.blending_conv == 'swe':
                do_swe_to_lake_conv(Sol, mpv, elem, node, ud, th, writer, label, debug)
                swe_to_lake = True
            else:
                Sol, mpv = do_comp_to_psinc_conv(Sol, mpv, bld, elem, node, th, ud, label, writer)

    ######################################################
    # Blending : Do full steps or transition steps?
    ######################################################
    if bld is not None:
        c_init = bld.criterion_init(window_step)
    else:
        c_init = False

    ######################################################
    # Blending : If full blending steps...
    ######################################################
    # check that blending switches are on
    if c_init and bld.cb and ud.blending_conv is not None:
        # distinguish between Euler and SWE blending
        if ud.blending_conv is not 'swe':
            Sol, mpv = do_psinc_to_comp_conv(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt)
    
    ######################################################
    # Initial Blending
    ######################################################
    # Is initial blending switch on, and if yes, are we in the 0th time-step?
    if ud.initial_blending == True and step < 1 and bld is not None:
        # Distinguish between SWE and Euler blendings
        if ud.blending_conv is not 'swe':
            if bld.psinc_init > 0:
                ud.is_compressible = 0
                ud.compressibility = 0.0
                Sol, mpv = do_comp_to_psinc_conv(Sol, mpv, bld, elem, node, th, ud, label, writer)
            elif bld.hydro_init > 0:
                Sol, mpv, t = do_nonhydro_to_hydro_conv(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt)
                ud.is_nonhydrostatic = 0
                ud.nonhydrostasy = 0.0
        else:
            do_swe_to_lake_conv(Sol, mpv, elem, node, ud, th, writer, label, debug)
            swe_to_lake = True
            ud.is_compressible = 0
            ud.compressibility = 0.0

    # Elif, is initial blending switch on and are we on the 1st time-step?
    # If we are on the first time-step, do we do comp-psinc blending?
    elif ud.initial_blending == True and step == ud.no_of_pi_initial and bld is not None:
        # Distinguish between SWE and Euler blendings
        if ud.blending_conv is not 'swe':
            Sol, mpv = do_psinc_to_comp_conv(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt)
            ud.is_compressible = 1
            ud.compressibility = 1.0
    # Else, do we do nonhydrostatic-hydrostatic blending?
    elif ud.initial_blending == True and step == ud.no_of_hy_initial and bld is not None:
        if ud.blending_conv is not 'swe':
            Sol, mpv = do_hydro_to_nonhydro_conv(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt)
            ud.is_nonhydrostatic = 1
            ud.nonhydrostasy = 1.0
    else:
        ud.is_compressible = is_compressible(ud,window_step)
        ud.compressibility = compressibility(ud,t,window_step)
        ud.is_nonhydrostatic = is_nonhydrostatic(ud,window_step)
        ud.nonhydrostasy = nonhydrostasy(ud,t,window_step)

    return swe_to_lake, Sol, mpv, t


def blending_after_timestep(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt, swe_to_lake, debug):
    ######################################################
    # Blending : Are we in the lake regime? And is this
    #            the window step where we go back to SWE?
    ######################################################
    if bld is not None and swe_to_lake and step > 0:
        initialise_lake_to_swe_conv = bld.criterion_init(window_step+1)

    ######################################################
    # Blending : If we are in the lake regime, is blending
    #            on? If yes, do lake-to-swe conversion.
    ######################################################
    if ud.blending_conv == 'swe' and swe_to_lake and initialise_lake_to_swe_conv and bld is not None:
        tmp_CFL = np.copy(ud.CFL)
        ud.CFL = 0.8
        Sol, mpv = do_lake_to_swe_conv(Sol, flux, mpv, elem, node, ud, th, writer, label, debug, step, window_step, t, dt)
        ud.CFL = tmp_CFL
        ud.is_compressible = 1
        ud.compressibility = 1.0

    return Sol, mpv