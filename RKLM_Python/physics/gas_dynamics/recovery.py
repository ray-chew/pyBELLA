import numpy as np

from management.variable import States, Characters
from management.enumerator import LimiterType

truefalse = True

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array.flat[idx], idx

def recovery(Sol, flux, lmbda, ud, th, elem):
    gamm = th.gamm
    
    order_two = 1 # always 1

    Sol.primitives(th)

    lefts_idx = (slice(None),slice(0,-1))
    rights_idx = (slice(None),slice(1,None))
    inner_idx = (slice(1,-1), slice(1,-1))

    # inner_idx here are where the interface fluxes are calculated with non-zero values.
    face_inner_idx = (slice(1,-1), slice(1,-1))
    u = np.zeros_like(Sol.rhoY)
    u[inner_idx] = 0.5 * (flux.rhoY[face_inner_idx][lefts_idx] + flux.rhoY[face_inner_idx][rights_idx]) / Sol.rhoY[inner_idx]

    Diffs = States(elem.sc,ud)
    Ampls = Characters(elem.sc)
    Lefts = States(elem.sc, ud)
    Rights = States(elem.sc, ud)

    Diffs.u[:,:-1] = Sol.u[rights_idx] - Sol.u[lefts_idx]
    Diffs.v[:,:-1] = Sol.v[rights_idx] - Sol.v[lefts_idx]
    Diffs.w[:,:-1] = Sol.w[rights_idx] - Sol.w[lefts_idx]
    Diffs.Y[:,:-1] = 1.0 / Sol.Y[rights_idx] - 1.0 / Sol.Y[lefts_idx]

    global truefalse
    if truefalse == True:
        # print(u[1])

        # print(Sol.rhou[-3])
        # print(flux.rhoY[1])

        idx = 0
        # print(Sol.u[idx])
        print(Sol.rhou[idx])
        print(Sol.rhov[idx])
        # print(Diffs.u[idx])
        # print(Diffs.v[idx])

        # val, idx = find_nearest(Sol.rhou,0.50000030887122227)
        # print("val = ", val)
        # print("idx = ", idx)
        truefalse = False

    Slopes = slopes(Sol, Diffs, ud, elem)

    Ampls.u[...] = 0.5 * Slopes.u * (1. - lmbda * u)
    Ampls.v[...] = 0.5 * Slopes.v * (1. - lmbda * u)
    Ampls.w[...] = 0.5 * Slopes.w * (1. - lmbda * u)
    Ampls.Y[...] = 0.5 * Slopes.Y * (1. - lmbda * u)
    
    Lefts.u[...] = Sol.u + order_two * Ampls.u
    Lefts.v[...] = Sol.v + order_two * Ampls.v
    Lefts.w[...] = Sol.w + order_two * Ampls.w
    Lefts.Y[...] = 1.0 / (1.0 / Sol.Y + order_two * Ampls.Y)
    
    Ampls.change_dir()

    Rights.u[...] = Sol.u + order_two * Ampls.u
    Rights.v[...] = Sol.v + order_two * Ampls.v
    Rights.w[...] = Sol.w + order_two * Ampls.w
    Rights.Y[...] = 1.0 / (1.0 / Sol.Y + order_two * Ampls.Y)

    Lefts.rhoY[lefts_idx] = Rights.rhoY[rights_idx] = 0.5 * (Sol.rhoY[lefts_idx] + Sol.rhoY[rights_idx]) \
        - order_two * 0.5 * lmbda * (Sol.u[rights_idx] * Sol.rhoY[rights_idx] - Sol.u[lefts_idx] * Sol.rhoY[lefts_idx]) 
    Lefts.p0[lefts_idx] = Rights.p0[rights_idx] = Lefts.rhoY[lefts_idx]**gamm

    get_conservatives(Rights, ud, th)
    get_conservatives(Lefts, ud, th)

    return Lefts, Rights

    # print(Slopes.u.shape)

def slopes(Sol, Diffs, ud, elem):
    limiter_type_velocity = ud.limiter_type_velocity
    limiter_type_scalar = ud.limiter_type_scalars

    lefts_idx = (slice(None),slice(0,-1))
    rights_idx = (slice(None),slice(1,None))

    # what are these?
    kp = ud.kp
    kz = ud.kz
    kY = ud.kY

    # amplitudes of the left state difference:
    aul = Diffs.u[lefts_idx][lefts_idx]
    avl = Diffs.v[lefts_idx][lefts_idx]
    awl = Diffs.w[lefts_idx][lefts_idx]
    aYl = Diffs.Y[lefts_idx][lefts_idx]

    aur = Diffs.u[lefts_idx][rights_idx]
    avr = Diffs.v[lefts_idx][rights_idx]
    awr = Diffs.w[lefts_idx][rights_idx]
    aYr = Diffs.Y[lefts_idx][rights_idx]

    Slopes = Characters(elem.sc)

    Slopes.u[:,1:-1] = limiters(limiter_type_velocity, aul, aur, kp)
    Slopes.v[:,1:-1] = limiters(limiter_type_velocity, avl, avr, kz)
    Slopes.w[:,1:-1] = limiters(limiter_type_velocity, awl, awr, kz)
    Slopes.Y[:,1:-1] = limiters(limiter_type_scalar, aYl, aYr, kY)

    return Slopes

def limiters(limiter_type, al, ar, kp):
    # write switch for limiter types
    # for now, just use LimiterType == None
    if limiter_type == LimiterType.NONE:
        return 0.5 * (al + ar)

def get_conservatives(U, ud, th):
    U.rho = U.rhoY / U.Y
    U.rhou = U.u * U.rho
    U.rhov = U.v * U.rho
    U.rhow = U.w * U.rho
    U.rhoY = U.Y * U.rho

    p = U.rhoY**th.gamminv
    U.rhoe = ud.rhoe(U.rho, U.u, U.v, U.w, p, ud, th)