import numpy as np

from management.variable import States, Characters
from management.enumerator import LimiterType
from physics.gas_dynamics.eos import rhoe

from debug import find_nearest

truefalse = True

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

    shape = Sol.u.shape
    Diffs = States(shape,ud)
    Ampls = Characters(shape)
    Lefts = States(shape, ud)
    Rights = States(shape, ud)

    Diffs.u[:,:-1] = Sol.u[rights_idx] - Sol.u[lefts_idx]
    Diffs.v[:,:-1] = Sol.v[rights_idx] - Sol.v[lefts_idx]
    Diffs.w[:,:-1] = Sol.w[rights_idx] - Sol.w[lefts_idx]
    Diffs.Y[:,:-1] = 1.0 / Sol.Y[rights_idx] - 1.0 / Sol.Y[lefts_idx]

    # global truefalse
    # if truefalse == True:
    #     # print(u[1])

    #     # print(Sol.rhou[-3])
    #     # print(flux.rhoY[1])

    #     idx = 130
    #     print("Sol.u[idx] = ", Sol.u[idx])
    #     print("Sol.rhou[idx] = ", Sol.rhou[idx])
    #     print("Sol.rhov[idx] = ", Sol.rhov[idx])
    #     print("Diffs.u[idx] = ", Diffs.u[idx])
    #     print("Diffs.v[idx] = ", Diffs.v[idx])

    #     # val, idx = find_nearest(Sol.rhou,0.50000030887122227)
    #     # print("val = ", val)
    #     # print("idx = ", idx)
    #     truefalse = False

    Slopes = slopes(Sol, Diffs, ud, elem)

    Ampls.u[...] = 0.5 * Slopes.u * (1. - lmbda * u)
    Ampls.v[...] = 0.5 * Slopes.v * (1. - lmbda * u)
    Ampls.w[...] = 0.5 * Slopes.w * (1. - lmbda * u)
    Ampls.Y[...] = 0.5 * Slopes.Y * (1. - lmbda * u)
    
    Lefts.u[...] = Sol.u + order_two * Ampls.u
    Lefts.v[...] = Sol.v + order_two * Ampls.v
    Lefts.w[...] = Sol.w + order_two * Ampls.w
    Lefts.Y[...] = 1.0 / (1.0 / Sol.Y + order_two * Ampls.Y)
    
    # Ampls.change_dir()

    # global truefalse
    # if truefalse == True:
    #     # print(u[1])

    #     # print(Sol.rhou[-3])
    #     # print(flux.rhoY[1])

    #     idx = 265
    #     # print("Lefts.u[idx] = ", Lefts.u.flatten()[idx])

    #     print("u[idx] = ", u.flatten()[idx])
    #     # print("Slopes.v[idx] = ", Slopes.v.flatten()[idx])
    #     print("Sol.rhou[idx] = ", Sol.u.flatten()[0])
    #     print("Sol.u[idx] = ", Sol.u.flatten()[idx])
    #     print("Sol.v[idx] = ", Sol.v.flatten()[idx])

    #     print("Ampls.u[idx] = ", Ampls.u.flatten()[idx])
    #     print("Ampls.v[idx] = ", Ampls.v.flatten()[idx])

    #     print("Slopes.u[idx] = ", Slopes.u.flatten()[idx])
    #     print("Slopes.v[idx] = ", Slopes.v.flatten()[idx])
    #     print("Diffs.u[idx] = ", Diffs.u.flatten()[idx])
    #     print("Diffs.v[idx] = ", Diffs.v.flatten()[idx])
        
    #     # print(Slopes.u.shape)

    #     val, idx = find_nearest(Sol.u,-0.66152648868480246)
    #     # print("val = ", val)
    #     # print("idx = ", idx)
        
    #     truefalse = False

    Ampls.u[...] = -0.5 * Slopes.u * (1. + lmbda * u)
    Ampls.v[...] = -0.5 * Slopes.v * (1. + lmbda * u)
    Ampls.w[...] = -0.5 * Slopes.w * (1. + lmbda * u)
    Ampls.Y[...] = -0.5 * Slopes.Y * (1. + lmbda * u)

    Rights.u[...] = Sol.u + order_two * Ampls.u
    Rights.v[...] = Sol.v + order_two * Ampls.v
    Rights.w[...] = Sol.w + order_two * Ampls.w
    Rights.Y[...] = 1.0 / (1.0 / Sol.Y + order_two * Ampls.Y)

    Lefts.rhoY[lefts_idx] = Rights.rhoY[rights_idx] = 0.5 * (Sol.rhoY[lefts_idx] + Sol.rhoY[rights_idx]) \
        - order_two * 0.5 * lmbda * (Sol.u[rights_idx] * Sol.rhoY[rights_idx] - Sol.u[lefts_idx] * Sol.rhoY[lefts_idx]) 
    Lefts.p0[lefts_idx] = Rights.p0[rights_idx] = Lefts.rhoY[lefts_idx]**gamm

    get_conservatives(Rights, ud, th)
    get_conservatives(Lefts, ud, th)

    # print(Rights.rhoe[0])

    return Lefts, Rights

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

    Slopes = Characters(Diffs.u.shape)

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
    # sgn = np.sign(U.rhoY)
    p = U.rhoY**th.gamminv
    
    U.rhoe = rhoe(U.rho, U.u, U.v, U.w, p, ud, th)