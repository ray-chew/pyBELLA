import numpy as np

from management.variable import States, Characters
from management.enumerator import LimiterType

def recovery(lefts, rights, Sol, flux, lmbda, ud, th, elem):
    gamm = th.gamm
    
    order_two = 1 # always 1

    Sol.primitives()

    flux_rhoY = flux.rhoY
    Sol_rhoY = Sol.rhoY[1:-1,1:-2]
    lefts_idx = (slice(None),slice(0,-1))
    rights_idx = (slice(None),slice(1,None))

    u = 0.5 * (flux_rhoY[lefts_idx] + flux_rhoY[rights_idx]) / Sol_rhoY

    diffs = States(elem.sc,ud)
    diffs.u[:,:-1] = Sol.u[rights_idx] - Sol.u[lefts_idx]
    diffs.v[:,:-1] = Sol.v[rights_idx] - Sol.v[lefts_idx]
    diffs.w[:,:-1] = Sol.w[rights_idx] - Sol.w[lefts_idx]
    diffs.Y[:,:-1] = 1.0 / Sol.Y[rights_idx] - 1.0 / Sol.Y[lefts_idx]

    slopes(Sol, diffs, ud, elem)

def slopes(Sol, diffs, ud, elem):
    limiter_type_velocity = ud.limiter_type_velocity
    limiter_type_scalar = ud.limiter_type_scalars

    lefts_idx = (slice(None),slice(0,-1))
    rights_idx = (slice(None),slice(1,None))

    # what are these?
    kp = ud.kp
    kz = ud.kz
    kY = ud.kY

    # amplitudes of the left state difference:
    aul = diffs.u[lefts_idx]
    avl = diffs.v[lefts_idx]
    awl = diffs.w[lefts_idx]
    aYl = diffs.Y[lefts_idx]

    aur = diffs.u[rights_idx]
    avr = diffs.v[rights_idx]
    awr = diffs.w[rights_idx]
    aYr = diffs.Y[rights_idx]

    Slopes = Characters(elem.sc)
    Slopes.u = limiters(limiter_type_velocity, aul, aur, kp)
    Slopes.v = limiters(limiter_type_velocity, avl, avr, kz)
    Slopes.w = limiters(limiter_type_velocity, awl, awr, kz)
    Slopes.Y = limiters(limiter_type_scalar, aYl, aYr, kY)

    return Slopes

def limiters(limiter_type, al, ar, kp):
    # write switch for limiter types
    # for now, just use LimiterType == None
    if limiter_type == LimiterType.NONE:
        return 0.5 * (al + ar)

    

