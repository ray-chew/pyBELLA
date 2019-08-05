import numpy as np

def recovery(lefts, rights, Sol, flux, lmbda, ud, th, elem):
    gamm = th.gamm
    
    order_two = 1 # always 1

    Sol.primitives()
    print(Sol.u)