import numpy as np
from scipy import signal

def recompute_advective_fluxes(flux, Sol, elem, dt):
    ################################################
    # 2D case for now - generalise in future
    ################################################

    left_idx = (slice(None),slice(0,-2))
    right_idx = (slice(None),slice(1,-1))

    # left_idx = (slice(None),slice(1,-1))
    # right_idx = (slice(None),slice(2,None))

    top_idx = (slice(0,-2),)
    bottom_idx = (slice(1,-1),)

    umid_idx = (slice(1,-1),)
    vmid_idx = (slice(None),slice(1,-1))
    
    BB = np.arange(25).reshape(5,5)
    print(BB)
    # print(BB[left_idx])
    # print(BB[right_idx])

    # print(BB[top_idx])
    # print(BB[bottom_idx])
    
    rhoYu = Sol.rhoY * Sol.rhou / Sol.rho
    rhoYu_interface = 0.5 * (rhoYu[left_idx] + rhoYu[right_idx])
    print(rhoYu_interface.shape)
    kernel_u = np.array([[1.],[2.],[1.]])

    print(flux[0].rhoY.shape)
    BBp = signal.convolve2d(BB,kernel_u,mode='valid')
    # print(BBp)
    print(signal.convolve2d(rhoYu_interface,kernel_u, mode='valid').shape)
    flux[0].rhoY[1:-2,1:-1] = signal.convolve2d(rhoYu_interface,kernel_u, mode='valid') / kernel_u.sum()

    rhoYv = Sol.rhoY * Sol.rhov / Sol.rho
    rhoYv_interface = 0.5 * (rhoYv[top_idx] + rhoYv[bottom_idx])
    print(rhoYv_interface.shape)

    kernel_v = np.array([[1.,2.,1.]])
    flux[1].rhoY[1:-1,1:-2] = signal.convolve2d(rhoYv_interface,kernel_v, mode='valid') / kernel_v.sum()

    # print(flux[0].rhoY)