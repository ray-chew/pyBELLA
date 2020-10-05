import numpy as np
from scipy import signal

def get_B(data,g, ud):
    # equation 3
    return 0.5*(data.u**2 + data.v**2) + g * (data.rho * ud.h_ref)
    
def get_Pi(data,ud):
    # equation 4
    # node-to-cell averaging, (nx+1)*(ny+1) -> nx*ny
    kernel = np.array([[1.0,1.0],[1.0,1.0]])
    kernel /= kernel.sum()
    vorty = signal.convolve(data.vorty, kernel, mode='valid')
    # get relative vorticity
    vorty /= (data.rho * ud.h_ref)
#     vorty *= 86400.0 / ud.t_ref # vorticity in days^(-1)
    vorty /= ud.t_ref # vorticity in s^(-1)
    
    f = ud.coriolis_strength[0] / ud.t_ref
    
    return (vorty + f) / (data.rho * ud.h_ref)
    
def grad(arr,dd,direction):
    # get partial derivatives
    if direction == 'x':
        axs = 0
    elif direction == 'y':
        axs = 1
    else:
        assert(0, 'direction unspported')
        
    return np.gradient(arr,dd,axis=axs)
    
def get_DSI_SW(data, g, ud, elem):
    # equation 1
#     print(ud.h_ref, ud.t_ref, ud.u_ref)
    
    u = data.rhou / data.rho * ud.u_ref
    w = data.rhow / data.rho * ud.u_ref
    setattr(data,'u',u)
    setattr(data,'v',w)
    
    B = get_B(data,g, ud) # icx * icy
    Pi = get_Pi(data,ud) # icx * icy
    
    dx = np.diff(elem.x)[0]
    dy = np.diff(elem.z)[0]
    
    return 1./data.rho * (grad(B,dx,'x') * grad(Pi,dy,'y') - grad(B,dy,'y') * grad(Pi,dx,'x')) # icx * icy
