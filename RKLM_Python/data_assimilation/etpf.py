import numpy as np
import pyemd

def da_interface(results, obs, obs_attributes, times, tout, N, loc=0):
    ig = 2
    inner = (slice(ig,-ig),slice(ig,-ig))


    # local_ens = np.array([getattr(results[:,loc,...][n],attr) for n in range(N)])
    # local_ens = np.array([mem[inner] for mem in local_ens])
    # local_ens = analysis(local_ens,attr)
    # print(results[:,0,...][0].rho)

    attributes = ['rho','rhou','rhov','rhow','rhoY','rhoX']
    tmp = np.array([getattr(results[:,loc,...][n],obs_attributes[0])[inner] for n in range(N)])
    tmp = tmp[:,np.newaxis,...]
    ensemble = np.array([getattr(results[:,loc,...][n],attributes[0])[inner] for n in range(N)])
    ensemble = ensemble[:,np.newaxis,...]
    obs_current = np.array(obs[np.where(np.isclose(times,tout))[0][0]][obs_attributes[0]])[inner]
    
    # obs_current = bin_func(obs_current,(Nx,Ny))
    obs_current = obs_current[np.newaxis,...]
    
    for attr in obs_attributes[1:]:
        tmp = np.hstack((tmp,np.array([getattr(results[:,loc,...][n],attr)[inner] for n in range(N)])[:,np.newaxis,...]))
        tmp01 = np.array(obs[np.where(np.isclose(times,tout))[0][0]][attr])[inner]
        # tmp01 = bin_func(tmp01,(Nx,Ny))
        # print(tmp01.shape)
        tmp01 = tmp01[np.newaxis,...]
        obs_current = np.vstack((obs_current,tmp01))

    for attr in attributes[1:]:
        ensemble = np.hstack((ensemble,np.array([getattr(results[:,loc,...][n],attr)[inner] for n in range(N)])[:,np.newaxis,...]))
        

    obs_current = obs_current[np.newaxis,...]
    # r = tmp - obs_current

    Hx = tmp.reshape(N,-1)
    obs_current = obs_current.reshape(-1)
    # ensemble = ensemble.reshape(N,-1)
    # print(ensemble.shape)
    
    etpf = analysis(ensemble)
    etpf.analyse(obs_current,0.01,Hx,N)

    analysis_ens = etpf.get_ensemble_from_X()

    
    for n in range(N):
        cnt = 0
        current = analysis_ens[n]
        for attr in attributes:
            data = current[cnt]
            data = np.pad(data,2,mode='wrap')
            
            setattr(results[:,loc,...][n],attr,data)
            cnt += 1
    return results

class analysis(object):
    def __init__(self,ensemble,identifier=None):
        self.ensemble = np.array(ensemble)
        self.ensemble_shape = self.ensemble.shape

        self.X = self.state_vector(ensemble)
        self.identifier = identifier

        # rejuvenation factor
        self.delta = 0.1

    def analyse(self,obs_current,obs_covar,Hx,N):
        print("starting ETPF analysis...")

        r = Hx - obs_current
        r = np.sum(r, axis=1)

        ww = np.exp(-r**2 / (2. * obs_covar))
        ww /= np.sum(ww)

        Co = self.X @ self.X.T
        diag = np.diag(Co)
        Co = diag * np.ones((1,N)) - 2. * Co + np.ones((N,1)) * diag.T

        _, T = pyemd.emd_with_flow(ww,np.ones(N)/N, Co, -1)
        T = np.array(T)
        T = T*N

        self.X = np.dot(self.X.T,T) + self.delta * np.random.randn(self.X.T.shape[0],self.X.T.shape[1])

    @staticmethod
    def state_vector(ensemble):
        return ensemble.reshape(ensemble.shape[0],-1)

    def get_ensemble_from_X(self):
        return self.X.reshape(self.ensemble_shape)

        