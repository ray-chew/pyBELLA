import numpy as np
from copy import deepcopy

class ensemble(object):
    def __init__(self):
        None
        
    def ensemble(self):
        for key,value in vars(self).items():
            setattr(self,key,value)

    def initialise_members(self,ic,N,sampler):
        for cnt in range(N):
            mem = deepcopy(ic)
            mem = sampler(mem)
            setattr(self,'mem_' + str(cnt),mem)

    def state_vector(self,ensemble):
        N = ensemble.members(ensemble).shape[0]
        return self.members(ensemble).reshape(N,-1)

    def set_members(self,analysis_ensemble):
        cnt = 0
        for xi in analysis_ensemble:
            setattr(self,'mem_' + str(cnt),xi.reshape(self.mem_0.shape))
            cnt += 1


    @staticmethod
    def members(ensemble):
        return np.array(list(ensemble.__dict__.values()))

class sampler(object):
    def __init__(self, ic):
        # ic is a class instance of type variable container
        None