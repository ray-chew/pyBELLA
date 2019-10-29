import numpy as np
from copy import deepcopy

class ensemble(object):
    def __init__(self, input_ensemble=None):
        if input_ensemble != None:
            cnt = 0
            for mem in input_ensemble:
                setattr(self,'mem_' + str(cnt),mem)
        else:
            None
        

    def initialise_members(self,ic,N):
        for cnt in range(N):
            mem = [deepcopy(arr) for arr in ic]
            # mem = sampler(mem)
            setattr(self,'mem_' + str(cnt),mem)

    # def state_vector(self,ensemble):
    #     N = self.members(ensemble).shape[0]
    #     return self.members(ensemble).reshape(N,-1)

    def set_members(self,analysis_ensemble):
        cnt = 0
        for xi in analysis_ensemble:
            setattr(self,'mem_' + str(cnt),np.array(xi))
            cnt += 1

    # rethink this eveutally....
    def ensemble_spreading(self, ens, sampler, attributes, loc=0):
        N = self.members(ens).shape[0]
        for attribute in attributes:
            for n in range(N):
                mem = getattr(self,'mem_' + str(n))
                value = getattr(mem[loc],attribute)
                setattr(mem[loc],attribute,sampler(value))

    @staticmethod
    def members(ensemble):
        return np.array(list(ensemble.__dict__.values()))

class sampler(object):
    def __init__(self, ic):
        # ic is a class instance of type variable container
        None