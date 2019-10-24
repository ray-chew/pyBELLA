import numpy as np

class ensemble(object):
    def __init__(self):
        None
        
    def ensemble(self):
        for key,value in vars(self).items():
            setattr(self,key,value)

    def set_members(self,ic,N,sampler):
        cnt = 0
        for _ in range(N):
            mem = deepcopy(ic)
            sampler(mem)
            setattr(self,'mem_' + str(cnt),mem)
            cnt += 1

class sampler(object):
    def __init__(self, ic):
        # ic is a class instance of type variable container
        None