# unit tester for DA infrastructure

import numpy as np
from data_assimilation.utils import ensemble
from data_assimilation.enkf import enkf

import numpy as np
from copy import deepcopy

input_array = np.arange(25).reshape(5,5)
N = 20
sampler = lambda ic : ic + np.random.randint(5)

# def sampler(ic):
#     # ic.A += np.random.randint(5)
#     ic += np.random.randint(5)

ens = ensemble()
ens.initialise_members(input_array,N,sampler)

print(ens.members(ens).reshape(N,-1).shape)
# print(sampler(input_array))

# print(sampler(5))