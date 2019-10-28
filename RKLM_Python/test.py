# unit tester for DA infrastructure

import numpy as np
from data_assimilation.utils import ensemble
from data_assimilation.letkf import letkf

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

localisation_matrix = np.eye(5,5)
forward_operator = np.eye(5,5)

obs = np.random.randn(5,5)
obs_covar = np.eye(25,25)

da = letkf(ens)
da.forward(forward_operator)
da.localisation(localisation_matrix)
analysis_ensemble = da.analyse(obs,obs_covar)
ens.set_members(analysis_ensemble)

print(ens.members(ens))

# print(sampler(input_array))

# print(sampler(5))