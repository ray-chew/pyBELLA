import numpy as np

# simulation parameters
debug = False
da_debug = False

random_seed = 888

print_precision = 18


output_path = './outputs'


# global constants
class global_constants(object):
    def __init__(self):
        self.nspec = 1
        self.buoy = 0

        self.grav = 9.81                 # [m s^{-2}]
        self.omega = 0.0                 # [s^{-1}]

        self.R_gas = 287.4               # [J kg^{-1} K^{-1}]
        self.R_vap = 461.0
        self.Q_vap = 2.53e+06
        self.gamm = 1.4

        self.p_ref = 8.61 * 1e4          # [N/m^2]
        self.T_ref = 300.00              # [K]

        self.h_ref = 10000.0              # [m]
        self.t_ref = 100.0               # [s]