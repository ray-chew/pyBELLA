import numpy as np

# simulation parameters
debug = False
da_debug = False

random_seed = 888

print_precision = 18


output_path = './outputs'


# global constants
NSPEC = 1
grav = 9.81                 # [m s^{-2}]
omega = 7.292 * 1e-5        # [s^{-1}]

R_gas = 287.4               # [J kg^{-1} K^{-1}]
R_vap = 461.0
Q_vap = 2.53e+06
gamma = 1.4
cp_gas = gamma * R_gas / (gamma-1.0)

p_ref = 1e+5
T_ref = 300.00              # [K]
rho_ref = p_ref / (R_gas * T_ref)
N_ref = grav / np.sqrt(cp_gas * T_ref)
Cs = np.sqrt(gamma * R_gas * T_ref)

h_ref = 10.0e3              # [m]
t_ref = 100.0               # [s]
u_ref = h_ref / t_ref