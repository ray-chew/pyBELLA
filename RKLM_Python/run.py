import sys
import subprocess
import numpy as np

import os
os.chdir('../')

N = 10

restart = False
#########################################################
# simulation restart parameters
#########################################################
path = '/home/ray/git-projects/RKLM_Reference/output_swe/'
fn = 'output_swe_ensemble=1_64_1_64_86400.0_dvortex_3D.h5'
path += fn
name = '_ensemble_mem=0_86400.000_after_full_step'
# time = np.arange(432000.0,864000.0+1200,1200)
time_start = 86400.0
time_end = 86400.0*2.0+1200
time_int = 1200.0

time_start, time_end, time_int = str(time_start), str(time_end), str(time_int)

#########################################################

if restart == False:
    subprocess.call([sys.executable, './RKLM_Python/__main__.py', '-ic', 'swe_bal_vortex', '-N', '%i' %N])

elif restart == True:
    subprocess.call([sys.executable, './RKLM_Python/__main__.py', '-ic', 'swe_test', 'restart', '-p', path, '-n', name, '-t', time_start, time_end, time_int])