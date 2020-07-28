import sys
import subprocess
import numpy as np

N = 1



restart = True
path = '/home/ray/git-projects/RKLM_Reference/output_swe/'
fn = 'output_swe_ensemble=1_64_1_64_432000.0_dvortex_3D.h5'
path += fn
name = '_ensemble_mem=0_432000.000_after_full_step'
# time = np.arange(432000.0,864000.0+1200,1200)
time_start = 432000.0
time_end = 864000.0+1200
time_int = 1200.0

time_start, time_end, time_int = str(time_start), str(time_end), str(time_int)

if restart == False:
    subprocess.call([sys.executable, './RKLM_Python/__main__.py', '-ic', 'swe_test', 'argument2'])

elif restart == True:
    subprocess.call([sys.executable, './RKLM_Python/__main__.py', '-ic', 'swe_test', 'restart', '-p', path, '-n', name, '-t', time_start, time_end, time_int])


