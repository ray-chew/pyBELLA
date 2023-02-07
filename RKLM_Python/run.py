import sys
import subprocess
import numpy as np

import os
os.chdir('../')


class run_params(object):
    N = 1
    tc = 'mark'

    def __init__(self):
        self.N = self.N
        self.tc = self.tc

        self.ud = "None"
        self.dap = "None"
        self.restart = False

    def single_run(self):
        subprocess.call([sys.executable, './RKLM_Python/__main__.py', '-ic', self.tc, '-N', '%i' %self.N])

    def queue_run(self):
        if self.ud is None and self.dap is None:
            assert 0, "ud or params must be defined"
        subprocess.call([sys.executable, './RKLM_Python/__main__.py', '-ic', self.tc, '-N', '%i' %self.N, 'queue', '-w', self.ud, self.dap])


    def restart_set(self, path, fn, name, ts, te, ti):
        path += fn
        ts, te, ti = str(ts), str(te), str(ti)

        self.path = path
        self.fn = fn
        self.name = name
        self.time_start = ts
        self.time_end = te
        self.time_int = ti
        self.restart = True

    def restart_run(self):
        if self.restart:
            subprocess.call([sys.executable, './RKLM_Python/__main__.py', '-ic', self.tc, '-N', '%i' %self.N, 'restart', '-p', self.path, '-n', self.name, '-t', self.time_start, self.time_end, self.time_int])
        else:
            print("restart parameters not found.")


if __name__ == '__main__':
    rp = run_params()
    rp.single_run()

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

    # if restart == False:
    #     subprocess.call([sys.executable, './RKLM_Python/__main__.py', '-ic', 'tv', '-N', '%i' %N])

    # elif restart == True:
    #     subprocess.call([sys.executable, './RKLM_Python/__main__.py', '-ic', 'swe_test', 'restart', '-p', path, '-n', name, '-t', time_start, time_end, time_int])