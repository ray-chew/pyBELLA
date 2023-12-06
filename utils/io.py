# input/output
import h5py
import os
import numpy as np
import shutil # for copying of simulation restart file
import dill as pickle # pickle jar to debug classes
import yaml # for parsing of dict-style arguments

import argparse

from utils.sim_params import output_path

class io(object):
    """
    HDF5 writer class. Contains methods to create HDF5 file, create data sets and populate them with output variables.

    """
    def __init__(self,ud,restart=False):
        """
        Creates HDF5 file based on filename given attribute `OUTPUT_FILENAME`.

        Parameters
        ----------
        ud : :class:`inputs.user_data.UserDataInit`
            Data container for the initial conditions

        """
        self.ud = ud

        self.FORMAT = ".h5"
        self.BASE_NAME = self.ud.output_base_name
        self.OUTPUT_FILENAME = self.ud.output_type + self.BASE_NAME
        self.OUTPUT_FOLDER = output_path + '/' + self.OUTPUT_FILENAME
        self.OUTPUT_FILENAME = self.OUTPUT_FOLDER + "/" + self.ud.output_type

        self.SUFFIX = self.ud.output_suffix
        if restart: self.OLD_SUFFIX = self.ud.old_suffix
                
        self.PATHS = [  #'buoy',
                        # 'dp2_c',
                        # 'dp2_nodes',
                        # 'dpdim',
                        # 'drhoY',
                        # 'dT',
                        # 'dY',
                        # 'p',
                        # 'p2_c',
                        'p2_nodes',
                        'rho',
                        'rhoY',
                        # 'S',
                        # 'T',
                        # 'u',
                        # 'v',
                        # 'w',
                        # 'vortz',
                        # 'vorty',
                        # 'Y',
                        # 'rhs'
                    ]
 
        self.io_create_file(self.PATHS,restart)
        self.time = 'None'

    def io_create_file(self,paths,restart):
        """
        Helper function to create file.

        Parameters
        ----------
        paths : list
            List of strings containing the name of the data sets. For now,

                PATH = ['dp2_nodes', 'drhoY', 'dT', 'dY', 'p2_nodes', 'rho', 'rhoY', 'Y']

        Notes
        -----
        Currently, if the filename of the HDF5 file already exists, this function will append the existing filename with '_old' and create an empty HDF5 file with the same filename in its place.

        """
        # If directory does not exist, create it.
        if not os.path.exists(self.OUTPUT_FOLDER):
            os.makedirs(self.OUTPUT_FOLDER)

        # If file exists, rename it with old.
        if os.path.exists(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT):
            os.rename(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + '_old' + self.FORMAT)
            
        # create a new output file for each rerun - old output will be overwritten.
        if restart:
            src = self.OUTPUT_FILENAME + self.BASE_NAME + self.OLD_SUFFIX + self.FORMAT
            dest = self.OUTPUT_FILENAME + self.BASE_NAME + self.OLD_SUFFIX + '_old' + self.FORMAT
            shutil.copy2(src, dest)
        else:
            file = h5py.File(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, 'a')
            for path in paths:
                # check if groups have been created
                # if not created, create empty groups
                if not (path in file):
                    file.create_group(path,track_order=True)
            
            file.close()


    def write_all(self,Sol,mpv,elem,node,th,name):
        """
        At a given time, write output from `Sol` and `mpv` to the HDF5 file.

        Parameters
        ----------
        Sol : :class:`management.variable.Vars`
            Solution data container
        mpv : :class:`physics.low_mach.mpv.MPV`
            Variables relating to the elliptic solver
        elem : :class:`discretization.kgrid.ElemSpaceDiscr`
            Cells grid
        node : :class:`discretization.kgrid.NodeSpaceDiscr`
            Nodes grid
        th : :class:`physics.gas_dynamics.thermodynamic.ThermodynamicInit`
            Thermodynamic variables of the system
        name: str
            The time and additional suffix label for the dataset, e.g. "_10.0_after_full_step", where 10.0 is the time and "after_full_step" denotes when the output was made.

        """
        print("writing hdf output..." + name)
        # rho
        self.populate(name,'rho',Sol.rho)
        # rhoY
        self.populate(name,'rhoY',Sol.rhoY)

        # rho u ,v w
        self.populate(name,'rhou',Sol.rhou)
        self.populate(name,'rhov',Sol.rhov)
        self.populate(name,'rhow',Sol.rhow)
        self.populate(name,'rhoX',Sol.rhoX)

        # if hasattr(Sol,'rhov_half'):
        #     self.populate(name,'rhov_half',Sol.rhov_half)
        # if hasattr(Sol,'rhou_half'):
        #     self.populate(name,'rhou_half',Sol.rhou_half)
        # if hasattr(Sol,'rhow_half'):
        #     self.populate(name,'rhow_half',Sol.rhow_half)
        # if hasattr(Sol,'rhoX_half'):
        #     self.populate(name,'rhoX_half',Sol.rhoX_half)
        # if hasattr(Sol,'rhoY_half'):
        #     self.populate(name,'rhoY_half',Sol.rhoY_half)
        # if hasattr(Sol,'rho_half'):
        #     self.populate(name,'rho_half',Sol.rho_half)
        # if hasattr(mpv,'p2_nodes_half'):
        #     self.populate(name,'p2_nodes_half',mpv.p2_nodes_half)
        # self.populate(name,'buoy',Sol.rhoX / Sol.rho)


        self.populate(name,'p2_nodes',mpv.p2_nodes)


        # vorticity in (x,z)
        # self.populate(name,'vortz', self.vortz(Sol,elem,node))
        # self.vorty(Sol,elem,node)
        # self.populate(name,'vorty', self.vorty(Sol,elem,node))

    def vortz(self,Sol,elem,node):
        """
        Calculate the vorticity of the solution in the x-y plane.

        Parameters
        ----------
        Sol : :class:`management.variable.Vars`
            Solution data container
        elem : :class:`discretization.kgrid.ElemSpaceDiscr`
            Cells grid
        node : :class:`discretization.kgrid.NodeSpaceDiscr`
            Nodes grid

        Returns
        -------
        ndarray
            An ndarray with the vorticity of the solution and with the shape of `node`.

        """
        if elem.ndim != 2:
            return np.zeros_like(Sol.rho)
        # 2d-case
        igs = elem.igs
        dx = elem.dx
        dy = elem.dy

        inner_domain_rho = Sol.rho[igs[0]-1:-igs[0], igs[1]-1:-igs[1]]
        inner_domain_rhou = Sol.rhou[igs[0]-1:-igs[0], igs[1]-1:-igs[1]]
        inner_domain_rhov = Sol.rhov[igs[0]-1:-igs[0], igs[1]-1:-igs[1]]

        top_left_idx     = (slice(0,-1)    , slice(0,-1))
        top_right_idx    = (slice(1, None) , slice(0,-1))
        bottom_left_idx  = (slice(0,-1)    , slice(1, None))
        bottom_right_idx = (slice(1, None) , slice(1,None))

        # print(inner_domain_rhov[bottom_right_idx])
        dvdx = 0.5 * ((inner_domain_rhov[bottom_right_idx] / inner_domain_rho[bottom_right_idx] - inner_domain_rhov[bottom_left_idx] / inner_domain_rho[bottom_left_idx]) + (inner_domain_rhov[top_right_idx] / inner_domain_rho[top_right_idx] - inner_domain_rhov[top_left_idx] / inner_domain_rho[top_left_idx])) / dx

        dudy = 0.5 * ((inner_domain_rhou[bottom_right_idx] / inner_domain_rho[bottom_right_idx] - inner_domain_rhou[top_right_idx] / inner_domain_rho[top_right_idx]) + (inner_domain_rhou[bottom_left_idx] / inner_domain_rho[bottom_left_idx] - inner_domain_rhou[top_left_idx] / inner_domain_rho[top_left_idx])) / dy

        vortz = np.zeros((node.sc)).squeeze()

        vortz[igs[0]:-igs[0]-1, igs[1]:-igs[1]-1] = dvdx - dudy
        return vortz

    def vorty(self,Sol,elem,node):
        """
        Calculate the vorticity of the solution in the x-z plane.

        Parameters
        ----------
        Sol : :class:`management.variable.Vars`
            Solution data container
        elem : :class:`discretization.kgrid.ElemSpaceDiscr`
            Cells grid
        node : :class:`discretization.kgrid.NodeSpaceDiscr`
            Nodes grid

        Returns
        -------
        ndarray
            An ndarray with the vorticity of the solution and with the shape of `node`.

        """
        if elem.ndim != 3 or (elem.ndim == 3 and elem.iicy > 1):
            return np.zeros_like(Sol.rho)
        # 2d-case
        igs = elem.igs
        dx = elem.dx
        dy = elem.dz

        rho = Sol.rho[:,igs[1],:]
        rhou = Sol.rhou[:,igs[1],:]
        rhov = Sol.rhow[:,igs[1],:]

        inner_domain_rho = rho
        inner_domain_rhou = rhou
        inner_domain_rhov = rhov

        top_left_idx     = (slice(0,-1)    , slice(0,-1))
        top_right_idx    = (slice(1, None) , slice(0,-1))
        bottom_left_idx  = (slice(0,-1)    , slice(1, None))
        bottom_right_idx = (slice(1, None) , slice(1,None))

        dvdx = 0.5 * ((inner_domain_rhov[bottom_right_idx] / inner_domain_rho[bottom_right_idx] - inner_domain_rhov[bottom_left_idx] / inner_domain_rho[bottom_left_idx]) + (inner_domain_rhov[top_right_idx] / inner_domain_rho[top_right_idx] - inner_domain_rhov[top_left_idx] / inner_domain_rho[top_left_idx])) / dx

        dudy = 0.5 * ((inner_domain_rhou[bottom_right_idx] / inner_domain_rho[bottom_right_idx] - inner_domain_rhou[top_right_idx] / inner_domain_rho[top_right_idx]) + (inner_domain_rhou[bottom_left_idx] / inner_domain_rho[bottom_left_idx] - inner_domain_rhou[top_left_idx] / inner_domain_rho[top_left_idx])) / dy

        vortz = np.zeros((node.sc)).squeeze()
        tmp = np.expand_dims((dvdx - dudy) * 1.0, axis=1)
        tmp = np.repeat(tmp, node.icy, axis=1)
        vortz[1:-1, :, 1:-1] = tmp
        return vortz


    def dpress_dim(self,mpv,ud,th):
        p0 = (th.Gamma * ud.Msq * mpv.p2_cells)
        p = np.power(p0,th.Gammainv, dtype=np.complex)
        p = p.real
        return (p - mpv.HydroState.p0[0,:]) * self.ud.p_ref


    def populate(self,name,path,data,options=None):
        """
        Helper function to write data into HDF5 dataset.

        Parameters
        ----------
        name : str
            The time and additional suffix label for the dataset
        path : str
            Path of the dataset, e.g. `rhoY`.
        data : ndarray
            The output data to write to the dataset
        options : list
            `default == None`. Additional options to write to dataset, currently unused.

        """
        # name is the simulation time of the output array
        # path is the array type, e.g. U,V,H, and data is it's data.
        file = h5py.File(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, 'r+')
        file.create_dataset(str(path) + '/' + str(path) + '_' + str(name), data=data, chunks=True, compression='gzip', compression_opts=4, dtype=np.float32)
        # add attributes, i.e. the simulation parameters to each dataset.
        # for key in options:
            # file[str(path) + '/' + str(name)].attrs.create(key,options[key])
        try:
            file[str(path)][str(path) + '_' + str(name)].attrs.create('t',self.time)
        except:
            # file.attrs.create(key,repr(value),dtype='<S' + str(len(repr(value))))
            file[str(path)][str(path) + '_' + str(name)].attrs.create('t',repr(self.time), dtype='<S' + str(len(repr(self.time))))
        # print("writing time = %.1f for arrays %s" %(name,path))
        file.close()

    def write_attrs(self):
        """
        Method to write all attributes in the userdata initial condition to HDF5 file.
        """
        file = h5py.File(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, 'a')
        for key, value in vars(self.ud).items():
            try:
                file.attrs.create(key,value)
            except:
                file.attrs.create(key,repr(value),dtype='<S' + str(len(repr(value))))
        file.close()

    def write_da_attrs(self,params):
        """
        Method to write all data-assimilation attributes in the userdata initial condition to HDF5 file.
        """
        file = h5py.File(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, 'a')
        path = 'da_parameters'
        if not (path in file):
            file.create_group(path,track_order=True)

        for key, value in vars(params).items():
            try:
                file['da_parameters'].attrs.create(key,value)
            except:
                # print(str(repr(value)))
                file['da_parameters'].attrs.create(key,repr(value),dtype='<S' + str(len(repr(value))))
        file.close()

    def close_everything(self):
        """
        In parallel, some workers do not close file correctly, this function forces the program to close the HDF file before exiting.
        """
        file = h5py.File(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, 'r')
        if file.__bool__():
            file.close()


    def check_jar(self):
        fn = self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + '.dat'
        if os.path.exists(fn):
            os.rename(fn, self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + '_old.dat')

    def jar(self, content = None):
        fn = self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + '.dat'
        if os.path.exists(fn):
            file = open(fn,"ab")
        else:
            file = open(fn,"wb")
        
        # let's fill the pickle jar
        if content is not None:
            for each_pickle in content:
                pickle.dump(each_pickle,file)
        file.close()


class read_input(object):
    
    def __init__(self, fn, path):
        self.fn = fn
        self.path = path

    def get_data(self, Sol, mpv, time_tag, half=False):
        file = h5py.File(self.path + '/' + self.fn, 'r')

        if half:
            half_tag = '_half'
        else:
            half_tag = ''

        Sol.rho[...] = file['rho' + half_tag]['rho' + half_tag + '_' + time_tag][:]
        Sol.rhou[...] = file['rhou' + half_tag]['rhou' + half_tag + '_' + time_tag][:]
        Sol.rhov[...] = file['rhov' + half_tag]['rhov' + half_tag + '_' + time_tag][:]
        Sol.rhow[...] = file['rhow' + half_tag]['rhow' + half_tag + '_' + time_tag][:]
        Sol.rhoY[...] = file['rhoY' + half_tag]['rhoY' + half_tag + '_' + time_tag][:]

        mpv.p2_nodes[...] = file['p2_nodes' + half_tag]['p2_nodes' + half_tag + '_' + time_tag][:]

        file.close()
    

def get_args():
    """
    Argument parser for initial conditions and ensemble size.

    """

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-N',action='store',dest='N',help='<Optional> Set ensemble size, if none is given N=1 is used.',required=False,type=int)

    parser.add_argument('-ic', '--initial_conditions',
                        action='store',
                        dest='ic',
                        help='<Required> Set initial conditions',
                        required=True,
                        choices={'aw',
                                 'tv',
                                 'tv_neg',
                                 'tv_2d',
                                 'tv_3d',
                                 'tv_corr',
                                 'rb',
                                 'rbc',
                                 'igw',
                                 'igw_3d',
                                 'lbw',
                                 'skl',
                                 'mark',
                                 'lw_p',
                                 'igw_bb',
                                 'swe',
                                 'swe_bal_vortex',
                                 'swe_icshear',
                                 'swe_dvortex',
                                 'test_travelling_vortex',
                                 'test_internal_long_wave',
                                 'test_lamb_wave'
                                 }
                        )

    subparsers = parser.add_subparsers(dest='subcommand')

    restart = subparsers.add_parser('restart')
    restart.add_argument('-p', '--path', action='store', dest='path', help='path to data for simulation restart.', required=True, type=str)
    restart.add_argument('-n', '--name', action='store', dest='name', help='name of datasets for simulation restart.', required=True, type=str)
    restart.add_argument('-t', '--time', nargs="*", help='time outputs for simulation restart in format [start,stop,interval). Use None for ud.tout settings.', type=float, required=False, default=None)


    queue = subparsers.add_parser('queue')
    queue.add_argument('-w', '--rewrite', nargs="*", help='', required=True, type=yaml.safe_load)

    args = parser.parse_args() # collect cmd line args
    ic = args.ic

    if ic == 'bi':
        from inputs.baroclinic_instability_periodic import UserData, sol_init
    elif ic == 'tv' or ic == 'tv_2d':
        from inputs.travelling_vortex_2D import UserData, sol_init
    elif ic == 'tv_neg':
        from inputs.travelling_vortex_2D_neg import UserData, sol_init
    elif ic == 'tv_3d':
        from inputs.travelling_vortex_3D import UserData, sol_init
    elif ic == 'tv_corr':
        from inputs.travelling_vortex_3D_Coriolis import UserData, sol_init
    elif ic == 'aw':
        from inputs.acoustic_wave_high import UserData, sol_init
    elif ic == 'igw':
        from inputs.internal_long_wave import UserData, sol_init
    elif ic == 'igw_3d':
        from inputs.internal_long_wave_3D import UserData, sol_init
    elif ic == 'lbw':
        from inputs.lamb_waves import UserData, sol_init
    elif ic == 'skl':
        from inputs.sk_lamb_wave import UserData, sol_init
    elif ic == 'mark':
        from inputs.mark import UserData, sol_init
    elif ic == 'lw_p':
        from inputs.lamb_wave_perturb import UserData, sol_init
    elif ic == 'igw_bb':
        from inputs.igw_baldauf_brdar import UserData, sol_init
    elif ic == 'rb':
        from inputs.rising_bubble import UserData, sol_init
    elif ic == 'rbc':
        from inputs.rising_bubble_cold import UserData, sol_init
    elif ic == 'swe_bal_vortex':
        from inputs.swe_bal_vortex import UserData, sol_init
    elif ic == 'swe':
        from inputs.shallow_water_3D import UserData, sol_init
    elif ic == 'swe_icshear':
        from inputs.shallow_water_3D_icshear import UserData, sol_init
    elif ic == 'swe_dvortex':
        from inputs.shallow_water_3D_dvortex import UserData, sol_init
    elif ic == 'test_travelling_vortex':
        from tests.test_travelling_vortex import UserData, sol_init
    elif ic == 'test_internal_long_wave':
        from tests.test_internal_long_wave import UserData, sol_init
    elif ic == 'test_lamb_wave':
        from tests.test_lamb_wave import UserData, sol_init


    if UserData is None or sol_init is None:
        assert 0, "Initial condition file is not well defined."
    if args.N is None:
        N = 1
    else:
        N = args.N


    if args.subcommand == 'restart':
        rstrt = True
        if args.time is not None:
            t_vals = args.time
            time = np.arange(t_vals[0],t_vals[1],t_vals[2])
        params = [args.path,args.name,time]
    else:
        rstrt = False
        params = None

    if args.subcommand == 'queue':
        ud = args.rewrite[0]
        dap = args.rewrite[1]
    else:
        ud = None
        dap = None


    return N, UserData, sol_init, rstrt, ud, dap, params



def sim_restart(path, name, elem, node, ud, Sol, mpv, restart_touts):
    """
    Function to restart simulation from a saved file. Dataset has to be structured in the same way as the output format of this file.

    Parameters
    ----------
    path : str
        path to the hdf5 file for simulation restart.

    """
    file = h5py.File(str(path), 'r')

    Sol_data = ['rho','rhou','rhov','rhow','rhoX', 'rhoY']
    mpv_data = ['p2_nodes']

    for data in Sol_data:
        value = file[data][data+name][:]
        
        if hasattr(Sol,data):
            shp = getattr(Sol,data).shape
            setattr(Sol,data,value)
            assert getattr(Sol,data).shape == shp
        else:
            assert 0, "Sol attribute mismatch"

    for data in mpv_data:
        value = file[data][data+name][:]
        if hasattr(mpv,data):
            shp = getattr(mpv,data).shape
            setattr(mpv,data,value)
            assert getattr(mpv,data).shape == shp
        else:
            assert 0, "mpv attribute mismatch"

    t = restart_touts

    ud.output_suffix = "_%i_%i_%i_%.1f_%s" %(ud.inx-1,ud.iny-1,ud.inz-1,t[-1],ud.aux)


    file.close()
    return Sol, mpv, t


def fn_gen(ud, dap, N):

    suffix = ""
    suffix += "_%i" %(ud.inx-1)
    suffix += "_%i" %(ud.iny-1)
    if ud.iny == 2:
        suffix += "_%i" %(ud.inz-1)
    suffix += "_%.6f" %ud.tout[-1]
    suffix = '_ensemble=%i%s' %(N, suffix)

    if len(dap.da_times) > 0 and N >1:
        suffix += '_wda'
        if dap.da_type == 'rloc':
            suffix += 'wloc'
        if len(dap.obs_attributes) < 5:
            for attr in dap.obs_attributes:
                suffix += '_%s' %attr
        else:
            suffix += '_all'

    if ud.aux is not None:
        suffix += '_' + ud.aux
    
    if ud.initial_blending:
        bw = int(ud.blending_weight * 16)
        suffix += '_ib-%i' %bw 

    if ud.continuous_blending == True:
        suffix += "_cont_blend"
        suffix += "_fs=%i_ts=%i" %(ud.no_of_pi_initial, ud.no_of_pi_transition)

    return suffix
