# input/output
import h5py
import os
import numpy as np
import shutil # for copying of simulation restart file
import pickle

from management.variable import Vars
from physics.low_mach.mpv import MPV

import argparse

class io(object):
    """
    HDF5 writer class. Contains methods to create HDF5 file, create data sets and populate them with output variables.

    """
    def __init__(self,ud,restart=False):
        """
        Creates HDF5 file based on filename given in user input data.

        Parameters
        ----------
        ud : :class:`inputs.user_data.UserDataInit`
            Data container for the initial conditions

        """
        self.ud = ud

        self.FORMAT = ".h5"
        self.BASE_NAME = self.ud.output_base_name
        self.OUTPUT_FILENAME = "./output" + self.BASE_NAME
        self.OUTPUT_FILENAME = self.OUTPUT_FILENAME + "/output"

        self.SUFFIX = self.ud.output_suffix
        if restart: self.OLD_SUFFIX = self.ud.old_suffix

        if self.ud.continuous_blending == True:
            self.SUFFIX += "_cont_blend"
            self.SUFFIX += "_fs=%i_ts=%i" %(ud.no_of_pi_initial, ud.no_of_pi_transition)
                
        self.PATHS = [  'buoy',
                        # 'dp2_c',
                        'dp2_nodes',
                        # 'dpdim',
                        'drhoY',
                        'dT',
                        'dY',
                        # 'p',
                        # 'p2_c',
                        'p2_nodes',
                        'rho',
                        # 'rhoe',
                        'rhoY',
                        # 'S',
                        # 'T',
                        # 'u',
                        # 'v',
                        # 'w',
                        # 'vortz',
                        'vorty',
                        'Y',
                        'rhs'
                    ]
 
        self.io_create_file(self.PATHS,restart)

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
        # If file exists, rename it with old.
        if os.path.exists(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT):
            os.rename(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + '_old' + self.FORMAT)
            
        # create a new output file for each rerun - old output will be overwritten.
        if restart:
            src = self.OUTPUT_FILENAME + self.BASE_NAME + self.OLD_SUFFIX + self.FORMAT
            dest = self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT
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
        th : :class:`physics.gas_dynamics.thermodynamic.ThemodynamicInit`
            Thermodynamic variables of the system
        name: str
            The time and additional suffix label for the dataset, e.g. "_10.0_after_full_step", where 10.0 is the time and "after_full_step" denotes when the output was made.

        """
        print("writing hdf output..." + name)
        # rho
        self.populate(name,'rho',Sol.rho)
        # rhoe
        # self.populate(name,'rhoe',Sol.rhoe)
        # rhoY
        self.populate(name,'rhoY',Sol.rhoY)

        # rho u ,v w
        self.populate(name,'rhou',Sol.rhou)
        self.populate(name,'rhov',Sol.rhov)
        self.populate(name,'rhow',Sol.rhow)
        self.populate(name,'rhoX',Sol.rhoX)
        # self.populate(name,'buoy',Sol.rhoX / Sol.rho)

        # dp2_nodes
        self.populate(name,'dp2_nodes',mpv.dp2_nodes)
        # p2_nodes
        self.populate(name,'p2_nodes',mpv.p2_nodes)
        # p2_nodes0
        # self.populate(name,'p2_nodes0',mpv.p2_nodes0)
        # p2_cells
        self.populate(name,'p2_cells',mpv.p2_cells)
        # dp2_cells
        self.populate(name,'dp2_cells',mpv.dp2_cells)

        # pressure
        # self.populate(name,'p',Sol.rhoY**th.gamm)
        # pressure difference
        # self.populate(name,'dpdim', self.dpress_dim(mpv,self.ud,th))

        # velocity (u,v,w)
        # self.populate(name,'u',Sol.rhou / Sol.rho)
        # self.populate(name,'v',Sol.rhov / Sol.rho)
        # self.populate(name,'w',Sol.rhow / Sol.rho)

        # vorticity in (x,z)
        # self.populate(name,'vortz', self.vortz(Sol,elem,node))
        # self.vorty(Sol,elem,node)
        self.populate(name,'vorty', self.vorty(Sol,elem,node))

        # temperature
        # self.populate(name,'T', Sol.rhoY**th.gamm / Sol.rho)
        # temperature difference
        # print(mpv.HydroState.p0[0,:])
        # print(mpv.HydroState.rho0[0,:])
        # self.populate(name,'dT', Sol.rhoY**th.gamm / Sol.rho - mpv.HydroState.p0[:] / mpv.HydroState.rho0[:])

        # self.populate(name,'drhoY', Sol.rhoY - mpv.HydroState.rho0[:] * mpv.HydroState.Y0[:])

        # species mass fraction(?)
        # self.populate(name,'Y', Sol.rhoY / Sol.rho)
        # species mass fraction perturbation
        # self.populate(name,'dY', Sol.rhoY / Sol.rho - self.ud.stratification(elem.y))

        # self.populate(name,'wplusx',mpv.wplus[0])
        # self.populate(name,'wplusy',mpv.wplus[1])
        # self.populate(name,'hcenter',mpv.wcenter)
        # self.populate(name,'p2',mpv.dp2_nodes)
        # self.populate(name,'rhs',mpv.rhs)
        # self.populate(name,'X',Sol.rhoX/Sol.rho)

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
            return
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

        # inner_domain_rho = rho[igs[0]-1:-igs[0], igs[1]-1:-igs[1]]
        # inner_domain_rhou = rhou[igs[0]-1:-igs[0], igs[1]-1:-igs[1]]
        # inner_domain_rhov = rhov[igs[0]-1:-igs[0], igs[1]-1:-igs[1]]

        inner_domain_rho = rho
        inner_domain_rhou = rhou
        inner_domain_rhov = rhov

        top_left_idx     = (slice(0,-1)    , slice(0,-1))
        top_right_idx    = (slice(1, None) , slice(0,-1))
        bottom_left_idx  = (slice(0,-1)    , slice(1, None))
        bottom_right_idx = (slice(1, None) , slice(1,None))

        # print(inner_domain_rhov[bottom_right_idx])
        dvdx = 0.5 * ((inner_domain_rhov[bottom_right_idx] / inner_domain_rho[bottom_right_idx] - inner_domain_rhov[bottom_left_idx] / inner_domain_rho[bottom_left_idx]) + (inner_domain_rhov[top_right_idx] / inner_domain_rho[top_right_idx] - inner_domain_rhov[top_left_idx] / inner_domain_rho[top_left_idx])) / dx

        dudy = 0.5 * ((inner_domain_rhou[bottom_right_idx] / inner_domain_rho[bottom_right_idx] - inner_domain_rhou[top_right_idx] / inner_domain_rho[top_right_idx]) + (inner_domain_rhou[bottom_left_idx] / inner_domain_rho[bottom_left_idx] - inner_domain_rhou[top_left_idx] / inner_domain_rho[top_left_idx])) / dy

        vortz = np.zeros((node.sc)).squeeze()
        tmp = np.expand_dims((dvdx - dudy) * 1.0, axis=1)
        tmp = np.repeat(tmp, node.icy, axis=1)
        vortz[1:-1, :, 1:-1] = tmp
        # vortz = np.repeat(vortz[:,igs[1],:], node.icy, axis=1)
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
        #     file[str(path) + '/' + str(name)].attrs.create(key,options[key])
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


def get_args():
    """
    Argument parser for initial conditions and ensemble size.

    """

    parser = argparse.ArgumentParser(description='Python solver for unified numerical model based on (Benacchio and Klein, 2019) with data assimilation and blending. Written by Ray Chew and based on Prof. Rupert Klein\'s C code.')

    parser.add_argument('-N',action='store',dest='N',help='<Optional> Set ensemble size, if none is given N=1 is used.',required=False,type=int)

    parser.add_argument('-ic','--initial_conditions',action='store',dest='ic',help='<Required> Set initial conditions',required=True,choices={'aw','tv','tv_2d','tv_3d','tv_corr','rb', 'rbc', 'igw','swe','swe_bal_vortex','swe_icshear', 'swe_dvortex'})

    # parser.add_argument('-r','--restart_sim',action='store',dest='rstrt',help='<Optional> Restart simulation?.',required=False,default=False,type=bool)

    subparsers = parser.add_subparsers(dest='subcommand')
    restart = subparsers.add_parser('restart')
    restart.add_argument('-p', '--path', action='store', dest='path', help='path to data for simulation restart.', required=True, type=str)
    restart.add_argument('-n', '--name', action='store', dest='name', help='name of datasets for simulation restart.', required=True, type=str)
    restart.add_argument('-t', '--time', nargs="*", help='time outputs for simulation restart in format [start,stop,interval). Use None for ud.tout settings.', type=float, required=False, default=None)


    args = parser.parse_args() # collect cmd line args
    ic = args.ic

    if ic == 'bi':
        from inputs.baroclinic_instability_periodic import UserData, sol_init
    elif ic == 'tv' or ic == 'tv_2d':
        from inputs.travelling_vortex_2D import UserData, sol_init
    elif ic == 'tv_3d':
        from inputs.travelling_vortex_3D import UserData, sol_init
    elif ic == 'tv_corr':
        from inputs.travelling_vortex_3D_Coriolis import UserData, sol_init
    elif ic == 'aw':
        from inputs.acoustic_wave_high import UserData, sol_init
    elif ic == 'igw':
        from inputs.internal_long_wave import UserData, sol_init
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



    if UserData is None or sol_init is None:
        assert(0, "Initial condition file is not well defined.")
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


    return N, UserData, sol_init, rstrt, params



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
            assert(getattr(Sol,data).shape == shp)
        else:
            assert(0, "Sol attribute mismatch")

    for data in mpv_data:
        value = file[data][data+name][:]
        if hasattr(mpv,data):
            shp = getattr(mpv,data).shape
            setattr(mpv,data,value)
            assert(getattr(mpv,data).shape == shp)
        else:
            assert(0, "mpv attribute mismatch")

    t = restart_touts

    ud.output_suffix = "_%i_%i_%i_%.1f_%s" %(ud.inx-1,ud.iny-1,ud.inz-1,t[-1],ud.aux)


    file.close()
    return Sol, mpv, t