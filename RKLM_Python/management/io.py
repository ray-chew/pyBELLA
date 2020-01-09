# input/output
import h5py
import os
import numpy as np

class io(object):
    """
    HDF5 writer class. Contains methods to create HDF5 file, create data sets and populate them with output variables.

    """
    def __init__(self,ud):
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
        # if ud.is_ArakawaKonor:
        #     self.SUFFIX = self.ud.output_name_ak
        # else:
        #     if self.ud.is_nonhydrostatic == 0:
        #         if self.ud.is_compressible == 1:
        #             self.SUFFIX = self.ud.output_name_hydro
        #         else:
        #             assert("Not implemented")
        #     elif self.ud.is_nonhydrostatic == 1:
        #         if self.ud.is_compressible == 1:
        #             self.SUFFIX = self.ud.output_name_comp
        #         else:
        #             self.SUFFIX = self.ud.output_name_psinc

        if self.ud.continuous_blending == True:
            self.SUFFIX += "_cont_blend"
            self.SUFFIX += "_fs=%i_ts=%i" %(ud.no_of_pi_initial, ud.no_of_pi_transition)
                
        # else:
        #     self.SUFFIX = "_PIs1=" + str(self.ud.no_of_pi_initial) + \
        #         "_PIs2=" + str(self.ud.no_of_pi_transition) + \
        #         "_HYs1=" + str(self.ud.no_of_hy_initial) + \
        #         "_HYs2=" + str(self.ud.no_of_hy_transition) + \
        #         "_contblend"


        self.PATHS = [  #'bouy',
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
                        'Y',
                        'rhs'
                    ]
        
        self.io_create_file(self.PATHS)

    def io_create_file(self,paths):
        """
        Helper function to create file.

        Parameters
        ----------
        paths : list
            List of strings containing the name of the data sets. For now,

                PATH = ['dp2_nodes', 'drhoY', 'dT', 'dY', 'p2_nodes', 'rho', 'rhoY', 'Y']

        Notes
        -----
        Currently, if the filename of the HDF5 file already exists, this function will delete the file and create an empty HDF5 file with the same filename in its place. This was enabled to prevent certain parallel-writing issues.

        """
        # If file exists, delete it.. Nuclear option
        if os.path.exists(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT):
            os.remove(self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT)
        # create a new output file for each rerun - old output will be overwritten.
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

        # dp2_nodes
        self.populate(name,'dp2_nodes',mpv.dp2_nodes)
        # p2_nodes
        self.populate(name,'p2_nodes',mpv.p2_nodes)
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

        # temperature
        # self.populate(name,'T', Sol.rhoY**th.gamm / Sol.rho)
        # temperature difference
        # print(mpv.HydroState.p0[0,:])
        # print(mpv.HydroState.rho0[0,:])
        # self.populate(name,'dT', Sol.rhoY**th.gamm / Sol.rho - mpv.HydroState.p0[:] / mpv.HydroState.rho0[:])

        # self.populate(name,'drhoY', Sol.rhoY - mpv.HydroState.rho0[:] * mpv.HydroState.Y0[:])

        # species mass fraction(?)
        self.populate(name,'Y', Sol.rhoY / Sol.rho)
        # species mass fraction perturbation
        self.populate(name,'dY', Sol.rhoY / Sol.rho - self.ud.stratification(elem.y))

        # self.populate(name,'wplusx',mpv.wplus[0])
        # self.populate(name,'wplusy',mpv.wplus[1])
        # self.populate(name,'hcenter',mpv.wcenter)
        # self.populate(name,'p2',mpv.dp2_nodes)
        self.populate(name,'rhs',mpv.rhs)
        # self.populate(name,'X',Sol.rhoX/Sol.rho)

    def vortz(self,Sol,elem,node):
        """
        Calculate the vorticity of the solution.

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
            # print(key)
            # print(value)
            # print(value.__name__)
            try:
                file.attrs.create(key,value)
            except:
                # print(str(repr(value)))
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