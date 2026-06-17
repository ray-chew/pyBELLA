"""HDF5 output writer and the ensemble-output bootstrap for pyBELLA."""

import os
import shutil
import logging

import h5py
import numpy as np
import dill as pickle  # pickle jar to debug classes

from .. import sim_params as params


def initialise(sst):
    es = sst.ensemble_state
    dp = sst.da_params
    ######################################################
    # Initialise writer class for I/O operations
    ######################################################
    writer = hdf5(sst.ud, sst.restart)
    writer.check_jar()
    writer.jar([sst.ud, es[0].elem, es[0].node, dp.dap])
    # sys.exit("Let's just dill the stuff and quit!")

    writer.write_attrs()
    wrtr = None
    if sst.N > 1:
        writer.write_da_attrs(dp.dap)
    elif params.output_timesteps == True:
        wrtr = writer

    for n in range(sst.N):  # write initial ensemble
        if params.label_type == "STEP":
            label = "ensemble_mem=%i_%.3d" % (n, sst.step)
        else:
            label = "ensemble_mem=%i_%.3f" % (n, 0.0)
        if not sst.restart:
            writer.write_all(es[n], str(label) + "_ic")

    if params.da_debug:
        # writer.jar([obs,obs_noisy,obs_noisy_interp,obs_mask,obs_covar])
        # obs = obs_noisy_interp
        writer.jar([dp.obs, dp.obs_noisy, dp.obs_mask, dp.obs_covar])

    return writer, wrtr


class hdf5(object):
    """
    HDF5 writer class. Contains methods to create HDF5 file, create data sets and populate them with output variables.

    """

    def __init__(self, ud, restart=False):
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
        self.OUTPUT_FOLDER = params.output_path + "/" + self.OUTPUT_FILENAME
        self.OUTPUT_FILENAME = self.OUTPUT_FOLDER + "/" + self.ud.output_type

        self.SUFFIX = self.ud.output_suffix
        if restart:
            self.OLD_SUFFIX = self.ud.old_suffix

        self.PATHS = [  #'buoy',
            # 'dp2_c',
            # 'dp2_nodes',
            # 'dpdim',
            # 'drhoY',
            # 'dT',
            # 'dY',
            # 'p',
            # 'p2_c',
            "p2_nodes",
            "rho",
            "rhoY",
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

        self.io_create_file(self.PATHS, restart)
        self.time = "None"

    def io_create_file(self, paths, restart):
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
        if os.path.exists(
            self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT
        ):
            os.rename(
                self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT,
                self.OUTPUT_FILENAME
                + self.BASE_NAME
                + self.SUFFIX
                + "_old"
                + self.FORMAT,
            )

        # create a new output file for each rerun - old output will be overwritten.
        if restart:
            src = self.OUTPUT_FILENAME + self.BASE_NAME + self.OLD_SUFFIX + self.FORMAT
            dest = (
                self.OUTPUT_FILENAME
                + self.BASE_NAME
                + self.OLD_SUFFIX
                + "_old"
                + self.FORMAT
            )
            shutil.copy2(src, dest)
        else:
            file = h5py.File(
                self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, "a"
            )
            for path in paths:
                # check if groups have been created
                # if not created, create empty groups
                if not (path in file):
                    file.create_group(path, track_order=True)

            file.close()

    def write_all(self, model_state, name):
        """
        At a given time, write output from `Sol` and `npf` to the HDF5 file.

        Parameters
        ----------
        Sol : :class:`management.variable.Vars`
            Solution data container
        npf : :class:`physics.low_mach.npf.MPV`
            Variables relating to the elliptic solver
        elem : :class:`discretization.kgrid.ElemSpaceDiscr`
            Cells grid
        node : :class:`discretization.kgrid.NodeSpaceDiscr`
            Nodes grid
        th : :class:`physics.gas_dynamics.thermodynamic.init`
            Thermodynamic variables of the system
        name: str
            The time and additional suffix label for the dataset, e.g. "_10.0_after_full_step", where 10.0 is the time and "after_full_step" denotes when the output was made.

        """

        Sol = model_state.sol
        npf = model_state.npf

        logging.info("writing hdf output..." + name)
        # rho
        self.populate(name, "rho", Sol.rho)
        # rhoY
        self.populate(name, "rhoY", Sol.rhoY)

        # rho u ,v w
        self.populate(name, "rhou", Sol.rhou)
        self.populate(name, "rhov", Sol.rhov)
        self.populate(name, "rhow", Sol.rhow)
        self.populate(name, "rhoX", Sol.rhoX)

        self.populate(name, "p2_nodes", npf.p2_nodes)

    def populate(self, name, path, data, options=None):
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
        file = h5py.File(
            self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, "r+"
        )
        file.create_dataset(
            str(path) + "/" + str(path) + "_" + str(name),
            data=data,
            chunks=True,
            compression="gzip",
            compression_opts=4,
            dtype=np.float32,
        )
        # add attributes, i.e. the simulation parameters to each dataset.
        # for key in options:
        # file[str(path) + '/' + str(name)].attrs.create(key,options[key])
        try:
            file[str(path)][str(path) + "_" + str(name)].attrs.create("t", self.time)
        except:
            # file.attrs.create(key,repr(value),dtype='<S' + str(len(repr(value))))
            file[str(path)][str(path) + "_" + str(name)].attrs.create(
                "t", repr(self.time), dtype="<S" + str(len(repr(self.time)))
            )
        # print("writing time = %.1f for arrays %s" %(name,path))
        file.close()

    def write_attrs(self):
        """
        Method to write all attributes in the userdata initial condition to HDF5 file.
        """
        file = h5py.File(
            self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, "a"
        )
        for key, value in vars(self.ud).items():
            try:
                file.attrs.create(key, value)
            except:
                file.attrs.create(key, repr(value), dtype="<S" + str(len(repr(value))))
        file.close()

    def write_da_attrs(self, params):
        """
        Method to write all data-assimilation attributes in the userdata initial condition to HDF5 file.
        """
        file = h5py.File(
            self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, "a"
        )
        path = "da_parameters"
        if not (path in file):
            file.create_group(path, track_order=True)

        for key, value in vars(params).items():
            try:
                file["da_parameters"].attrs.create(key, value)
            except:
                # print(str(repr(value)))
                file["da_parameters"].attrs.create(
                    key, repr(value), dtype="<S" + str(len(repr(value)))
                )
        file.close()

    def close_everything(self):
        """
        In parallel, some workers do not close file correctly, this function forces the program to close the HDF file before exiting.
        """
        file = h5py.File(
            self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + self.FORMAT, "r"
        )
        if file.__bool__():
            file.close()

    def check_jar(self):
        fn = self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + ".dat"
        if os.path.exists(fn):
            os.rename(
                fn, self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + "_old.dat"
            )

    def jar(self, content=None):
        fn = self.OUTPUT_FILENAME + self.BASE_NAME + self.SUFFIX + ".dat"
        if os.path.exists(fn):
            file = open(fn, "ab")
        else:
            file = open(fn, "wb")

        # let's fill the pickle jar
        if content is not None:
            for each_pickle in content:
                pickle.dump(each_pickle, file)
        file.close()
