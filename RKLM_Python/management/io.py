# input/output
import h5py

class io(object):
    def __init__(self):
        self.FORMAT = ".h5"
        self.OUTPUT_FILENAME = "output"
        self.SUFFIX = "_low_mach_gravity_comp"

        self.PATHS = [   'bouy',
                        'dp2_c',
                        'dp2_nodes',
                        'dpdim',
                        'drhoY',
                        'dT',
                        'dY',
                        'p',
                        'p2_c',
                        'p2_nodes',
                        'rho',
                        'rhoe',
                        'rhoY',
                        'S',
                        'T',
                        'u',
                        'v',
                        'w',
                        'vortz',
                        'Y'
                    ]
        
        self.io_create_file(self.PATHS)

    def io_create_file(self,paths):
        # create a new output file for each rerun - old output will be overwritten.
        file = h5py.File(self.OUTPUT_FILENAME + self.SUFFIX + self.FORMAT, 'w')
        for path in paths:
            # check if groups have been created
            # if not created, create empty groups
            if not (path in file):
                file.create_group(path,track_order=True)
        file.close()

    def write_all(self,Sol,mpv,elem,node,th,name):
        # rho
        self.populate(name,'rho',Sol.rho)
        # rhoe
        self.populate(name,'rhoe',Sol.rhoe)
        # rhoY
        self.populate(name,'rhoY',Sol.rhoY)
        # dp2_nodes
        self.populate(name,'dp2_nodes',mpv.dp2_nodes)
        # p2_nodes
        self.populate(name,'p2_nodes',mpv.p2_nodes)
        # p2_cells
        self.populate(name,'p2_cells',mpv.p2_cells)

        # pressure
        self.populate(name,'p',Sol.rhoY**th.gamm)

        # velocity (u)
        self.populate(name,'u',Sol.rhou / Sol.rho)

    def populate(self,name,path,data,options=None):
        # name is the simulation time of the output array
        # path is the array type, e.g. U,V,H, and data is it's data.
        file = h5py.File(self.OUTPUT_FILENAME + self.SUFFIX + self.FORMAT, 'r+')
        file.create_dataset(str(path) + '/' + str(path) + '_' + str(name), data=data, chunks=True, compression='gzip', compression_opts=4)
        # add attributes, i.e. the simulation parameters to each dataset.
        # for key in options:
        #     file[str(path) + '/' + str(name)].attrs.create(key,options[key])
        # print("writing time = %.1f for arrays %s" %(name,path))
        file.close()