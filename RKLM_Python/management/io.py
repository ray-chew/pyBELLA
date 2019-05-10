# input/output
import h5py
import numpy as np

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

    def write_all(self,Sol,mpv,elem,node,th,name,ud):
        # rho
        self.populate(name,'rho',Sol.rho)
        # rhoe
        self.populate(name,'rhoe',Sol.rhoe)
        # rhoY
        self.populate(name,'rhoY',Sol.rhoY)

        # rho u ,v w
        self.populate(name,'rhou',Sol.rhou)
        self.populate(name,'rhov',Sol.rhov)
        self.populate(name,'rhow',Sol.rhow)

        # dp2_nodes
        self.populate(name,'dp2_nodes',mpv.dp2_nodes)
        # p2_nodes
        self.populate(name,'p2_nodes',mpv.p2_nodes)
        # p2_cells
        self.populate(name,'p2_cells',mpv.p2_cells)
        # dp2_cells
        self.populate(name,'dp2_cells',mpv.dp2_cells)

        # pressure
        self.populate(name,'p',Sol.rhoY**th.gamm)
        # pressure difference
        self.populate(name,'dpdim', self.dpress_dim(mpv,ud,th))

        # velocity (u,v,w)
        self.populate(name,'u',Sol.rhou / Sol.rho)
        self.populate(name,'v',Sol.rhov / Sol.rho)
        self.populate(name,'w',Sol.rhow / Sol.rho)

        # vorticity in (x,z)
        self.populate(name,'vortz', self.vortz(Sol,elem,node))

        # temperature
        self.populate(name,'T', Sol.rhoY**th.gamm / Sol.rho)
        # temperature difference
        self.populate(name,'dT', Sol.rhoY**th.gamm / Sol.rho - mpv.HydroState.p0[0,:] / mpv.HydroState.rho0[0,:])

        # species mass fraction(?)
        self.populate(name,'Y', Sol.rhoY / Sol.rho)
        # species mass fraction perturbation
        self.populate(name,'dY', Sol.rhoY / Sol.rho - ud.stratification(elem.y))

    def vortz(self,Sol,elem,node):
        # 2d-case
        igs = elem.igs
        dx = elem.dx
        dy = elem.dy
        inner_domain_rho = Sol.rho[igs[0]-1:-igs[0]+1, igs[1]-1:-igs[1]+1]
        inner_domain_rhou = Sol.rho[igs[0]-1:-igs[0]+1, igs[1]-1:-igs[1]+1]
        inner_domain_rhov = Sol.rho[igs[0]-1:-igs[0]+1, igs[1]-1:-igs[1]+1]

        top_left_idx     = (slice(0,-1)    , slice(0,-1))
        top_right_idx    = (slice(0,-1)    , slice(1, None))
        bottom_left_idx  = (slice(1, None) , slice(0,-1))
        bottom_right_idx = (slice(1, None) , slice(1,None))

        dvdx = 0.5 * ((inner_domain_rhov[bottom_right_idx] / inner_domain_rho[bottom_right_idx] - inner_domain_rhov[bottom_left_idx] / inner_domain_rho[bottom_left_idx]) + (inner_domain_rhov[top_right_idx] / inner_domain_rho[top_right_idx] - inner_domain_rhov[top_left_idx] / inner_domain_rho[top_left_idx])) / dx

        dudy = 0.5 * ((inner_domain_rhou[bottom_right_idx] / inner_domain_rho[bottom_right_idx] - inner_domain_rhou[bottom_left_idx] / inner_domain_rho[bottom_left_idx]) + inner_domain_rhou[top_right_idx] / inner_domain_rho[top_right_idx] - inner_domain_rhou[top_left_idx] / inner_domain_rho[top_left_idx]) / dy

        vortz = np.zeros((node.sc)).squeeze()

        vortz[igs[0]:-igs[0], igs[1]:-igs[1]] = dvdx - dudy
        return vortz

    def dpress_dim(self,mpv,ud,th):
        p = (th.Gamma * ud.Msq * mpv.p2_cells)**th.Gammainv
        return (p - mpv.HydroState.p0[0,:]) * ud.p_ref



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