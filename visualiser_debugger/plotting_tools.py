import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation
import numpy as np
import itertools

class plotter(object):
    def __init__(self,arr_lst, ncols=4, figsize=(12,8), sharex=False, sharey=False):
        self.arr_lst = np.array(arr_lst)
        N = self.arr_lst.shape[0]
        
        if N > ncols:
            self.nrows = int(np.ceil(N/ncols))
            self.ncols = ncols
        if N <= ncols:
            self.nrows = 1
            self.ncols = N
        self.N = N
        
        ridx = np.arange(self.nrows)
        cidx = np.arange(self.ncols)

        self.idx = []
        if self.nrows > 1:
            for pair in itertools.product(ridx, cidx):
                self.idx.append(pair)
        else:
            self.idx = cidx
         
        self.visualise = self.visualise
        self.fig, self.ax = plt.subplots(ncols=self.ncols,nrows=self.nrows,figsize=figsize,sharex=sharex,sharey=sharey)
        
        self.img = plt
            
    def set_axes(self,**kwargs):
        for key, value in kwargs.items():
            setattr(self,key,value)
        #self.x_locs= x_locs
        #self.x_axs = x_axs
        #self.y_locs = y_locs
        #self.y_axs = y_axs
        #self.x_label = x_label
        #self.y_label = y_label
        #self.axhline = axhline
        #self.axvline = axvline
        
    def set_cax_axes(self,cax):
        if hasattr(self, 'x_locs') : cax.set_xticks(self.x_locs)
        if hasattr(self, 'x_axs') : cax.set_xticklabels(self.x_axs)
        if hasattr(self, 'y_locs') : cax.set_yticks(self.y_locs)
        if hasattr(self, 'y_axs') : cax.set_yticklabels(self.y_axs)
        if hasattr(self, 'x_label') : cax.set_xlabel(self.x_label)
        if hasattr(self, 'y_label') : cax.set_ylabel(self.y_label)
        if hasattr(self, 'axhline'): cax.axhline(self.axhline,c='k',lw=0.5)
        if hasattr(self, 'axvline'): cax.axvline(self.axvline,c='k',lw=0.5)
        
        
    def plot(self,method='imshow',inner=False,suptitle="",rect=[0, 0.03, 1, 0.95],fontsize=14,aspect='auto',lvls=None):
        plt.rcParams.update({'font.size': fontsize})
        if method != 'imshow' or method != 'contour':
            assert(0, "Visualisation method not implemented!")
            
        if self.N > 1:
            ims, caxs = [], []
            for n, arr in enumerate(self.arr_lst):
                arr, title = arr[0], arr[1]
                if inner == True:
                    arr = arr[2:-2,2:-2]
                cax = self.ax[self.idx[n]]
                
                lvl = lvls[n] if lvls is not None else None
                
                im = self.visualise(method,cax,arr,aspect,lvl)
                cax.set_title(title)
                loc = cax.get_xticklabels()
                
#                 if hasattr(self, 'x_locs') : cax.set_xticks(self.x_locs)
#                 if hasattr(self, 'x_axs') : cax.set_xticklabels(self.x_axs)
#                 if hasattr(self, 'y_locs') : cax.set_yticks(self.y_locs)
#                 if hasattr(self, 'y_axs') : cax.set_yticklabels(self.y_axs)
#                 if hasattr(self, 'x_label') : cax.set_xlabel(self.x_label)
#                 if hasattr(self, 'y_label') : cax.set_ylabel(self.y_label)
#                 if hasattr(self, 'axhline'): cax.axhline(self.axhline,c='k',lw=0.5)
#                 if hasattr(self, 'axvline'): cax.axvline(self.axvline,c='k',lw=0.5)
                self.set_cax_axes(cax)
                caxs.append(cax)
                divider = make_axes_locatable(cax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(im, cax=cax, format='%.4f')
                
                ims.append(im)
                
            for i in range(n+1,self.nrows*self.ncols):
                self.fig.delaxes(self.ax[self.idx[i]])
            
        else:
            arr, title = self.arr_lst[0][0], self.arr_lst[0][1]
            if inner == True:
                arr = arr[2:-2,2:-2]
            cax = self.fig.gca()
            im = self.visualise(method,cax,arr,aspect,lvls)
            cax.set_title(title)
#             if hasattr(self, 'x_locs') : cax.set_xticks(self.x_locs)
#             if hasattr(self, 'x_axs') : cax.set_xticklabels(self.x_axs)
#             if hasattr(self, 'y_locs') : cax.set_yticks(self.y_locs)
#             if hasattr(self, 'y_axs') : cax.set_yticklabels(self.y_axs)
#             if hasattr(self, 'x_label') : cax.set_xlabel(self.x_label)
#             if hasattr(self, 'y_label') : cax.set_ylabel(self.y_label)
#             if hasattr(self, 'axhline'):
#                 cax.set_axhline(self.axhline,c='k',lw=0.5)
#             if hasattr(self, 'axvline'):
#                 cax.axvline(self.axvline,c='k',lw=0.5)
            self.set_cax_axes(cax)
            divider = make_axes_locatable(cax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)
            ims = [im]
            
        plt.suptitle(suptitle)
        plt.tight_layout(rect=rect)
        
        return ims, caxs
        
    def save_fig(self, fn, format='.pdf'):
        self.img.savefig(fn + format, bbox_inches = 'tight', pad_inches = 0)
        
        
    @staticmethod
    def visualise(method,cax,arr,aspect,lvls):
        if method == 'imshow':
            im = cax.imshow(arr,aspect=aspect,origin='lower')
        elif method == 'contour':
            if lvls is None:
                cax.set_aspect(aspect)
                im = cax.contour(arr,colors='k')
                im = cax.contourf(arr)
                cax.set_aspect(aspect)
            else:
                cax.set_aspect(aspect)
                im = cax.contour(arr,linewidths=0.5,levels=lvls,colors='k')
                im = cax.contourf(arr,levels=lvls,extend='both')
                cax.set_aspect(aspect)
        return im
    

class animator_2D(plotter):
    def __init__(self,time_series,ncols,figsize=(16,8)):
        self.time_series = time_series
        self.frns = time_series.shape[0]
        super().__init__(self.time_series[0], ncols, figsize)
        
        self.update_plot = self.update_plot
        self.fig.tight_layout()
        self.suptitle = None
        self.method = None 
        
    def animate(self, interval=100, **kwargs):
        self.ims, self.caxs = self.plot(**kwargs)
        anim = animation.FuncAnimation(self.fig, self.update_plot, self.frns, fargs=(self.time_series, self.ims, self.caxs, self.fig, self.suptitle, self.method), interval=interval) 
        return anim

    @staticmethod
    def update_plot(frame_number, time_series, ims, caxs, img, title, method):
        
        if method == 'imshow':
            for ii,im in enumerate(ims):
                arr = time_series[frame_number][ii][0]
                im.set_array(arr)
                im.set_clim(arr.min(),arr.max()) 
        elif method == 'contour':
            for ii, cax in enumerate(caxs):
                arr = time_series[frame_number][ii][0]
                im = ims[ii]
                for c in cax.collections:
                    cax.collections.remove(c)
                for c in cax.collections:
                    cax.collections.remove(c)
                for c in cax.collections:
                    cax.collections.remove(c)
                im = cax.contourf(arr)
                im = cax.contour(arr,colors='k')
        if title is not None:
            stt = title(frame_number)
            img.suptitle(stt)
        else:
            img.suptitle(frame_number)
    
    
class plotter_1d(object):
    def __init__(self,ncols=3,nrows=2,figsize=(12,12),fontsize=16):
        plt.rcParams.update({'font.size': fontsize})
        self.fig, self.ax = plt.subplots(ncols=ncols,nrows=nrows, sharex=False, sharey=False, figsize=figsize)
        self.nrows = nrows
        self.ncols = ncols
        
        self.img = plt
        
    def set_suptitle(self,suptitle):
        self.img.suptitle(suptitle)
        
    def set_x(self,x_axs):
        self.x = x_axs
        
    def get_ax(self,i):
        
        #if not hasattr(self,'x'):
            #assert 0, "x-axis has not been set, use set_x(x_axs)."
        if self.ncols == 1 and self.nrows == 1:
            return self.ax
        elif self.nrows == 1 or self.ncols == 1:
            return self.ax[i]
        else:
            row = int(np.floor(i/self.ncols))
            col = int(i%self.ncols)
            
            return self.ax[row,col]

    def save_fig(self, fn, format='.pdf'):
        self.img.savefig(fn + format, bbox_inches = 'tight', pad_inches = 0)

    
def labels():
    labels_dict = {
        'rho'       : r'$\rho$, density',
        'rhou'      : r'$\rho u$, horizontal momentum',
        'rhov'      : r'$\rho v$, vertical momentum',
        'rhow'      : r'$\rho w$, horizontal momentum',
        'buoy'      : r'buoyancy',
        'rhoX'      : r'$\rho / \Theta$, mass-weighted inverse pot. temp.',
        'rhoY'      : r'$\rho \Theta$, mass-weighted potential temperature',
        'p2_nodes'  : r'$\pi$, nodal Exner pressure'
        }
    return labels_dict

def labels_increment():
    labels_dict = labels()
    labels_dict['p2_nodes'] = r'$\delta \pi$, nodal Exner pressure increment'
    return labels_dict

def short_labels():
    labels_dict = {
        'rho'       : r'$\rho$',
        'rhou'      : r'$\rho u$',
        'rhov'      : r'$\rho v$',
        'rhow'      : r'$\rho w$',
        'buoy'      : r'buoyancy',
        'rhoX'      : r'$\rho / \Theta$',
        'rhoY'      : r'$\rho \Theta$',
        'p2_nodes'  : r'$\pi$'
        }

    return labels_dict
