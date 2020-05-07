import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import itertools

class plotter(object):
    def __init__(self,arr_lst, ncols=4, figsize=(12,8)):
        self.arr_lst = np.array(arr_lst)
        N = self.arr_lst.shape[0]
        
        if N > ncols:
            self.nrows = int(int(N/ncols)+1)
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
        self.fig, self.ax = plt.subplots(ncols=self.ncols,nrows=self.nrows,figsize=figsize)
        
        self.img = plt
            
    def set_axes(self,x_locs=None,x_axs=None,y_locs=None,y_axs=None):
       self.x_locs= x_locs
       self.x_axs = x_axs
       self.y_locs = y_locs
       self.y_axs = y_axs
       
    def plot(self,method='imshow',inner=False,suptitle="",rect=[0, 0.03, 1, 0.95],fontsize=14):
        plt.rcParams.update({'font.size': fontsize})
        if method != 'imshow' or method != 'contour':
            assert(0, "Visualisation method not implemented!")
            
        if self.N > 1:
            for n, arr in enumerate(self.arr_lst):
                arr, title = arr[0], arr[1]
                if inner == True:
                    arr = arr[2:-2,2:-2]
                cax = self.ax[self.idx[n]]
                
                im = self.visualise(method,cax,arr)
                cax.set_title(title)
                loc = cax.get_xticklabels()
                
                if hasattr(self, 'x_locs') : cax.set_xticks(self.x_locs)
                if hasattr(self, 'x_axs') : cax.set_xticklabels(self.x_axs)
                if hasattr(self, 'y_locs') : cax.set_yticks(self.y_locs)
                if hasattr(self, 'y_axs') : cax.set_yticklabels(self.y_axs)
                
                divider = make_axes_locatable(cax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(im, cax=cax)
                
            for i in range(n+1,self.nrows*self.ncols):
                self.fig.delaxes(self.ax[self.idx[i]])
        else:
            arr, title = self.arr_lst[0][0], self.arr_lst[0][1]
            if inner == True:
                arr = arr[2:-2,2:-2]
            cax = self.fig.gca()
            im = self.visualise(method,cax,arr)
            cax.set_title(title)
            divider = make_axes_locatable(cax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)
            
        plt.suptitle(suptitle)
        plt.tight_layout(rect=rect)
        
    def save_fig(self, fn, format='.pdf'):
        self.img.savefig(fn + format, bbox_inches = 'tight', pad_inches = 0)
        
        
    @staticmethod
    def visualise(method,cax,arr,lvls=[None]):
        if method == 'imshow':
            im = cax.imshow(arr,origin='lower')
        elif method == 'contour':
            if lvls[0] == None:
                im = cax.contour(arr)
            else:
                im = cax.contour(arr,levels=lvls)
        return im
    
class plotter_1d(object):
    def __init__(self,ncols=3,nrows=2,figsize=(12,12),fontsize=16):
        plt.rcParams.update({'font.size': fontsize})
        self.fig, self.ax = plt.subplots(ncols=ncols,nrows=nrows, figsize=figsize)
        self.nrows = nrows
        self.ncols = ncols
        
        self.img = plt
        
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
