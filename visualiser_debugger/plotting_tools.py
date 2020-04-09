import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import itertools

class plotter(object):
    def __init__(self,arr_lst, ncols=4):
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
            
    def plot(self,method='imshow',figsize=[None],inner=False,suptitle=""):
        if method != 'imshow' or method != 'contour':
            assert(0, "Visualisation method not implemented!")
            
        if figsize[0] != None:
            figsize = figsize
        else:
            figsize = (12,8)
        fig, ax = plt.subplots(ncols=self.ncols,nrows=self.nrows,figsize=figsize)
        
        if self.N > 1:
            for n, arr in enumerate(self.arr_lst):
                arr, title = arr[0], arr[1]
                if inner == True:
                    arr = arr[2:-2,2:-2]
                cax = ax[self.idx[n]]
                im = self.visualise(method,cax,arr)
                cax.set_title(title)
                divider = make_axes_locatable(cax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(im, cax=cax)
                
            for i in range(n+1,self.nrows*self.ncols):
                fig.delaxes(ax[self.idx[i]])
        else:
            arr, title = self.arr_lst[0][0], self.arr_lst[0][1]
            if inner == True:
                arr = arr[2:-2,2:-2]
            cax = fig.gca()
            im = self.visualise(method,cax,arr)
            cax.set_title(title)
            divider = make_axes_locatable(cax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)
            
        plt.suptitle(suptitle)
        plt.tight_layout()
        plt.show()
        
        
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
    def __init__(self,ncols=2,nrows=3,figsize=(12,12)):
        self.fig, self.ax = plt.subplots(ncols=2,nrows=3, figsize=(12,12))
        
    def set_x(self,x_axs):
        self.x = x_axs
        
    def get_ax(self,i):
        
        #if not hasattr(self,'x'):
            #assert 0, "x-axis has not been set, use set_x(x_axs)."
            
        row = int(np.floor(i/2))
        col = int(i%2)
        
        return self.ax[row,col]
