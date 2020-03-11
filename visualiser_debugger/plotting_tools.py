import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import itertools

class plotter(object):
    def __init__(self,arr_lst):
        self.arr_lst = np.array(arr_lst)
        N = self.arr_lst.shape[0]
        
        if N > 4:
            self.nrows = int(int(N/4)+1)
            self.ncols = 4
        if N <= 4:
            self.nrows = 1
            self.ncols = N
        self.N = N
            
        ridx = np.arange(self.nrows)
        cidx = np.arange(self.ncols)

        self.idx = []
        for pair in itertools.product(ridx, cidx):
            self.idx.append(pair)
            
        self.visualise = self.visualise
    
            
    def plot(self,method='imshow',figsize=[None]):
        if method != 'imshow' or method != 'contour':
            assert(0, "Visualisation method not implemented!")
            
        if figsize[0] != None:
            figsize = figsize
        else:
            figsize = (12,8)
        fig, ax = plt.subplots(ncols=self.ncols,nrows=self.nrows,figsize=figsize)
        
        if self.N > 1:
            for n, arr in enumerate(self.arr_lst):
                cax = ax[self.idx[n]] 
                im = self.visualise(method,cax,arr)
                divider = make_axes_locatable(cax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(im, cax=cax)
                
            for i in range(n+1,self.nrows*self.ncols):
                fig.delaxes(ax[self.idx[i]])
        else:
            cax = fig.gca()
            im = self.visualise(method,cax,self.arr_lst[0])
            divider = make_axes_locatable(cax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)
        
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
            
            
            
            
        
