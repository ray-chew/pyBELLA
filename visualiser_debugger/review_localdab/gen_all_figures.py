import sys
# Notebook needs to see the RKLM_Python module to load pickled class instances
sys.path.append('../../RKLM_Python/')
# Path to output files
sys.path.append('../')

import numpy as np
import utils
import plotting_tools as pt

from importlib import reload
import utils
utils = reload(utils)
pt = reload(pt)

from scipy import signal
from numpy import linalg as la
import pickle

cm = 1/2.54
dt = './output/tmp'
dd = './output'

def create_output_dir():
    import os
    os.chdir(os.getcwd())
    # If directory does not exist, create it.
    if not os.path.exists('output'):
        os.mkdir('output')
    if not os.path.exists('output/tmp'):
        os.mkdir('output/tmp')
        
create_output_dir()

#######################################
#
# Helper function to merge figures
#
#######################################
def merge_plots(ps, out_name, shape):
    import subprocess
    import os

    os.chdir(os.getcwd())
    
    for p in ps:
        proc0 = subprocess.Popen(['pdfcrop', '%s' %(p), '%s' %(p)], stdout=subprocess.PIPE)
        proc0a = subprocess.check_output(('grep', '-v', 'Heiko'), stdin=proc0.stdout)
        proc0.wait()

    call = ['pdfjam', '-q', '--nup', shape]
    for p in ps:
        call += [p]
    call += ['--outfile', '%s/%s.pdf' %(dd,out_name)]
    call += ['--delta', '0.2cm 0.2cm']

    proc1 = subprocess.call(call)
    if not proc1:
        proc2 = subprocess.Popen(['pdfcrop', '%s/%s.pdf' %(dd,out_name), '%s/%s.pdf' %(dd,out_name)], stdout=subprocess.PIPE)
        proc2a = subprocess.check_output(('grep', '-v', 'Heiko'), stdin=proc2.stdout)
        proc2.wait()

        
#######################################
#
# Helper function to generate figures
# 4 and 6
#
#######################################

def figures_4and6(switches, attribute, output_end_time):
    euler = switches[0]
    rb = switches[1]
    end = output_end_time
    rhou = attribute[0]
    vortz = attribute[1]
    pot_temp = attribute[2]
    exner = attribute[3]
    
    # initial time-step
    times = [0]
    tag = 'ic'
    if rb and end:
        times = [163]
    #     times = [215]
        tag = 'after_full_step'

    N = 1 # ensemble size
    l_typ = 'WINDOW_STEP'

    def ic_loader(tc,sfx1, swe):
        # load pickled instances of data used in simulation
        fn_pickle = tc.get_filename(N,sfx1,format='dat')
        path_pickle = tc.get_path(fn_pickle)

        file = open(path_pickle,'rb')
        ud = pickle.load(file)
        file.close()
        return ud

    def get_ens(tc, sfx , diff, attribute, swe=False, rhou=False):
        ens = tc.get_ensemble(times, N, attribute, sfx, label_type=l_typ, avg=True, diff=diff, tag=tag)[1]
        ud = ic_loader(tc,sfx,swe)
        if euler:
            if attribute == 'p2_nodes': # convert to pressure p
                ens -= ens.mean() # pi is correct up to a mean, dimensionless
                ens *= ud.Msq # pi is the perturbation variable
                ens *= 1e4 # scale colourbar
            if attribute == 'rhou':
                ens = ens * ud.rho_ref * ud.u_ref # units in kg m^{-1} s^{-1}
            if attribute == 'vortz':
                ens = ens # units in 0.01 s^{-1}
            if attribute == 'rhoY':
                rho = tc.get_ensemble(times, N, 'rho', sfx, label_type=l_typ, avg=True, diff=diff, tag=tag)[1]
                ens /= rho # units in 300 K

        label = sfx + '_' + attribute
        return label, ens.T

    diff = False

    if euler:
        if rhou:
            attribute = 'rhou'
        elif vortz:
            attribute = 'vortz'
        elif pot_temp:
            attribute = 'rhoY'
        elif exner:
            attribute = 'p2_nodes'
        et = 1.0
        base_fn = "output_travelling_vortex"
        directory = "output_travelling_vortex"
        py_directory = "../../%s/" %directory

        Nx, Ny = 64, 64
        tc = utils.test_case(base_fn,py_directory,Nx,Ny,et)

        sfx1 = 'comp_bal_noib'
        l2, a2 = get_ens(tc, sfx1, diff, attribute, rhou=rhou)
    #     la = 'Euler vortex'
        la = r'$\rho u$' if rhou else r'$\pi$'
        la = ''
        aa = a2

        lbl = 'euler'
        if rhou:
            lvls = np.arange(50,120,5)
            lbl += '_rhou'
    #         clbl = r'kg m$^{-1}$ s$^{-1}$'
            clbl = '   '
        elif vortz:
            lvls = np.arange(-0.022,0.038,0.010) * 100
            lbl += '_vortz'
            clbl = r'$\times 10^{-2}$ s$^{-1}$'
            clbl = r'$\times 10^{-2}$'
        elif pot_temp:
            lvls = np.arange(1.10,2.0,0.10) #/ 300.0
            lbl += '_pt'
            clbl = r'$\times 300$'
    #         lvls = None
        elif exner:
            lvls = np.arange(-0.0005,0.0001, 0.0001) * 1e4
            lbl += '_exner'
            clbl = r'$\times 10^{-4}$'

    elif rb:
        et = 1.0
        base_fn = "output_rising_bubble"
        directory = "output_rising_bubble"
        py_directory = "../../%s/" %directory

        Nx, Ny = 160, 80
        tc = utils.test_case(base_fn,py_directory,Nx,Ny,et)

        sfx1 = 'comp_imbal_half_CFLfixed_ib-0'
    #     sfx1 = 'truth_CFLfixed_ib-0'
        _, a1 = get_ens(tc, sfx1, diff, 'rho')
        _, a2 = get_ens(tc, sfx1, diff, 'rhoY')
        la = 'Rising bubble'
        la = ''
        aa = a2/a1 * 300
        if end:
            lvls = np.linspace(299.75,302.0,10)
            lvls = np.linspace(300.25,302.0,8)
    #         lvls = np.linspace(299.5,302.0,11)
    #         lvls = np.arange(300.25,301.75,.25)
            lvls = np.arange(300.25,302.0,.25)
        else:
    #         lvls = np.linspace(300.0,302.0,11)
            lvls = np.linspace(299.75,302.0,10)
            lvls = np.linspace(300.25,302.0,8)
            lvls = np.arange(300.25,302,0.25)
        lbl = 'rb'
        clbl = 'K'

    ll = [aa, la]
    pl_lst = [ll]

    if euler:
        x_axs = [-0.5,-0.25,0.0,0.25,0.5]
        y_axs = [-0.5,-0.25,0.0,0.25,0.5]
        pl = pt.plotter(pl_lst,ncols=1,figsize=(5,5),sharey=False)
    else:
        x_axs = [-1.0,-0.5,0.0,0.5,1.0]
        y_axs = [0.0,0.25,0.5,0.75,1.0]
        pl = pt.plotter(pl_lst,ncols=1,figsize=(8,4),sharey=False)

    x_label = r'x [$\times 10$ km]'
    y_label = r'z [$\times 10$ km]'
    if exner or vortz:
        x_loc = np.linspace(0,Nx,5)
        y_loc = np.linspace(0,Ny,5)
    else:
        x_loc = np.linspace(0,Nx-1,5)
        y_loc = np.linspace(0,Ny-1,5)
    pl.set_axes(x_locs=x_loc, y_locs=y_loc, x_axs=x_axs, y_axs=y_axs, x_label=x_label, y_label=y_label,cbar_label=clbl)
    _ = pl.plot(aspect='equal',method='contour',lvls=lvls)

    out_t = 'final' if end else 'initial'
    # pl.fig.tight_layout()
    pl.save_fig('./output/tmp/%s_%s' %(out_t,lbl))
    
    
#######################################
#
# Helper function to generate figure 5
#
#######################################
    
def figure_5():
    # Top panel in figure 5
    attributes = ['p2_nodes']

    et = 1.0
    N = 1

    base_fn = "output_travelling_vortex"
    directory = "output_travelling_vortex"
    # path to output directory
    py_directory = "../../%s/" %directory

    Nx, Ny = 64, 64
    tc = utils.test_case(base_fn,py_directory,Nx,Ny,et)

    times = np.arange(0,176)
    t_axs = np.linspace(0.0,1.0,5)
    t_pos = np.linspace(0,175,5)
    probe_loc = [32,32]
    l_typ = 'WINDOW_STEP'

    attr_labels = pt.labels_increment()

    sfx1 = 'psinc_noib'
    sfx2a = 'comp_imbal_noib'
    sfx2b = 'comp_imbal_full_ib-16'
    sfx2c = 'comp_imbal_half_ib-0'
    sfx3 = 'comp_imbal_noib'

    # load pickled instances of data used in simulation
    fn_pickle = tc.get_filename(N,sfx1,format='dat')
    path_pickle = tc.get_path(fn_pickle)

    file = open(path_pickle,'rb')
    ud = pickle.load(file)
    file.close()

    p_ref = ud.p_ref
    Msq = ud.Msq


    def p_converter(ens,rhoY,ud):
        dp2n = np.array([ (at_t - at_t.mean()) * ud.Msq for at_t in ens ])
        kernel = np.ones((2,2))
        dp2c = np.array([signal.fftconvolve(item, kernel, mode='valid') / kernel.sum() for item in dp2n])

        P0 = (rhoY**(ud.gamm-1.0) - dp2c)**(1.0/(ud.gamm-1.0))
        p = rhoY**(ud.gamm) - P0**(ud.gamm)
        p *= ud.p_ref
        return p

    def get_ens(sfx , diff, attribute):
        ens = tc.get_ensemble(times, N, attribute, sfx, label_type=l_typ, avg=True, diff=diff)
        rhoY = tc.get_ensemble(times, N, 'rhoY', sfx, label_type=l_typ, avg=True, diff=diff)
        
        ens = ens[1:]
        rhoY = rhoY[1:]
        p = p_converter(ens,rhoY,ud)

        probe = p[:,probe_loc[0],probe_loc[1]]
        probe = probe[1:] - probe[:-1]

        label = sfx + '_' + attribute
        return probe

    pls = []

    for i,attribute in enumerate(attributes):
        diff = False
    #     diff = True if attribute == 'p2_nodes' else False

        p1 = get_ens(sfx1, diff, attribute)
        p2b = get_ens(sfx2b, diff, attribute)
        p2c = get_ens(sfx2c, diff, attribute)
        p3 = get_ens(sfx3, diff, attribute)

        fs = (8,4) # fs used in draft
        pl = pt.plotter_1d(figsize=fs,fontsize=14,ncols=1,nrows=1)
        ax = pl.get_ax(i)

        times = times[1:]

        ref = 'psinc'
        full = 'comp'
        ic = 'imbal.' if 'imbal' in sfx3 else 'bal.'
        l1 = '%s, imbal. IC' %(ref)
        l2c = r'%s, imbal. IC, w/ blending ($\pi_{half}$)' %full
        l2b = r'%s, imbal. IC, w/ blending ($\pi_{full}$)' %full
        l3 = '%s, %s IC w/o blending' %(full, ic)
        ic = ic[:-1]

        ax.plot(times, p1, 'k', label=l1)
    #     ax.plot(times, p2b, 'C0o', ms=4, markevery=4, label=l2b)
        ax.plot(times, p2c, 'C1o', ms=4, markevery=4, label=l2c)
        ax.plot(times, p3, 'C0', label=l3)

        ax.set_xlim([times[0],times[-1]])
        ax.set_xlabel(r'time [$\times 100$ s]')
        ax.set_ylabel(r'$\delta p^\prime$ [Pa]')

        ax.set_xticks(t_pos)
        ax.set_xticklabels(t_axs)
        ax.grid()
        ax.legend()
        pl.img.tight_layout()

        fn = 'euler'
        pl.save_fig('./output/tmp/%s_w_%s' %(fn, ic))
        
    # Bottom panel in figure 5
    attributes = ['p2_nodes']

    et = 1.0
    N = 1
 
    base_fn = "output_travelling_vortex"
    directory = "output_travelling_vortex"
    # path to output directory
    py_directory = "../../%s/" %directory

    Nx, Ny = 64, 64
    tc = utils.test_case(base_fn,py_directory,Nx,Ny,et)

    times = np.arange(0,176)
    t_axs = np.linspace(0.0,1.0,5)
    t_pos = np.linspace(0,175,5)
    probe_loc = [32,32]
    l_typ = 'WINDOW_STEP'

    attr_labels = pt.labels_increment()

    sfx1 = 'psinc_noib'
    sfx2a = 'comp_imbal_noib'
    sfx2b = 'comp_imbal_full_ib-16'
    sfx2c = 'comp_imbal_half_ib-0'
    sfx3 = 'comp_bal_noib'

    # load pickled instances of data used in simulation
    fn_pickle = tc.get_filename(N,sfx1,format='dat')
    path_pickle = tc.get_path(fn_pickle)

    file = open(path_pickle,'rb')
    ud = pickle.load(file)
    file.close()

    p_ref = ud.p_ref
    Msq = ud.Msq

    def get_ens(sfx , diff, attribute):
        ens = tc.get_ensemble(times, N, attribute, sfx, label_type=l_typ, avg=True, diff=diff)
        rhoY = tc.get_ensemble(times, N, 'rhoY', sfx, label_type=l_typ, avg=True, diff=diff)
    
        ens = ens[1:]
        rhoY = rhoY[1:]
        p = p_converter(ens,rhoY,ud)

        probe = p[:,probe_loc[0],probe_loc[1]]
        probe = probe[1:] - probe[:-1]

        label = sfx + '_' + attribute
        return probe

    pls = []

    for i,attribute in enumerate(attributes):
        diff = False

        p1 = get_ens(sfx1, diff, attribute)
        p2b = get_ens(sfx2b, diff, attribute)
        p2c = get_ens(sfx2c, diff, attribute)
        p3 = get_ens(sfx3, diff, attribute)

        fs = (8,4) # fs used in draft
        pl = pt.plotter_1d(figsize=fs,fontsize=14,ncols=1,nrows=1)
        ax = pl.get_ax(i)

        times = times[1:]

        ref = 'psinc'
        full = 'comp'
        ic = 'imbal.' if 'imbal' in sfx3 else 'bal.'
        l1 = '%s, imbal. IC' %(ref)
        l2c = r'%s, imbal. IC, w/ blending ($\pi_{half}$)' %full
        l2b = r'%s, imbal. IC, w/ blending ($\pi_{full}$)' %full
        l3 = '%s, %s IC w/o blending' %(full, ic)
        ic = ic[:-1]

        ax.plot(times, p2c, 'C1o', ms=4, markevery=4, label=l2c)
        ax.plot(times, p2b, 'C4o', ms=4, markevery=4, label=l2b)
        ax.plot(times, p3, 'C2', label=l3)

        ax.set_xlim([times[0],times[-1]])
        ax.set_xlabel(r'time [$\times 100$ s]')
        ax.set_ylabel(r'$\delta p^\prime$ [Pa]')

        ax.set_xticks(t_pos)
        ax.set_xticklabels(t_axs)
        ax.grid()
        ax.legend()
        pl.img.tight_layout()

        fn = 'euler'
        pl.save_fig('./output/tmp/%s_w_%s' %(fn, ic))

        print("figure 5b: errors in the vortex blended solution...")
        print("probe_loc = (%i, %i)" %(probe_loc[0],probe_loc[1]))
        print("blending (half) - psinc: %.6f" %(la.norm(p2c-p1)/la.norm(p1)))
        print("blending (full) - psinc: %.6f" %(la.norm(p2b-p1)/la.norm(p1))) 
        

#######################################
#
# Helper function to generate figure 7
#
#######################################

def figure_7():
    # Generate the small time-step plots
    attributes = ['p2_nodes']

    base_fn = "output_rising_bubble"
    directory = "output_rising_bubble"

    py_directory = "../../%s/" %directory

    et = 1.0
    N = 1
    Nx, Ny = 160, 80
    tc = utils.test_case(base_fn,py_directory,Nx,Ny,et)
    times = np.arange(0,185)[1:]

    probe_loc = [20,40]
    l_typ = 'WINDOW_STEP'

    rhoY = False
    CFL = False

    attr_labels = pt.labels_increment()

    def ic_loader(tc,sfx1, swe):
        # load pickled instances of data used in simulation
        fn_pickle = tc.get_filename(N,sfx1,format='dat')
        path_pickle = tc.get_path(fn_pickle)

        file = open(path_pickle,'rb')
        ud = pickle.load(file)
        file.close()

        return ud

    def p_converter(ens,rhoY,ud):
        dp2n = np.array([ (at_t - at_t.mean()) * ud.Msq for at_t in ens ])
        kernel = np.ones((2,2))
        dp2c = np.array([signal.fftconvolve(item, kernel, mode='valid') / kernel.sum() for item in dp2n])

        P0 = (rhoY**(ud.gamm-1.0) - dp2c)**(1.0/(ud.gamm-1.0))
        p = rhoY**(ud.gamm) - P0**(ud.gamm)
        p *= ud.p_ref
        return p

    def get_ens(sfx , diff, attribute):
        ens = tc.get_ensemble(times, N, attribute, sfx, label_type=l_typ, avg=True, diff=False)
        rhoY = tc.get_ensemble(times, N, 'rhoY', sfx, label_type=l_typ, avg=True, diff=False)

        tt = tc.t_arr
        ens = ens[1:]
        rhoY = rhoY[1:]
        tt = tt[1:]
        ud = ic_loader(tc,sfx,False)

        p = p_converter(ens,rhoY,ud)

        probe = p[:,probe_loc[0],probe_loc[1]]
        probe = probe[1:] - probe[:-1]

        label = sfx + '_' + attribute
        return label, probe, tt

    pls = []
    ss = 'idt=2'

    for i,attribute in enumerate(attributes):
        diff = True if attribute == 'p2_nodes' or attribute == 'p2_half' else False

        sfx1 = 'psinc_noib'
        if CFL : sfx1 = 'psinc_noib_CFLfixed'
        attribute1 = 'p2_nodes'
        if rhoY : attribute1 = 'rhoY'
        l1, p1, tt = get_ens(sfx1, diff, attribute1)
        l1 = 'psinc'

        sfx2 = 'comp_imbal_noib'
        if CFL : sfx2 = 'comp_imbal_noib_CFLfixed'
        attribute2 = 'p2_nodes'
        if rhoY : attribute2 = 'rhoY'
        l2_attr = np.copy(attribute2)
        l2, p2, tt = get_ens(sfx2, diff, attribute2)
        l2 = 'comp w/o blending'

        sfx3 = 'comp_imbal_half_ib-0'
        if CFL : sfx3 = 'comp_imbal_CFLfixed_ib-0'
        attribute3 = 'p2_nodes'
        if rhoY : attribute3 = 'rhoY'
        l3, p3, tt = get_ens(sfx3, diff, attribute3)
        l3 = 'comp w/ blending'

        ########################################
        #
        ########################################
#         print("====================")

        fs = (20*cm,10*cm) # used in the publication
        pl = pt.plotter_1d(figsize=fs,fontsize=14,ncols=1,nrows=1)
        ax = pl.get_ax(i)

        if diff == True:
            tt = tt[1:]

        ax.plot(tt, p1, 'k', ms=4, label=l1)
        ax.plot(tt, p2, 'C0', label=l2)
        ax.plot(tt, p3, 'C1o', ms=6, markevery=8, label=l3)

        ax.set_ylabel(r'$\delta p^\prime$ [Pa]')
        ax.set_xlabel(r'time [$\times 1000$ s]')
        ax.set_xlim([0.0,0.35*et])
        ax.grid()
        ax.legend()

        pl.img.tight_layout()

        fn = 'rb_%i_%i_w_imbal' %(probe_loc[0],probe_loc[1])
        pl.save_fig('./output/tmp/%s' %fn)

#         print("probe_loc = (%i, %i)" %(probe_loc[0],probe_loc[1]))
#         print("w/o blending - psinc: %.6f" %(la.norm(p2-p1)))
#         print("blending - psinc: %.6f" %(la.norm(p3-p1)))


        ########################################
        #
        ########################################
        print("====================")

        _, p1, _ = get_ens(sfx1, diff, attribute1)
        _, p2, _ = get_ens(sfx2, diff, attribute2)
        _, p3, _ = get_ens(sfx3, diff, attribute3)

        pl = pt.plotter_1d(figsize=fs,fontsize=14,ncols=1,nrows=1)
        ax = pl.get_ax(i)

        dist = la.norm(p3-p1)
    #     l3a = 'dist=%.4f' %dist
        l3a = ''
        l3a += r'$\delta_{\rm{i}}$=%.4f' %((la.norm(p2-p1))/la.norm(p1))
        l3a += '\n'
        l3a += r'$\delta_{\rm{b}}$=%.4f' %((la.norm(p3-p1))/la.norm(p1))
    #     ax.annotate(l3a, (0.025,0.0031))

        ax.plot(tt, p1, 'k', ms=4, label=l1)
        ax.plot(tt, p3, 'C1', ms=6, markevery=1, label=l3)

        ax.set_ylabel(r'$\delta p^\prime$ [Pa]')
        ax.set_xlabel(r'time [$\times 1000$ s]')
        ax.set_xlim([0.0,0.35])
        ax.grid()
        ax.legend()
    #     ax.text(0.025,0.0033,l3a)

        pl.img.tight_layout()

        fn = 'rb_%i_%i_w_ib' %(probe_loc[0],probe_loc[1])
        pl.save_fig('./output/tmp/%s' %fn)

        print("figure 7: small time-step error at...")
        print("probe_loc = (%i, %i)" %(probe_loc[0],probe_loc[1]))
        print("blending - psinc: %.6f" %(la.norm(p3-p1)))

        print(l3a)


        ########################################
        #
        ########################################
        print("====================")

        probe_loc = [80,40]
        _, p1, _ = get_ens(sfx1, diff, attribute1)
        _, p2, _ = get_ens(sfx2, diff, attribute2)
        _, p3, _ = get_ens(sfx3, diff, attribute3)

        pl = pt.plotter_1d(figsize=fs,fontsize=14,ncols=1,nrows=1)
        ax = pl.get_ax(i)

        dist = la.norm(p3-p1)
    #     l3b = 'dist=%.4f' %dist
        l3b = ''
        l3b += r'$\delta_{\rm{i}}$=%.4f' %((la.norm(p2-p1))/la.norm(p1))
        l3b += '\n'
        l3b += r'$\delta_{\rm{b}}$=%.4f' %((la.norm(p3-p1))/la.norm(p1))
    #     ax.annotate(l3b, (0.025,0.05))
    #     ax.annotate(l3b, (0.82,-0.225))

        ax.plot(tt, p1, 'k', ms=4, label=l1)
        ax.plot(tt, p3, 'C3', ms=6, markevery=1, label=l3)

        ax.set_ylabel(r'$\delta p^\prime$ [Pa]')
        ax.set_xlabel(r'time [$\times 1000$ s]')
        ax.set_xlim([0.0,0.35])
        ax.grid()
        ax.legend()
    #     ax.text(0.025,0.05,l3b)

        pl.img.tight_layout()

        fn = 'rb_%i_%i_w_ib' %(probe_loc[0],probe_loc[1])
        pl.save_fig('./output/tmp/%s' %fn)

        print("figure 7: small time-step error at...")
        print("probe_loc = (%i, %i)" %(probe_loc[0],probe_loc[1]))
        print("blending - psinc: %.6f" %(la.norm(p3-p1)))

        print(l3b)

        
    # Generate contour plot in figure 7
    attribute = 'p2_nodes'
    et = 1.0
    N = 1

    base_fn = "output_rising_bubble"
    directory = "output_rising_bubble"
    py_directory = "../../%s/" %directory

    Nx, Ny = 160, 80
    euler_tc = utils.test_case(base_fn,py_directory,Nx,Ny,et)
    tags = euler_tc.get_tag_dict()

    times = [14,15]
    l_typ = 'WINDOW_STEP'

    prt = utils.prt_time(debug=False)
    attr_labels = pt.labels_increment()

    def get_ens(tc, sfx , diff, attribute, swe=False):
        ens = tc.get_ensemble(times, N, attribute, sfx, label_type=l_typ, avg=True, diff=diff, tag=tags[9])
        rhoY = tc.get_ensemble(times, N, 'rhoY', sfx, label_type=l_typ, avg=True, diff=diff, tag=tags[9])
        ens = ens[1:]
        rhoY = rhoY[1:]
        ud = ic_loader(tc,sfx,swe)

        ens = p_converter(ens, rhoY, ud)
        ens = ens[1] - ens[0]

        label = sfx + '_' + attribute
        return label, ens.T

    diff = False
    sfx1 = 'comp_imbal_noib'
    attribute = 'p2_nodes'
    l2, a2 = get_ens(euler_tc, sfx1, diff, attribute)

#     print(a2.min(),a2.max())
    
    # pl_lst = [[a1/a2 * 300.0, '']]
    pl_lst = [[a2, '']]

    pl = pt.plotter(pl_lst,ncols=2,figsize=(30*cm,10*cm),sharey=False,fontsize=13)
    x_axs = [-10.0,-5.0,0.0,5.0,10.0]
    y_axs = [0.0,5.0,10.0]
    if attribute is not 'p2_nodes':
        x_loc = np.linspace(0,Nx-1,5)
        y_loc = np.linspace(0,Ny-1,3)
    else:
        x_loc = np.linspace(0,Nx-1,5)
        y_loc = np.linspace(0,Ny-1,3)  
    x_label = r'x [km]'
    y_label = r'y [km]'
    marker = [(20,40,'C1'),(80,40,'C3')]
    # marker = [(20,40,'r')]
    lvls = np.arange(-1.5,2.5,0.5)
    lvls = np.arange(-2.8,3.2,0.8)
    # clbl = 'Pa'

    pl.set_axes(x_locs=x_loc, y_locs=y_loc, x_axs=x_axs, y_axs=y_axs,marker=marker, x_label=x_label, y_label=y_label)
    _ = pl.plot(aspect='equal',method='contour',lvls=lvls)
    pl.save_fig('./output/tmp/rb_deltap_contour')

    # Generate plots with CFL-determined time-steps
    attributes = ['p2_nodes']

    base_fn = "output_rising_bubble"
    directory = "output_rising_bubble"

    py_directory = "../../%s/" %directory

    et = 1.0
    N = 1
    Nx, Ny = 160, 80
    tc = utils.test_case(base_fn,py_directory,Nx,Ny,et)

    p_ref = 8.61 * 1e4
    times = np.arange(0,164)#[:19]
    times = np.arange(0,30)

    probe_loc = [20,40]
    l_typ = 'WINDOW_STEP'

    rhoY = False
    CFL = True

    attr_labels = pt.labels_increment()

    def get_ens(sfx , diff, attribute):
        ens = tc.get_ensemble(times, N, attribute, sfx, label_type=l_typ, avg=True, diff=False)
        rhoY = tc.get_ensemble(times, N, 'rhoY', sfx, label_type=l_typ, avg=True, diff=False)
        tt = tc.t_arr
        ens = ens[1:]
        rhoY = rhoY[1:]
        tt = tt[1:]

        ud = ic_loader(tc,sfx,False)
        p = p_converter(ens,rhoY,ud)
        probe = p[:,probe_loc[0],probe_loc[1]]
    #     probe *= p_ref * Msq
        probe = probe[1:] - probe[:-1]

        label = sfx + '_' + attribute
        return label, probe, tt

    pls = []
    ss = 'idt=2'

    for i,attribute in enumerate(attributes):
        diff = True if attribute == 'p2_nodes' or attribute == 'p2_half' else False

        sfx1 = 'psinc_noib'
        if CFL : sfx1 = 'psinc_noib_CFLfixed'
        attribute1 = 'p2_nodes'
        if rhoY : attribute1 = 'rhoY'
        l1, p1, tt = get_ens(sfx1, diff, attribute1)
        l1 = 'psinc'

        sfx2 = 'comp_imbal_noib'
        if CFL : sfx2 = 'comp_imbal_noib_CFLfixed'
        attribute2 = 'p2_nodes'
        if rhoY : attribute2 = 'rhoY'
        l2_attr = np.copy(attribute2)
        l2, p2, tt = get_ens(sfx2, diff, attribute2)
        l2 = 'comp w/o blending'

        sfx3 = 'comp_imbal_half_ib-0'
        if CFL : sfx3 = 'comp_imbal_half_CFLfixed_ib-0'
        attribute3 = 'p2_nodes'
        if rhoY : attribute3 = 'rhoY'
        l3, p3, tt = get_ens(sfx3, diff, attribute3)
        l3 = 'comp w/ blending'

        fs = (20*cm,10*cm) # used in the publication
        pl = pt.plotter_1d(figsize=fs,fontsize=14,ncols=1,nrows=1)
        ax = pl.get_ax(i)

        if diff == True:
            tt = tt[1:]

        ax.plot(tt, p1, 'k', ms=4, label=l1)
        ax.plot(tt, p3, 'C1', ms=6, markevery=1, label=l3)

        fn = 'rb_slides'

        ax.set_ylabel(r'$\delta p^\prime$ [Pa]')
        ax.set_xlabel(r'time [$\times 1000$ s]')
        ax.set_xlim([0.0,0.35])
        ax.grid()
        ax.legend()

        l3a = ''
        l3a += r'$\delta_{\rm{i}}$=%.4f' %((la.norm(p2[1:]-p1[1:]))/la.norm(p1[1:]))
        l3a += '\n'
        l3a += r'$\delta_{\rm{b}}$=%.4f' %((la.norm(p3[1:]-p1[1:]))/la.norm(p1[1:]))

        pl.img.tight_layout()

        fn = 'rb_cfl_%i_%i' %(probe_loc[0],probe_loc[1])
        pl.save_fig('./output/tmp/%s' %fn)

        print("====================")
        print("figure 7: CFL-determined time-step error at...")
        print("probe_loc = (%i, %i)" %(probe_loc[0],probe_loc[1]))
        print("w/o blending - psinc: %.6f" %(la.norm(p2-p1)))
        print("blending - psinc: %.6f" %(la.norm(p3[1:]-p1[1:])))

        print(l3a)

        ########################################
        #
        ########################################
        print("====================")

        probe_loc = [80,40]
        _, p1, _ = get_ens(sfx1, diff, attribute1)
        _, p2, _ = get_ens(sfx2, diff, attribute2)
        _, p3, _ = get_ens(sfx3, diff, attribute3)

        pl = pt.plotter_1d(figsize=fs,fontsize=14,ncols=1,nrows=1)
        ax = pl.get_ax(i)

        ax.plot(tt, p1, 'k', ms=4, label=l1)
    #     ax.plot(tt, p2, 'C0', label=l2)
        ax.plot(tt, p3, 'C3', ms=6, markevery=1, label=l3)

        fn = 'rb_slides'

        ax.set_ylabel(r'$\delta p^\prime$ [Pa]')
        ax.set_xlabel(r'time [$\times 1000$ s]')
        ax.set_xlim([0.0,0.35])
        ax.grid()
        ax.legend()

        l3b = ''
        l3b += r'$\delta_{\rm{i}}$=%.4f' %((la.norm(p2[1:]-p1[1:]))/la.norm(p1[1:]))
        l3b += '\n'
        l3b += r'$\delta_{\rm{b}}$=%.4f' %((la.norm(p3[1:]-p1[1:]))/la.norm(p1[1:]))
    #     ax.annotate(l3b, (0.82,-0.68))

        pl.img.tight_layout()

        fn = 'rb_cfl_%i_%i' %(probe_loc[0],probe_loc[1])
        pl.save_fig('./output/tmp/%s' %fn)

        print("====================")
        print("figure 7: CFL-determined time-step error at...")
        print("probe_loc = (%i, %i)" %(probe_loc[0],probe_loc[1]))
        print("w/o blending - psinc: %.6f" %(la.norm(p2-p1)))
        print("blending - psinc: %.6f" %(la.norm(p3[1:]-p1[1:])))
        
        print(l3b)
    
    
#######################################
#
# Helper function to generate figures
# 8 and 9
#
#######################################

def figures_8and9(switches,attribute):
    euler = switches[0]
    rb = switches[1]
    rhou = attribute[0]
    p2_nodes = attribute[1]
    
    # Truth plotter    
    l_typ = 'TIME'
    N = 1

    def ic_loader(tc, sfx1, swe, lbl):
        # load pickled instances of data used in simulation
        fn_pickle = tc.get_filename(N,sfx1,format='dat')
        path_pickle = tc.get_path(fn_pickle)

        file = open(path_pickle,'rb')
        ud = pickle.load(file)
        file.close()

        if lbl == 'euler':
            if rhou:
                p_ref = ud.rho_ref * ud.u_ref
            else:
                p_ref = ud.p_ref * ud.Msq
        elif lbl == 'rb':
            p_ref = ud.rho_ref * ud.u_ref
        return p_ref

    def get_ens(tc, times, sfx , diff, attribute, lbl, swe=False):
        ens = tc.get_ensemble(times, N, attribute, sfx, label_type=l_typ, avg=True, diff=diff, tag='after_full_step')[1]
        p_ref = ic_loader(tc,sfx,swe,lbl)
        if attribute == 'p2_nodes':
            ens = ens * p_ref
        if attribute == 'rhou':
            ens = ens * p_ref

        label = sfx + '_' + attribute
        return label, ens.T

    diff = False
    if euler:
        base_fn = "output_travelling_vortex"
        directory = "output_travelling_vortex"
        py_directory = "../../%s/" %directory

        et = 3.0
        times = [3.0]
        attr = 'rhou' if rhou else 'p2_nodes'
        lbl = 'euler'
        Nx, Ny = 64, 64
        euler_tc = utils.test_case(base_fn,py_directory,Nx,Ny,et)

        sfx1 = 'truth_ib-0'
        l2, a2 = get_ens(euler_tc, times, sfx1, diff, attr, lbl)
        la = 'Euler vortex'
        aa = a2
        lvls = np.arange(-0.0065,0.0015,0.0010) * 100.0
    elif rb:
        base_fn = "output_rising_bubble"
        directory = "output_rising_bubble"
        py_directory = "../../%s/" %directory

        et = 1.0
        times = [et]
        times = [1.0]
        attr = 'rhou'
        lbl = 'rb'
        Nx, Ny = 160, 80
        rb_tc = utils.test_case(base_fn,py_directory,Nx,Ny,et)

        sfx1 = 'truth_CFLfixed_ib-0'
    #     sfx1 = 'comp_imbal_half_CFLfixed_ib-0'
        _, a1 = get_ens(rb_tc, times, sfx1, diff, attr, lbl)
        la = 'Rising bubble'
        aa = a1
        lvls = np.linspace(300.0,302.0,11)

    ll = [aa, '']
    pl_lst = [ll]

    if euler:
        x_label = r'x [$\times 10$ km]'
        y_label = r'y [$\times 10$ km]'
        if rhou:
            x_loc = np.linspace(0,Nx-1,5)
            y_loc = np.linspace(0,Ny-1,5)
            axvline = 31.5
            axhline = 31.5
            lvls = np.arange(60,120,10)
            clbl = r'$\times 100$'
            clbl = r'   '
        else:
            x_loc = np.linspace(0,Nx,5)
            y_loc = np.linspace(0,Ny,5)
            axvline = 32
            axhline = 32
            lvls = np.arange(-60,10,10)
            clbl = r'$\times 10^{-5}$'
        x_axs = [-0.5,-0.25,0.0,0.25,0.5]
        y_axs = [-0.5,-0.25,0.0,0.25,0.5]

        pl = pt.plotter(pl_lst,ncols=1,figsize=(5,5),sharey=False)
    elif rb:
        x_label = r'x [$\times 10$ km]'
        y_label = r'y [$\times 10$ km]'
        x_loc = np.linspace(0,Nx-1,5)
        y_loc = np.linspace(0,Ny-1,5) 
        x_axs = [-1.0,-0.5,0.0,0.5,1.0]
        y_axs = [0.0,0.25,0.5,0.75,1.0]
        pl = pt.plotter(pl_lst,ncols=1,figsize=(8,4),sharey=False)
        lvls = np.arange(-8,10,2)

    if rb:
        pl.set_axes(x_locs=x_loc, y_locs=y_loc, x_axs=x_axs, y_axs=y_axs, x_label=x_label, y_label=y_label)
    else:
        pl.set_axes(x_locs=x_loc, y_locs=y_loc, x_axs=x_axs, y_axs=y_axs, x_label=x_label, y_label=y_label, axhline=axhline, axvline=axvline, cbar_label=clbl)
    _ = pl.plot(aspect='equal',method='contour',lvls=lvls)

    rhou_lbl = '_rhou' if rhou else ''
    pl.save_fig('./output/tmp/%s_truth%s' %(lbl,rhou_lbl))
    
    # Observations plotter
    import matplotlib.patches as patches

    if euler:
        Nx, Ny = 64, 64
        et = 3.0
        base_fn = 'output_travelling_vortex'
        pydir = '../../%s/' %base_fn
        tc = utils.test_case(base_fn, pydir, Nx, Ny, et)
    elif rb:
        Nx, Ny = 160, 80
        et = 1.0
        base_fn = 'output_rising_bubble'
        pydir = '../../%s/' %base_fn
        tc = utils.test_case(base_fn, pydir, Nx, Ny, et)

    tags = tc.get_tag_dict()

    N = 10
    if euler:
        sfx0 = 'wdawloc_all_wda_ib-0'
        lbl = 'euler'
    elif rb:
        sfx0 = 'wdawloc_rhou_rhov_wda_CFLfixed_ib-0'
        lbl = 'rb'
    sfx0 = tc.cb_suffix(1,0, '%s' %(sfx0))

    fn_pickle = tc.get_filename(N,sfx0,format='dat')
    path_pickle = tc.get_path(fn_pickle)

    i2 = (slice(2,-2),slice(2,-2))

    file = open(path_pickle,'rb')
    ud = pickle.load(file)
    elem = pickle.load(file)
    node = pickle.load(file)
    obs = pickle.load(file)
    obs_noisy = pickle.load(file)
    obs_mask = pickle.load(file)
    obs_covar = pickle.load(file)
    file.close()

    if rb:
        attr = 'rhou'
    elif euler:
        if rhou:
            attr = 'rhou'
        else:
            attr = 'p2_nodes'

    if euler:
        time_index = 11 # last observation time
    elif rb:
        time_index = 10 # last observation time

    attribute = attr
    obs_arr = obs[time_index][attribute][i2].T
    obs_arr = [obs_arr, 'observation']
    obs_n_arr = obs_noisy[time_index][attribute][i2].T
    obs_mask_arr = obs_mask[time_index][attribute]
    obs_mask_arr = obs_mask_arr[i2].T
    obs_noisy_masked = np.ma.array(obs_n_arr,mask=obs_mask_arr).filled(fill_value=np.nan)
    if euler:
        if rhou:
            pl_lst = [[obs_noisy_masked * ud.rho_ref * ud.u_ref,'']]
        else:
            data = obs_noisy_masked * ud.p_ref * ud.Msq
            data = data[~np.isnan(data)]
            pl_lst = [[obs_noisy_masked * ud.p_ref * ud.Msq,'']]  # dimensionless, multiplied by p_ref = 10^5.
        pl = pt.plotter(pl_lst,ncols=1,figsize=(5,5),sharey=False)
    elif rb:
        data = obs_noisy_masked * ud.rho_ref * ud.u_ref
        data = data[~np.isnan(data)]
        pl_lst = [[obs_noisy_masked * ud.rho_ref * ud.u_ref,'']] # units in kgms^{-1}
        pl = pt.plotter(pl_lst,ncols=1,figsize=(8,4),sharey=False)

    if euler:
        x_label = r'x [$\times 10$ km]'
        y_label = r'y [$\times 10$ km]'
        if rhou:
            x_loc = np.linspace(-0.5,Nx-1+0.5,5)
            y_loc = np.linspace(-0.5,Ny-1+0.5,5)
            x_axs = [-0.5,-0.25,0.0,0.25,0.5]
            y_axs = [-0.5,-0.25,0.0,0.25,0.5]
            axvline = 31.5
            axhline = 31.5
            shft = 15.0
            rect0 = patches.Rectangle((25.5+shft,25.5+shft),1.0,1.0,linewidth=1,edgecolor='none',facecolor='red')
            rect = patches.Rectangle((20.5+shft,20.5+shft),11,11,linewidth=1,edgecolor='r',facecolor='none')
            lvls = np.arange(50,120,10)
    #         lvls = np.arange(60,120,10)
            clbl = r'$\times 100$'
            clbl = r'   '
        else:
            x_loc = np.linspace(-0.5,Nx+0.5,5)
            y_loc = np.linspace(-0.5,Ny+0.5,5)
            x_axs = [-0.5,-0.25,0.0,0.25,0.5]
            y_axs = [-0.5,-0.25,0.0,0.25,0.5]
            axvline = 32
            axhline = 32
            rect0 = patches.Rectangle((25.5,25.5),1.1,1.1,linewidth=1,edgecolor='none',facecolor='red')
            rect = patches.Rectangle((20.5,20.5),11,11,linewidth=1,edgecolor='r',facecolor='none')
            lvls = np.arange(-70,20,10)
            lvls = np.arange(-60,20,10)
            clbl = r'$\times 10^{-5}$'
        rects = [rect0,rect]
        pl.set_axes(x_locs=x_loc, y_locs=y_loc, x_axs=x_axs, y_axs=y_axs, x_label=x_label, y_label=y_label, axhline=axhline, axvline=axvline, rects=rects, cbar_label=clbl)
    elif rb:
        x_label = r'x [$\times 10$ km]'
        y_label = r'y [$\times 10$ km]'
        x_loc = np.linspace(-0.5,Nx-0.5,5)
        y_loc = np.linspace(-0.5,Ny-0.5,5)
        x_axs = [-1.0,-0.5,0.0,0.5,1.0]
        y_axs = [0.0,0.25,0.5,0.75,1.0]
        axvline = 79.5
        axhline = 39.5
        lvls = np.arange(-10,12.5,2.5)
        lvls = np.arange(-8,10,2)
        pl.set_axes(x_locs=x_loc, y_locs=y_loc, x_axs=x_axs, y_axs=y_axs, x_label=x_label, y_label=y_label, axhline=axhline, axvline=axvline)

    _ = pl.plot(aspect='equal',method='imshow',lvls=lvls)

    pl.save_fig('./output/tmp/obs_%s_%s' %(lbl,attr))
    

#######################################
#
# Helper function to generate figures
# 10 to 13
#
#######################################

def figures_10to13(switches):
    euler = switches[0]
    rb = switches[1]
    MWR = True
    arXiv = False
    
    # Load test case specific parameters...
    def load(base_fn,Nx,Ny,et,Nz=None):
        pydir = '../../%s/' %base_fn

        tc = utils.test_case(base_fn, pydir, Nx, Ny, et, Nz=Nz)
        tags = tc.get_tag_dict()
        return tc, tags

    if euler:
        base_fn = "output_travelling_vortex"
        Nx, Ny = 64, 64
        et = 3.0
        tc, tags = load(base_fn,Nx,Ny,et)

        attr_labels = pt.labels()
        attributes = ['rho','rhoY','rhou','rhov']
        obs_attrs = 'rhou_rhov'
        aux = 'wda'
        aux_ref = ''
        lbl = 'euler'
        times = np.arange(0.0,3.01,0.01)[1:]
    #     times = np.arange(0.0,3.01,0.01)[1:10]

        x_axs = [-0.5,0.0,0.5]
        y_axs = [-0.5,0.0,0.5]
        axs0, axs1 = Nx, Ny

    elif rb:
        base_fn = "output_rising_bubble"
        Nx, Ny = 160,80
        et = 1.0
        tc, tags = load(base_fn,Nx,Ny,et)

        attr_labels = pt.labels()
        attributes = ['rho','rhoY','rhou','rhov']
        obs_attrs = 'rhou_rhov'
        aux = 'wda_CFLfixed'
        aux_ref = '_CFLfixed'
        lbl = 'rb'
        times = np.arange(1.0,10.1,0.1) / 10.0

        x_axs = [-10,-5,0.0,5,10]
        y_axs = [0.0,5.0,10.0]
        axs0, axs1 = Nx, Ny

    ens_noda_suffix = 'noda%s_ib-0' %aux_ref
    ens_noda_label = 'EnNoDA'

    ens_1_suffix = 'wdawloc_%s_%s_ib-0' %(obs_attrs,aux)
    ens_1_label = r'EnDA'

    ens_2_suffix = tc.cb_suffix(1,0, '%s' %(ens_1_suffix))
    ens_2_label = r'EnDAB'

    if euler:
        ens_1_label += r', $\{ \rho u, \rho v \}$'
        ens_2_label += r', $\{ \rho u, \rho v \}$'

        obs_attrs_34 = 'all'
        ens_3_suffix = 'wdawloc_%s_%s_ib-0' %(obs_attrs_34,aux)
        ens_3_label = r'EnDA, $\{ \rho, \rho u, \rho v, \rho \Theta, \pi \}$'
        ens_4_suffix = tc.cb_suffix(1,0, '%s' %(ens_3_suffix))
        ens_4_label = r'EnDAB, $\{ \rho, \rho u, \rho v, \rho \Theta, \pi \}$'

    # Ensemble size is the same across all test cases
    N = 10
    
    # Plot RMSES
    import matplotlib.pyplot as plt

    plt.style.use('default')

    fs = (16,8) # MWR and arXiv
    nrows = 2
    if MWR:
        pl = pt.plotter_1d(figsize=fs,fontsize=16,ncols=2,nrows=nrows)
    elif arXiv:
        pl = pt.plotter_1d(figsize=fs,fontsize=16,ncols=2,nrows=nrows)

    def ic_loader(tc, N, sfx, lbl):
        # load pickled instances of data used in simulation
        fn_pickle = tc.get_filename(N,sfx,format='dat')
        path_pickle = tc.get_path(fn_pickle)

        file = open(path_pickle,'rb')
        ud = pickle.load(file)
        file.close()

        return ud

    def unitify(arr,attribute,tc,N,sfx,lbl):
        ud = ic_loader(tc,N,sfx,lbl)
        if attribute == 'rho':
            arr *= ud.rho_ref
        if attribute == 'rhou' or attribute == 'rhov':
            arr *= ud.rho_ref * ud.u_ref
        if attribute == 'rhoY' or attribute == 'p2_nodes':
            arr *= ud.p_ref
            if lbl == 'euler':
                arr /= 1000.0 # units in kPa for the Euler RMSEs in rhoY and p2_nodes

        if attribute == 'p2_nodes':
            arr *= ud.Msq
        return arr


    def get_ens(sfx,attribute,lbl,diff=False):
        ens = tc.get_ensemble(times, N, attribute, sfx, label_type='TIME', avg=True, diff=diff)[1:]
        ens = unitify(ens, attribute, tc, N, sfx, lbl)

        return ens

    for i,attribute in enumerate(attributes):
        gt = 'n' if attribute == 'p2_nodes' else 'c'

        ens_noda = get_ens(ens_noda_suffix,attribute,lbl)
        ens_1 = get_ens(ens_1_suffix,attribute,lbl)
        ens_2 = get_ens(ens_2_suffix,attribute,lbl)

        if euler:
            ens_3 = get_ens(ens_3_suffix,attribute,lbl)
            ens_4 = get_ens(ens_4_suffix,attribute,lbl)

        truth = tc.get_ensemble(times, 1, attribute, 'truth%s_ib-0' %aux_ref, label_type='TIME',avg=True)
        truth = unitify(truth, attribute, tc, 1, 'truth%s_ib-0' %aux_ref,lbl)

        ax = pl.get_ax(i)
        avg = False

        diff_noda = tc.spatially_averaged_rmse(ens_noda,truth,avg=avg,grid_type=gt)
        diff_ens_1 = tc.spatially_averaged_rmse(ens_1,truth,avg=avg,grid_type=gt)
        diff_ens_2 = tc.spatially_averaged_rmse(ens_2,truth,avg=avg,grid_type=gt)
        if euler:
            diff_ens_3 = tc.spatially_averaged_rmse(ens_3,truth,avg=avg,grid_type=gt)
            diff_ens_4 = tc.spatially_averaged_rmse(ens_4,truth,avg=avg,grid_type=gt)

        l1 = ax.plot(times,diff_noda, 'k--', label=ens_noda_label)
        l2 = ax.plot(times,diff_ens_1, 'C1', label=ens_1_label)
        l3 = ax.plot(times,diff_ens_2, 'C2', label=ens_2_label)
        if euler:
            l4 = ax.plot(times,diff_ens_3, 'C1--', label=ens_3_label)
            l5 = ax.plot(times,diff_ens_4, 'C2--', label=ens_4_label)

        ax.set_title("%s" %attr_labels[attribute])
        ax.set_xlim([0.0,times[-1]])
        ax.grid()
        if not euler:
            ax.legend()

        plt.tight_layout(rect=[-0.0, -0.0, 1.0, 1.0])

        if MWR:
            if attribute == 'rho':
                ax.set_ylabel("spatially and ens. avg. RMSE")
            elif attribute == 'rhou':
                ax.set_ylabel("spatially and ens. avg. RMSE")
        elif arXiv:
            if attribute == 'rho':
                ax.set_ylabel("RMSE")
            elif attribute == 'rhou':
                ax.set_ylabel("RMSE")
        if euler:
            if attribute == 'rhou':
                ax.set_xlabel(r'time [$\times 100$ s]')
            elif attribute == 'rhov':
                ax.set_xlabel(r'time [$\times 100$ s]')
        elif rb:
            if attribute == 'rhou':
                ax.set_xlabel(r'time [$\times 1000$ s]')
            elif attribute == 'rhov':
                ax.set_xlabel(r'time [$\times 1000$ s]')

    if euler:
        pl.fig.tight_layout()
        if MWR:
            pl.fig.subplots_adjust(bottom=0.22)
        elif arXiv:
            pl.fig.subplots_adjust(bottom=0.22)
        pl.fig.legend((l1, l2, l3, l4, l5), 
                      labels=('EnNoDA', r'EnDA $\{ \rho u, \rho w \}$', r'EnDAB $\{ \rho u, \rho w \}$', 
                              r'EnDA $\{ \rho, \rho u, \rho w, P, \pi \}$', r'EnDAB $\{ \rho, \rho u, \rho w, P, \pi \}$'),
    #                   ncol=5,
                      ncol=3,
                      loc='lower center',
                      bbox_to_anchor=(0.5, -0.01, 0.0, 0.0),
                      fontsize='large'
                     )

    # plt.tight_layout()
    plt.savefig('./output/%sRmse.pdf' %(lbl), bbox_inches="tight")
    
    print("")
    print("RMSE figure plotted...")
    print("")
    
    ################################################################
    # Plot ensemble members
    ################################################################
    attributes = ['p2_nodes']

    # plot the very last output
    ens_time = [times[-1]] 
    tag = tags[9]

    class oo(object): pass
    ens0_oo = oo()
    ens1_oo = oo()
    ens_ref_oo = oo()

    if euler:
        ens_1_suffix = ens_3_suffix
        ens_2_suffix = ens_4_suffix

    if len(attributes) > 1:
        # for plotting of derived quantities, e.g. u,v,w, and Theta.
        for attribute in attributes:
            ens0 = tc.get_ensemble(times, N, attribute, ens_1_suffix, tag=tag, inner=True)[-1]
            ens1 = tc.get_ensemble(times, N, attribute, ens_2_suffix, tag=tag, inner=True)[-1]
            ens_ref = tc.get_ensemble(times, N, attribute, ens_noda_suffix, tag=tag, inner=True)[-1]
            setattr(ens0_oo,attribute,ens0)
            setattr(ens1_oo,attribute,ens1)
            setattr(ens_ref_oo,attribute,ens_ref)

        ens0 = getattr(ens0_oo,attributes[1]) / getattr(ens0_oo,attributes[0])
        ens1 = getattr(ens1_oo,attributes[1]) / getattr(ens1_oo,attributes[0])
        ens_ref = getattr(ens_ref_oo,attributes[1]) / getattr(ens_ref_oo,attributes[0])
    else:
        ens0 = tc.get_ensemble(times, N, attributes[0], ens_1_suffix, tag=tag, inner=True)[-1]
        ens1 = tc.get_ensemble(times, N, attributes[0], ens_2_suffix, tag=tag, inner=True)[-1]
        ens_ref = tc.get_ensemble(times, N, attributes[0], ens_noda_suffix, tag=tag, inner=True)[-1]

        ens0_rhoY = tc.get_ensemble(times, N, 'rhoY', ens_1_suffix, tag=tag, inner=True)[-1]
        ens1_rhoY = tc.get_ensemble(times, N, 'rhoY', ens_2_suffix, tag=tag, inner=True)[-1]
        ens_ref_rhoY = tc.get_ensemble(times, N, 'rhoY', ens_noda_suffix, tag=tag, inner=True)[-1]

    def p_converter(ens,rhoY,ud):
        dp2n = np.array([ (mem - mem.mean()) * ud.Msq for mem in ens ])
        kernel = np.ones((2,2))
        dp2c = np.array([signal.fftconvolve(mem, kernel, mode='valid') / kernel.sum() for mem in dp2n])

        P0 = (rhoY**(ud.gamm-1.0) - dp2c)**(1.0/(ud.gamm-1.0))
        p = rhoY**(ud.gamm) - P0**(ud.gamm)
        p *= ud.p_ref
        return p


    def ic_loader(tc, sfx, lbl):
        # load pickled instances of data used in simulation
        fn_pickle = tc.get_filename(N,sfx,format='dat')
        path_pickle = tc.get_path(fn_pickle)

        file = open(path_pickle,'rb')
        ud = pickle.load(file)
        file.close()

        return ud

    ud = ic_loader(tc,ens_1_suffix,lbl)

    ens_ref = p_converter(ens_ref, ens_ref_rhoY, ud)
    ens0 = p_converter(ens0, ens0_rhoY, ud)
    ens1 = p_converter(ens1, ens1_rhoY, ud)

    enses = [ens_ref, ens0, ens1]
    sfxes = [ens_noda_suffix, ens_1_suffix, ens_2_suffix]

    mean0 = np.mean(ens0,axis=0)
    mean1 = np.mean(ens1,axis=0)
    mean_ref = np.mean(ens_ref,axis=0)

    diff = False
    if diff == True:
        mean0 -= mean_ref
        mean1 -= mean_ref

    if MWR:
        mean0 = [mean0.T,'ensemble mean']
        mean1 = [mean1.T,'ensemble mean']
        mean_ref = [mean_ref.T,'ensemble mean']
    elif arXiv:
        mean0 = [mean0.T,'']
        mean1 = [mean1.T,'']
        mean_ref = [mean_ref.T,'ensemble mean']
    means = [mean_ref, mean0, mean1]

    arr_lst = []
    lvls = []
    ens_store, ens_diff = [], []
    for bb,ens in enumerate(enses):
        sfx = sfxes[bb]
        fs = (14,8)

        ens_arr = []
        for n,arr in enumerate(ens):
            arr_ref = ens_ref[n]
            arr_ref = arr_ref.T
            arr = arr.T

            if diff == True:
                arr -= arr_ref

            if rb:
                if bb == 1:
                    ens_store.append(arr)
                if bb == 2:
                    if MWR:
                        ens_diff.append([ens_store[n] - (arr), 'member index %i' %n])
                    elif arXiv:
                        ens_diff.append([ens_store[n] - (arr), ''])

            if bb == 0:
                if MWR:
                    ens_arr.append([arr, 'member index %i' %n])
                elif arXiv:
                    ens_arr.append([arr, 'member index %i' %n])
            else:
                ens_arr.append([arr, ''])


        # which members to plot?
        arr_lst += [ens_arr[1],ens_arr[5],ens_arr[7]]
        # append mean to plot
        arr_lst.append(means[bb])

        if rb:
            # pressure levels
            if bb == 0:
                lvls += [np.arange(-100,60,20)]*4
            elif bb == 1:
    #             lvls += [np.arange(-70,60,20)]*4
                lvls += [np.arange(-160,120,40)]*4
            elif bb == 2:
    #             lvls += [np.arange(-50,40,10)]*4
                lvls += [np.arange(-100,60,20)]*4


        if euler:
            # pressure levels
            if bb == 0:
                lvls += [np.arange(-200,40,40)]*4
            elif bb == 1:
    #             lvls += [np.arange(-550.0,650.0,150.0)]*4
                lvls += [np.arange(-1900,1300+640,640)]*4
            elif bb == 2:
                lvls += [np.arange(-200,40,40)]*4

    if rb:
        # which member-diffs to plot?
        arr_lst += [ens_diff[1],ens_diff[5],ens_diff[7]]
        # append mean to plot
        if MWR:
            mean_diff = [means[1][0] - means[2][0], 'ensemble mean']
        elif arXiv:
            mean_diff = [means[1][0] - means[2][0], '']
        arr_lst.append(mean_diff)
        lvls += [np.arange(-80,130,30)]*4 # pressure levels

    for idx, arr in enumerate(arr_lst):
        if idx % 4 == 0:
            print("===========")
            print("ensemble member index, minimum value, maximum value")
        print(idx, arr[0].min(),arr[0].max())

    nnx, nnz = axs0-1,axs1-1

    narr_lst = np.array(arr_lst)
    ncolslen = 4
    if MWR:
        fs = (20,10) if rb else (16,10)
    elif arXiv:
        fs = (16,9) if rb else (13,7)
    pl = pt.plotter(arr_lst,ncols=ncolslen,figsize=fs,sharexlabel=True,shareylabel=True)

    x_label = r'x [$\times 10$ km]'
    y_label = r'z [$\times 10$ km]'

    if rb:
        x_loc = np.linspace(0,nnx,5)
        y_loc = np.linspace(0,nnz,3)
        axh, axv = 39.5, 79.5
        pl.set_axes(x_locs=x_loc, y_locs=y_loc, x_axs=x_axs, y_axs=y_axs, x_label=x_label, y_label=y_label)
    #     pl.set_axes(x_locs=x_loc, y_locs=y_loc, x_axs=x_axs, y_axs=y_axs, axhline = axh, axvline = axv, x_label=x_label, y_label=y_label)
    else:
        x_loc = np.linspace(0,nnx,3)
        y_loc = np.linspace(0,nnz,3)
        axh, axv = 31.5, 31.5
        pl.set_axes(x_locs=x_loc, y_locs=y_loc, x_axs=x_axs, y_axs=y_axs, axhline = axh, axvline = axv, x_label=x_label, y_label=y_label)

    # mthd = 'imshow' if rb else 'contour'
    mthd = 'contour'
    _ = pl.plot(aspect='equal',method=mthd, lvls=lvls)

    pl.save_fig('./output/%sEnses' %lbl)
    
    print("")
    print("ensemble members plotted...")
    print("")
    
#######################################
#
# Helper function to generate figure B1
#
#######################################

def figure_B1(switches):
    import matplotlib.pyplot as plt
    euler = switches[0]
    rb = switches[1]
    
    output_anim = False
    
    # load experiment-specific parameters
    def load(base_fn,Nx,Ny,et,Nz=None):
        pydir = '../../%s/' %base_fn

        tc = utils.test_case(base_fn, pydir, Nx, Ny, et, Nz=Nz)
        tags = tc.get_tag_dict()
        return tc, tags

    if euler:
        base_fn = "output_travelling_vortex"
        Nx, Ny = 64, 64
        et = 3.0
        tc, tags = load(base_fn,Nx,Ny,et)

        attr_labels = pt.labels()
        attributes = ['rho','rhou','rhov']
        obs_attrs = 'rhou_rhov'
        aux = 'wda'
        aux_ref = ''
        lbl = 'euler'
        times = np.arange(0.0,3.01,0.01)[1:]
        da_times = np.arange(0.0,3.25,0.25)[1:]
        times, da_times = np.around(times, 2), np.around(da_times, 2)

        p_ref = 1e+2 # unit of contours in [x kPa]

    elif rb:
        base_fn = "output_rising_bubble"
        Nx, Ny = 160,80
        et = 1.0
        tc, tags = load(base_fn,Nx,Ny,et)

        attr_labels = pt.labels()
        attributes = ['rho','rhou','rhov']
        obs_attrs = 'rhou_rhov'
        aux = 'wda_CFLfixed'
        aux_ref = '_CFLfixed'
        lbl = 'rb'
        times = np.arange(1.0,10.1,0.1)[1:]/10.0
        da_times = np.arange(5.0,10.5,0.5)/10.0
        times, da_times = np.around(times, 2), np.around(da_times, 2)

        p_ref = 8.61 * 1e1 # unit of contours in [kPa]

    sfx_ref = 'noda%s_ib-0' %aux_ref
    lbl_ref = 'EnNoDA'

    sfx0 = 'wdawloc_%s_%s_ib-0' %(obs_attrs,aux)
    lbl0 = r'EnDA'

    sfx1 = tc.cb_suffix(1,0, '%s' %(sfx0))
    lbl1 = r'EnDAB'

    # Ensemble size is the same across all test cases
    N = 10
    tag = tags[9]
    
    ################################################################
    # Do and plot scale analysis
    ################################################################
    # load pickled instances of data used in simulation
    fn_pickle = tc.get_filename(N,sfx0,format='dat')
    path_pickle = tc.get_path(fn_pickle)

    file = open(path_pickle,'rb')
    ud = pickle.load(file)
    elem = pickle.load(file)
    node = pickle.load(file)
    file.close()

    # get dimensionless grid-size
    dx = np.diff(elem.x)[0]
    dy = np.diff(elem.y)[0]

    # remove all outputs that are right after the assimilation of data,
    # since we want to capture the velocities that go into the flow fields.
    times = np.setdiff1d(times,da_times)

    class oo(object): pass
    ens0_oo = oo()
    ens1_oo = oo()
    ens_ref_oo = oo()

    # load ensembles
    def get_ens_by_attr(times, N, attributes, sfx, tag=tag, inner=True):
        obj = oo()
        for attribute in attributes:
            ens = tc.get_ensemble(times, N, attribute, sfx, tag=tag, inner=False)#[-1]

            mean = np.mean(ens,axis=1)
            ens = np.append(ens,mean[:,np.newaxis,...],axis=1)
            setattr(obj,attribute,ens)

        # recover u and v, the velocity fields
        ens_u = getattr(obj,attributes[1]) / getattr(obj,'rho')
        ens_v = getattr(obj,attributes[2]) / getattr(obj,'rho')
        setattr(obj, 'u', ens_u)
        setattr(obj, 'v', ens_v)

        return obj

    # helper function to get partial derivatives
    def grad(arr,dd,direction,ud):
        dd *= ud.h_ref
        arr *= ud.u_ref
        # get partial derivatives
        if direction == 'x':
            axs = 0
        elif direction == 'y':
            axs = 1
        else:
            assert(0, 'direction unspported')

        return np.gradient(arr,dd,axis=axs)

    ens0 = get_ens_by_attr(times,N,attributes,sfx0)
    ens1 = get_ens_by_attr(times,N,attributes,sfx1)
    ens_ref = get_ens_by_attr(times,N,attributes,sfx_ref)

    enses = [ens_ref, ens0, ens1]
    sfxes = [sfx_ref, sfx0, sfx1]

    # recover magic numbers
    Ma = np.sqrt(ud.Msq)
    c = ud.u_ref / Ma
#     print("p_ref / rho_ref =", ud.p_ref/ ud.rho_ref)
#     print("c^2 =", c**2)

    arr_lst = []
    barPu = np.zeros((3,len(times)))


    def title_gen(frn):
        return "time = %.2f" %times[frn]

    # for each array, over all time, we calculate the scale analysis
    for bb,ens in enumerate(enses):
#         print("ens %i" %bb)
        arr_plt = np.empty_like(times, dtype='object')

        for tt,time in enumerate(times):
            us, vs = ens.u[tt], ens.v[tt]
            rho = ens.rho[tt][:,2:-2,2:-2] # remove ghost cells

            divu = np.zeros_like(us)
            divu = divu[:,2:-2,2:-2] # empty array without ghost cells
            for ii,uv in enumerate(zip(us,vs)):
                u, v = np.copy(uv[0]), np.copy(uv[1])
                # calculate divergence of the velocity fields with ghost cells
                div = grad(u,dx,'x',ud) + grad(v,dy,'y',ud)
                # take only the inner domain, i.e. without ghost cells
                divu[ii] = div[2:-2,2:-2]

            # put array for plotting if animation is enabled
            # this outputs plot for the 0th ensemble member
    #         if output_anim: arr_plt[tt] = [[divu[0].T,'ens %i, div(u)' %bb]]
            # this outputs plot for the ensemble mean
            if output_anim: arr_plt[tt] = [[divu.mean(axis=0).T,'ens %i, div(u) mean' %bb]]

            # calculate scale analysis
            Ptu = divu * rho * c**2 * ud.rho_ref
            Ptu *=  0.5*11.0*dx * ud.h_ref * Ma / ud.u_ref
#             Ptu /= ud.p_ref
            if euler:
                Ptu /= 1000.0

            # get norm for the contribution from the velocity fields at each time-point
            Ptu_norm = np.array([ np.sqrt((mem**2).mean()) for mem in Ptu[:-1]]).mean()
            Ptu_norm *= 2.0/np.pi
            barPu[bb,tt] = Ptu_norm * 1.0/1.4

        # output an animation of the div(u) evolution if output_anim is True
        if output_anim:
            fs = (8,6) if rb else (8,8)
            a2d = pt.animator_2D(arr_plt,ncols=1,figsize=fs)
            a2d.suptitle = title_gen
            a2d.method = 'contour'
            anim = a2d.animate(interval=500, aspect='equal', method='contour')

            import matplotlib.animation as animation

            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=10, metadata=dict(artist=''), bitrate=6000)
            anim.save('./output/%s_div_ens_%i.mp4' %(lbl,bb), writer=writer)

    # plotting...
    t_axs = times
    barPu00 = 0.0

    plt.figure(figsize=(5,3))

    plt.plot(t_axs, barPu[0] - barPu00, 'k', label=lbl_ref)
    plt.plot(t_axs, barPu[1] - barPu00, 'C1', label=lbl0)
    plt.plot(t_axs, barPu[2] - barPu00, 'C2', label=lbl1)

    if euler:
        plt.xlabel(r'time $[\times 100$ s$]$')
        plt.ylabel(r'$\hat{P}$ [kPa]')
        plt.axvline(0.25, c='k', lw=0.5)
    elif rb:
        plt.xlabel(r'time $[\times 100$ s$]$')
        plt.ylabel(r'$\hat{P}$ [Pa]')
        plt.axvline(5.0, c='k', lw=0.5)

    if rb:
        xet = 1.0
    elif euler:
        xet = 3.0
    plt.xlim([0, xet])
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig('./output/tmp/scale_analysis_%s.pdf' %lbl)
    
    
import matplotlib.pyplot as plt
#######################################
#
# Generate figure 4:
# initial plots for the vortex
#
#######################################
    
print("====================")
print("generating figure 4...")

# The inputs are:
# [euler, rb], [rhou, vortz, pot_temp, exner], end
figures_4and6([True,False],[1,0,0,0],False)
figures_4and6([True,False],[0,1,0,0],False)
figures_4and6([True,False],[0,0,1,0],False)
figures_4and6([True,False],[0,0,0,1],False)

print("merging figure 4...")

p1 = '%s/initial_euler_exner.pdf' %dt
p2 = '%s/initial_euler_rhou.pdf' %dt
p3 = '%s/initial_euler_vortz.pdf' %dt
p4 = '%s/initial_euler_pt.pdf' %dt

out_name = 'initialVortices'
ps = [p1,p2,p3,p4]
shape = '2x2'
merge_plots(ps,out_name,shape)

print("figure 4 done")
print("====================")
print("")
plt.close("all")

#######################################
#
# Generate figure 5:
# Initial blending plots for the vortex
#
#######################################

print("====================")
print("generating figure 5...")
figure_5()

print("merging figure 5...")

p1 = '%s/euler_w_imbal.pdf' %dt
p2 = '%s/euler_w_bal.pdf' %dt

out_name = 'comparisonswBalImbal'
ps = [p1,p2]
shape = '1x2'
merge_plots(ps,out_name,shape)

print("figure 5 done")
print("====================")
print("")
plt.close("all")

#######################################
#
# Generate figure 6:
# initial and end time plots for the 
# rising bubble
#
#######################################

print("====================")
print("generating figure 6...")
figures_4and6([False,True],[0,0,1,0],False)
figures_4and6([False,True],[0,0,1,0],True)

print("merging figure 6...")

p1 = '%s/initial_rb.pdf' %dt
p2 = '%s/final_rb.pdf' %dt

out_name = 'initialRb'
ps = [p1,p2]
shape = '1x2'
merge_plots(ps,out_name,shape)

print("figure 6 done")
print("====================")
print("")
plt.close("all")

#######################################
#
# Generate figure 7:
# initial blending plots for the rising
# bubble
#
#######################################

print("====================")
print("generating figure 7...")
figure_7()

print("merging figure 7...")

p1 = '%s/rb_deltap_contour.pdf' %dt
p2 = '%s/rb_20_40_w_imbal.pdf' %dt
p3 = '%s/rb_20_40_w_ib.pdf' %dt
p4 = '%s/rb_cfl_20_40.pdf' %dt
p5 = '%s/rb_80_40_w_ib.pdf' %dt
p6 = '%s/rb_cfl_80_40.pdf' %dt

out_name = 'risingBubble'
ps = [p1,p2,p3,p4,p5,p6]
shape = '2x3'
merge_plots(ps,out_name,shape)

print("figure 7 done")
print("====================")
print("")
plt.close("all")

#######################################
#
# Generate figure 8:
# initial truth and observation plots
# for the travelling vortex
#
#######################################

print("====================")
print("generating figure 8...")
figures_8and9([True,False],[1,0])
figures_8and9([True,False],[0,1])

print("merging figure 8...")
p1 = '%s/obs_euler_rhou.pdf' %dt
p2 = '%s/euler_truth_rhou.pdf' %dt
p3 = '%s/obs_euler_p2_nodes.pdf' %dt
p4 = '%s/euler_truth.pdf' %dt

out_name = 'obsTruth'
ps = [p1,p2,p3,p4]
shape ='2x2'
merge_plots(ps,out_name,shape)

print("figure 8 done")
print("====================")
print("")
plt.close("all")

#######################################
#
# Generate figure 9:
# initial truth and observation plots
# for the rising bubble
#
#######################################

print("====================")
print("generating figure 9...")
figures_8and9([False,True],[1,0])
figures_8and9([False,True],[1,0])

print("merging figure 9...")
p1 = '%s/obs_rb_rhou.pdf' %dt
p2 = '%s/rb_truth_rhou.pdf' %dt

out_name = 'rbObsTruth'
ps = [p1,p2]
shape ='1x2'
merge_plots(ps,out_name,shape)

print("figure 9 done")
print("====================")
print("")
plt.close("all")

#######################################
#
# Generate figures 10 and 11:
# RMSE and ensemble member plots for
# the travelling vortex ensemble
# experiments
#
#######################################

print("====================")
print("Subsequent figures take some time to generate...")
print("generating figures 10 and 11...")
figures_10to13([True,False])

print("figures 10 and 11 done")
print("====================")
print("")
plt.close("all")

#######################################
#
# Generate figures 12 and 13:
# RMSE and ensemble member plots for
# the rising bubble ensemble
# experiments
#
#######################################

print("====================")
print("generating figures 12 and 13...")
figures_10to13([False,True])

print("figures 12 and 13 done")
print("====================")
print("")

#######################################
#
# Generate figure B1:
# Scale analysis
#
#######################################

print("====================")
print("generating figure B1...")
figure_B1([True,False])
figure_B1([False,True])

print("merging figure B1...")
p1 = '%s/scale_analysis_euler.pdf' %dt
p2 = '%s/scale_analysis_rb.pdf' %dt

out_name = 'scaleAnalysis'
ps = [p1,p2]
shape = '2x1'
merge_plots(ps,out_name,shape)

print("figure B1 done")
print("====================")
print("")
plt.close("all")

print("All figures generated")