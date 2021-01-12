Quickstart
==========
This is a brief introduction on how to use the RKLM-Py code. By the end of this document, you should

1. have installed the dependencies necessary to run RKLM-Py,
2. have ran a test simulation, and
3. be able to visualise the output.

Requirements
------------

These are the Python packages you have to install before RKLM-Py runs. The requirements can be found in `requirements.txt`::
   
   python==3.7.3
   numba==0.45.1
   numpy==1.17.3
   h5py==2.8.0
   scipy==1.3.1
   PyYAML==5.3.1
   termcolor==1.1.0
   dask==2.6.0
   matplotlib==3.1.3
    
While RKLM-Py code has been tested to work with these package versions, other versions may still work. 

If you have Anaconda installed, you may setup a virtual enviroment::

   conda create --name <env name> --file requirement.txt
    
Getting RKLM-Py
---------------
Having installed the required packages, clone the RKLM-Py repository. If you have SSH access::

   git@git.imp.fu-berlin.de:raychew/RKLM_Reference.git

Note that the latest stable branch is `develop`.

    
Running the flow solver
-----------------------
By default, the simulation outputs are saved to folders `./output_<test_case>`. Therefore, at the top directory of the cloned repository::

    mkdir output_travelling_vortex
    
This creates a directory to store the simulation outputs for the 2D Euler vortex experiment, which is used as an illustration in this documentation.

Now, `cd` to the `RKLM_Python` directory. From here, there are three ways to start a simulation.

From the command line
~~~~~~~~~~~~~~~~~~~~~
To run the code from the command line::

   python3 ./__main__.py -ic tv
   
The argument `-ic tv` runs the 2D Euler vortex experiment and outputs the results to the folder created above.

With run.py
~~~~~~~~~~~
`run.py` is a driver script used to run the program. Edit the attributes of the class `run_params` then run the file with::

   python3 ./run.py

With queue_run.py
~~~~~~~~~~~~~~~~~
`queue_run.py` allows you to set a queue. The following code shows a queue of two Euler vortex simulations that will be run sequentially, one with the compressible flow equations and the other with the pseudo-incompressible equations.

.. code-block:: python

    rp.N = 1
    rp.tc = 'tv'
    ud = {
        'is_compressible' : 1,
        'aux' : 'compressible'
    }

    dap = {
        'None' : None
    }

    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ud = {
        'is_compressible' : 0,
        'aux' : 'pseudo_compressible'
    }

    rp.ud = json.dumps(ud)
    rp.queue_run()


In this case, all initial parameters will remain the same between the two simulations, apart from `is_compressible` and `aux`. The former toggles compressibility and the latter specifies the tag used in generating the output filename. Running the script is then::

   python3 ./queue_run.py


Further examples
----------------
Brief introduction to more specify use cases of RKLM-Py is in this section.


Changing the initial conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Input files specifying the initial conditions are in the `inputs` module. For the 2D Euler vortex, the input file is `travelling_vortex_2D.py`. Each input file has a `UserData` class with the method `sol_init`. The class is initialised with the simulation parameters, e.g. grid-size, while the `sol_init` method populates the initial data containers :py:class:`management.variable.Vars` and :py:class:`physics.low_mach.mpv.MPV`. The former for the cell-based quantities :math:`\rho, \rho u, \rho v, \rho w, P, \rho \chi` and the latter for the node-based quantity :math:`\pi`. 

.. todo::

   Link to article on UserData.

Ensemble simulation
~~~~~~~~~~~~~~~~~~~
The optional argument `[-N <ensemble size>]` toggles ensemble simulation. This command runs an ensemble simulation for the 2D Euler vortex experiment with an ensemble size of 10 members::

   python3 ./__main__.py -ic tv -N 10
   
The attribute `N` in the driver files sets the ensemble size.

Data assimilation
~~~~~~~~~~~~~~~~~
Data assimilation only works for ensemble simulations, `N>1`. :py:class:`da_params` in :py:mod:`data_assimilation.params` defines the data assimilation parameters. To run an experiment with data assimilation, you will minimally need to specify the following attributes in :py:class:`da_params`, 

1. :py:attr:`data_assimilation.params.da_params.da_times` - time points to do data assimilation
2. :py:attr:`data_assimilation.params.da_params.obs_attributes` - the quantities to assimilate
3. :py:attr:`data_assimilation.params.da_params.obs_path` - the path to output file containing the observation fields

.. todo::

   Link to article on data assimilation.

Visualisation
-------------
RKLM-Py comes with some tools that aid in the visualisation and analysis of the output. Below is a detailed working example on how the scripts `utils.py` and `plotting_tools.py` in `visualiser_debugger` can be used.

.. code-block:: python

    import plotting_tools as pt
    import utils

    import numpy as np

    # quantities to read
    attributes = ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes']

    # end time and grid-size of the simulation
    et = 1.0
    Nx, Ny = 64, 64

    # the base filename of the output file
    base_fn = "output_travelling_vortex" 

    # path to the output file
    directory = "output_travelling_vortex"
    py_directory = "../%s/" %directory

    # load the output arrays
    tc = utils.test_case(base_fn,py_directory,Nx,Ny,et)

    # time label time, 'TIME by default'
    l_typ = 'TIME'
    # tag == 'after_full_step' by default
    tag = tc.get_tag_dict()[9]

    # get output at this time
    times = [0.01]

    # load plot titles
    attr_labels = pt.labels_increment()

    # helper function to load the ensemble, 
    def get_ens(tc, sfx, attribute):
        # ensemble size
        N = 1
        ens = tc.get_ensemble(times, attribute, sfx, label_type=l_typ, avg=True, tag=tag)[1]
        return ens

    # `aux` tag of the filename
    sfx1 = 'comp_bal'

    ll = []
    # loop through all the attributes and store them for plotting
    for acnt, attribute in enumerate(attributes):
        a2 = get_ens(tc, sfx1, attribute)
        
        # plotting_tools reads a list, each element of the list has size two.
        # the first element contains the array and the second element is the
        # plot title.
        ll.append([a2,attribute])
        
        # recover velocity fields from momenta and density fields.
        if attribute == 'rho':
            rho = np.copy(a2)
        if attribute == 'rhou' or attribute == 'rhov':
            vel = a2/rho
            ll.append([vel,attribute[-1]])

    # Setup plotter
    pl = pt.plotter(ll,ncols=3,figsize=(15,12),sharey=False)

    # plot settings
    x_axs = [-0.5,-0.25,0.0,0.25,0.5]
    y_axs = [-0.5,-0.25,0.0,0.25,0.5]
    x_loc = np.linspace(0,Nx-1,5)
    y_loc = np.linspace(0,Ny-1,5)
    x_label = r'x [$\times 10$ km]'
    y_label = r'y [$\times 10$ km]'

    pl.set_axes(x_locs=x_loc, y_locs=y_loc, x_axs=x_axs, y_axs=y_axs, x_label=x_label, y_label=y_label)

    # plot arrays of each attribute
    _ = pl.plot(aspect='equal',method='contour')

    # save a pdf output
    pl.save_fig('first_results')
