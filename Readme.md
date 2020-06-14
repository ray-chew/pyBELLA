# RKLM_Python

## Running the atmospheric flow solver

In the top directory, run command
```code
python ./RKLM_Python/__main__.py -ic {tv_corr,rb,aw,tv_3d,igw,tv_2d,tv} [-N N]
```

* `-ic` specifies the initial condition file used. Choices are
    1. `tv` or `tv_2d` - 2D travelling vortex horizontal slice
    2. `tv_3d` - 3D travelling vortex
    3. `tv_corr` - 2D stack (3D) travelling vortex with coriolis
    4. `rb` - Rising bubble, 2D vertical slice
    5. `aw` - Acoustic wave, 2D slice
    6. `igw` - Skamarock / Klemp internal gravity waves, 2D vertical slice

* `-N` (Optional, int) defines the ensemble size. If not given, `N=1`. Setting `N=1` will disable the data assimilation part of the code.

* The output filenames and directories are defined in the respective initial condition data files in `./inputs/`.

## Example

Running the solver for the ensemble travelling vortex test case, `python ./RKLM_Python/__main__.py -ic tv -N 10`

1. In `./inputs/travelling_vortex_2D.py`, adjust
   a. `self.tout` to determine end time and time to make outputs.
   b. `self.output_base_name` for the output directory.
   c. `self.inx` and `self.iny` for the grid resolution
2. Adjust data assimilation parameters at `./data_assimilation/params.py`
5. Outputs are stored in a HDF5 file with the filename, e.g. `output_travelling_vortex_ensemble=10_32_32_10.0.h5` (ensemble size, grid_x, grid_y,end_time), in the directory `output_travelling_vortex`.

## Documentation

A draft documentation can be compiled by the code
```code
cd ./docs/
make html
```
or
```code
make latexpdf
```

This requires:
1. `Sphinx 2.1.2`,
2. the [Napolean extension](https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html) for Numpy docstrings and
3. [the wild theme of Sphinx](https://pypi.org/project/wild_sphinx_theme/), since I had problems using the default themes.

The documentation and docstrings are far from complete.

## Requirements

1. `python 3.6`
2. `numpy 1.16`
3. `scipy 1.3.1`
4. `h5py 2.9.0` - HDF5 output support
5. `numba 0.45.1` - JIT compilation of BiCGStab bottleneck
6. `matplotlib 3.1.1` - Plot output
7. `dask 2.6.0` - (Optional) parallelisation support



