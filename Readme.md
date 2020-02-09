# RKLM_Python

## Running the atmospheric flow solver

To run the solve,
`python ./RKLM_Python/__main__.py`

* Initial condition files are used according to lines 21-24 of `__main__.py`.
* Set line line 64, `N=1` to run a non-Ensemble run, i.e. a run of only 1 member.
* Setting `N=1` will disable the data assimilation part of the code.
* The output filenames and directories are given in the initial condition data files in `./inputs/`.

## Example

Running the solver for the travelling vortex test case.

1. Leave only the line 21 uncommented of the initial conditions,
`from inputs.travelling_vortex_3D_48 import UserData, sol_init`
2. In `./inputs/travelling_vortex_3D_48.py`, adjust
    a. `self.tout` to determine end time and time to make outputs.
    b. `self.output_base_name` for the output directory.
    c. `self.inx` and `self.iny` for the grid resolution
3. The filename is currently determined by line 137 of `__main__.py`:
    `ud.output_name_comp = ("_ensemble=%i_%i_%i_%.1f" %(N,elem.icx-2*elem.igx,elem.icy-2*elem.igy,ud.tout[-1]))`
4. Run `python ./RKLM_Python/__main__.py`.
5. A file with the filename, e.g. `output_travelling_vortex_ensemble=1_32_32_10.0`, will be created in the directory `output_travelling_vortex` in the root directory of the git repository.

## Documentation

A draft documentation can be compiled by the code
`./docs/make html`
or
`./docs/make latexpdf`.

This requires:
1. `Sphinx 2.1.2`,
2. the [Napolean extension](https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html) for Numpy docstrings and
3. [the wild theme of Sphinx](https://pypi.org/project/wild_sphinx_theme/), since I had problems using the default themes.

The documentation will be improved in the coming days...

## Requirements

0. `python 3.6`
1. `numpy 1.16`
2. `scipy 1.3.1`
3. `h5py 2.9.0` - HDF5 output support
4. `numba 0.45.1` - JIT compilation of BiCGStab bottleneck
5. `dask 2.6.0` - Parallelisation support
6. `matplotlib 3.1.1` - Plot output



