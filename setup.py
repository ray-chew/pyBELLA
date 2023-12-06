from setuptools import setup, find_packages

setup(name='RKLM', 
      version='0.5', 
      packages=find_packages(),
      install_requires=[
        'numba',
        'numpy==1.22.1',
        'h5py',
        'scipy',
        'PyYAML',
        'termcolor',
        'dask',
        'dask[distributed]',
        'matplotlib',
        'dill',
      ]
      )