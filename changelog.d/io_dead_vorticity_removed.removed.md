Removed dead code from the HDF5 writer: the never-called ``vortz`` / ``vorty``
vorticity methods and ``dpress_dim`` (the last also referenced ``np.complex``,
removed in NumPy >= 1.24), plus their commented-out dataset entries.
