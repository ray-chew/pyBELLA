

def get_inner_slice(ndim):
    """Get slice tuple for inner cells (excluding outermost ghost cells)."""
    return tuple([slice(1, -1)] * ndim)