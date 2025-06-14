
def get_neighbor_indices(ndim):
    """Create left and right neighbor indices for n-dimensional arrays."""
    lefts_idx = [slice(None)] * ndim
    rights_idx = [slice(None)] * ndim
    
    lefts_idx[-1] = slice(0, -1)
    rights_idx[-1] = slice(1, None)
    
    return tuple(lefts_idx), tuple(rights_idx)


def get_inner_slice(ndim):
    """Get slice tuple for inner cells (excluding outermost ghost cells)."""
    return tuple([slice(1, -1)] * ndim)