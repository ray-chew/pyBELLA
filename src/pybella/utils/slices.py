
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


def get_interface_indices(ndim):
    """
    Get complete set of indices for interface calculations.
    
    Parameters
    ----------
    ndim : int
        Number of dimensions
        
    Returns
    -------
    tuple
        (lefts_idx, rights_idx, inner_idx) where:
        - lefts_idx: indices for left neighbors
        - rights_idx: indices for right neighbors  
        - inner_idx: indices for inner cells (face_inner_idx)
    """
    lefts_idx, rights_idx = get_neighbor_indices(ndim)
    inner_idx = get_inner_slice(ndim)
    
    return lefts_idx, rights_idx, inner_idx


def get_all_slice_indices(ndim):
    """
    Get all commonly used slice indices for finite volume calculations.
    
    Parameters
    ----------
    ndim : int
        Number of dimensions
        
    Returns
    -------
    dict
        Dictionary containing all slice indices:
        - 'lefts': left neighbor indices
        - 'rights': right neighbor indices
        - 'inner': inner cell indices
        - 'face_inner': face inner indices (alias for inner)
    """
    lefts_idx, rights_idx, inner_idx = get_interface_indices(ndim)
    
    return {
        'lefts': lefts_idx,
        'rights': rights_idx, 
        'inner': inner_idx,
        'face_inner': inner_idx  # alias for clarity in interface flux calculations
    }