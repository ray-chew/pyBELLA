
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

def get_last_dim_inner_slice(ndim):
    """Get slice for inner faces in the last dimension only."""
    idx = [slice(None)] * ndim
    idx[-1] = slice(1, -1)
    return tuple(idx)


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



def get_averaging_indices(ndim, axis):
    """
    Get indices for averaging along a specific axis.
    
    Parameters
    ----------
    ndim : int
        Number of dimensions
    axis : int
        Axis along which to average
        
    Returns
    -------
    tuple
        (left_idx, right_idx) for averaging operation
    """
    left_idx = [slice(None)] * ndim
    right_idx = [slice(None)] * ndim
    
    left_idx[axis] = slice(0, -1)
    right_idx[axis] = slice(1, None)
    
    return tuple(left_idx), tuple(right_idx)


def get_periodic_inner_slice(igs, is_periodic, ndim):
    """
    Get inner slice accounting for periodic boundary conditions.
    
    Parameters
    ----------
    igs : array-like
        Number of ghost cells in each dimension
    is_periodic : array-like
        Boolean array indicating periodic boundaries
    ndim : int
        Number of dimensions
        
    Returns
    -------
    tuple
        Slice tuple for inner region with periodic boundaries
    """
    inner_idx = []
    for dim in range(ndim):
        start = igs[dim] - int(is_periodic[dim])
        end = -igs[dim] + int(is_periodic[dim]) if igs[dim] > int(is_periodic[dim]) else None
        inner_idx.append(slice(start, end))
    
    return tuple(inner_idx)


def get_face_center_averaging_indices(ndim):
    """
    Get indices for averaging finite difference results to face centers.
    
    Parameters
    ----------
    ndim : int
        Number of dimensions
        
    Returns
    -------
    dict
        Dictionary with 'y_avg' and 'x_avg' index tuples for 2D,
        or 'y_avg', 'x_avg', 'z_avg' for 3D averaging operations
    """
    indices = {}
    
    if ndim >= 2:
        # For averaging in y-direction (axis=1)
        indices['y_avg'] = (slice(None, -1), slice(None))
        indices['y_avg_right'] = (slice(1, None), slice(None))
        
        # For averaging in x-direction (axis=0)  
        indices['x_avg'] = (slice(None), slice(None, -1))
        indices['x_avg_right'] = (slice(None), slice(1, None))
        
    if ndim == 3:
        # For averaging in z-direction (axis=2)
        indices['z_avg'] = (slice(None), slice(None), slice(None, -1))
        indices['z_avg_right'] = (slice(None), slice(None), slice(1, None))
        
        # 3D specific averaging combinations
        indices['xy_avg'] = (slice(None, -1), slice(None, -1), slice(None))
        indices['xy_avg_right'] = (slice(1, None), slice(1, None), slice(None))
        
    return indices