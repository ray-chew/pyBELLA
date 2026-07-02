"""Member-field access for the data-assimilation layer.

The DA algorithms operate on stacked ``(N, ...)`` arrays; the ensemble members
(``ModelState``) own the storage. This module is the single place that knows on
which container a DA attribute lives (``CellSolField`` vs ``NodePressureField``),
replacing the pre-refactor ``results[:, loc, ...]`` container-index convention.

Rebinding semantics are intentional: ``set_field`` re-binds the attribute to the
analysis array, exactly as the reference implementation did
(``setattr(results[n][loc], attr, data)``).
"""

import numpy as np

# Node-grid attributes; every other DA attribute lives on member.sol.
# Face-based (flux) attributes were never supported by the DA layer.
NODE_ATTRS = frozenset({"p2_nodes"})


def is_node_attr(attr):
    """True if `attr` lives on the node-pressure container."""
    return attr in NODE_ATTRS


def get_field(member, attr):
    """Return the array bound to `attr` on a ModelState member."""
    container = member.npf if attr in NODE_ATTRS else member.sol
    return getattr(container, attr)


def set_field(member, attr, data):
    """Bind `data` as `attr` on the right container of a ModelState member."""
    container = member.npf if attr in NODE_ATTRS else member.sol
    setattr(container, attr, data)


def stack_fields(members, attr, slc=None, pads=None, pad_mode="wrap"):
    """Stack `attr` across members into an (N, ...) array.

    Parameters
    ----------
    members : iterable of ModelState
        Ensemble members (works with EnsembleState.members and with the
        object ndarray the driver passes around).
    slc : tuple of slice, optional
        Slice applied per member before stacking (e.g. the inner domain).
    pads : sequence, optional
        ``np.pad`` pad-width spec applied per member after slicing.
    pad_mode : str
        Mode forwarded to ``np.pad`` when `pads` is given.
    """
    fields = [get_field(mem, attr) for mem in members]
    if slc is not None:
        fields = [field[slc] for field in fields]
    if pads is not None:
        fields = [np.pad(field, pads, mode=pad_mode) for field in fields]
    return np.array(fields)


def scatter_fields(members, attr, data):
    """Write the rows of an (N, ...) array back onto the members."""
    for member, row in zip(members, data):
        set_field(member, attr, row)
