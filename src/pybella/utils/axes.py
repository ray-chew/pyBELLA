"""Axis geometry — the single source of truth for the role/array-axis map.

The solver's dynamics is written in **role space** ``(h1, v, h2)`` — first
horizontal, vertical, second horizontal. The cyclic permutation

    role_perm(v) = ((v - 1) % 3, v, (v + 1) % 3)

maps roles onto array axes; ``v = 1`` (the default y-vertical convention)
gives the identity, so y-vertical configurations remain unchanged.

Only *cyclic* (even) permutations are allowed: the rotation vector is a
pseudovector, so odd axis swaps would flip every Coriolis cross-term sign.
This is why the vertical choice selects a cyclic role layout instead of an
arbitrary one.

Everything that needs to know "which axis is vertical" must consume this
module rather than hardcoding axis 1.
"""

import numpy as np

VERTICAL_DEFAULT = 1

# axis-indexed (NOT role-indexed) component names
MOMENTA = ("rhou", "rhov", "rhow")
VELOCITIES = ("u", "v", "w")


def vertical_axis(ud):
    """The vertical/gravity array axis, validated in {0, 1, 2}."""
    v = int(getattr(ud, "gravity_direction", VERTICAL_DEFAULT))
    if v not in (0, 1, 2):
        raise ValueError(f"gravity_direction must be 0, 1 or 2, got {v}")
    return v


def role_perm(v):
    """Roles -> axes: (axis of h1, axis of v, axis of h2). v=1 -> (0, 1, 2)."""
    return ((v - 1) % 3, v, (v + 1) % 3)


def role_of_axis(v):
    """Axes -> roles: inverse of role_perm. role_of_axis(v)[axis] = role index."""
    perm = role_perm(v)
    inv = [0, 0, 0]
    for role, axis in enumerate(perm):
        inv[axis] = role
    return tuple(inv)


def horizontal_axes(v):
    """The (h1, h2) axis pair for vertical v."""
    perm = role_perm(v)
    return (perm[0], perm[2])


def role_attrs(attrs, v):
    """Reorder an axis-indexed 3-tuple (e.g. MOMENTA) into role order."""
    perm = role_perm(v)
    return tuple(attrs[axis] for axis in perm)


def vertical_momentum(ud):
    return MOMENTA[vertical_axis(ud)]


def vertical_velocity(ud):
    return VELOCITIES[vertical_axis(ud)]


def coords_along(grid_obj, axis):
    """Coordinate array of a SpaceDiscr-like object along an axis."""
    return (grid_obj.x, grid_obj.y, grid_obj.z)[axis]


def extent_along(grid_obj, axis):
    """(cell count incl. ghosts, ghost count, spacing) along an axis."""
    return (grid_obj.sc[axis], grid_obj.igs[axis], grid_obj.dxyz[axis])


def wall_slabs(ndim, axis, depth=2):
    """Index tuples selecting the low/high boundary slabs along an axis.

    wall_slabs(3, 1) -> ((:, :2, :), (:, -2:, :)) as slice tuples — the
    any-axis form of the axis-1-specific ``[:, :2, ...]`` / ``[:, -2:, ...]``
    slabs.
    """
    lo = [slice(None)] * ndim
    hi = [slice(None)] * ndim
    lo[axis] = slice(None, depth)
    hi[axis] = slice(-depth, None)
    return tuple(lo), tuple(hi)


def expand_profile(profile_1d, ndim, vaxis, counts):
    """Broadcast a 1D vertical profile to the full grid by repetition.

    Equivalent to ``for dim in range(0, ndim, 2): expand_dims + repeat``
    when vaxis == 1 (ascending non-vertical dims), generalised to any
    vertical axis. ``counts[dim]`` is the target size along each
    non-vertical dim (e.g. ``elem.sc``).
    """
    out = profile_1d
    for dim in range(ndim):
        if dim == vaxis:
            continue
        out = np.expand_dims(out, dim)
        out = np.repeat(out, counts[dim], axis=dim)
    return out


def degenerate_axes(node):
    """Axes with a single interior cell layer (quasi-2D broadcast targets)."""
    return [dim for dim in range(node.ndim) if node.iisc[dim] == 2]


def permute_axes(arr, sigma):
    """Move reference axis i to twin axis sigma[i] (for permutation checks)."""
    n = arr.ndim
    return np.moveaxis(arr, list(range(n)), list(sigma[:n]))


def validate(ud, ndim):
    """2D runs are x-y by convention: the vertical must be axis 1."""
    v = vertical_axis(ud)
    if ndim == 2 and v != 1:
        raise ValueError(
            f"2D runs require gravity_direction == 1 (x-y plane, y vertical); got {v}"
        )
    return v
