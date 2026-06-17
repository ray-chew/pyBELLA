"""
For more details on this module, refer to the write-up :ref:`boundary_handling`.
"""

import numpy as np
from ....utils import axes
from ....utils import options as opts
from .common import get_ghost_padding
from ....backends import is_jax_backend


class CellBoundaryHandler:
    """Handles different types of boundary conditions for ghost cells."""

    def __init__(self, mem, ud):
        self.mem = mem
        self.ud = ud
        self.igs = mem.elem.igs
        self.ndim = mem.elem.ndim
        # physical vertical axis and the (axis-named) momentum components
        self.v_phys = axes.vertical_axis(ud)
        self.vert_mom = axes.MOMENTA[self.v_phys]
        self.hor_moms = tuple(m for i, m in enumerate(axes.MOMENTA) if i != self.v_phys)
        # terrain metric (None on uniform-Cartesian runs); during advection
        # sweeps it is flipped alongside the solution arrays
        self.metric = mem.elem.metric
        if self.metric is not None:
            a_h1, a_h2 = axes.horizontal_axes(self.v_phys)
            # physical-component momentum names matching G1/G2
            self.slope_moms = (axes.MOMENTA[a_h1], axes.MOMENTA[a_h2])

    def _slope_terms(self, sol, idx):
        """G1*mom_h1 + G2*mom_h2 at the given index (terrain only)."""
        m = self.metric
        out = m.G1[idx] * getattr(sol, self.slope_moms[0])[idx]
        if m.G2 is not None:
            out = out + m.G2[idx] * getattr(sol, self.slope_moms[1])[idx]
        return out

    def apply_no_gravity_boundary(self, sol, current_step, ghost_padding, idx):
        """Apply boundary conditions for axes without gravity."""
        bdry_type = self.ud.bdry_type[current_step]

        if bdry_type == opts.BdryType.PERIODIC:
            _set_boundary(sol, ghost_padding, "wrap", idx)
        elif bdry_type == opts.BdryType.WALL:
            # the wall-normal momentum is the component along the wall axis
            _set_boundary(
                sol,
                ghost_padding,
                "symmetric",
                idx,
                normal_mom=axes.MOMENTA[current_step],
            )
        elif bdry_type == opts.BdryType.RAYLEIGH:
            raise AssertionError("Rayleigh boundary only defined on the gravity axis.")

    def apply_gravity_boundary(self, sol, dim, ghost_padding, step):
        """Apply boundary conditions for axes with gravity.

        ``dim`` is the ARRAY axis of the boundary (during advection sweeps
        the data is flipped, so it differs from the physical vertical);
        gravity_strength is indexed by the PHYSICAL axis.
        """
        gravity_axis = dim
        g = self.ud.gravity_strength[self.v_phys]
        direction = -1.0
        offset = 0

        for side in ghost_padding[gravity_axis]:
            direction *= -1
            self._process_ghost_cells_side(sol, side, dim, direction, offset, step, g)
            offset += 1

    def _process_ghost_cells_side(self, sol, side, dim, direction, offset, step, g):
        """Process ghost cells for one side of the boundary."""
        y_axs = self.ndim - 1 if step is not None else self.v_phys

        for current_idx in np.arange(side)[::-1]:
            indices = self._get_gravity_indices(current_idx, direction, offset, y_axs)
            ghost_values = self._calculate_ghost_values(
                sol, indices, direction, g, y_axs
            )
            self._assign_ghost_values(sol, indices["image"], ghost_values)

    def _get_gravity_indices(self, current_idx, direction, offset, y_axs):
        """Get the indices for last, source, and image cells."""
        nlast, nsource, nimage = _get_gravity_padding(
            self.ndim,
            current_idx,
            direction,
            offset,
            self.mem.elem.sc[self.v_phys],
            self.mem.elem.igs[self.v_phys],
            y_axs=y_axs,
        )
        return {"last": nlast, "source": nsource, "image": nimage}

    def _calculate_ghost_values(self, sol, indices, direction, g, y_axs):
        """Calculate values for ghost cells with gravity."""
        nlast, nsource, nimage = indices["last"], indices["source"], indices["image"]

        # Calculate basic quantities
        Y_last = sol.rhoY[nlast] / sol.rho[nlast]
        Y_source = sol.rhoY[nsource] / sol.rho[nsource]

        vert = getattr(sol, self.vert_mom)
        if self.metric is not None:
            # the metric must be oriented like the (possibly sweep-flipped)
            # solution arrays — compute_advection flips them together
            assert self.metric.vaxis == y_axs, "metric not sweep-oriented"
            # free slip through the terrain surface: reflect the
            # CONTRAVARIANT momentum (mom_v - G.mom_h), not the Cartesian one
            contra_source = vert[nsource] - self._slope_terms(sol, nsource)
            rhoYv_image = -contra_source * sol.rhoY[nsource] / sol.rho[nsource]
            # stratification at the PHYSICAL height of the image cell
            S = 1.0 / self.ud.stratification(self.metric.z[nimage])
        else:
            rhoYv_image = -vert[nsource] * sol.rhoY[nsource] / sol.rho[nsource]
            y_coords = axes.coords_along(self.mem.elem, self.v_phys)
            S = 1.0 / self.ud.stratification(y_coords[nimage[y_axs]])

        # Calculate pressure difference
        dpi = self._calculate_pressure_difference(
            nlast, nimage, direction, g, Y_last, S, y_axs
        )

        # Calculate density and mass fraction
        rho, rhoY = self._calculate_density_and_mass_fraction(
            sol, nlast, nimage, dpi, S, y_axs
        )
        Y_image = rhoY / rho

        # Calculate velocity components
        velocities = self._calculate_velocities(
            sol, nsource, rhoYv_image, rhoY, Y_source, Y_image, direction
        )

        return {
            "rho": rho,
            "rhoY": rhoY,
            "hor": {
                m: getattr(sol, m)[nsource] / sol.rho[nsource] for m in self.hor_moms
            },
            "v": velocities["v"],
            "X": sol.rhoX[nsource] / sol.rho[nsource],
            "Th_slc": velocities.get("Th_slc", 1.0),
        }

    def _calculate_pressure_difference(
        self, nlast, nimage, direction, g, Y_last, S, y_axs
    ):
        """Calculate pressure difference for ghost cells."""
        if hasattr(self.ud, "ATMOSPHERIC_EXTENSION"):
            return (
                self.mem.npf.HydroState.p20[nimage[y_axs]]
                - self.mem.npf.HydroState.p20[nlast[y_axs]]
            ) * self.ud.Msq
        else:
            deta = self.mem.elem.dxyz[self.v_phys]
            if self.metric is not None:
                # local vertical cell extent dz = z_eta * deta across the
                # last -> image interval; z_eta = J / (N_v)_v (== J for
                # vertical-line maps, bit-exactly — on stretched grids J
                # carries the horizontal stretch factors too)
                m = self.metric
                nv = m.N[m.vaxis][m.cart_v]
                dz = 0.5 * (m.J[nimage] / nv[nimage] + m.J[nlast] / nv[nlast]) * deta
            else:
                dz = deta
            return direction * (self.mem.th.Gamma * g) * 0.5 * dz * (1.0 / Y_last + S)

    def _calculate_density_and_mass_fraction(self, sol, nlast, nimage, dpi, S, y_axs):
        """Calculate density and mass fraction for ghost cells."""
        if self.ud.is_compressible == 1:
            rhoY = ((sol.rhoY[nlast] ** self.mem.th.gm1) + dpi) ** self.mem.th.gm1inv
        else:
            rhoY = self.mem.npf.HydroState.rhoY0[nimage[y_axs]]

        rho = rhoY * S
        return rho, rhoY

    def _calculate_velocities(
        self, sol, nsource, rhoYv_image, rhoY, Y_source, Y_image, direction
    ):
        """Calculate velocity components for ghost cells."""
        result = {}

        vert = getattr(sol, self.vert_mom)
        if hasattr(self.ud, "ATMOSPHERIC_EXTENSION"):
            if direction > 0:  # bottom boundary
                result["v"] = (
                    vert[nsource] * Y_source / sol.rho[nsource] * rhoY / Y_image
                )
            else:  # top boundary
                result["v"] = vert[nsource] * Y_source
            result["Th_slc"] = (
                rhoY / (rhoY / Y_image) / (sol.rhoY[nsource] / sol.rho[nsource])
            )
        else:
            result["v"] = rhoYv_image / rhoY
            result["Th_slc"] = 1.0

        return result

    def _assign_ghost_values(self, sol, nimage, ghost_values):
        """Assign calculated values to ghost cells."""
        sol.rho[nimage] = ghost_values["rho"]
        for m, val in ghost_values["hor"].items():
            getattr(sol, m)[nimage] = ghost_values["rho"] * val * ghost_values["Th_slc"]
        sol.rhoY[nimage] = ghost_values["rhoY"]
        sol.rhoX[nimage] = ghost_values["rho"] * ghost_values["X"]

        # Handle the vertical component differently for atmospheric extension
        vert = getattr(sol, self.vert_mom)
        if hasattr(self.ud, "ATMOSPHERIC_EXTENSION"):
            vert[nimage] = -ghost_values["v"] / (
                ghost_values["rhoY"] / ghost_values["rho"]
            )
        elif self.metric is not None:
            # rho*v carries the reflected CONTRAVARIANT momentum; rebuild the
            # Cartesian vertical momentum with the ghost cell's slope terms
            # (the horizontal momenta were assigned just above)
            vert[nimage] = ghost_values["rho"] * ghost_values["v"] + self._slope_terms(
                sol, nimage
            )
        else:
            vert[nimage] = ghost_values["rho"] * ghost_values["v"]


def set_ghost_cells(mem, ud, step=None, sol=None):
    """
    In-place update of the ghost cells in :class:`management.variable.Vars`
    given the boundary conditions specified by :class:`inputs.user_data.UserDataInit`.

    Parameters
    ----------
    mem : object
        Memory container with sol, elem, th, and npf attributes
    ud : :class:`inputs.user_data.UserDataInit`
        Data container for the initial conditions
    step : int, optional
        Current step for advection directional Strang-splitting
    sol : object, optional
        Solution object, defaults to mem.sol
    """
    if is_jax_backend(ud):
        from ....backends.jax_ops import boundary as jax_boundary

        return jax_boundary.set_ghost_cells(mem, ud, step=step, sol=sol)

    if sol is None:
        sol = mem.sol

    # this is inefficient but whatever for now.
    handler = CellBoundaryHandler(mem, ud)
    dims = _get_dimensions_to_process(handler.ndim, step)

    for dim in dims:
        current_step = step if step is not None else dim
        ghost_padding, idx = get_ghost_padding(handler.ndim, dim, handler.igs)

        if ud.gravity_strength[current_step] == 0.0:
            handler.apply_no_gravity_boundary(sol, current_step, ghost_padding, idx)
        else:
            handler.apply_gravity_boundary(sol, dim, ghost_padding, step)


def _get_dimensions_to_process(ndim, step):
    """Determine which dimensions to process based on step parameter."""
    if step is None:
        return np.arange(ndim)
    else:
        return [ndim - 1]


# Functional approach with helper functions
def _pad_field(sol, field_name, idx, pads, mode):
    """Helper function to pad a single field"""
    if hasattr(sol, field_name):
        field = getattr(sol, field_name)
        if mode == "negative_symmetric":
            field[...] = np.pad(field[idx], pads, _negative_symmetric)
        else:
            field[...] = np.pad(field[idx], pads, mode)


def _set_boundary(sol, pads, btype, idx, normal_mom="rhov"):
    """
    Functional approach to setting the boundary. ``normal_mom`` names the
    wall-normal momentum component (mirrored with a sign flip); historically
    this was hardcoded to rhov, which broke walls on non-vertical axes.
    """
    tangential = ["rho", "rhoY", "rhoX"] + [
        m for m in ("rhou", "rhov", "rhow") if m != normal_mom
    ]

    # Define field groupings for each boundary type
    boundary_specs = {
        "symmetric": [
            (tangential, "symmetric"),
            ([normal_mom], "negative_symmetric"),
        ],
        "constant": [
            (tangential, "symmetric"),
            ([normal_mom], "constant"),
        ],
        "wrap": [(["rho", "rhou", "rhov", "rhow", "rhoY", "rhoX"], "wrap")],
    }

    if btype not in boundary_specs:
        raise ValueError(f"Unsupported boundary type: {btype}")

    # Apply padding to each group of fields
    for field_names, padding_mode in boundary_specs[btype]:
        for field_name in field_names:
            _pad_field(sol, field_name, idx, pads, padding_mode)


def _negative_symmetric(vector, pad_width, iaxis, kwargs=None):
    """
    Taken from the reference:

    Parameters
    ----------
    vector : ndarray
        A rank 1 array already padded with zeros. Padded values are vector `[:iaxis_pad_width[0]] and vector[-iaxis_pad_width[1]:]`.
    iaxis_pad_width : tuple
        A 2-tuple of ints, `iaxis_pad_width[0]` represents the number of values padded at the beginning of vector where `iaxis_pad_width[1]` represents the number of values padded at the end of vector.
    iaxis : int
        The axis currently being calculated.
    kwargs : dict
        Any keyword arguments the function requires.

    References
    ----------
    https://docs.scipy.org/doc/numpy/reference/generated/numpy.pad.html

    """
    if pad_width[1] > 0:
        sign = -1
        vector[: pad_width[0]] = sign * vector[pad_width[0] : 2 * pad_width[0]][::-1]
        vector[-pad_width[1] :] = sign * vector[-2 * pad_width[1] : -pad_width[1]][::-1]
        return vector
    else:  # axis must have length > 0 for padding
        return vector


def _get_gravity_padding(ndim, cur_idx, direction, offset, icv, igv, y_axs=None):
    """
    Parameters
    ----------
    ndim : int
        Number of dimensions.
    cur_idx : int
        The current index of the ghost cell in the gravity direction to be updated.
    direction : int
        Top of the domain, `direction=+1`, bottom of the domain, `direction=-1`.
    offset : int
        `offset=0`, index starts counting from 0,1.... `offset=1`, index starts counting from -1,-2,..., i.e. end-selection of the array.
    icv, igv : int
        Cell count (incl. ghosts) and ghost count along the gravity axis.
    y_axs : int, optional
        `Default == None`. Specifies the direction of the gravity axis. If `None`, then direction is the the y-axis.

    """
    cur_i = np.copy(cur_idx)
    cur_idx += offset * ((icv - 1) - 2 * cur_idx)
    gravity_padding = [slice(None)] * ndim
    if y_axs == None:
        y_axs = 1

    nlast = np.copy(gravity_padding)
    nlast[y_axs] = int(cur_idx + direction)

    nsource = np.copy(gravity_padding)
    nsource[y_axs] = int(offset * (icv) + direction * (2 * igv - (1 - offset) - cur_i))

    nimage = np.copy(gravity_padding)
    nimage[y_axs] = int(cur_idx)
    return tuple(nlast), tuple(nsource), tuple(nimage)
