"""
For more details on this module, refer to the write-up :ref:`boundary_handling`.
"""

import numpy as np
from ....utils import axes
from ....utils import options as opts
from . import common as bdry
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
            if self.metric is not None and not self.metric.vertical_line:
                # curved wall (sphere): mirror the CONTRAVARIANT momentum
                # triple with the local normals — flipping one Cartesian
                # component is wrong where the wall is not a Cartesian plane
                _set_boundary(sol, ghost_padding, "symmetric", idx)
                _mirror_momenta_general(
                    sol, self.metric, current_step, self.ndim, self.igs[current_step]
                )
            else:
                # the wall-normal momentum is the component along the wall
                # axis (Cartesian walls; vertical-line maps keep Cartesian
                # wall planes on the non-gravity axes)
                _set_boundary(
                    sol,
                    ghost_padding,
                    "symmetric",
                    idx,
                    normal_mom=axes.MOMENTA[current_step],
                )
        elif bdry_type == opts.BdryType.POLE:
            self._apply_pole_boundary(sol)
        elif bdry_type == opts.BdryType.RAYLEIGH:
            raise AssertionError("Rayleigh boundary only defined on the gravity axis.")

    def _apply_pole_boundary(self, sol):
        """Pole ghost exchange on the phi axis.

        The pole is a coordinate singularity, not a wall: a ghost cell past
        |phi| = pi/2 IS the interior cell on the far side of the pole
        (longitude lambda + pi, mirror latitude). Because momenta are global
        Cartesian, every field — scalars and momenta alike — copies with no
        sign flip. lambda and phi array axes are read from the (possibly
        sweep-flipped) metric orientation.
        """
        m = self.metric
        assert m is not None and not m.vertical_line, "POLE needs a spherical metric"
        lam_axis, phi_axis = m.haxes[0], m.haxes[1]
        ig = self.igs[0]
        ncx = sol.rho.shape[lam_axis]
        ncz = sol.rho.shape[phi_axis]
        src_lam, src_phi, slabs = bdry.pole_source_indices(ncx, ncz, ig, nodal=False)
        for name in ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX"):
            bdry.pole_exchange_field(
                getattr(sol, name), lam_axis, phi_axis, src_lam, src_phi, slabs
            )

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
        tang = None
        v_coord = None
        if self.metric is not None:
            # the metric must be oriented like the (possibly sweep-flipped)
            # solution arrays — compute_advection flips them together
            assert self.metric.vaxis == y_axs, "metric not sweep-oriented"
            if self.metric.vertical_line:
                # free slip through the terrain surface: reflect the
                # CONTRAVARIANT momentum (mom_v - G.mom_h), not the
                # Cartesian one
                contra_source = vert[nsource] - self._slope_terms(sol, nsource)
            else:
                # general (sphere): the up-momentum is (N_v.m)/(N_v.e_up);
                # the tangential VELOCITY (u - (u.e)e at the source) is
                # copied per Cartesian component (generalizes copying the
                # horizontal components when e_up is a Cartesian axis)
                m = self.metric
                Nv = m.N[m.vaxis]
                moms = [getattr(sol, axes.MOMENTA[k]) for k in range(self.ndim)]
                Nv_dot_m = sum(
                    Nv[k][nsource] * moms[k][nsource] for k in range(self.ndim)
                )
                Nv_dot_e = sum(
                    Nv[k][nsource] * m.e_up[k][nsource] for k in range(self.ndim)
                )
                contra_source = Nv_dot_m / Nv_dot_e
                u_dot_e = (
                    sum(moms[k][nsource] * m.e_up[k][nsource] for k in range(self.ndim))
                    / sol.rho[nsource]
                )
                tang = [
                    moms[k][nsource] / sol.rho[nsource] - u_dot_e * m.e_up[k][nsource]
                    for k in range(self.ndim)
                ]
                # Well-balanced free slip: the wall mass flux vanishes iff the
                # rhoY flux Y*(N_v.m) (Y = rhoY/rho) is ODD about the wall, so
                # its convolution -> 0 gives a zero coordinate velocity at the
                # wall face. The plain up-velocity reflection below builds
                # (N_v.m)_g = rho_g v E_g with v routed through the SOURCE cell's
                # N_v.e_up, leaving Y_g (N_v.m)_g = -Y_s (N_v.m)_s * (E_g/E_s):
                # an O(dz) antisymmetry defect from the one-cell metric offset
                # -> an O(dz^2) mass leak under dynamics. Using the IMAGE cell's
                # E_g = N_v(img).e_up(img) instead makes it exact. Both are the
                # Jacobian J for the sphere; identical at rest (N_v.m = 0).
                E_image = sum(
                    Nv[k][nimage] * m.e_up[k][nimage] for k in range(self.ndim)
                )
                Y_src = sol.rhoY[nsource] / sol.rho[nsource]
                v_coord = -Y_src * Nv_dot_m / E_image  # * 1/rhoY_image (below)
            rhoYv_image = -contra_source * sol.rhoY[nsource] / sol.rho[nsource]
            # stratification at the PHYSICAL height of the image cell
            # (generalized altitude: == z for vertical-line maps)
            S = 1.0 / self.ud.stratification(self.metric.height[nimage])
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
        if v_coord is not None:
            # general free-slip: set v so ``_assign_ghost_values`` builds
            # (N_v.m)_g = rho_g v E_g with Y_g (N_v.m)_g = -Y_s (N_v.m)_s
            # exactly (rhoY here is the image rhoY, known now)
            velocities["v"] = v_coord / rhoY

        return {
            "rho": rho,
            "rhoY": rhoY,
            "hor": {
                m: getattr(sol, m)[nsource] / sol.rho[nsource] for m in self.hor_moms
            },
            "tang": tang,
            "v": velocities["v"],
            "X": sol.rhoX[nsource] / sol.rho[nsource],
            "Th_slc": velocities.get("Th_slc", 1.0),
        }

    def _hydro_at(self, arr, idx, y_axs):
        """Index a HydroState array at a gravity-ghost index.

        Profile-mode hydrostates are 1D vertical profiles: index by the
        vertical component of the ghost slice only. Field-mode hydrostates
        (terrain / sphere runs) are full grid-shaped fields that vary per
        column: index by the WHOLE ghost slice tuple, exactly like
        ``metric.height[idx]`` / ``sol.rhoY[idx]`` in the sibling branches.
        Indexing a field-mode array with the scalar ``idx[y_axs]`` slices
        the wrong axis and mis-shapes the result.
        """
        return arr[idx] if self.mem.npf.HydroState.field_mode else arr[idx[y_axs]]

    def _calculate_pressure_difference(
        self, nlast, nimage, direction, g, Y_last, S, y_axs
    ):
        """Calculate pressure difference for ghost cells."""
        if hasattr(self.ud, "ATMOSPHERIC_EXTENSION"):
            p20 = self.mem.npf.HydroState.p20
            return (
                self._hydro_at(p20, nimage, y_axs) - self._hydro_at(p20, nlast, y_axs)
            ) * self.ud.Msq
        else:
            deta = self.mem.elem.dxyz[self.v_phys]
            if self.metric is not None and not self.metric.vertical_line:
                # general map: vertical arc length per unit eta is |t_v|
                m = self.metric
                dz = 0.5 * (m.h_v[nimage] + m.h_v[nlast]) * deta
            elif self.metric is not None:
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
            rhoY = self._hydro_at(self.mem.npf.HydroState.rhoY0, nimage, y_axs)

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
        general = ghost_values["tang"] is not None
        if not general:
            for m, val in ghost_values["hor"].items():
                getattr(sol, m)[nimage] = (
                    ghost_values["rho"] * val * ghost_values["Th_slc"]
                )
        sol.rhoY[nimage] = ghost_values["rhoY"]
        sol.rhoX[nimage] = ghost_values["rho"] * ghost_values["X"]

        # Handle the vertical component differently for atmospheric extension
        vert = getattr(sol, self.vert_mom)
        if hasattr(self.ud, "ATMOSPHERIC_EXTENSION"):
            vert[nimage] = -ghost_values["v"] / (
                ghost_values["rhoY"] / ghost_values["rho"]
            )
        elif general:
            # general map: momenta = tangential part + beta e_up, with beta
            # enforcing the reflected up-momentum (N_v.m)/(N_v.e) = rho v
            # at the image cell's own metric
            m = self.metric
            Nv = m.N[m.vaxis]
            moms = [getattr(sol, axes.MOMENTA[k]) for k in range(self.ndim)]
            for k in range(self.ndim):
                moms[k][nimage] = (
                    ghost_values["rho"]
                    * ghost_values["tang"][k]
                    * ghost_values["Th_slc"]
                )
            Nv_dot_mt = sum(Nv[k][nimage] * moms[k][nimage] for k in range(self.ndim))
            Nv_dot_e = sum(Nv[k][nimage] * m.e_up[k][nimage] for k in range(self.ndim))
            beta = ghost_values["rho"] * ghost_values["v"] - Nv_dot_mt / Nv_dot_e
            for k in range(self.ndim):
                moms[k][nimage] += beta * m.e_up[k][nimage]
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


def _solve_normal_system(N_at, c, ndim):
    """Solve sum_k N[a][k] m_k = c_a per point (Cramer; det N = J^2 > 0).

    ``N_at[a][k]`` are the (already indexed) normal components at the
    target cells, ``c`` the contravariant triple to realize there.
    """
    if ndim == 2:
        det = N_at[0][0] * N_at[1][1] - N_at[0][1] * N_at[1][0]
        m0 = (c[0] * N_at[1][1] - c[1] * N_at[0][1]) / det
        m1 = (N_at[0][0] * c[1] - N_at[1][0] * c[0]) / det
        return [m0, m1]
    det = (
        N_at[0][0] * (N_at[1][1] * N_at[2][2] - N_at[1][2] * N_at[2][1])
        - N_at[0][1] * (N_at[1][0] * N_at[2][2] - N_at[1][2] * N_at[2][0])
        + N_at[0][2] * (N_at[1][0] * N_at[2][1] - N_at[1][1] * N_at[2][0])
    )
    out = []
    for k in range(3):
        M = [[c[a] if kk == k else N_at[a][kk] for kk in range(3)] for a in range(3)]
        det_k = (
            M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
        )
        out.append(det_k / det)
    return out


def _mirror_momenta_general(sol, metric, dim, ndim, igs_d):
    """Contravariant free-slip mirror of the momenta at both walls of
    ``dim`` (general curvilinear walls, e.g. the sphere's phi/r walls).

    For each ghost/source pair mirrored about the wall: the wall-normal
    contravariant component flips, the tangential contravariant
    components copy, each with its own LOCAL normals —
    N_b(ghost).m(ghost) = -N_b(src).m(src) exactly, so the (collocated,
    face-averaged) wall flux vanishes to roundoff. With N = identity
    this is exactly the Cartesian component flip.

    Scalars must already be padded (symmetric) so ghost rho/rhoY match.
    """
    moms = [getattr(sol, axes.MOMENTA[k]) for k in range(ndim)]
    N = metric.N
    n_cells = sol.rho.shape[dim]
    for k_layer in range(igs_d):
        for low in (True, False):
            if low:
                i_ghost = k_layer
                i_src = 2 * igs_d - 1 - k_layer
            else:
                i_ghost = n_cells - 1 - k_layer
                i_src = n_cells - 2 * igs_d + k_layer
            sl_g = [slice(None)] * ndim
            sl_s = [slice(None)] * ndim
            sl_g[dim] = i_ghost
            sl_s[dim] = i_src
            sl_g, sl_s = tuple(sl_g), tuple(sl_s)

            c = [
                sum(N[a][kk][sl_s] * moms[kk][sl_s] for kk in range(ndim))
                for a in range(ndim)
            ]
            c[dim] = -c[dim]
            N_ghost = [[N[a][kk][sl_g] for kk in range(ndim)] for a in range(ndim)]
            m_new = _solve_normal_system(N_ghost, c, ndim)
            for kk in range(ndim):
                moms[kk][sl_g] = m_new[kk]


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
    wall-normal momentum component (mirrored with a sign flip); it must
    follow the wall's axis — pinning it to rhov breaks walls on any
    non-vertical axis.
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
    ``np.pad`` callback mirroring with a sign flip; the signature is the
    one numpy.pad requires (see References).

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
