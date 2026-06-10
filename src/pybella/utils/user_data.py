import numpy as np
from collections import defaultdict

from . import sim_params
from . import options as opts


class DependencyManager:
    """Manages computational dependencies between attributes."""

    def __init__(self):
        # dependency_graph[computed_attr] = [list of required attributes]
        self.dependency_graph = {
            "u_ref": ["h_ref", "t_ref"],
            "Msq": ["u_ref", "R_gas", "T_ref"],
            "gravity_strength": [
                "grav",
                "h_ref",
                "R_gas",
                "T_ref",
                "gravity_direction",
            ],
            "i_gravity": ["grav", "h_ref", "R_gas", "T_ref", "gravity_direction"],
            "coriolis_strength": ["omega", "t_ref", "gravity_direction"],
            "cp_gas": ["gamm", "R_gas"],
            "N_ref": ["grav", "cp_gas", "T_ref"],
            "Nsq_ref": ["grav", "cp_gas", "T_ref"],
            "rho_ref": ["p_ref", "R_gas", "T_ref"],
            "Cs": ["gamm", "R_gas", "T_ref"],
        }

        # Reverse mapping: which computations depend on each attribute
        self.reverse_deps = defaultdict(list)
        for computed, deps in self.dependency_graph.items():
            for dep in deps:
                self.reverse_deps[dep].append(computed)

    def get_computations_to_update(self, changed_attr):
        """Get list of computations that need to be updated when an attribute changes."""
        return self.reverse_deps.get(changed_attr, [])

    def can_compute(self, obj, computation):
        """Check if all dependencies are available for a computation."""
        required_attrs = self.dependency_graph[computation]
        return all(hasattr(obj, attr) for attr in required_attrs)


class UserDataInit:
    """
    Loads user defined initial conditions with automatic dependency management.
    """

    def __init__(self, **kwargs):
        # Initialise dependency manager
        self._dep_manager = DependencyManager()
        self._updating = False  # Prevent infinite recursion

        # Load global constants first
        gconsts = sim_params.global_constants()
        for key, value in vars(gconsts).items():
            setattr(self, key, value)

        # Initialise with default values
        self._init_defaults()

        # Apply any user-provided kwargs
        if kwargs:
            for key, value in kwargs.items():
                setattr(self, key, value)

        # Initialise all computable attributes that can be computed
        self._initialise_computed_attributes()

    def _initialise_computed_attributes(self):
        """Initialise all computed attributes that have their dependencies available."""
        # Get all computed attributes from the dependency graph
        computed_attrs = list(self._dep_manager.dependency_graph.keys())

        for computation in computed_attrs:
            # Only compute if not already set and dependencies are available
            if not hasattr(self, computation) and self._dep_manager.can_compute(
                self, computation
            ):
                method_name = f"compute_{computation}"
                if hasattr(self, method_name):
                    getattr(self, method_name)()

    def _set_example_global_constants(self):
        """Set example global constants for demonstration."""
        self.nspec = 1
        self.buoy = 0
        self.grav = 9.81  # [m s^{-2}]
        self.omega = 0.0  # [s^{-1}]
        self.R_gas = 287.4  # [J kg^{-1} K^{-1}]
        self.R_vap = 461.0
        self.Q_vap = 2.53e06
        self.gamm = 1.4
        self.p_ref = 8.61 * 1e4  # [N/m^2]
        self.T_ref = 300.00  # [K]
        self.h_ref = 10000.0  # [m]
        self.t_ref = 100.0  # [s]

    def _init_defaults(self):
        """Initialise default values."""
        # Spatial grid
        self.inx = 64 + 1
        self.iny = 64 + 1
        self.inz = 1

        self.xmin = -1.0
        self.xmax = 1.0
        self.ymin = 0.0
        self.ymax = 1.0
        self.zmin = -1.0
        self.zmax = 1.0

        # Blending choices
        self.initial_blending = False

        self.continuous_blending = False
        self.no_of_pi_initial = 1
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.perturb_type = "pos_perturb"
        self.blending_mean = "rhoY"  # 1.0, rhoY
        self.blending_conv = "rho"  # theta, rho
        self.blending_type = "half"  # half, full
        self.blending_weight = 0.0 / 16

        # Vertical/gravity axis (array axis index; see utils/axes.py).
        # 2D runs are x-y by convention and require 1.
        self.gravity_direction = 1

        # Boundary conditions
        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = opts.BdryType.PERIODIC
        self.bdry_type[1] = opts.BdryType.WALL
        self.bdry_type[2] = opts.BdryType.WALL

        # Temporal
        self.CFL = 0.5
        self.dtfixed0 = 100.0
        self.dtfixed = 100.0
        self.acoustic_timestep = 0
        self.tout = np.arange(0.0, 1.01, 0.01)[10:]
        self.stepmax = 10000

        # Model regimes
        self.is_compressible = 1
        self.is_nonhydrostatic = 1
        self.is_ArakawaKonor = 0
        self.compressibility = 1.0

        # Physics and background wind
        self.u_wind_speed = 0.0
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0
        self.stratification = self.stratification_function

        # Explicit diffusion (off by default; see flow_solver/numerics/diffusion.py)
        self.diffusion = False
        self.diffusion_coeff = 0.0

        # Numerics
        self.do_advection = True
        self.limiter_type_scalars = opts.LimiterType.NONE
        self.limiter_type_velocity = opts.LimiterType.NONE
        self.tol = 1.0e-8
        self.max_iterations = 6000

        # Other attributes
        self.diag = False
        self.diag_state = None
        self.autogen_fn = False
        self.output_timesteps = False
        self.output_type = "output"
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

    def __setattr__(self, name, value):
        """Override setattr to handle dependency updates."""
        # Always set the attribute first
        super().__setattr__(name, value)

        # Skip dependency updates during initialisation or recursive updates
        if name.startswith("_") or not hasattr(self, "_dep_manager") or self._updating:
            return

        # Update dependent computations
        self._update_dependencies(name)

    def _update_dependencies(self, changed_attr):
        """Update all computations that depend on the changed attribute."""
        if self._updating:  # Prevent infinite recursion
            return

        self._updating = True
        try:
            computations_to_update = self._dep_manager.get_computations_to_update(
                changed_attr
            )

            for computation in computations_to_update:
                if self._dep_manager.can_compute(self, computation):
                    method_name = f"compute_{computation}"
                    if hasattr(self, method_name):
                        getattr(self, method_name)()
        finally:
            self._updating = False

    # Computation methods
    def compute_u_ref(self):
        """Compute reference velocity."""
        self.u_ref = self.h_ref / self.t_ref
        # u_ref change triggers Msq update automatically

    def compute_Msq(self):
        """Compute Mach number squared."""
        if hasattr(self, "u_ref"):
            self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

    def compute_gravity_strength(self):
        """Compute gravity-related parameters along the configured vertical axis."""
        from . import axes

        v = axes.vertical_axis(self)
        self.i_gravity = np.zeros(3)
        self.gravity_strength = np.zeros(3)

        self.gravity_strength[v] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        self.i_gravity[v] = 1

    # Alias for backward compatibility
    compute_i_gravity = compute_gravity_strength

    def compute_coriolis_strength(self):
        """Compute Coriolis parameters on the two horizontal axes."""
        from . import axes

        h1, h2 = axes.horizontal_axes(axes.vertical_axis(self))
        self.i_coriolis = np.zeros(3)
        self.coriolis_strength = np.zeros(3)

        self.coriolis_strength[h1] = self.omega * self.t_ref
        self.coriolis_strength[h2] = self.omega * self.t_ref

    def compute_cp_gas(self):
        """Compute specific heat at constant pressure."""
        self.cp_gas = self.gamm * self.R_gas / (self.gamm - 1.0)

    def compute_rho_ref(self):
        """Compute reference density."""
        self.rho_ref = self.p_ref / (self.R_gas * self.T_ref)

    def compute_N_ref(self):
        """Compute Brunt-Väisälä frequency."""
        if hasattr(self, "cp_gas"):
            self.N_ref = self.grav / np.sqrt(self.cp_gas * self.T_ref)
            self.Nsq_ref = self.N_ref * self.N_ref

    # Alias for backward compatibility
    compute_Nsq_ref = compute_N_ref

    def compute_Cs(self):
        """Compute sound speed."""
        self.Cs = np.sqrt(self.gamm * self.R_gas * self.T_ref)

    @staticmethod
    def stratification_function(y):
        """Default stratification function."""
        return 1.0

    def update_ud(self, obj):
        """Update multiple attributes at once."""
        # Temporarily disable dependency updates
        old_updating = self._updating
        self._updating = True

        try:
            # Set all attributes first
            for key, value in obj.items():
                super(UserDataInit, self).__setattr__(key, value)
        finally:
            self._updating = old_updating

        # Now update all dependencies at once
        if not self._updating:
            for key in obj.keys():
                self._update_dependencies(key)
