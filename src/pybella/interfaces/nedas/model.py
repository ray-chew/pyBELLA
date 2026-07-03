"""PyBellaModel — NEDAS Model adapter wrapping the pyBELLA solver.

Contract verified against NEDAS==1.2.0 source (NEDAS/core/model.py); design
rationale and the risk register live in dev_notes/nedas_interface.md
("Phase N0/N1 findings"). Conventions pinned here:

- NEDAS fields are ghost-free 2D arrays in (y, x) order. pyBELLA 2D fields
  are (x, y) with 2-cell ghost frames. read_var/write_var own ALL
  transposition + ghost stripping/refill — nowhere else.
- io_mode is 'online' (single process): ``self.members`` (member id ->
  live ModelState) is the authority for tag 'current'; snapshot tags
  ('prior', 'post') use the plain ndarray ``self.memory`` dict of the base
  class, written by preprocess/postprocess.
- Time: NEDAS uses datetimes with hour-based arithmetic. Convention:
  1 NEDAS hour == 1 pyBELLA nondimensional time unit, so the MWR TV OSSE
  (da_times 0.25..3.0) runs with cycle_period 0.25 from time_start.
- The frozen ``data_assimilation`` layer is NOT imported; its member-field
  conventions (which attr lives on which container, ghost refill after an
  analysis write) are reimplemented here.
"""

import importlib

import numpy as np
from NEDAS.core import Model
from NEDAS.core.types import VarDesc
from NEDAS.grid import RegularGrid

from ...flow_solver.discretisation import grid as dis_grid
from ...flow_solver.discretisation import time_update as dis_time_update
from ...flow_solver.physics import thermodynamics as gd_thermodynamics
from ...flow_solver.utils import cache as fs_cache
from ...flow_solver.utils import fields
from ...flow_solver.utils.boundary import cell_boundary as bdry_c
from ...utils import axes, data_structures, options as opts, user_data
from ...utils.io.debug import NullDebugWriter
from ..dynamics_blending import schemes as blending_schemes
from ..ic_config import IC_MODULES

# DA-visible cell fields. p2_nodes is exposed for snapshots/diagnostics only:
# it lives on the node grid and must never enter state_def (the analysis grid
# is the cell grid); the MWR-parity experiments update {rhou, rhov} only.
CELL_VARS = ("rho", "rhou", "rhov", "rhoY")
NODE_VARS = ("p2_nodes",)


class PyBellaModel(Model[RegularGrid]):
    """pyBELLA as a NEDAS model: one live ModelState per ensemble member.

    Config keys (default.yml next to this file; override via model_def.pybella):
        ic: IC_MODULES key of the case (e.g. 'test_travelling_vortex')
        random_seed: member-seed chain root (archive convention: 888)
        ud_overrides: dict of UserData attribute overrides (stepmax, aux,
            initial_blending, continuous_blending, ...)
    """

    ic: str
    random_seed: int
    ud_overrides: dict | None
    restart_dt: float
    memory: dict

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # the base class derives model_name from the directory name ('nedas')
        self.model_name = "pybella"
        self.memory = {}
        self.members = {}
        self._seeds = None

        case = importlib.import_module(IC_MODULES[self.ic])
        self._sol_init = case.sol_init

        ud = user_data.UserDataInit(**vars(case.UserData()))
        for key, value in (self.ud_overrides or {}).items():
            setattr(ud, key, value)
        ud.coriolis_strength = np.array(ud.coriolis_strength)
        ud.diag = False  # the CompareSol regression gate is not for DA members
        self.ud = ud

        self.elem, self.node = dis_grid.grid_init(ud)
        axes.validate(ud, self.elem.ndim)
        self.th = gd_thermodynamics.ThermodynamicalQuantities(ud)
        self.bld = blending_schemes.Blend(ud)

        if self.elem.ndim != 2:
            raise NotImplementedError("PyBellaModel supports 2D x-y cases only")

        # analysis-grid geometry: inner cell centres, nondimensional units
        xc = self.elem.x[self.elem.igx : -self.elem.igx]
        yc = self.elem.y[self.elem.igy : -self.elem.igy]
        xx, yy = np.meshgrid(xc, yc)  # (ny, nx) — NEDAS rows-are-y convention
        cyclic = "".join(
            dim
            for dim, bt in zip("xy", ud.bdry_type[:2])
            if bt == opts.BdryType.PERIODIC
        )
        self.grid = RegularGrid(None, xx, yy, cyclic_dim=cyclic or None)
        self.grid.mask = np.full(xx.shape, False)

        levels = np.array([0])
        self.variables = {
            name: VarDesc(
                name=name,
                dtype="float",
                is_vector=False,
                dt=self.restart_dt,
                levels=levels,
                units="nondim",
                z_units="nondim",
            )
            for name in CELL_VARS + NODE_VARS
        }

    # --- helpers ----------------------------------------------------------

    def nondim_time(self, time) -> float:
        """NEDAS datetime -> nondim solver time (1 hour == 1 time unit)."""
        return round((time - self.c.config.time_start).total_seconds() / 3600.0, 3)

    def _inner(self, member: int, name: str) -> np.ndarray:
        """Ghost-free (y, x) view-copy of a member field."""
        mem = self.members[member]
        if name in NODE_VARS:
            arr, grid = mem.npf.p2_nodes, self.node
        else:
            arr, grid = getattr(mem.sol, name), self.elem
        inner = np.asarray(arr[grid.i2])
        return inner.reshape(grid.iicx, grid.iicy).T.copy()

    # --- grid / geometry ----------------------------------------------------

    def read_grid(self, **kwargs) -> None:
        pass  # built in __init__ from UserData

    def z_coords(self, **kwargs) -> np.ndarray:
        return np.zeros(self.grid.x.shape)  # 2D x-y: no vertical

    def filename(self, **kwargs):
        raise NotImplementedError(
            "PyBellaModel is online (in-memory) only; offline file cycling "
            "is the MPI fallback and is not implemented"
        )

    # --- state I/O (tag 'current' <-> live ModelState) -----------------------

    def read_var_from_memory(self, **kwargs) -> np.ndarray:
        kwargs = self.parse_kwargs(kwargs)
        if kwargs.get("tag", "current") != "current":
            return super().read_var_from_memory(**kwargs)
        member = kwargs["member"] if kwargs["member"] is not None else 0
        return self._inner(member, kwargs["name"])

    def write_var_to_memory(self, var, **kwargs) -> None:
        kwargs = self.parse_kwargs(kwargs)
        if kwargs.get("tag", "current") != "current":
            super().write_var_to_memory(var, **kwargs)
            return
        name = kwargs["name"]
        if name in NODE_VARS:
            raise NotImplementedError(
                "node-based fields are diagnostics-only; they cannot be "
                "written through the (cell-grid) analysis state"
            )
        member = kwargs["member"] if kwargs["member"] is not None else 0
        mem = self.members[member]
        arr = getattr(mem.sol, name)
        arr[self.elem.i2] = np.asarray(var).T.reshape(arr[self.elem.i2].shape)
        # analysis writes inner values only; refill the ghost frame, exactly
        # as the native analysis step does after set_field
        bdry_c.set_ghost_cells(mem, self.ud)

    # --- scheme hooks ---------------------------------------------------------

    def _snapshot(self, snap_tag: str, **kwargs) -> None:
        kwargs.pop("tag", None)  # the io backend injects tag='current'
        member = kwargs["member"]
        for name in self.variables:
            var = self._inner(member, name)
            super().write_var_to_memory(
                var, **{**kwargs, "tag": snap_tag, "name": name}
            )

    def preprocess(self, **kwargs) -> None:
        """Snapshot 'current' -> 'prior' (per member) for O-B statistics."""
        kwargs = self.parse_kwargs(kwargs)
        self._snapshot("prior", **kwargs)

    def postprocess(self, **kwargs) -> None:
        """Snapshot 'current' (analysis, post-updator) -> 'post' (per member)."""
        kwargs = self.parse_kwargs(kwargs)
        self._snapshot("post", **kwargs)

    def run(self, **kwargs) -> None:
        """Advance one member by forecast_period (ens_run_strategy 'scheduler').

        Mirrors the native queue driver (interfaces.dynamics_blending.prepare.
        init_da_window + __main__ window_step handling): the Blend object is
        handed to the stepper whenever the window STARTS at an analysis time
        (whether it acts is governed by ud.continuous_blending), plus the
        first two windows under ud.initial_blending; window_step resets to 0
        at every analysed window start.
        """
        kwargs = self.parse_kwargs(kwargs)
        self.run_status = "running"

        member = kwargs["member"]
        mem = self.members[member]
        fp = kwargs["forecast_period"]
        t0 = self.nondim_time(kwargs["time"])
        tout = round(t0 + fp, 3)

        cfg = self.c.config
        analysis_at_t0 = (
            cfg.run_analysis
            and kwargs["time"] >= cfg.time_analysis_start
            and kwargs["time"] <= cfg.time_analysis_end
        )
        outer_step = int(round(t0 / fp))

        blend = self.bld if analysis_at_t0 else None
        if self.ud.initial_blending and outer_step in (0, 1):
            blend = self.bld
        if analysis_at_t0:
            mem.time.window_step = 0

        self.members[member] = dis_time_update.do(
            mem, self.ud, tout, blend, debug_writer=NullDebugWriter()
        )
        self.run_status = "complete"

    # --- OSSE hooks -------------------------------------------------------------

    def generate_init_ensemble(self, **kwargs) -> None:
        """Build one member from FRESH containers with the archive seed chain.

        sol_init must run exactly once per member on fresh fields — re-running
        it on an initialised member double-applies the += initialisations
        (dev_notes/da_reinstatement.md, phase D1).
        """
        kwargs = self.parse_kwargs(kwargs)
        member = kwargs["member"]
        if self._seeds is None:
            rng = np.random.RandomState(self.random_seed)
            self._seeds = rng.randint(10000, size=self.c.nens)

        sol = fields.CellSolField(self.elem.sc)
        npf = fields.NodePressureField(self.elem, self.node, self.ud)
        sol = self._sol_init(
            sol, npf, self.elem, self.node, self.th, self.ud, seed=self._seeds[member]
        )
        mem = data_structures.ModelState(
            elem=self.elem,
            node=self.node,
            sol=sol,
            npf=npf,
            th=self.th,
            cache=fs_cache.FlowSolverCache(),
        )
        bdry_c.set_ghost_cells(mem, self.ud)
        self.members[member] = mem

    def generate_truth(self, **kwargs) -> None:
        raise NotImplementedError(
            "truth + obs come from the native pipeline "
            "(run_scripts/osse_mwr2022.py --runs obs truth); PyBellaObs "
            "serves that obs file, so NEDAS never generates its own truth"
        )
