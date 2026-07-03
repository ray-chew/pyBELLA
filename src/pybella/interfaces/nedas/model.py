"""PyBellaModel — NEDAS Model adapter wrapping the pyBELLA solver (N1 skeleton).

Contract verified against NEDAS==1.2.0 source (NEDAS/core/model.py); design
rationale and the risk register live in dev_notes/nedas_interface.md
("Phase N0 findings"). All array/axis conventions are pinned here so the
implementation has a single reference:

- NEDAS fields are ghost-free 2D arrays in (y, x) order; vector variables
  carry a leading axis of length 2. pyBELLA arrays are (x, y[, z]) with
  2-cell ghost frames. read_var/write_var own ALL transposition + ghost
  stripping/refill — nowhere else.
- io_mode is 'online' (single process): the dict ``self.members`` of live
  ModelState objects is the authority for tag 'current'; snapshot tags
  ('prior', 'post', 'truth') use the plain ndarray ``self.memory`` dict of
  the base class.
- Time: NEDAS uses datetimes with hour-based arithmetic. Convention:
  1 NEDAS hour == 1 pyBELLA nondimensional time unit, so the MWR TV OSSE
  (da_times 0.25..3.0) runs with cycle_period 0.25 from time_start.
"""

import numpy as np
from NEDAS.core import Model
from NEDAS.core.types import VarDesc
from NEDAS.grid import RegularGrid

# Node-based fields (p2_nodes) are deliberately NOT part of the DA state for
# the MWR-parity experiments: the native rloc LETKF updates obs_attributes =
# {rhou, rhov} only, and blending — not the filter — acts on the pressure.
# See "PyBellaModel design" in dev_notes/nedas_interface.md.
CELL_SCALARS = ("rho", "rhoY", "rhoX")
MOMENTUM = ("rhou", "rhov")


class PyBellaModel(Model[RegularGrid]):
    """pyBELLA as a NEDAS model: one ModelState per ensemble member, in memory.

    Config keys (default.yml next to this file; override via model_def.pybella):
        ic: IC_MODULES key of the case (e.g. 'test_travelling_vortex')
        random_seed: member-seed chain root (archive convention: 888)
        truth_seed / obs settings are owned by the obs dataset + UserData aux.
    """

    ic: str
    random_seed: int
    memory: dict = {}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # base class sets model_name from the directory basename ('nedas')
        self.model_name = "pybella"

        # TODO(N1): resolve self.ic via interfaces.ic_config.IC_MODULES,
        # build UserData / elem / node / th once (shared geometry), then:
        #   self.grid = RegularGrid(None, x_cells, y_cells, cyclic_dim='xy')
        # from INNER cell-centre coords (nondim units; hroi must be consistent).
        self.members: dict[int, object] = {}  # member id -> ModelState

        levels = np.array([0])
        self.variables = {
            "momentum": VarDesc(
                name=MOMENTUM,
                dtype="float",
                is_vector=True,
                dt=self.restart_dt,
                levels=levels,
                units="nondim",
                z_units="nondim",
            ),
            **{
                attr: VarDesc(
                    name=attr,
                    dtype="float",
                    is_vector=False,
                    dt=self.restart_dt,
                    levels=levels,
                    units="nondim",
                    z_units="nondim",
                )
                for attr in CELL_SCALARS
            },
        }

    # --- grid / geometry -------------------------------------------------

    def read_grid(self, **kwargs) -> None:
        pass  # grid is built in __init__ from UserData; nothing to read

    def z_coords(self, **kwargs) -> np.ndarray:
        return np.zeros(self.grid.x.shape)  # 2D x-y cases: no vertical

    def filename(self, **kwargs):
        # online mode only; offline file cycling is the MPI fallback (not N1)
        raise NotImplementedError("PyBellaModel is online (in-memory) only")

    # --- state I/O (tag 'current' <-> live ModelState) --------------------

    def read_var_from_memory(self, **kwargs) -> np.ndarray:
        """Slice attr off members[m] (tag 'current') or the snapshot dict.

        TODO(N1): for 'current', use data_assimilation-frozen conventions
        REIMPLEMENTED here (that layer stays untouched): attr -> mem.sol.<attr>,
        strip ghosts with elem.i2, transpose (x,y)->(y,x); stack (rhou,rhov)
        for 'momentum'. Other tags -> super().read_var_from_memory(**kwargs).
        """
        raise NotImplementedError

    def write_var_to_memory(self, var, **kwargs) -> None:
        """Write inner-domain values back; refill ghost cells.

        TODO(N1): for 'current', transpose (y,x)->(x,y), write into
        mem.sol.<attr>[elem.i2], then bdry set_ghost_cells(mem, ud).
        Other tags -> super().write_var_to_memory(var, **kwargs).
        """
        raise NotImplementedError

    # --- scheme hooks ------------------------------------------------------

    def preprocess(self, **kwargs) -> None:
        """Snapshot 'current' -> 'prior' (per member) for O-B statistics."""
        raise NotImplementedError

    def postprocess(self, **kwargs) -> None:
        """Snapshot 'current' (now the analysis) -> 'post' (per member)."""
        raise NotImplementedError

    def run(self, **kwargs) -> None:
        """Advance one member by forecast_period (ens_run_strategy 'scheduler').

        TODO(N1): t_next = hours(kwargs['time'] + forecast_period - time_start);
        reset mem.time.window_step = 0 (blending re-enters its schedule each
        assimilation window, matching the native driver); then
            mem = dis_time_update.do(mem, self.ud, tout=t_next, bld=self.bld)
        with bld from InterfaceParameters(ud) built once in __init__.
        """
        raise NotImplementedError

    # --- OSSE hooks ---------------------------------------------------------

    def generate_init_ensemble(self, **kwargs) -> None:
        """Build one member from FRESH containers, seeded.

        TODO(N1): fresh CellSolField(elem.sc) + NodePressureField(elem, node,
        ud); sol_init(sol, npf, elem, node, th, ud, seed=seeds[member]); seeds
        from np.random.seed(self.random_seed) -> randint(10000, size=nens)
        (archive chain). NEVER re-run sol_init on an initialised member —
        the += initialisations double-apply (da_reinstatement.md, D1).
        """
        raise NotImplementedError

    def generate_truth(self, **kwargs) -> None:
        """N=1 truth trajectory (TV IC seed 2233), stored per cycle as 'truth'."""
        raise NotImplementedError
