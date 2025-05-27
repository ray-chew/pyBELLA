from typing import Optional, Callable, List, Any
from dataclasses import dataclass, field, fields

from ..flow_solver.discretisation.grid import Grid
from ..flow_solver.utils.variable import Vars
from ..flow_solver.physics.low_mach.mpv import MPV
from ..flow_solver.physics.gas_dynamics.thermodynamics import ThermodynamicalQuantities


@dataclass
class IntegrationTime:
    step: int = 0
    t: float = 0.0
    window_step: int = 0


@dataclass
class ModelState:
    elem: Grid
    node: Grid
    sol: Vars
    flux: List[Vars]
    mpv: MPV
    th: ThermodynamicalQuantities
    time: IntegrationTime = field(init=False)

    def __post_init__(self):
        self.time = IntegrationTime()

    def __iter__(self):
        return iter(getattr(self, field.name) for field in fields(self))


@dataclass
class InterfaceParameters:
    bld: Optional[Any] = None


@dataclass
class DataAssimilationParameters:
    # DA user-input parameters
    dap: Any
    # r-localisation function
    rloc: Any
    # solution ensemble
    # sol_ens : Any

    # observation related attributes
    obs: Any
    obs_noisy: Any
    obs_mask: Any
    obs_covar: Any


@dataclass
class RestartParameters:
    ud_rewrite: Optional[object] = None
    dap_rewrite: Optional[object] = None
    r_params: Optional[object] = None


@dataclass
class EnsembleState:
    members: List[ModelState] = field(default_factory=list)

    def update_member(
        self,
        elem: Grid,
        node: Grid,
        sol: Vars,
        mpv: MPV,
        flux: List[Vars],
        th: ThermodynamicalQuantities,
    ):
        new_state = ModelState(elem, node, sol, flux, mpv, th)
        self.members.append(new_state)

    def set_members(self, members: List[ModelState]):
        assert len(self.set_members == members)
        self.members = members

    def get_member(self, index: int) -> ModelState:
        return self.members[index]

    def get_all_members(self) -> List[ModelState]:
        return self.members

    def get_grid(self) -> tuple[Grid, Grid]:
        # Assuming identical underlying grid for all ensemble memebers
        elem = self.memebers[0].elem
        node = self.memebers[0].node
        return elem, node

    def __getitem__(self, index):
        return self.members[index]


@dataclass
class DiagnosticState:
    """
    Initialise diagnostic state for tests

    Consider removing run-specific parameters, e.g., (Nx, Ny), in future.
    """

    # the name to look up in test_targets.yml
    test_name: str
    # filename of the reference (if updt_target = True or plot_compare = True)
    file_name: str
    # details related to loading the reference fields
    Nx: int
    Ny: int
    steps: list
    path: str = "./outputs/"

    # plot the comparison?
    plot_compare: bool = True


@dataclass
class SimulationState:
    N: int
    restart: bool

    ud: object
    sol_init: Callable

    ensemble_state: EnsembleState
    restart_params: RestartParameters
    interface_params: Optional[InterfaceParameters] = None
    da_params: Optional[DataAssimilationParameters] = None

    diag_comparison: Optional[object] = None
