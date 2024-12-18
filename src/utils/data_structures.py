from typing import Optional, Callable, Any
from dataclasses import dataclass

@dataclass
class ModelParameters:
    elem: Any
    node: Any
    Sol: Any
    flux: Any
    mpv: Any
    th: Any
    bld: Any

@dataclass
class DataAssimilationParameters:
    dap : Any
    rloc: Any
    sol_ens : Any


@dataclass
class RestartParameters:
    ud_rewrite: Optional[object] = None
    dap_rewrite: Optional[object] = None
    r_params: Optional[object] = None


@dataclass
class SimulationState:
    step: int
    t: float
    N: int
    restart: bool

    ud: object
    sol_init: Callable

    model_params: ModelParameters
    restart_params: RestartParameters
    da_params: Optional[DataAssimilationParameters] = None

    diag_comparison: Optional[object] = None