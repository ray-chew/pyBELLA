"""Input/output for pyBELLA.

Split from the former monolithic ``io.py`` into cohesive submodules; this
package re-exports the same public names so existing ``from ..utils import io``
call sites are unchanged:

* :mod:`.writer`  — the HDF5 ``hdf5`` writer and the ensemble-output bootstrap
  ``initialise``.
* :mod:`.restart` — ``read_input``, ``sim_restart``, ``fn_gen``.
* :mod:`.debug`   — ``NullDebugWriter`` / ``DebugWriter`` / ``create_debug_writer``.
* :mod:`.cli`     — ``get_args`` (argument parsing), ``init_logger``, ``mkdir_p``.
"""

from .writer import hdf5, initialise
from .restart import read_input, sim_restart, fn_gen
from .debug import NullDebugWriter, DebugWriter, create_debug_writer
from .cli import get_args, init_logger, mkdir_p

__all__ = [
    "hdf5",
    "initialise",
    "read_input",
    "sim_restart",
    "fn_gen",
    "NullDebugWriter",
    "DebugWriter",
    "create_debug_writer",
    "get_args",
    "init_logger",
    "mkdir_p",
]
