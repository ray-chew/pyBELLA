Split the 800-line monolithic ``utils/io.py`` into an ``utils/io`` package with
cohesive submodules — ``writer`` (the ``hdf5`` writer + ``initialise`` bootstrap),
``restart`` (``read_input`` / ``sim_restart`` / ``fn_gen``), ``debug`` (the debug
writers), and ``cli`` (``get_args`` / ``init_logger`` / ``mkdir_p``). The public
names are re-exported from the package ``__init__`` so every ``from ..utils import
io`` call site is unchanged; output is bit-identical.
