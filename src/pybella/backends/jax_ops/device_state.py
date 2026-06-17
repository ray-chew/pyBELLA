"""Device-resident step: host<->device state marshalling (7-leaf pytree)."""

import numpy as np
import jax.numpy as jnp

from .device_config import _SOL_FIELDS


def to_device(mem):
    s = {name: jnp.asarray(getattr(mem.sol, name)) for name in _SOL_FIELDS}
    s["p2_nodes"] = jnp.asarray(mem.npf.p2_nodes)
    return s


def write_back(s, mem):
    for name in _SOL_FIELDS:
        getattr(mem.sol, name)[...] = np.asarray(s[name])
    mem.npf.p2_nodes[...] = np.asarray(s["p2_nodes"])
