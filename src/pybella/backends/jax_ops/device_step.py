"""Device-resident jitted time step (`ud.backend = "jax-device"`).

One full semi-implicit step — RK half-advection, the two explicit/implicit
pairs, Strang advection, Rayleigh damping/forcing, optional diffusion, and
every ghost fill — traced into a single ``step(state, dt, ...)`` jit.
The host loop syncs one dt scalar per step and pulls the full state only at
output times (see :func:`run_window`).

State is a plain dict pytree with exactly 7 leaves
(rho, rhou, rhov, rhow, rhoY, rhoX, p2_nodes). Everything else —
geometry, thermodynamic constants, hydrostatic profiles, metric arrays in
every sweep orientation, boundary fill plans, laplacian gather plans,
sponge profiles — is frozen into a :class:`DeviceConfig` at window start
and closure-captured by ``make_step`` (XLA constants). Static jit
structure: Strang parity (2 compiled variants) and the regime ints, which
are constant per blending window (blending is guarded out).

Faithfulness notes:
- the numpy step's ``deepcopy`` saves (sol0, sol_half) become free
  reference holds on immutable arrays;
- ``npf`` scratch (p2_nodes0/rhs/wcenter/wplus/u,v,w) are trace locals;
- ``sol.pwchi`` (RK advection) has no live consumer and is dropped;
- the laplacian operators are rebuilt per solve from traced coefficients
  using static gather/mask plans — same arithmetic as the hybrid twins;
- bicgstab runs inside the jit (jax.scipy, scipy-matching semantics).
"""

# Split into device_config (host config + lap plans), device_state (state

# marshalling), device_kernels (traced substeps + compiled step) and

# device_loop (host window driver).

# Public entry points are re-exported so `device_step.run_window` is unchanged.

from .device_config import build_device_config

from .device_state import to_device, write_back

from .device_kernels import make_step

from .device_loop import run_window

__all__ = ["run_window", "make_step", "build_device_config", "to_device", "write_back"]
