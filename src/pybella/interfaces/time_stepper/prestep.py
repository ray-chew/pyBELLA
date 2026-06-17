"""
Pre-step modifications for the time stepper.
"""

# Legacy warm-bubble fixed timestep, in seconds (nondimensionalised by ud.t_ref
# at use). Pins dt for the first two steps of "CFLfixed" cases — e.g.
# test_blending_warm_bubble — to reproduce the reference run's startup.
_CFLFIXED_DT_SECONDS = 21.69


def apply_modifcations(dt, ud, step):
    """Apply modifications to the timestep based on user-defined parameters."""

    dt = _apply_cfl_override(dt, ud, step)

    return dt


def _apply_cfl_override(dt, ud, step):
    """Apply CFL override for fixed timestep cases."""
    if "CFLfixed" in ud.aux and step < 2:
        return _CFLFIXED_DT_SECONDS / ud.t_ref
    return dt
