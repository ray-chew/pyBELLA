"""
Pre-step modifications for the time stepper.
"""

def apply_modifcations(dt, ud, step):
    """Apply modifications to the timestep based on user-defined parameters."""

    dt = _apply_cfl_override(dt, ud, step)

    return dt


def _apply_cfl_override(dt, ud, step):
    """Apply CFL override for fixed timestep cases."""
    if "CFLfixed" in ud.aux and step < 2:
        return 21.69 / ud.t_ref
    return dt