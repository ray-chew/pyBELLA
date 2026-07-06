"""Tangent-plane momentum constraint for thin-shell (SWE) sphere runs.

With fixed Cartesian momentum components, advecting flow that follows a
curved surface generates the physical centripetal acceleration
-|v|^2/a e_r. In genuine shallow-water dynamics that acceleration is
absorbed by the constraint normal force (the free surface / rigid lid),
not by radial motion — so thin-shell runs project it out after every
momentum-modifying substep:

    m <- m - (m . e_up) e_up

Active iff ``ud.constrain_to_surface`` (set only by SWE-on-sphere cases;
3D compressible shells never set it — their radial dynamics is real).
Elementwise with the local up direction, ghost cells included; canonical
orientation (all call sites are outside the advection sweep flips).
"""


def apply(mem, ud):
    if not getattr(ud, "constrain_to_surface", False):
        return
    metric = mem.elem.metric
    e = metric.e_up
    sol = mem.sol
    moms = (sol.rhou, sol.rhov, sol.rhow)
    m_dot_e = moms[0] * e[0] + moms[1] * e[1] + moms[2] * e[2]
    for k in range(3):
        moms[k][...] -= m_dot_e * e[k]
