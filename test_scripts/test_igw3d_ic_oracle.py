"""3D IC oracle for the internal-long-wave case (Phase E 3D DA OSSE).

The 3D branch of ``test_internal_long_wave.sol_init`` is a transcription of
the 2D path with z as a trailing broadcast axis. For an UNPERTURBED IC
(seed None) every operation is within-slab, so every z-slab of every 3D
field must be BITWISE identical to the 2D IC on the same (x, y) grid —
max-abs == 0.0, no tolerance. Also checks the seeded 'igw_theta'
perturbation is deterministic and actually perturbs.

Run: python test_scripts/test_igw3d_ic_oracle.py
"""

import importlib
import sys

import numpy as np

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.physics import thermodynamics as gd_thermodynamics
from pybella.flow_solver.utils import fields
from pybella.interfaces.ic_config import IC_MODULES
from pybella.utils import axes, user_data

INX, INY, INZ = 66, 17, 17  # 65 x-cells (initial_pressure needs odd), 16, 16
FIELDS = ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX")


def build(inz, seed=None, perturb_type=None):
    case = importlib.import_module(IC_MODULES["test_internal_long_wave"])
    ud = user_data.UserDataInit(
        **{**vars(case.UserData()), "inx": INX, "iny": INY, "inz": inz}
    )
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    ud.diag = False
    ud.aux = "nedas_wda"
    if perturb_type is not None:
        ud.perturb_type = perturb_type
    elem, node = dis_grid.grid_init(ud)
    axes.validate(ud, elem.ndim)
    th = gd_thermodynamics.ThermodynamicalQuantities(ud)
    sol = fields.CellSolField(elem.sc)
    npf = fields.NodePressureField(elem, node, ud)
    sol = case.sol_init(sol, npf, elem, node, th, ud, seed=seed)
    return sol, npf, elem


def main() -> int:
    sol2, npf2, elem2 = build(inz=1)
    sol3, npf3, elem3 = build(inz=INZ)
    assert elem3.ndim == 3, "3D grid did not build"

    worst = 0.0
    for name in FIELDS:
        f2 = np.asarray(getattr(sol2, name))            # (nx, ny)
        f3 = np.asarray(getattr(sol3, name))            # (nx, ny, nz)
        d = float(np.max(np.abs(f3 - f2[:, :, None])))
        worst = max(worst, d)
        assert d == 0.0, f"{name}: 3D z-slabs differ from 2D by {d:.3e}"
    # p2_nodes: node z axis; every node z-plane must equal the 2D field
    p2_2 = np.asarray(npf2.p2_nodes)
    p2_3 = np.asarray(npf3.p2_nodes)
    d = float(np.max(np.abs(p2_3 - p2_2[:, :, None])))
    assert d == 0.0, f"p2_nodes: 3D z-planes differ from 2D by {d:.3e}"

    # seeded perturbation: deterministic, non-trivial, z-modulated
    sa, _, _ = build(inz=INZ, seed=1234, perturb_type="igw_theta")
    sb, _, _ = build(inz=INZ, seed=1234, perturb_type="igw_theta")
    for name in FIELDS:
        assert np.array_equal(
            np.asarray(getattr(sa, name)), np.asarray(getattr(sb, name))
        ), f"seeded IC not deterministic for {name}"
    dpert = float(np.max(np.abs(np.asarray(sa.rho) - np.asarray(sol3.rho))))
    assert dpert > 0.0, "igw_theta perturbation is a no-op"
    zvar = float(np.std(np.asarray(sa.rho)[10, 8, :]))
    assert zvar > 0.0, "seeded member has no z-structure"

    print(
        "PASS: unperturbed 3D IC bitwise slab-identical to 2D "
        f"(max-abs {worst:.1e}); seeded member deterministic, "
        f"perturbation {dpert:.2e}, z-std {zvar:.2e}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
