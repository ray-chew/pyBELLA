"""Unit tests for the axis-geometry module (axial-agnosticity Phase 0)."""

import numpy as np
import pytest

from pybella.utils import axes


class _UD:
    def __init__(self, v=None):
        if v is not None:
            self.gravity_direction = v


def _perm_sign(p):
    sign = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                sign = -sign
    return sign


@pytest.mark.parametrize("v", [0, 1, 2])
def test_role_perm_is_cyclic_and_invertible(v):
    perm = axes.role_perm(v)
    assert sorted(perm) == [0, 1, 2]
    assert perm[1] == v  # vertical sits in the middle role slot
    assert _perm_sign(perm) == +1  # even (cyclic) only — pseudovector safety
    inv = axes.role_of_axis(v)
    for role, axis in enumerate(perm):
        assert inv[axis] == role


def test_v1_is_identity():
    assert axes.role_perm(1) == (0, 1, 2)
    assert axes.role_of_axis(1) == (0, 1, 2)
    assert axes.horizontal_axes(1) == (0, 2)
    assert axes.role_attrs(axes.MOMENTA, 1) == ("rhou", "rhov", "rhow")
    assert axes.vertical_momentum(_UD()) == "rhov"  # default v = 1
    assert axes.vertical_axis(_UD()) == 1


def test_role_attrs_rotates():
    assert axes.role_attrs(axes.MOMENTA, 2) == ("rhov", "rhow", "rhou")
    assert axes.role_attrs(axes.MOMENTA, 0) == ("rhow", "rhou", "rhov")


def test_vertical_axis_validation():
    with pytest.raises(ValueError):
        axes.vertical_axis(_UD(3))
    with pytest.raises(ValueError):
        axes.validate(_UD(2), ndim=2)
    assert axes.validate(_UD(1), ndim=2) == 1
    assert axes.validate(_UD(0), ndim=3) == 0


def test_wall_slabs_match_legacy_pattern():
    lo, hi = axes.wall_slabs(3, 1)
    a = np.arange(4 * 6 * 5).reshape(4, 6, 5)
    assert np.array_equal(a[lo], a[:, :2, :])
    assert np.array_equal(a[hi], a[:, -2:, :])
    lo0, hi0 = axes.wall_slabs(2, 0)
    b = np.arange(20).reshape(4, 5)
    assert np.array_equal(b[lo0], b[:2, :])
    assert np.array_equal(b[hi0], b[-2:, :])


@pytest.mark.parametrize("ndim", [2, 3])
def test_expand_profile_reproduces_legacy_for_v1(ndim):
    rng = np.random.default_rng(0)
    counts = (4, 7, 5)[:ndim]
    prof = rng.standard_normal(counts[1])

    # legacy construction (fields.py get_dSdy / get_S0c)
    legacy = prof
    for dim in range(0, ndim, 2):
        legacy = np.expand_dims(legacy, dim)
        legacy = np.repeat(legacy, counts[dim], axis=dim)

    new = axes.expand_profile(prof, ndim, 1, counts)
    assert new.shape == legacy.shape
    assert np.array_equal(new, legacy)  # bit-identical


def test_expand_profile_other_axes():
    prof = np.arange(5.0)
    out = axes.expand_profile(prof, 3, 2, (3, 4, 5))
    assert out.shape == (3, 4, 5)
    assert np.array_equal(out[1, 2, :], prof)
    out0 = axes.expand_profile(prof, 3, 0, (5, 3, 4))
    assert out0.shape == (5, 3, 4)
    assert np.array_equal(out0[:, 1, 2], prof)


def test_permute_axes_roundtrip():
    rng = np.random.default_rng(1)
    a = rng.standard_normal((3, 4, 5))
    sigma = (1, 2, 0)  # cyclic: ref axis i -> twin axis sigma[i]
    t = axes.permute_axes(a, sigma)
    assert t.shape == (5, 3, 4)
    inv = tuple(np.argsort(sigma))
    assert np.array_equal(axes.permute_axes(t, inv), a)
    # spot value check: a[i,j,k] should equal t[k,i,j] for sigma=(1,2,0)
    assert a[2, 1, 3] == t[3, 2, 1]
