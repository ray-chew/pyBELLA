# -*- coding: utf-8 -*-
"""Example NumPy style docstrings.

This module demonstrates documentation as specified by the `NumPy
Documentation HOWTO`_. Docstrings may extend over multiple lines. Sections
are created with a section header followed by an underline of equal length.

Example
-------
Examples can be given using either the ``Example`` or ``Examples``
sections. Sections support any reStructuredText formatting, including
literal blocks::

    $ python example_numpy.py

Section breaks are created with two blank lines. Section breaks are also
implicitly created anytime a new section starts. Section bodies *may* be
indented:

.. code-block::

    rhoYv = Sol.rhoY * Sol.rhov / Sol.rho
    kernel_v = np.array([[0.5,1.,0.5],[0.5,1.,0.5]]).T
    flux[1].rhoY[1:-1,1:-1] = signal.convolve2d(rhoYv, kernel_v, mode='valid') /kernel_v.sum()
    flux[1].rhoY[:,-2] *= 0.
    ...

Notes
-----
    This is an example of an indented section. It's like any other section,
    but the body is indented to help it stand out from surrounding text.

If a section is indented, then a section break is created by
resuming unindented text.

Attributes
----------
module_level_variable1 : int
    Module level variables may be documented in either the ``Attributes``
    section of the module docstring, or in an inline docstring immediately
    following the variable.

    Either form is acceptable, but the two should not be mixed. Choose
    one convention to document module level variables and be consistent
    with it.


.. _NumPy Documentation HOWTO:
   https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

"""

import numpy as np
from scipy import signal
from management.debug import find_nearest

def recompute_advective_fluxes(flux, Sol):
    """
    Todo
    ----
    * 2D case for now - generalise in future

    """
    rhoYu = Sol.rhoY * Sol.rhou / Sol.rho

    kernel_u = np.array([[0.5, 0.5],[1., 1.],[0.5, 0.5]]).T
    flux[0].rhoY[1:-1,1:-1] = signal.convolve2d(rhoYu, kernel_u, mode='valid').T / kernel_u.sum()
    # flux[0].rhoY[:,-2] *= 0.

    rhoYv = Sol.rhoY * Sol.rhov / Sol.rho
    kernel_v = np.array([[0.5,1.,0.5],[0.5,1.,0.5]]).T
    flux[1].rhoY[1:-1,1:-1] = signal.convolve2d(rhoYv, kernel_v, mode='valid') /kernel_v.sum()
    # flux[1].rhoY[:,-2] *= 0.


def hll_solver(flux, Lefts, Rights, Sol, lmbda, ud, th):
    # flux: index 1 to end = Left[inner_idx]: index 0 to -1 = Right[inner_idx]: index 1 to end
    remove_cols_idx = (slice(None),slice(1,-1))
    left_idx = (slice(None),slice(0,-1))
    right_idx = (slice(None),slice(1,None))

    Lefts.primitives(th)
    Rights.primitives(th)
    
    upwind = 0.5 * (1.0 + np.sign(flux.rhoY))
    upl = upwind[right_idx]
    upr = (1.0 - upwind[left_idx]) 

    flux.rhou[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.u[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.u[right_idx])
    flux.rho[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * 1.0 + upr[right_idx] / Rights.Y[right_idx] * 1.0)

    Hl = Lefts.rhoe[left_idx] + Lefts.p[left_idx]
    Hr = Rights.rhoe[right_idx] + Rights.p[right_idx]
    flux.rhoe[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Hl + upr[right_idx] / Rights.Y[right_idx] * Hr)

    flux.rhov[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.v[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.v[right_idx])
    flux.rhow[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.w[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.w[right_idx])
    flux.rhoX[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.X[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.X[right_idx])