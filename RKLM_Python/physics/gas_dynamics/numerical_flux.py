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

    # keep for debugging purpose... last checked, function is correct
    # will need to think of a better method to generalise from 2d

    # idx = 874
    # print("rhoYu_cm = ", rhoYu.flatten()[idx-52])
    # print("rhoYu_mm = ", rhoYu.flatten()[idx-52-1])
    # print("rhoYu_cc = ", rhoYu.flatten()[idx])
    # print("rhoYu_mc = ", rhoYu.flatten()[idx-1])
    # print("rhoYu_cp = ", rhoYu.flatten()[idx+52])
    # print("rhoYu_mp = ", rhoYu.flatten()[idx+52-1])

    kernel_u = np.array([[0.5, 0.5],[1., 1.],[0.5, 0.5]]).T
    flux[0].rhoY[1:-1,1:-1] = signal.convolve2d(rhoYu, kernel_u, mode='valid').T / kernel_u.sum()
    flux[0].rhoY[:,-2] *= 0.
    # print(flux[0].rhoY.shape)
    # flux[0].rhoY[:,-2] = 0.
    # flux[0].rhoY = flux[0].rhoY.T

    # print("rhoYu = ", flux[0].rhoY.flatten()[idx-5:idx+5])
    # find_nearest(flux[0].rhoY, 1.2420080775988327)

    rhoYv = Sol.rhoY * Sol.rhov / Sol.rho
    kernel_v = np.array([[0.5,1.,0.5],[0.5,1.,0.5]]).T
    flux[1].rhoY[1:-1,1:-1] = signal.convolve2d(rhoYv, kernel_v, mode='valid') /kernel_v.sum()
    flux[1].rhoY[:,-2] *= 0.
    # flux[1].rhoY = flux[1].rhoY.T


def hll_solver(flux, Lefts, Rights, Sol, lmbda, ud, th):
    # flux: index 1 to end = Left[inner_idx]: index 0 to -1 = Right[inner_idx]: index 1 to end
    remove_cols_idx = (slice(None),slice(1,-1))
    left_idx = (slice(None),slice(0,-1))
    right_idx = (slice(None),slice(1,None))

    Lefts.primitives(th)
    Rights.primitives(th)
    # print(Lefts.rhoe[0][:10])
    # print(Lefts.rhoe[0][:10])
    
    upwind = 0.5 * (1.0 + np.sign(flux.rhoY))
    upl = upwind[right_idx]
    upr = (1.0 - upwind[left_idx]) 

    flux.rhou[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.u[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.u[right_idx])
    flux.rho[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * 1.0 + upr[right_idx] / Rights.Y[right_idx] * 1.0)

    Hl = Lefts.rhoe[left_idx] + Lefts.p[left_idx]
    Hr = Rights.rhoe[right_idx] + Rights.p[right_idx]
    flux.rhoe[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Hl + upr[right_idx] / Rights.Y[right_idx] * Hr)

    # print("flux.rhoY = ", flux.rhoY[1][:10])
    # print("upl = ", upl[1][:10] / Lefts.Y[1][:10])
    # # print("Lefts.Y = ", Lefts.Y[1][:10])
    # print("Hl = ", Hl[1][:10])
    # print("upr = ", upr[1][:10] / Rights.Y[1][:10])
    # # print("Rights.Y = ", Rights.Y[1][:10])
    # print("Hr = ", Hr[1][:10])
    # print("flux.rhoe = ", flux.rhoe[1][:10])

    flux.rhov[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.v[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.v[right_idx])
    flux.rhow[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.w[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.w[right_idx])
    flux.rhoX[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.X[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.X[right_idx])