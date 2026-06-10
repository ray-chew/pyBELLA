`lap2D_manual` now treats RAYLEIGH boundaries as walls, matching `lap3D`,
`lap2D_numba` and the divergence slab-zeroing. Previously a sponged top in a
native-2D run got a periodic-in-y elliptic stencil. No-op for the existing 2D
golden masters (their vertical handling goes through the atmospheric-extension
branch — verified bit-identical, and the stable lamb wave still propagates at
0.99 Cs with stable amplitude at twice the regression horizon).
