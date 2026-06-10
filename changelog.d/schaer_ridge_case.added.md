Schär (2002) ridge golden-master case (`test_schaer_ridge`): two-scale
orography (h0 250 m, a 5 km, lambda 4 km) under SLEVE coordinates, native 2D,
with analytic gradients and the exact cos^2 smooth/residual split. New linear
FFT mountain-wave oracle (`schaer_linear_analytic.py`, self-tested against
Smith's closed-form drag) gates the wave field, drag and flux constancy, and
the Gal-Chen-vs-SLEVE discriminator: spurious small-scale w aloft collapses
~17x under SLEVE (E_ss 0.004 vs 0.070) where the true lambda-scale response
is evanescent-dead. Smoke + analytic tests wired into CI.
