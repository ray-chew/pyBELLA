Routed backend detection through the existing ``backends.is_jax_backend``
helper. The literal ``getattr(ud, "backend", "numpy") in ("jax", "jax-device")``
check was copy-pasted inline across eight hot/parity-path modules (cell/node/
common/rayleigh boundaries, advective flux, diffusion, coriolis, implicit
elliptic); they now call ``is_jax_backend(ud)`` — the single place the
backend-name set is defined. Behaviour-identical (the helper *is* that
expression): numpy bit-identical, jax/jax-device green.
