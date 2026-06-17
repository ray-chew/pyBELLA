Split the 1273-line ``backends/jax_ops/device_step.py`` into ``device_config``
(the host-side ``DeviceConfig`` + laplacian gather/mask-plan builders),
``device_state`` (host↔device 7-leaf-pytree marshalling), and ``device_kernels``
(the traced substeps, the compiled ``make_step``, and the host ``run_window``
loop). ``device_step`` keeps its design docstring and becomes a re-export shim,
so ``device_step.run_window`` is unchanged. No arithmetic moved: jax-device
output is bit-identical to the pre-split run and the inline gate stays green.
