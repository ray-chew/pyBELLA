NEDAS Phase D2: ensemble-vmapped GPU forecasts. `device_batch.run_window_batch`
advances the whole ensemble in lockstep on device — batch-min host-controlled
dt, the full compiled step (bicgstab included) vmapped over the member axis,
optional member-axis sharding across GPUs (`ens_devices`) — wired into the
`ens_run_strategy: batch` branch of `PyBellaModel.run()` for JAX device
backends (`ens_batch_mode: vmap|loop`; `loop` is the same-dt gate comparator).
Gate: `test_scripts/test_nedas_vmap_gate.py`. Design + validation:
dev_notes/nedas_interface.md Phase D.
