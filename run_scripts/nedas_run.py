"""Launch a NEDAS scheme with the pyBELLA adapter registered.

Usage (from the repository root, env with `pip install "pybella[nedas]"`):

    python run_scripts/nedas_run.py -c run_scripts/nedas_tv_osse.yml

Equivalent to `python -m NEDAS -c <config>` plus (a) registering
PyBellaModel/PyBellaObs in the NEDAS registries first and (b) dumping ALL
in-memory prior/post snapshots to npy checkpoints at the end (the filter
scheme itself only checkpoints up to the second-to-last cycle).
"""

import sys

import pybella.interfaces.nedas as pybella_nedas


def main() -> None:
    pybella_nedas.register()

    from NEDAS import get_scheme

    scheme = get_scheme(parse_args=True)

    step = scheme.config.step
    if step:
        scheme.run_step(step)
        return

    scheme()

    # checkpoint every cycle's snapshots (time=None -> all times in memory)
    for _, model in scheme.c.models.items():
        for tag in ("prior", "post"):
            model.save_memory(tag)


if __name__ == "__main__":
    sys.exit(main())
