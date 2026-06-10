Terrain-following coordinates, phase 0: `VerticalTransform`/`GalChenTransform` +
`MetricFields` scaffolding (`flow_solver/discretisation/terrain.py`), built in
`grid_init` and attached as `elem.metric`/`node.metric` (`None` without
`ud.orography` — uniform-Cartesian path untouched). Transform unit tests, h≡0
identity oracle, and a quasi-2D mountain-wave smoke case (`smoke_agnesi`)
through the full-tensor 3D elliptic path.
