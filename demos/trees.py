import dtcc_core as dtcc


bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

trees = dtcc.datasets.trees(
    bounds=bounds,
    tree_type="urban",
)
trees.info()
trees.plot()
