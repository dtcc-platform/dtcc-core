"""Plot trees inferred from the point cloud."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)

trees = dtcc.datasets.trees(
    bounds=bounds,
    tree_type="urban",
)
trees.info()
trees.plot()
