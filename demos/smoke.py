import dtcc_core as dtcc


bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

smoke_slice = dtcc.datasets.smoke(
    bounds=bounds,
    product="slice",
    resolution=64,
)
smoke_slice.info()
smoke_slice.plot()

smoke_streamlines = dtcc.datasets.smoke(
    bounds=bounds,
    product="streamlines",
    streamline_count=32,
    streamline_steps=90,
)
smoke_streamlines.info()
smoke_streamlines.plot()
