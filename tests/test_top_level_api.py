import dtcc_core


def test_top_level_bounds_and_datasets_exports():
    assert dtcc_core.Bounds is dtcc_core.model.Bounds
    assert callable(dtcc_core.datasets.city_flat_mesh)
    assert callable(dtcc_core.datasets.city_surface_mesh)
    assert callable(dtcc_core.datasets.city_volume_mesh)
