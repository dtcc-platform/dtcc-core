import dtcc_core


def test_top_level_bounds_datasets_and_city_loader_exports(tmp_path):
    assert dtcc_core.Bounds is dtcc_core.model.Bounds
    assert callable(dtcc_core.datasets.city_flat_mesh)
    assert callable(dtcc_core.datasets.city_surface_mesh)
    assert callable(dtcc_core.datasets.city_volume_mesh)
    loaders = {name for name in vars(dtcc_core.io) if name.startswith("load_")}
    for name in loaders:
        assert name in dtcc_core.__all__
        assert getattr(dtcc_core, name) is getattr(dtcc_core.io, name)
    city = dtcc_core.model.City(id="top-level-city")
    path = tmp_path / "city.dtcc"
    dtcc_core.io.save_model(city, path)
    restored = dtcc_core.load_city(path)
    assert isinstance(restored, dtcc_core.model.City)
    assert restored.id == city.id


def test_datasets_publish_api_exports():
    assert hasattr(dtcc_core.datasets, "DatasetPackageError")
    assert hasattr(dtcc_core.datasets, "DatasetPublication")
    assert hasattr(dtcc_core.datasets, "DatasetPublishConfigurationError")
    assert hasattr(dtcc_core.datasets, "DatasetPublishError")
    assert hasattr(dtcc_core.datasets, "DatasetUploadClient")
    assert hasattr(dtcc_core.datasets, "DatasetUploadConflictError")
    assert hasattr(dtcc_core.datasets, "DatasetUploadError")
    assert hasattr(dtcc_core.datasets, "DatasetUploadInProgressError")
    assert hasattr(dtcc_core.datasets, "DatasetUploadRateLimitError")
    assert hasattr(dtcc_core.datasets, "PublishedFile")
    assert hasattr(dtcc_core.datasets.smoke, "publish")
