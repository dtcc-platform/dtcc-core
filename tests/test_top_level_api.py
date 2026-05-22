import dtcc_core


def test_top_level_bounds_and_datasets_exports():
    assert dtcc_core.Bounds is dtcc_core.model.Bounds
    assert callable(dtcc_core.datasets.city_flat_mesh)
    assert callable(dtcc_core.datasets.city_surface_mesh)
    assert callable(dtcc_core.datasets.city_volume_mesh)


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
