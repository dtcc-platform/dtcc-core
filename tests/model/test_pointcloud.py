import numpy as np
import pytest

from dtcc_core.model import Bounds, PointCloud, TimeSeries, exchange


@pytest.fixture
def pc():
    pc = PointCloud()
    pc.points = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    pc.calculate_bounds()
    pc.crs = "EPSG:3857"
    pc.origin = (0, 0)
    pc.classification = np.array([1, 2, 3])
    pc.intensity = np.array([1, 2, 3])
    pc.return_number = np.array([1, 1, 1])
    pc.num_returns = np.array([1, 1, 1])

    return pc


def test_empty():
    pc = PointCloud()
    assert pc.bounds.tuple == (0, 0, 0, 0)


def test_calc_bounds(pc):
    pc.calculate_bounds()
    assert pc.bounds.tuple == (0, 0, 2, 2)


def test_to_proto(pc):

    proto_pc = pc.to_proto()
    assert list(proto_pc.geometry.point_cloud.points.shape) == [3, 3]
    assert list(proto_pc.geometry.point_cloud.classification.shape) == [3]


def test_to_protobuf_missing_fields():
    pc = PointCloud()
    pc.points = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    pc.calculate_bounds()
    proto_pc = pc.to_proto()
    assert list(proto_pc.geometry.point_cloud.points.shape) == [3, 3]
    assert proto_pc.geometry.point_cloud.classification.data == b""


def test_from_proto(pc):
    proto_pc = pc.to_proto()
    pc2 = PointCloud()
    pc2.from_proto(proto_pc)
    assert pc2.points.tolist() == pc.points.tolist()
    assert pc2.classification.tolist() == pc.classification.tolist()
    assert pc2.intensity.tolist() == pc.intensity.tolist()
    assert pc2.return_number.tolist() == pc.return_number.tolist()
    assert pc2.num_returns.tolist() == pc.num_returns.tolist()


if __name__ == "__main__":
    pytest.main()


def _tile(offset, classes):
    count = len(classes)
    tile = PointCloud()
    tile.points = np.column_stack(
        [np.arange(count, dtype=float) + offset, np.zeros(count), np.zeros(count)]
    )
    tile.classification = np.asarray(classes, dtype=np.uint8)
    tile.intensity = np.full(count, 100, dtype=np.uint16)
    tile.return_number = np.ones(count, dtype=np.uint8)
    tile.num_returns = np.ones(count, dtype=np.uint8)
    return tile


@pytest.mark.parametrize(
    "name", ["classification", "intensity", "return_number", "num_returns"]
)
def test_merge_into_empty_point_cloud_keeps_integer_attributes(name):
    tile = _tile(0.0, [1, 2, 9])
    merged = PointCloud()

    merged.merge(tile)

    values = getattr(merged, name)
    assert values.dtype == getattr(tile, name).dtype
    assert values.tolist() == getattr(tile, name).tolist()


@pytest.mark.parametrize(
    "name", ["classification", "intensity", "return_number", "num_returns"]
)
@pytest.mark.parametrize("mutate_source", [True, False])
def test_merge_keeps_attributes_independent(name, mutate_source):
    """Merging into an empty cloud copies attributes instead of sharing them."""
    tile = _tile(0.0, [1, 2, 9])
    merged = PointCloud().merge(tile)

    changed, unchanged = (tile, merged) if mutate_source else (merged, tile)
    expected = getattr(unchanged, name).copy()

    getattr(changed, name)[0] = 99

    np.testing.assert_array_equal(getattr(unchanged, name), expected)


def test_merging_tiles_keeps_classification_integer():
    merged = PointCloud()
    merged.merge(_tile(0.0, [1, 2, 9]))
    merged.merge(_tile(10.0, [2, 2, 7]))

    assert merged.classification.dtype.kind == "u"
    assert merged.classification.tolist() == [1, 2, 9, 2, 2, 7]
    assert len(merged.classification) == len(merged.points)


@pytest.mark.parametrize("operation", ["remove", "keep", "crop", "remove_all"])
def test_point_selection_preserves_time_series(operation):
    cloud = _tile(0, [1, 2, 9])
    times = ["2026-01-01", "2026-01-02"]
    temperatures = np.arange(6, dtype=np.int32).reshape(2, 3)
    vectors = np.arange(18, dtype=np.float32).reshape(2, 3, 3)
    cloud.time_series = [
        TimeSeries("temp", times=times, values=temperatures, association="sample"),
        TimeSeries("wind", times=times, values=vectors, dim=3, association="vertex"),
        TimeSeries("mean", times=times, values=[10, 20], association="geometry"),
    ]
    original_sample = cloud.time_series[0]
    kept = [1, 2]
    if operation == "remove":
        selected = cloud.remove_points([0])
    elif operation == "keep":
        selected = cloud.keep_points(kept)
    elif operation == "crop":
        selected = cloud.crop(Bounds(1, -1, 2, 1))
        assert len(cloud) == 3
    else:
        kept = []
        selected = cloud.remove_points([0, 1, 2])

    restored = exchange.loads(exchange.dumps(selected))
    np.testing.assert_array_equal(restored.points[:, 0], kept)
    sample, vertex, geometry = restored.time_series
    np.testing.assert_array_equal(sample.times, original_sample.times)
    np.testing.assert_array_equal(sample.values, temperatures[:, kept])
    np.testing.assert_array_equal(vertex.values, vectors[:, kept])
    np.testing.assert_array_equal(geometry.values, [[10], [20]])
    np.testing.assert_array_equal(original_sample.values, temperatures)


def test_remove_points_rejects_misaligned_series_before_mutating_cloud():
    cloud = _tile(0, [1, 2, 9])
    cloud.time_series = [
        TimeSeries(
            "temp", times=["2026-01-01"], values=[[1, 2, 3]], association="sample"
        ),
        TimeSeries("wind", times=["2026-01-01"], values=[[4, 5]], association="vertex"),
    ]
    with pytest.raises(ValueError, match="2 elements for 3 points"):
        cloud.remove_points([0])

    np.testing.assert_array_equal(cloud.points[:, 0], [0, 1, 2])
    np.testing.assert_array_equal(cloud.classification, [1, 2, 9])
    np.testing.assert_array_equal(cloud.time_series[0].values, [[1, 2, 3]])
