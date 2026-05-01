import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg", force=True)
plt = pytest.importorskip("matplotlib.pyplot")

from dtcc_core.model import GeometryType, LineString, MultiLineString, RoadNetwork


def _graph_roadnetwork():
    roadnetwork = RoadNetwork()
    roadnetwork.vertices = np.array([(0, 0), (1, 0), (1, 1)], dtype=float)
    roadnetwork.edges = np.array([(0, 1), (1, 2)], dtype=np.int64)
    roadnetwork.length = np.array([1.0, 2.0])
    roadnetwork.attributes = {
        "maxspeed_kmh": [30, 50],
        "highway": ["residential", "primary"],
    }
    return roadnetwork


def _geometry_roadnetwork():
    roadnetwork = RoadNetwork()
    roadnetwork.vertices = np.array([(0, 0), (1, 0)], dtype=float)
    roadnetwork.edges = np.array([(0, 1)], dtype=np.int64)
    roadnetwork.length = np.array([1.0])
    roadnetwork.attributes = {"maxspeed_kmh": [30]}

    line = LineString()
    line.vertices = np.array([(0, 0), (1, 0), (2, 1)], dtype=float)
    multilinestring = MultiLineString()
    multilinestring.linestrings.append(line)
    roadnetwork.geometry[GeometryType.MULTILINESTRING] = multilinestring
    return roadnetwork


def test_roadnetwork_plot_graph_only():
    roadnetwork = _graph_roadnetwork()

    ax = roadnetwork.plot(show=False)

    assert len(ax.collections) == 1
    assert len(ax.collections[0].get_segments()) == 2
    plt.close(ax.figure)


def test_roadnetwork_plot_geometry_numeric_column():
    roadnetwork = _geometry_roadnetwork()

    ax = roadnetwork.plot(column="maxspeed_kmh", show=False)

    assert len(ax.collections) == 1
    assert len(ax.collections[0].get_segments()) == 1
    assert len(ax.figure.axes) == 2
    plt.close(ax.figure)


def test_roadnetwork_plot_categorical_column():
    roadnetwork = _graph_roadnetwork()

    ax = roadnetwork.plot(column="highway", show=False)

    assert len(ax.collections) == 1
    assert ax.get_legend() is not None
    plt.close(ax.figure)


def test_roadnetwork_plot_missing_column_raises():
    roadnetwork = _graph_roadnetwork()

    with pytest.raises(KeyError, match="missing"):
        roadnetwork.plot(column="missing", show=False)


def test_roadnetwork_plot_misaligned_column_raises():
    roadnetwork = _graph_roadnetwork()
    roadnetwork.attributes["bad"] = [1]

    with pytest.raises(ValueError, match="1 value"):
        roadnetwork.plot(column="bad", show=False)
