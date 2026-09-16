import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg", force=True)
plt = pytest.importorskip("matplotlib.pyplot")

from dtcc_core.model.geometry import Point
from dtcc_core.model.object import Object, VehicleCollection
from dtcc_core.plotting.live import LiveVehiclePlot, vehicle_positions


def _vehicles(*items):
    collection = VehicleCollection()
    collection.attributes = {"retrieval_time": "2026-05-08T10:00:00+00:00"}
    for vehicle_id, x, y in items:
        vehicle = Object()
        vehicle.attributes = {"vehicle_id": vehicle_id, "mode": "bus"}
        vehicle.add_geometry(Point(x=x, y=y, z=0.0), "location")
        collection.add_vehicle(vehicle)
    return collection


def test_vehicle_positions_extracts_keyed_points():
    positions = vehicle_positions(_vehicles(("bus-1", 1.0, 2.0)))

    assert positions == {"bus-1": (1.0, 2.0)}


def test_live_vehicle_plot_updates_current_markers_and_trails():
    plotter = LiveVehiclePlot(bounds=(0.0, 0.0, 10.0, 10.0), stale_after_s=5.0)

    plotter.update(_vehicles(("bus-1", 1.0, 2.0)), now=0.0)
    plotter.update(
        _vehicles(("bus-1", 2.0, 3.0), ("bus-2", 4.0, 5.0)),
        now=1.0,
    )

    assert set(plotter.trails) == {"bus-1", "bus-2"}
    assert list(plotter.trails["bus-1"]) == [(1.0, 2.0), (2.0, 3.0)]
    assert plotter.current.get_offsets().shape == (2, 2)
    plt.close(plotter.fig)


def test_live_vehicle_plot_prunes_stale_trails():
    plotter = LiveVehiclePlot(bounds=(0.0, 0.0, 10.0, 10.0), stale_after_s=1.0)

    plotter.update(_vehicles(("bus-1", 1.0, 2.0)), now=0.0)
    plotter.update(_vehicles(), now=2.0)

    assert plotter.trails == {}
    assert plotter.trail_lines == {}
    plt.close(plotter.fig)


def test_animate_vehicle_collection_polls_and_returns_updated_plot(monkeypatch):
    from dtcc_core.plotting.live import animate_vehicle_collection

    fetched = []
    pauses = []

    def fetch():
        fetched.append(True)
        return _vehicles(("bus-1", 1.0, 2.0))

    def pause(interval):
        pauses.append(interval)
        plt.close(plt.gcf())

    monkeypatch.setattr(plt, "pause", pause)
    plotter = animate_vehicle_collection(
        fetch, interval_s=0.125, bounds=(0, 0, 10, 10), trail_length=3
    )
    try:
        assert fetched == [True]
        assert pauses == [0.125]
        assert list(plotter.trails["bus-1"]) == [(1.0, 2.0)]
        assert plotter.current.get_offsets().tolist() == [[1.0, 2.0]]
        assert not plt.fignum_exists(plotter.fig.number)
    finally:
        plt.close(plotter.fig)
