import numpy as np
import pytest
from pyproj import Transformer

from dtcc_core.model import (
    City,
    PeopleSimulation,
    Person,
    Transform,
    exchange,
)
from dtcc_core.model.object.people import TimedPath


def daily_person(second=10):
    return Person(
        id="person",
        household_id="household",
        paths=(
            TimedPath(
                "walk",
                "work",
                [[0, 0, 5, 0], [second, 5, 5, 0]],
                journey_id="out",
                edge_ids=("west",),
            ),
            TimedPath(
                "wait", "work", [[second, 5, 5, 0], [20, 5, 5, 0]], journey_id="out"
            ),
            TimedPath(
                "walk",
                "work",
                [[20, 5, 5, 0], [30, 10, 5, 0]],
                journey_id="out",
                edge_ids=("east",),
            ),
        ),
    )


def test_person_native_city_round_trip():
    person = daily_person()
    city = City(id="city")
    city.add_child(person)
    reconstructed = exchange.loads(exchange.dumps(city)).get_children(Person)[0]
    assert reconstructed.id == person.id
    assert reconstructed.paths == person.paths
    assert reconstructed.position_at(10).state == "wait"
    assert reconstructed.position_at(20).heading == pytest.approx(np.pi / 2)
    assert reconstructed.position_at(25).x == 7.5
    assert reconstructed.position_at(30) is None


def test_passage_at_signal_is_counted_once_and_filtered_by_grade():
    model = PeopleSimulation(
        id="people",
        transform=Transform(srs="EPSG:3006"),
        attributes={
            "view_bounds": [0, 0, 20, 20],
            "households": {"household": ["person"]},
        },
    )
    model.add_child(daily_person())
    model = exchange.loads(exchange.dumps(model))
    assert model.household_members("household") == ("person",)
    start, end = (5, 0), (5, 10)
    assert model.passages(start, end, interval=(0, 10)) == []
    (event,) = model.passages(start, end, interval=(10, 30))
    assert (event.person_id, event.second, event.direction, event.edge_id) == (
        "person",
        10,
        1,
        "west",
    )
    assert model.passages(start, end, height_range=(3, 6)) == []
    assert model.passages(start, end, edge_ids={"bridge"}) == []
    counts, _ = model.segment_passages(10)
    assert counts == {(0, "west", 1): 1, (2, "east", 1): 1}


def test_invalid_timed_path_rejected():
    with pytest.raises(ValueError, match="strictly timed"):
        TimedPath("walk", "work", [[4, 0, 0, 0], [4, 1, 0, 0]])
    with pytest.raises(ValueError, match="stationary"):
        TimedPath("wait", "work", [[0, 0, 0, 0], [2, 1, 0, 0]])
    with pytest.raises(ValueError, match="journey_id"):
        TimedPath("walk", "work", [[0, 0, 0, 0], [2, 1, 0, 0]], journey_id=1)


@pytest.mark.parametrize("srs", ["EPSG:4326", "EPSG:2263"])
def test_canonical_simulation_requires_metric_projected_coordinates(srs):
    model = PeopleSimulation(
        id="people",
        transform=Transform(srs=srs),
        attributes={"view_bounds": [0, 0, 20, 20]},
    )
    with pytest.raises(ValueError, match="projected CRS in metres"):
        exchange.dumps(model)


def test_canonical_simulation_rejects_untimed_vehicle_tracks():
    model = PeopleSimulation(
        id="people",
        transform=Transform(srs="EPSG:3006"),
        attributes={
            "view_bounds": [0, 0, 20, 20],
            "vehicles": [
                {
                    "id": "route",
                    "kind": "bus",
                    "line": "19",
                    "samples": [[0, 0, 0, 0], [0, 5, 5, 0]],
                }
            ],
        },
    )
    with pytest.raises(ValueError, match="Scheduled vehicle.*strictly timed"):
        exchange.dumps(model)


def test_animated_plot_updates_pedestrian_positions(monkeypatch):
    import matplotlib.animation
    import matplotlib.pyplot as plt

    model = PeopleSimulation(attributes={"view_bounds": [0, 0, 20, 20]})
    model.add_child(daily_person())

    def animate(figure, update, **options):
        markers = figure.axes[0].collections[-1]
        np.testing.assert_allclose(markers.get_offsets(), [[0, 5]])
        update(25)
        np.testing.assert_allclose(markers.get_offsets(), [[7.5, 5]])
        update(30)
        assert len(markers.get_offsets()) == 0

    monkeypatch.setattr(matplotlib.animation, "FuncAnimation", animate)
    ax = model.plot(second=0, animate=True, show=False)
    plt.close(ax.figure)


def test_playback_derivative_carries_presentation_facts(tmp_path):
    import json

    from dtcc_core import datasets
    from dtcc_core.datasets import load_model_package

    child = Person(
        id="child",
        appearance_key=np.uint32(7),
        attributes={"synthetic_age": 7},
        paths=(
            TimedPath(
                "indoor",
                "home",
                [[0, 319900, 6397900, 0], [60, 319900, 6397900, 0]],
                "osm-way-1",
            ),
            TimedPath(
                "walk",
                "education",
                [[60, 319900, 6397900, 0], [120, 319950, 6397900, 0]],
            ),
        ),
    )
    simulation = PeopleSimulation(
        id="people",
        transform=Transform(srs="EPSG:3006"),
        attributes={
            "view_bounds": [319720, 6397660, 320220, 6398160],
            "screenline": {
                "name": "test line",
                "start": [319925, 6397890],
                "end": [319925, 6397910],
                "bin_seconds": 300,
            },
            "buildings": [
                {
                    "id": "osm-way-1",
                    "ring": [[319890, 6397890], [319910, 6397890], [319910, 6397910]],
                    "height": 12.0,
                    "height_source": "mapped levels",
                    "kind": "home",
                }
            ],
            "vehicles": [
                {
                    "id": "trip-1",
                    "kind": "tram",
                    "line": "7",
                    "samples": [[0, 319800, 6397800, 0], [90, 320100, 6397800, 0]],
                }
            ],
            "daylight": {"date": "2019-10-15", "basis": "transit service date"},
        },
    )
    simulation.add_child(child)
    derived = simulation.to_geojson()
    indoor, walk = derived["features"]
    assert indoor["properties"]["place"] == "osm-way-1"
    assert "place" not in walk["properties"]
    assert walk["properties"]["child"] is True
    building = derived["people"]["buildings"][0]
    assert (building["height"], building["height_source"]) == (12.0, "mapped levels")
    vehicle = derived["people"]["vehicles"][0]
    assert vehicle["seconds"] == [0, 90]
    assert vehicle["coordinates"][0][0] == pytest.approx(11.98, abs=0.02)
    assert vehicle["coordinates"][0][1] == pytest.approx(57.69, abs=0.02)
    assert derived["people"]["daylight"]["date"] == "2019-10-15"
    line = derived["people"]["screenline"]
    assert line["start"] == pytest.approx(
        Transformer.from_crs(3006, 4326, always_xy=True).transform(319925, 6397890)
    )
    assert line["counts"] == {"0": 1}
    with pytest.raises(ValueError, match="EPSG:4326"):
        simulation.to_geojson(crs="EPSG:3006")
    restored = exchange.loads(exchange.dumps(simulation))
    assert restored.attributes["vehicles"] == simulation.attributes["vehicles"]
    descriptor = datasets.calibration_grid
    simulation.dataset_context = descriptor.create_context(
        descriptor.validate({"bounds": simulation.attributes["view_bounds"]})
    )
    native = simulation.export(tmp_path / "native", canonical=True)
    (canonical,) = native.manifest.artifacts
    assert canonical.role == "canonical_model"
    package = simulation.export(tmp_path / "playback", canonical=True, format="geojson")
    artifacts = {item.role: item for item in package.manifest.artifacts}
    primary, playback = artifacts["canonical_model"], artifacts["derived"]
    assert primary.sha256 == canonical.sha256
    assert playback.role == "derived" and playback.derived_from == primary.path
    assert playback.crs == "EPSG:4326"
    restored = load_model_package(package.path)
    assert (
        json.loads((package.path / playback.path).read_text()) == restored.to_geojson()
    )


def test_edge_counts_preserve_directed_traversals_across_intermediate_samples():
    model = PeopleSimulation(attributes={"view_bounds": [0, 0, 20, 20]})
    path = TimedPath(
        "walk",
        "work",
        [[0, 1, 2, 0], [1, 2, 2, 0], [2, 3, 2, 0], [4, 1, 2, 0]],
        edge_ids=("street",) * 3,
    )
    model.add_child(Person(paths=(path,)))
    counts, endpoints = model.segment_passages()
    assert counts == {(0, "street", 1): 1, (0, "street", -1): 1}
    assert endpoints["street", 1] == ((1, 2), (3, 2))
    assert endpoints["street", -1] == ((3, 2), (1, 2))


def test_simulation_requires_people_in_its_coordinate_frame():
    simulation = PeopleSimulation(
        transform=Transform(srs="EPSG:3006"), attributes={"view_bounds": [0, 0, 20, 20]}
    )
    person = daily_person()
    person.transform.affine[0, 3] = 100
    simulation.add_child(person)
    with pytest.raises(ValueError, match="simulation coordinate frame"):
        exchange.dumps(simulation)


def test_simulation_rejects_string_vehicle_coordinates():
    track = {
        "id": "route",
        "kind": "bus",
        "line": "1",
        "samples": [["0", "0", "0", "0"], ["1", "1", "1", "0"]],
    }
    simulation = PeopleSimulation(
        transform=Transform(srs="EPSG:3006"),
        attributes={"view_bounds": [0, 0, 20, 20], "vehicles": [track]},
    )
    with pytest.raises(ValueError, match="samples must be numeric"):
        exchange.dumps(simulation)
    with pytest.raises(ValueError, match="samples must be numeric"):
        simulation.to_geojson()
