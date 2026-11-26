import json

import numpy as np
import pytest
from pyproj import Transformer

from dtcc_core.model import (
    City,
    GeometryType,
    LineString,
    MultiLineString,
    PeopleSimulation,
    Person,
    RoadNetwork,
    Transform,
    exchange,
)
from dtcc_core.model.object import TimedPath


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
        [[0, 1, 2, 0], [1, 2, 2, 0], [2, 3, 2, 0], [3, 3, 2, 0], [4, 1, 2, 0]],
        edge_ids=("street",) * 4,
    )
    model.add_child(Person(paths=(path,)))
    counts, endpoints = model.segment_passages()
    assert counts == {(0, "street", 1): 1, (0, "street", -1): 1}
    assert endpoints["street", 1] == ((1, 2), (3, 2))
    assert endpoints["street", -1] == ((3, 2), (1, 2))


def test_edge_counts_join_a_wait_on_the_same_edge_in_wall_clock_time():
    model = PeopleSimulation(attributes={"view_bounds": [-10, -10, 20, 20]})
    model.add_child(
        Person(
            paths=(
                TimedPath(
                    "walk", "work", [[0, 0, 0, 0], [50, 5, 0, 0]], edge_ids=("main",)
                ),
                TimedPath("wait", "work", [[50, 5, 0, 0], [140, 5, 0, 0]]),
                TimedPath(
                    "walk",
                    "work",
                    [[140, 5, 0, 0], [160, 10, 0, 0]],
                    edge_ids=("main",),
                ),
            )
        )
    )
    assert model.segment_passages(60)[0] == {(1, "main", 1): 1}


@pytest.mark.parametrize("network", [False, True])
def test_curved_edge_counts_have_shared_orientation_and_stable_endpoints(network):
    model = PeopleSimulation(attributes={"view_bounds": [-100, -100, 100, 100]})
    xy = np.array([[0, 0], [10, 0], [10, -10], [5, -10]])
    if network:
        roads = RoadNetwork(
            vertices=xy[[0, -1]],
            edges=np.array([[0, 1]]),
            attributes={"id": ["bend"]},
        )
        roads.add_geometry(
            MultiLineString(linestrings=[LineString(vertices=xy)]),
            GeometryType.MULTILINESTRING,
        )
        model.add_child(roads)
    for ident, trace in (("full", xy), ("partial", xy[1:]), ("reverse", xy[-1:0:-1])):
        samples = np.column_stack((np.arange(len(trace)), trace, np.zeros(len(trace))))
        path = TimedPath("walk", "work", samples, edge_ids=("bend",) * (len(trace) - 1))
        model.add_child(Person(id=ident, paths=(path,)))
    counts, endpoints = model.segment_passages()
    assert counts == {(0, "bend", 1): 2, (0, "bend", -1): 1}
    assert endpoints == {
        ("bend", 1): ((0, 0), (5, -10)),
        ("bend", -1): ((5, -10), (0, 0)),
    }


def test_loop_counts_follow_network_order_across_the_seam():
    model = PeopleSimulation(attributes={"view_bounds": [-10, -10, 20, 20]})
    xy = np.array([[0, 0], [0, 5], [5, 5], [5, 0], [0, 0]])
    roads = RoadNetwork(
        vertices=xy[:1],
        edges=np.array([[0, 0]]),
        attributes={"id": 42, "osm_way_id": [1234], "segment_index": [8]},
    )
    roads.add_geometry(
        MultiLineString(linestrings=[LineString(vertices=xy)]),
        GeometryType.MULTILINESTRING,
    )
    model.add_child(roads)
    for ident, trace in (
        ("forward", xy),
        ("reverse", xy[::-1]),
        ("seam", xy[[3, 4, 1]]),
    ):
        samples = np.column_stack((np.arange(len(trace)), trace, np.zeros(len(trace))))
        path = TimedPath(
            "walk", "work", samples, edge_ids=("1234:8",) * (len(trace) - 1)
        )
        model.add_child(Person(id=ident, paths=(path,)))
    assert model.segment_passages()[0] == {
        (0, "1234:8", 1): 2,
        (0, "1234:8", -1): 1,
    }


def test_edge_counts_detect_reversal_with_roundoff_at_projected_northings():
    y = 6_400_000.0
    model = PeopleSimulation(attributes={"view_bounds": [-1, y - 1, 2000, y + 1]})
    path = TimedPath(
        "walk",
        "work",
        [[0, 0, y, 0], [100, 1000, y, 0], [200, 0, np.nextafter(y, np.inf), 0]],
        edge_ids=("street", "street"),
    )
    model.add_child(Person(paths=(path,)))
    assert model.segment_passages()[0] == {(0, "street", 1): 1, (0, "street", -1): 1}


@pytest.mark.parametrize(("first", "last", "second"), [(15, 5, 15), (5, 15, 10)])
def test_passage_along_screenline_counts_the_first_entry(first, last, second):
    model = PeopleSimulation()
    path = TimedPath(
        "walk",
        "work",
        [[0, -1, first, 0], [10, 0, first, 0], [20, 0, last, 0], [30, 1, last, 0]],
    )
    model.add_child(Person(paths=(path,)))
    (passage,) = model.passages((0, 0), (0, 10))
    assert passage.second == second


@pytest.mark.parametrize(
    ("attributes", "error"),
    [
        (
            {
                "screenline": {
                    "name": "line",
                    "start": [0, 0],
                    "end": [0, 10],
                    "bin_seconds": 0,
                }
            },
            "Screenline bin_seconds",
        ),
        (
            {
                "screenline": {
                    "name": "line",
                    "start": [0, 0],
                    "end": [0, 0],
                    "bin_seconds": 300,
                }
            },
            "Screenline requires distinct",
        ),
        (
            {"buildings": [{"id": "building", "ring": [[0, 0], [1, 1]]}]},
            "finite XY ring",
        ),
        ({"buildings": [{"ring": [[0, 0], [1, 0], [1, 1]]}]}, "string id"),
    ],
)
def test_exchange_and_geojson_reject_invalid_playback_geometry(attributes, error):
    model = PeopleSimulation(
        transform=Transform(srs="EPSG:3006"),
        attributes={"view_bounds": [0, 0, 20, 20], **attributes},
    )
    with pytest.raises(ValueError, match=error):
        exchange.dumps(model)
    with pytest.raises(ValueError, match=error):
        model.to_geojson()


@pytest.mark.parametrize("dtype", [np.float64, np.int64])
def test_vehicle_numpy_samples_export_as_json_numbers(dtype):
    samples = np.array([[0, 0, 0, 0], [10, 5, 5, 0]], dtype=dtype)
    model = PeopleSimulation(
        transform=Transform(srs="EPSG:3006"),
        attributes={
            "view_bounds": [0, 0, 20, 20],
            "vehicles": [
                {"id": "trip", "kind": "bus", "line": "1", "samples": samples}
            ],
        },
    )
    expected = json.loads(json.dumps(model.to_geojson()))
    assert expected["people"]["vehicles"][0]["seconds"] == [0, 10]
    if dtype is np.float64:
        model.attributes["vehicles"][0]["samples"] = [list(row) for row in samples]
        assert exchange.loads(exchange.dumps(model)).to_geojson() == expected


def test_simulation_requires_people_in_its_coordinate_frame():
    simulation = PeopleSimulation(
        transform=Transform(srs="EPSG:3006"), attributes={"view_bounds": [0, 0, 20, 20]}
    )
    person = daily_person()
    person.transform.affine[0, 3] = 100
    simulation.add_child(person)
    with pytest.raises(ValueError, match="simulation coordinate frame"):
        exchange.dumps(simulation)


@pytest.mark.parametrize(
    "samples",
    [
        [["0", "0", "0", "0"], ["1", "1", "1", "0"]],
        [[False, 2, 2, 0], [True, 5, 2, 0]],
    ],
)
def test_simulation_rejects_nonnumeric_vehicle_samples(samples):
    track = {
        "id": "route",
        "kind": "bus",
        "line": "1",
        "samples": samples,
    }
    simulation = PeopleSimulation(
        transform=Transform(srs="EPSG:3006"),
        attributes={"view_bounds": [0, 0, 20, 20], "vehicles": [track]},
    )
    with pytest.raises(ValueError, match="samples must be numeric"):
        exchange.dumps(simulation)
    with pytest.raises(ValueError, match="samples must be numeric"):
        simulation.to_geojson()
