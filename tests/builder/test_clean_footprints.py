"""The single footprint cleaning entry point for buildings and cities."""

import numpy as np
import pytest
from shapely.geometry import Polygon, box

from dtcc_core import io
from dtcc_core.builder import UnresolvedFootprintCleaningError, clean_footprints
from dtcc_core.builder.cleaning import clean as clean_module
from dtcc_core.builder.cleaning.blocks import block_walls
from dtcc_core.builder.cleaning.clean import (
    MERGED_FROM_KEY,
    RECORD_KEY,
    cleaning_record,
    footprint_fingerprint,
)
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.model import (
    Bounds,
    Building,
    BuildingPart,
    City,
    GeometryType,
    MultiSurface,
    PointCloud,
    Raster,
    Surface,
    Terrain,
)


def building(id, polygon, height, *, measured=None, ground=0.0, **attributes):
    surface = Surface()
    surface.from_polygon(polygon, ground + height)
    value = Building(id=id)
    value.add_geometry(surface, GeometryType.LOD0)
    value.estimated_height = float(height)
    value.attributes["ground_height"] = float(ground)
    if measured is not None:
        value.measured_height = measured
    value.attributes.update(attributes)
    return value


def city_of(*buildings, terrain=True):
    city = City(id="city")
    if terrain:
        raster = Raster()
        raster.data = np.zeros((4, 4))
        raster.set_bounds(Bounds(-20, -20, 80, 60))
        value = Terrain()
        value.add_geometry(raster, GeometryType.RASTER)
        city.add_terrain(value)
    city.add_buildings(list(buildings))
    return city


def footprint(value):
    return value.get_geometry(GeometryType.LOD0).to_polygon(simplify=0.0)


def test_touching_buildings_become_a_block_whose_parts_are_the_originals():
    low = building("low", box(10, 10, 30, 25), 6, measured=6.5, name="A")
    high = building("high", box(30, 10, 40, 25), 30, name="B")
    city = city_of(low, high)

    cleaned = city.clean_footprints()

    assert [b.id for b in city.buildings] == ["low", "high"]  # input unchanged
    assert len(cleaned.buildings) == 1
    block = cleaned.buildings[0]
    assert block.id not in {"low", "high"}
    assert block.attributes[MERGED_FROM_KEY] == ["low", "high"]
    parts = block.building_parts
    assert [type(p) for p in parts] == [BuildingPart, BuildingPart]
    assert [p.id for p in parts] == ["low", "high"]
    assert [p.attributes["name"] for p in parts] == ["A", "B"]
    assert [p.measured_height for p in parts] == [6.5, None]
    assert [p.estimated_height for p in parts] == [6, 30]
    # Each part keeps its own footprint, at its own roof.
    assert footprint(parts[0]).equals(box(10, 10, 30, 25))
    assert footprint(parts[1]).equals(box(30, 10, 40, 25))
    assert [p.get_geometry(GeometryType.LOD0).zmax for p in parts] == [6, 30]
    # The block's roof and height are those of its highest part.
    assert block.estimated_height == 30
    assert block.measured_height is None
    assert block.get_geometry(GeometryType.LOD0).zmax == 30
    assert footprint(block).equals(box(10, 10, 40, 25))
    assert block.footprint().to_polygon(simplify=0.0).equals(box(10, 10, 40, 25))
    record = cleaned.attributes[RECORD_KEY]
    assert record["merged_count"] == 1
    assert record["walls"]["strict"] == 1
    assert cleaned.terrain is not city.terrain
    assert cleaned.terrain.raster is not None


def test_block_has_the_roof_of_its_highest_part_and_the_lowest_ground():
    city = city_of(
        building("a", box(0, 0, 20, 10), 10, ground=1.0),
        building("b", box(20.2, 0, 30.2, 10), 11),
    )
    block = city.clean_footprints().buildings[0]
    assert block.estimated_height == 11
    assert block.get_geometry(GeometryType.LOD0).zmax == 11
    assert block.attributes["ground_height"] == 0


def test_separate_buildings_keep_identity_order_and_non_footprint_geometry():
    first = building("first", box(20, 0, 30, 12), 12, kind="office")
    second = building("second", box(0, 0, 10, 10), 9)
    roof = PointCloud(points=np.array([[5.0, 5.0, 9.0]]))
    second.add_geometry(roof, GeometryType.POINT_CLOUD)
    lod1 = Surface()
    lod1.from_polygon(box(0, 0, 10, 10), 9)
    second.add_geometry(lod1, GeometryType.LOD1)

    cleaned = clean_footprints([first, second])

    assert [b.id for b in cleaned] == ["first", "second"]
    assert cleaned[0].attributes["kind"] == "office"
    assert cleaned[1].get_geometry(GeometryType.POINT_CLOUD) is not None
    assert cleaned[1].get_geometry(GeometryType.LOD1) is None
    assert footprint(cleaned[1]).equals(box(0, 0, 10, 10))
    assert cleaned[0] is not first and cleaned[1] is not second
    assert second.get_geometry(GeometryType.LOD1) is not None


def test_isolated_cleaning_never_merges():
    a = building("a", box(0, 0, 10, 10), 10)
    b = building("b", box(9.8, 0, 20, 10), 20)  # overlaps a
    city = city_of(a, b)

    cleaned = city.clean_footprints(isolated=True)

    assert [x.id for x in cleaned.buildings] == ["a", "b"]
    assert all(not x.building_parts for x in cleaned.buildings)
    assert footprint(cleaned.buildings[0]).intersects(footprint(cleaned.buildings[1]))
    assert cleaned.attributes[RECORD_KEY]["isolated"] is True
    # Isolated results are no valid input for meshing.
    assert cleaning_record(cleaned, min_building_detail=0.5) is None


def test_zero_detail_keeps_overlaps_and_warns_of_them(caplog):
    tower = building("tower", box(0, 0, 10, 10), 60)
    house = building("house", box(5, 0, 40, 10), 6)  # overlaps the tower
    shed = building("shed", box(40, 0, 45, 10), 3)  # shares an edge with the house
    city = city_of(tower, house, shed, terrain=False)

    cleaned = city.clean_footprints(min_building_detail=0)

    assert [x.id for x in cleaned.buildings] == ["tower", "house", "shed"]
    assert footprint(cleaned.buildings[0]).overlaps(footprint(cleaned.buildings[1]))
    overlaps = [r.getMessage() for r in caplog.records if "overlap" in r.getMessage()]
    assert len(overlaps) == 1
    assert "1 pair(s)" in overlaps[0] and "tower/house" in overlaps[0]
    assert "shed" not in overlaps[0]

    # Buildings cleaned in isolation may overlap; nothing to warn of.
    caplog.clear()
    city.clean_footprints(min_building_detail=0, isolated=True)
    assert "overlap" not in caplog.text


def test_strict_zero_detail_cleaning_does_not_accept_overlaps():
    tower = building("tower", box(0, 0, 10, 10), 60)
    house = building("house", box(5, 0, 40, 10), 6)  # overlaps the tower
    shed = building("shed", box(40, 0, 45, 10), 3)  # shares an edge with the house

    with pytest.raises(UnresolvedFootprintCleaningError, match="tower/house") as error:
        clean_footprints([tower, house, shed], min_building_detail=0, strict=True)

    (group,) = error.value.diagnostics["unresolved_groups"]
    assert (group["reason"], group["source_indices"]) == ("overlap", [0, 1])
    # Touching is no overlap, and buildings cleaned in isolation may overlap.
    assert len(clean_footprints([house, shed], min_building_detail=0, strict=True)) == 2
    cleaned = clean_footprints(
        [tower, house], min_building_detail=0, strict=True, isolated=True
    )
    assert len(cleaned) == 2


def test_min_area_and_missing_footprints_are_reported_by_id():
    city = city_of(
        building("big", box(0, 0, 10, 10), 10),
        building("tiny", box(20, 0, 22, 2), 5),
        Building(id="empty"),
    )
    cleaned = city.clean_footprints(min_area=10)
    record = cleaned.attributes[RECORD_KEY]
    assert [b.id for b in cleaned.buildings] == ["big"]
    assert record["missing_footprint_ids"] == ["empty"]
    assert record["excluded"] == [
        {"ids": ["tiny"], "area": pytest.approx(4.0), "reason": "minimum_area"}
    ]


def test_a_building_without_an_occupied_core_is_reported_when_it_vanishes():
    # A 0.3 m triangle next to a house has no core at epsilon = 0.25 m.
    from shapely.geometry import Polygon

    city = city_of(
        building("house", box(0, 0, 10, 10), 10),
        building("speck", Polygon([(10.3, 2), (10.6, 2), (10.3, 2.3)]), 5),
    )
    cleaned = city.clean_footprints()
    assert [b.id for b in cleaned.buildings] == ["house"]
    assert cleaned.attributes[RECORD_KEY]["excluded"] == [
        {"ids": ["speck"], "area": pytest.approx(0.045), "reason": "no_protected_core"}
    ]


def test_cleaning_a_cleaned_city_changes_nothing():
    city = city_of(
        building("a", box(0, 0, 10, 10), 10),
        building("b", box(10, 0, 20, 10), 20),
        building("c", box(30, 0, 40, 10), 15),
    )
    once = city.clean_footprints()
    twice = once.clean_footprints()
    assert [b.id for b in twice.buildings] == [b.id for b in once.buildings]
    assert all(
        footprint(x).equals_exact(footprint(y), 0.0)
        for x, y in zip(once.buildings, twice.buildings)
    )
    assert [b.estimated_height for b in twice.buildings] == [
        b.estimated_height for b in once.buildings
    ]
    assert twice.attributes[RECORD_KEY]["outcome"] == "unchanged"
    parts = [(p.id, footprint(p)) for b in once.buildings for p in b.building_parts]
    again = [(p.id, footprint(p)) for b in twice.buildings for p in b.building_parts]
    assert [i for i, _ in again] == [i for i, _ in parts] == ["a", "b"]
    assert all(x.equals_exact(y, 0.0) for (_, x), (_, y) in zip(parts, again))


def test_merging_blocks_again_lifts_their_parts():
    city = city_of(
        building("a", box(0, 0, 10, 10), 10),
        building("b", box(10, 0, 20, 10), 10),
        building("c", box(20.8, 0, 30, 10), 10),
    )
    first = city.clean_footprints(min_building_detail=0.5)
    assert [len(b.building_parts) for b in first.buildings] == [2, 0]
    second = first.clean_footprints(min_building_detail=1.0)
    assert len(second.buildings) == 1
    block = second.buildings[0]
    assert [p.id for p in block.building_parts] == ["a", "b", "c"]
    assert block.attributes[MERGED_FROM_KEY] == ["a", "b", "c"]
    assert all(not p.building_parts for p in block.building_parts)
    # The new block is divided by the old parts' footprints and c; the gap
    # between b and c is split in the middle.
    pieces = [footprint(p) for p in block.building_parts]
    expected = [box(0, 0, 10, 10), box(10, 0, 20.4, 10), box(20.4, 0, 30, 10)]
    assert all(p.symmetric_difference(q).area < 1e-9 for p, q in zip(pieces, expected))


def test_building_shortcut_cleans_one_building_in_isolation():
    value = building("one", box(0, 0, 10, 10), 10)
    cleaned = value.clean_footprint()
    assert isinstance(cleaned, Building)
    assert cleaned.id == "one" and cleaned is not value
    assert footprint(cleaned).equals(box(0, 0, 10, 10))
    assert clean_footprints(value, min_area=200) is None


def test_holes_and_courtyard_buildings_are_preserved_by_default():
    courtyard = box(0, 0, 20, 20).difference(box(5, 5, 15, 15))
    values = [
        building("outer", courtyard, 10),
        building("inner", box(8, 8, 12, 12), 30),
    ]

    cleaned = clean_footprints(values)

    assert [value.id for value in cleaned] == ["outer", "inner"]
    assert footprint(cleaned[0]).equals(courtyard)


@pytest.mark.parametrize("delta", [0.0, 0.5])
@pytest.mark.parametrize("isolated", [False, True])
@pytest.mark.parametrize("as_city", [False, True])
def test_remove_holes_drops_nested_polygons_before_merging(delta, isolated, as_city):
    outer = Polygon(
        box(0, 0, 30, 20).exterior,
        [box(5, 5, 15, 15).exterior, box(20, 5, 25, 15).exterior],
    )
    middle = box(5.1, 5.1, 14.9, 14.9).difference(box(8, 8, 12, 12))
    values = [
        building("nested", box(9, 9, 11, 11), 50),
        building("other_hole", box(21, 6, 24, 14), 40),
        building("outside", box(40, 0, 50, 10), 20),
        building("middle", middle, 30),
        building("outer", outer, 10),
    ]
    target = city_of(*values) if as_city else values
    options = dict(min_building_detail=delta, isolated=isolated, remove_holes=True)

    result = (
        target.clean_footprints(**options)
        if as_city
        else clean_footprints(target, **options)
    )
    cleaned = result.buildings if as_city else result

    assert [value.id for value in cleaned] == ["outside", "outer"]
    assert footprint(cleaned[1]).equals(box(0, 0, 30, 20))
    assert cleaned[1].estimated_height == 10
    assert not cleaned[1].building_parts
    assert footprint(values[4]).equals(outer)
    assert footprint(values[3]).equals(middle)
    if as_city:
        record = result.attributes[RECORD_KEY]
        assert record["remove_holes"] is True
        assert {row["ids"][0] for row in record["excluded"]} == {
            "nested",
            "other_hole",
            "middle",
        }
        assert all(row["reason"] == "inside_hole" for row in record["excluded"])


def test_remove_holes_building_shortcut_uses_filled_area_for_selection():
    polygon = box(0, 0, 20, 20).difference(box(2, 2, 18, 18))
    value = building("courtyard", polygon, 10)

    cleaned = value.clean_footprint(remove_holes=True)
    selected = clean_footprints(value, remove_holes=True, min_area=300)

    assert footprint(cleaned).equals(box(0, 0, 20, 20))
    assert footprint(selected).equals(box(0, 0, 20, 20))
    assert clean_footprints(value, min_area=300) is None
    assert footprint(value).equals(polygon)


def test_remove_holes_preserves_components_outside_a_hole():
    polygon = box(0, 0, 20, 20).difference(box(5, 5, 15, 15))
    value = Building(id="multipart")
    surfaces = []
    for piece in [polygon, box(7, 7, 10, 10), box(30, 0, 40, 10)]:
        surface = Surface()
        surface.from_polygon(piece, 10)
        surfaces.append(surface)
    value.add_geometry(MultiSurface(surfaces=surfaces), GeometryType.LOD0)

    cleaned = clean_footprints([value], remove_holes=True, min_building_detail=0)

    assert len(cleaned) == 2
    assert all(not footprint(item).interiors for item in cleaned)
    assert sum(footprint(item).area for item in cleaned) == pytest.approx(500)
    assert len(value.get_geometry(GeometryType.LOD0).surfaces) == 3


def test_remove_holes_keeps_polygon_crossing_hole_boundary_in_isolation():
    polygon = box(0, 0, 20, 20).difference(box(5, 5, 15, 15))
    values = [
        building("outer", polygon, 10),
        building("crossing", box(10, 7, 17, 13), 20),
    ]

    cleaned = clean_footprints(values, remove_holes=True, isolated=True)

    assert [value.id for value in cleaned] == ["outer", "crossing"]
    assert footprint(cleaned[1]).equals(box(10, 7, 17, 13))


@pytest.mark.parametrize("preclean", [False, True])
def test_remove_holes_fills_holes_formed_by_merging_and_drops_their_islands(preclean):
    values = [
        building("bottom", box(0, 0, 20, 5), 10),
        building("top", box(0, 15, 20, 20), 10),
        building("left", box(0, 5, 5, 15), 10),
        building("right", box(15, 5, 20, 15), 10),
        building("island", box(8, 8, 12, 12), 40),
    ]

    if preclean:
        values = clean_footprints(values)
    result = city_of(*values).clean_footprints(remove_holes=True)
    cleaned = result.buildings

    assert len(cleaned) == 1
    assert footprint(cleaned[0]).equals(box(0, 0, 20, 20))
    assert "island" not in cleaned[0].attributes[MERGED_FROM_KEY]
    assert result.attributes[RECORD_KEY]["excluded"] == [
        {"ids": ["island"], "area": 16.0, "reason": "inside_hole"}
    ]
    assert all(
        not footprint(part).interiors
        for part in cleaned[0].building_parts
        if part.get_geometry(GeometryType.LOD0) is not None
    )


@pytest.mark.parametrize("value", [None, 0, 1, "yes"])
def test_remove_holes_requires_a_boolean(value):
    with pytest.raises(TypeError, match="remove_holes must be a boolean"):
        clean_footprints([], remove_holes=value)


def test_record_tracks_footprints_but_not_heights(tmp_path):
    cleaned = city_of(
        building("a", box(0, 0, 10, 10), 10),
        building("b", box(10, 0, 20, 10), 20),
    ).clean_footprints()
    record = cleaning_record(cleaned, min_building_detail=0.5)
    assert record is not None
    assert record["fingerprint"] == footprint_fingerprint(cleaned.buildings)
    assert cleaning_record(cleaned, min_building_detail=1.0) is None
    assert (
        cleaning_record(cleaned, min_building_detail=0.5, lod=GeometryType.LOD1) is None
    )

    cleaned.buildings[0].get_geometry(GeometryType.LOD0).set_z(42.0)
    assert cleaning_record(cleaned, min_building_detail=0.5) is not None

    path = tmp_path / "cleaned.dtcc"
    io.save_model(cleaned, path)
    loaded = io.load_model(path)
    assert cleaning_record(loaded, min_building_detail=0.5) == record
    assert [p.id for p in loaded.buildings[0].building_parts] == ["a", "b"]

    loaded.add_building(building("new", box(50, 0, 60, 10), 10))
    assert cleaning_record(loaded, min_building_detail=0.5) is None

    # The footprints of parts are checked by the walls record of their block:
    # changed ones are recomputed from the sources.
    moved = io.load_model(path)
    block = moved.buildings[0]
    block.building_parts[0].get_geometry(GeometryType.LOD0).translate(x=0.1)
    assert cleaning_record(moved, min_building_detail=0.5) is not None
    walls = block_walls(block)
    assert not walls.current
    assert [
        cell.equals(box(0, 0, 10, 10) if part == 0 else box(10, 0, 20, 10))
        for cell, part in walls.cells
    ] == [True, True]


def test_strict_cleaning_does_not_accept_residual_separation(monkeypatch):
    seen = []
    original = clean_module.condition_polygon_coverage

    def spy(raw, *, source_map, options):
        seen.append(options.allow_residual_separation)
        return original(raw, source_map=source_map, options=options)

    monkeypatch.setattr(clean_module, "condition_polygon_coverage", spy)
    clean_footprints([building("a", box(0, 0, 10, 10), 10)])
    clean_footprints([building("a", box(0, 0, 10, 10), 10)], strict=True)
    assert seen == [True, False]


def test_mesh_builders_reuse_a_recorded_cleaning(monkeypatch):
    city = city_of(
        building("a", box(0, 0, 10, 10), 10),
        building("b", box(10, 0, 20, 10), 20),
        building("c", box(30, 0, 40, 10), 15),
        building("tiny", box(50, 0, 52, 2), 5),
    )
    cleaned = city.clean_footprints()
    expected = meshes.build_city_flat_mesh(city, report_mesh_quality=False)

    def forbidden(*args, **kwargs):
        raise AssertionError("cleaned twice")

    monkeypatch.setattr(meshes, "clean_regions", forbidden)
    mesh = meshes.build_city_flat_mesh(cleaned, report_mesh_quality=False)
    assert np.array_equal(mesh.vertices, expected.vertices)
    assert np.array_equal(mesh.markers, expected.markers)
    # Markers are indices of buildings in the cleaned city; the tiny one is
    # dropped by min_building_area.
    assert sorted(set(mesh.markers[mesh.markers >= 0].tolist())) == [0, 1]
    surface = meshes.build_city_surface_mesh(cleaned, report_mesh_quality=False)
    assert surface.vertices[:, 2].max() == pytest.approx(20)
