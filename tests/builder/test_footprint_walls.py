"""Walls between merged buildings: every part keeps a cleaned footprint."""

import json
from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import box, shape
from shapely.ops import unary_union

from dtcc_core.builder.cleaning import clean as clean_module
from dtcc_core.builder.cleaning.clean import MERGED_FROM_KEY, RECORD_KEY, clean_regions
from dtcc_core.builder.cleaning.contract import check_labelled_contract
from dtcc_core.builder.cleaning.walls import LabelledRegion, label_region
from dtcc_core.builder.geometry_builders.buildings import (
    compute_building_heights,
    extract_roof_points,
    set_building_heights_from_attribute,
)
from dtcc_core.model import (
    Bounds,
    Building,
    City,
    GeometryType,
    MultiSurface,
    PointCloud,
    Raster,
    Surface,
)

DATA = Path(__file__).resolve().parents[1] / "data" / "cleaning"


def building(id, polygon, height, **attributes):
    surface = Surface()
    surface.from_polygon(polygon, height)
    value = Building(id=id)
    value.add_geometry(surface, GeometryType.LOD0)
    value.estimated_height = float(height)
    value.attributes["ground_height"] = 0.0
    value.attributes.update(attributes)
    return value


def city_of(*buildings):
    city = City(id="city")
    city.add_buildings(list(buildings))
    return city


def footprint(obj):
    geometry = obj.get_geometry(GeometryType.LOD0)
    if isinstance(geometry, MultiSurface):
        return unary_union([s.to_polygon(simplify=0.0) for s in geometry.surfaces])
    return geometry.to_polygon(simplify=0.0)


def parts_by_id(block):
    return {part.id: part for part in block.building_parts}


def assert_parts_tile_block(block):
    pieces = [footprint(p) for p in block.building_parts if p.get_geometry(GeometryType.LOD0)]
    outline = footprint(block)
    assert sum(p.area for p in pieces) == pytest.approx(outline.area, rel=1e-12)
    assert unary_union(pieces).symmetric_difference(outline).area < 1e-9


def test_overlap_goes_to_the_taller_building_then_the_larger():
    low = building("low", box(0, 0, 10, 10), 10)
    high = building("high", box(9.8, 0, 20, 10), 20)
    block = city_of(low, high).clean_footprints().buildings[0]
    parts = parts_by_id(block)
    assert footprint(parts["high"]).equals(box(9.8, 0, 20, 10))
    assert footprint(parts["low"]).equals(box(0, 0, 9.8, 10))
    assert_parts_tile_block(block)

    # Equal heights: the larger footprint wins.
    small = building("small", box(0, 0, 10, 10), 10)
    large = building("large", box(9.8, 0, 22, 10), 10)
    block = city_of(small, large).clean_footprints().buildings[0]
    parts = parts_by_id(block)
    assert footprint(parts["large"]).equals(box(9.8, 0, 22, 10))
    assert footprint(parts["small"]).equals(box(0, 0, 9.8, 10))


def test_gap_between_merged_buildings_is_divided_within_the_budget():
    a, b = box(0, 0, 20, 15), box(20.3, 0, 30, 15)
    cleaned = city_of(building("a", a, 6), building("b", b, 30)).clean_footprints()
    block = cleaned.buildings[0]
    parts = parts_by_id(block)
    assert_parts_tile_block(block)
    # Each part stays within delta / 2 of its input footprint.
    for part, raw in ((parts["a"], a), (parts["b"], b)):
        cell = footprint(part)
        assert cell.difference(raw.buffer(0.25)).area < 1e-9
        assert raw.buffer(-0.25).difference(cell).area < 1e-9
    report = check_labelled_contract(
        {0: a, 1: b}, [1, 0], footprint(block),
        [(footprint(parts["a"]), 0), (footprint(parts["b"]), 1)],
        delta=0.5, epsilon=0.25,
    )
    assert report["status"] == "strict", report
    assert cleaned.attributes[RECORD_KEY]["walls"]["strict"] == 1


def test_building_inside_a_taller_one_has_no_footprint_but_stays_a_part():
    tower = building("tower", box(0, 0, 20, 20), 40)
    shed = building("shed", box(5, 5, 10, 10), 4, kind="shed")
    cleaned = city_of(tower, shed).clean_footprints()
    block = cleaned.buildings[0]
    parts = parts_by_id(block)
    assert set(parts) == {"tower", "shed"}
    assert footprint(parts["tower"]).equals(box(0, 0, 20, 20))
    assert parts["shed"].get_geometry(GeometryType.LOD0) is None
    assert parts["shed"].attributes["kind"] == "shed"

    again = cleaned.clean_footprints()
    assert [p.id for p in again.buildings[0].building_parts] == ["tower", "shed"]
    assert again.attributes[RECORD_KEY]["outcome"] == "unchanged"


def test_a_part_split_by_a_taller_building_has_a_multisurface_footprint(tmp_path):
    from dtcc_core import io
    from dtcc_core.builder.cleaning.clean import cleaning_record

    cleaned = city_of(
        building("long", box(0, 0, 40, 10), 6), building("tower", box(15, -5, 25, 15), 30)
    ).clean_footprints()
    block = cleaned.buildings[0]
    parts = parts_by_id(block)
    long_footprint = parts["long"].get_geometry(GeometryType.LOD0)
    assert isinstance(long_footprint, MultiSurface)
    assert footprint(parts["long"]).equals(box(0, 0, 15, 10).union(box(25, 0, 40, 10)))
    assert_parts_tile_block(block)

    again = cleaned.clean_footprints()
    assert again.attributes[RECORD_KEY]["outcome"] == "unchanged"
    assert footprint(parts_by_id(again.buildings[0])["long"]).equals(footprint(parts["long"]))

    io.save_model(cleaned, tmp_path / "block.dtcc")
    loaded = io.load_model(tmp_path / "block.dtcc")
    assert isinstance(
        parts_by_id(loaded.buildings[0])["long"].get_geometry(GeometryType.LOD0), MultiSurface
    )
    assert cleaning_record(loaded, min_building_detail=0.5) is not None


def test_cleaning_a_cleaned_block_with_walls_changes_nothing():
    city = city_of(
        building("a", box(0, 0, 20, 15), 6),
        building("b", box(20.3, 0, 30, 15), 30),
        building("c", box(29.6, 0, 40, 15), 12),
    )
    once = city.clean_footprints()
    twice = once.clean_footprints()
    assert twice.attributes[RECORD_KEY]["outcome"] == "unchanged"
    for x, y in zip(once.buildings, twice.buildings):
        assert footprint(x).equals_exact(footprint(y), 0.0)
        assert [p.id for p in x.building_parts] == [p.id for p in y.building_parts]
        for p, q in zip(x.building_parts, y.building_parts):
            assert footprint(p).equals_exact(footprint(q), 0.0)
            assert p.get_geometry(GeometryType.LOD0).zmax == q.get_geometry(GeometryType.LOD0).zmax
        assert x.estimated_height == y.estimated_height


def test_a_block_cleaned_on_its_own_keeps_its_parts_footprints():
    block = city_of(
        building("a", box(0, 0, 10, 10), 10), building("b", box(10, 0, 20, 10), 20)
    ).clean_footprints().buildings[0]
    again = block.clean_footprint()
    assert again.id == block.id
    assert [footprint(p) for p in again.building_parts] == [
        footprint(p) for p in block.building_parts
    ]


def fallback(polygon, sources, order, **kwargs):
    return LabelledRegion(polygon, [], "fallback", {"reason": "test"})


def test_failed_walls_fall_back_to_one_footprint(monkeypatch, caplog):
    monkeypatch.setattr(clean_module, "label_region", fallback)
    cleaned = city_of(
        building("a", box(0, 0, 10, 10), 10), building("b", box(10, 0, 20, 10), 11)
    ).clean_footprints()
    block = cleaned.buildings[0]
    assert all(p.get_geometry(GeometryType.LOD0) is None for p in block.building_parts)
    assert footprint(block).equals(box(0, 0, 20, 10))
    assert block.estimated_height == pytest.approx(10.5)
    walls = cleaned.attributes[RECORD_KEY]["walls"]
    assert walls["fallback"] == 1 and walls["fallback_ids"] == [["a", "b"]]
    assert "keep one footprint" in caplog.text

    # A block without part footprints is cleaned as one footprint again.
    monkeypatch.undo()
    again = cleaned.clean_footprints()
    assert again.attributes[RECORD_KEY]["outcome"] == "unchanged"
    assert again.attributes[RECORD_KEY]["walls"]["fallback"] == 0
    assert again.buildings[0].attributes[MERGED_FROM_KEY] == ["a", "b"]


def test_strict_cleaning_falls_back_instead_of_wall_warnings(monkeypatch):
    def warn(polygon, sources, order, **kwargs):
        labelled = label_region(polygon, sources, order, **kwargs)
        return LabelledRegion(labelled.region, labelled.cells, "warning", labelled.report)

    monkeypatch.setattr(clean_module, "label_region", warn)
    city = city_of(building("a", box(0, 0, 10, 10), 10), building("b", box(10, 0, 20, 10), 20))
    assert city.clean_footprints().attributes[RECORD_KEY]["walls"]["warning"] == 1
    walls = city.clean_footprints(strict=True).attributes[RECORD_KEY]["walls"]
    assert walls["warning"] == 0 and walls["fallback"] == 1


def survey_group(key):
    payload = json.loads((DATA / "survey-merge-groups.geojson").read_text())
    return [
        shape(f["geometry"])
        for f in payload["features"]
        if f"{f['properties']['case_id']}#{f['properties']['group']}" == key
    ]


@pytest.mark.parametrize(
    "key",
    [
        "city_grid:helsingborg:015#2",
        "city_grid:norrkoping:039#23",
        "city_grid:lund:016#8",
        "helsingborg_courtyard",
    ],
)
def test_walls_of_real_blocks_satisfy_the_identity_contract(key):
    if key == "helsingborg_courtyard":
        payload = json.loads((DATA / "helsingborg-courtyard.geojson").read_text())
        raw = [shape(f["geometry"]) for f in payload["features"]]
    else:
        raw = survey_group(key)
    rng = np.random.default_rng(0)
    heights = [float(rng.choice([6, 9, 12, 21, 30])) for _ in raw]
    buildings = [building(str(i), p, h) for i, (p, h) in enumerate(zip(raw, heights))]
    regions = clean_regions(buildings, enable_logging=False)
    reports = {tuple(row["labels"]): row for row in regions.diagnostics["walls"]}
    assert reports
    for output, indices in enumerate(regions.output_sources):
        if len(indices) < 2:
            continue
        region = regions.polygons[regions.region_outputs.index(output)]
        report = reports[tuple(indices)]
        assert report["outcome"] in ("strict", "warning")
        order = sorted(indices, key=lambda i: (-heights[i], -raw[i].area, i))
        contract = check_labelled_contract(
            {i: raw[i] for i in indices},
            order,
            region,
            regions.output_cells[output],
            delta=0.5,
            epsilon=report["fidelity_budget"],
        )
        assert contract["status"] == report["outcome"], contract


def test_heights_per_part_from_roof_points():
    city = city_of(
        building("low", box(0, 0, 20, 10), 0), building("high", box(20, 0, 30, 10), 0)
    )
    block = city.clean_footprints().buildings[0]
    xs, ys = np.meshgrid(np.arange(0.5, 30, 1.0), np.arange(0.5, 10, 1.0))
    zs = np.where(xs < 20, 6.0, 30.0)
    points = np.column_stack([xs.ravel(), ys.ravel(), zs.ravel()])
    cloud = PointCloud(points=points, classification=np.full(len(points), 6))
    terrain = Raster()
    terrain.data = np.zeros((4, 4))
    terrain.set_bounds(Bounds(-10, -10, 40, 20))

    extract_roof_points([block], cloud, statistical_outlier_remover=False)
    compute_building_heights([block], terrain, overwrite=True)

    parts = parts_by_id(block)
    assert parts["low"].estimated_height == pytest.approx(6.0)
    assert parts["high"].estimated_height == pytest.approx(30.0)
    assert parts["high"].get_geometry(GeometryType.LOD0).zmax == pytest.approx(30.0)
    # The block has the roof of its highest part.
    assert block.estimated_height == pytest.approx(30.0)
    assert block.get_geometry(GeometryType.LOD0).zmax == pytest.approx(30.0)


def test_heights_per_part_from_attribute():
    city = city_of(
        building("a", box(0, 0, 20, 10), 0, storey_height=9.0),
        building("b", box(20, 0, 30, 10), 0, storey_height=10.0),
    )
    block = city.clean_footprints().buildings[0]
    terrain = Raster()
    terrain.data = np.zeros((4, 4))
    terrain.set_bounds(Bounds(-10, -10, 40, 20))
    set_building_heights_from_attribute([block], terrain, height_attribute="storey_height")
    parts = parts_by_id(block)
    assert [parts["a"].estimated_height, parts["b"].estimated_height] == [9.0, 10.0]
    assert block.estimated_height == 10.0
    assert block.get_geometry(GeometryType.LOD0).zmax == 10.0
