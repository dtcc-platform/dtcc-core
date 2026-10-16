"""Cleaned blocks through heights, LOD1, saving and cleaning again.

Regression tests for the review of issue #147: the order of cleaning and
height estimation, the lifecycle of parts, unique IDs with provenance, and
traversal of part geometry from a city.
"""

from dtcc_core.builder.cleaning import walls as walls_module

import json
from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import box, shape

from dtcc_core import io
from dtcc_core.builder import build_lod1_buildings
from dtcc_core.builder.cleaning.blocks import (
    MERGED_FROM_KEY,
    SOURCE_FOOTPRINT_KEY,
    SOURCE_ID_KEY,
    WALLS_KEY,
    block_walls,
    set_part_roof,
    set_source_footprint,
    source_footprint,
    source_ids,
    update_walls,
)
from dtcc_core.builder.cleaning.clean import RECORD_KEY, cleaning_record, footprint_cells
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.builder.geometry_builders.buildings import (
    compute_building_heights,
    extrude_building,
    extract_roof_points,
    set_building_heights_from_attribute,
)
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

DATA = Path(__file__).resolve().parents[1] / "data" / "cleaning"


_original_xy = walls_module._WallRepair._xy

def building(id, polygon, roof=0.0, *, height=None, ground=None):
    surface = Surface()
    surface.from_polygon(polygon, roof)
    value = Building(id=id)
    value.add_geometry(surface, GeometryType.LOD0)
    if height is not None:
        value.estimated_height = float(height)
    if ground is not None:
        value.attributes["ground_height"] = float(ground)
    return value


def city_of(*buildings):
    city = City(id="city")
    raster = Raster()
    raster.data = np.zeros((4, 4))
    raster.set_bounds(Bounds(-20, -20, 60, 60))
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    city.add_buildings(list(buildings))
    return city


def polygon_of(geometry):
    if geometry is None:
        return None
    if isinstance(geometry, MultiSurface):
        from shapely.ops import unary_union

        return unary_union([s.to_polygon(simplify=0.0) for s in geometry.surfaces])
    return geometry.to_polygon(simplify=0.0)


def parts_by_id(block):
    return {part.id: part for part in block.building_parts}


def cells_by_id(block):
    """The roofs and cells of a block's parts with cells, by part ID."""
    from shapely.ops import unary_union

    cells, roofs = footprint_cells(block)
    parts = block.building_parts
    pieces = {}
    for polygon, index in cells:
        pieces.setdefault(parts[index].id, []).append(polygon)
    return (
        {parts[index].id: roofs[index] for _, index in cells},
        {id: unary_union(polygons) for id, polygons in pieces.items()},
    )


TOWER = box(5, 5, 10, 10)
LOW = box(0, 0, 20, 20)


def tower_heights(buildings):
    """Roof points of a 30 m tower inside a 6 m building, on flat ground."""
    xs, ys = np.meshgrid(np.arange(0.25, 20, 0.5), np.arange(0.25, 20, 0.5))
    zs = np.where((xs > 5) & (xs < 10) & (ys > 5) & (ys < 10), 30.0, 6.0)
    points = np.column_stack([xs.ravel(), ys.ravel(), zs.ravel()])
    cloud = PointCloud(points=points, classification=np.full(len(points), 6))
    terrain = Raster()
    terrain.data = np.zeros((4, 4))
    terrain.set_bounds(Bounds(-20, -20, 60, 60))
    extract_roof_points(buildings, cloud, statistical_outlier_remover=False)
    compute_building_heights(buildings, terrain, overwrite=True)
    return buildings


def raw_tower_city():
    return city_of(building("low", LOW), building("tower", TOWER))


def test_a_tower_inside_a_low_building_survives_cleaning_before_heights():
    # Height first: the tower wins its overlap.
    height_first = raw_tower_city()
    tower_heights(height_first.buildings)
    expected_roofs, expected_cells = cells_by_id(height_first.clean_footprints().buildings[0])
    assert expected_roofs == {"low": pytest.approx(6.0), "tower": pytest.approx(30.0)}
    assert expected_cells["tower"].equals(TOWER)

    # Geometry first: without heights the larger footprint wins, but the
    # tower keeps the footprint of its source, and its height.
    cleaned = raw_tower_city().clean_footprints()
    block = cleaned.buildings[0]
    assert block.attributes[WALLS_KEY]["precedence"] == "area"
    assert parts_by_id(block)["tower"].get_geometry(GeometryType.LOD0) is None
    assert source_footprint(parts_by_id(block)["tower"]).equals(TOWER)
    tower_heights(cleaned.buildings)
    assert parts_by_id(block)["tower"].estimated_height == pytest.approx(30.0)
    assert block.estimated_height == pytest.approx(30.0)

    # The walls follow the heights: the height builder recomputes them, so
    # meshing, LOD1 and cleaning again see the same cells as when the
    # heights came first.
    assert block.attributes[WALLS_KEY]["precedence"] == "roof"
    assert block_walls(block).current
    roofs, cells = cells_by_id(block)
    assert roofs == expected_roofs
    assert all(cells[id].equals(expected_cells[id]) for id in cells)
    assert cleaning_record(cleaned, min_building_detail=0.5) is not None
    mesh = meshes.build_city_surface_mesh(cleaned, report_mesh_quality=False)
    assert mesh.vertices[:, 2].max() == pytest.approx(30.0)
    again = cleaned.clean_footprints()
    assert again.buildings[0].attributes[WALLS_KEY]["precedence"] == "roof"
    assert polygon_of(parts_by_id(again.buildings[0])["tower"].lod0).equals(TOWER)
    build_lod1_buildings([block])
    assert block_walls(block).current
    tower = parts_by_id(block)["tower"]
    assert polygon_of(tower.lod0).equals(TOWER)
    assert tower.lod1.bounds.zmax == pytest.approx(30.0)


def test_roof_levels_decide_overlaps_without_heights():
    block = city_of(building("low", LOW, 6.0), building("tower", TOWER, 30.0)).clean_footprints()
    block = block.buildings[0]
    assert block.attributes[WALLS_KEY]["precedence"] == "roof"
    assert block.attributes[WALLS_KEY]["order"] == ["tower", "low"]
    assert polygon_of(parts_by_id(block)["tower"].lod0).equals(TOWER)


def test_on_sloping_ground_the_higher_roof_wins_over_the_taller_building():
    uphill = building("uphill", box(0, 0, 10, 10), 18.0, height=8, ground=10)
    downhill = building("downhill", box(8, 0, 20, 10), 12.0, height=12, ground=0)
    block = city_of(uphill, downhill).clean_footprints().buildings[0]
    parts = parts_by_id(block)
    assert block.attributes[WALLS_KEY]["precedence"] == "roof"
    assert polygon_of(parts["uphill"].lod0).equals(box(0, 0, 10, 10))
    assert polygon_of(parts["downhill"].lod0).equals(box(10, 0, 20, 10))

    # Heights without roof levels: the taller building wins.
    uphill = building("uphill", box(0, 0, 10, 10), height=8, ground=10)
    downhill = building("downhill", box(8, 0, 20, 10), height=12, ground=0)
    block = city_of(uphill, downhill).clean_footprints().buildings[0]
    assert block.attributes[WALLS_KEY]["precedence"] == "height"
    assert polygon_of(parts_by_id(block)["downhill"].lod0).equals(box(8, 0, 20, 10))


def test_parts_survive_lod1_saving_and_cleaning_again(tmp_path):
    city = city_of(
        building("low", box(0, 0, 10, 10), 10.0, height=10, ground=0),
        building("high", box(10, 0, 20, 10), 20.0, height=20, ground=0),
    )
    cleaned = city.clean_footprints()
    build_lod1_buildings(cleaned.buildings)

    def roofs(city):
        (block,) = city.buildings
        return sorted(part.lod0.zmax for part in block.building_parts if part.lod0 is not None)

    # After LOD1 the default footprint choice would be LOD1: a block is still
    # cleaned from its LOD0 footprint and its parts.
    assert roofs(cleaned.clean_footprints()) == [10.0, 20.0]
    assert cleaned.clean_footprints().attributes[RECORD_KEY]["outcome"] == "unchanged"
    # Zero detail leaves the footprints, and the walls, as they are.
    zero = cleaned.clean_footprints(min_building_detail=0)
    assert roofs(zero) == [10.0, 20.0]
    assert all(
        polygon_of(p.lod0).equals_exact(polygon_of(q.lod0), 0.0)
        for p, q in zip(zero.buildings[0].building_parts, cleaned.buildings[0].building_parts)
    )
    # Another detail cleans again, from the source footprints.
    coarse = cleaned.clean_footprints(min_building_detail=1.0)
    assert roofs(coarse) == [10.0, 20.0]
    assert coarse.buildings[0].attributes[WALLS_KEY]["min_building_detail"] == 1.0
    # Saved and loaded, the walls are current and the record valid.
    io.save_model(cleaned, tmp_path / "cleaned.dtcc")
    loaded = io.load_model(tmp_path / "cleaned.dtcc")
    assert cleaning_record(loaded, min_building_detail=0.5) is not None
    assert block_walls(loaded.buildings[0]).current
    assert roofs(loaded.clean_footprints()) == [10.0, 20.0]


def test_city_traversal_finds_part_geometry_once():
    city = city_of(
        building("a", box(0, 0, 10, 10), 10.0, height=10, ground=0),
        building("b", box(10, 0, 20, 10), 20.0, height=20, ground=0),
        building("solo", box(30, 0, 40, 10), 9.0, height=9, ground=0),
    ).clean_footprints()
    block, solo = city.buildings
    build_lod1_buildings(city.buildings)
    assert block.lod1 is None
    prisms = sum(len(part.lod1.surfaces) for part in block.building_parts)
    assert len(block.flatten_geometry(GeometryType.LOD1).surfaces) == prisms
    lod1 = city.flatten_geometry(GeometryType.LOD1)
    assert len(lod1.surfaces) == prisms + len(solo.lod1.surfaces)
    # The footprint of a block covers its parts' footprints: no double count.
    assert isinstance(block.flatten_geometry(GeometryType.LOD0), Surface)
    footprints = city.flatten_geometry(GeometryType.LOD0)
    assert len(footprints.surfaces) == 2


def split_source():
    polygon = shape(json.loads((DATA / "split-source.geojson").read_text()))
    return building("source", polygon, 10.0, height=10, ground=0), polygon


def all_ids(city):
    return [obj.id for b in city.buildings for obj in [b, *b.building_parts]]


def test_a_source_split_over_two_blocks_has_unique_ids_and_provenance(tmp_path):
    source, polygon = split_source()
    cleaned = city_of(source).clean_footprints()
    assert len(cleaned.buildings) == 2
    ids = all_ids(cleaned)
    assert len(ids) == len(set(ids))
    parts = [p for b in cleaned.buildings for p in b.building_parts]
    assert [p.id for p in parts] == ["source#1", "source#2"]
    assert all(p.attributes[SOURCE_ID_KEY] == "source" for p in parts)
    assert all(source_footprint(p).equals(polygon) for p in parts)
    assert [source_ids(b) for b in cleaned.buildings] == [["source"], ["source"]]
    io.save_model(cleaned, tmp_path / "split.dtcc")
    loaded = io.load_model(tmp_path / "split.dtcc")
    assert all_ids(loaded) == ids
    # Cleaning is deterministic, IDs included.
    assert all_ids(city_of(split_source()[0]).clean_footprints()) == ids
    assert all_ids(cleaned.clean_footprints()) == ids


def test_a_split_source_merged_with_a_neighbour():
    source, polygon = split_source()
    neighbour = building("neighbour", box(536780, 6474050, 536790, 6474056), 12.0, height=12, ground=0)
    cleaned = city_of(source, neighbour).clean_footprints()
    ids = all_ids(cleaned)
    assert len(ids) == len(set(ids))
    merged = [b for b in cleaned.buildings if len(b.building_parts) == 2]
    assert len(merged) == 1
    assert [p.id for p in merged[0].building_parts] == ["source#1", "neighbour"]
    assert merged[0].attributes[MERGED_FROM_KEY] == ["source", "neighbour"]
    assert sorted(source_ids(b)[0] for b in cleaned.buildings) == ["source", "source"]


def test_new_block_ids_are_derived_from_their_parts():
    def block_id(*ids):
        city = city_of(*(building(id, box(10 * i, 0, 10 * i + 10, 10), 10.0) for i, id in enumerate(ids)))
        return city.clean_footprints().buildings[0].id

    assert block_id("a", "b") == block_id("a", "b")
    assert block_id("a", "b") != block_id("a", "c")


def test_a_block_read_back_without_part_attributes_is_divided_by_its_parts(tmp_path):
    # The CityJSON reader keeps neither part IDs nor part attributes: parts
    # then have their own footprints as their sources.
    cleaned = city_of(
        building("a", box(0, 0, 10, 10), 10.0, height=10, ground=0),
        building("b", box(10, 0, 20, 10), 20.0, height=20, ground=0),
    ).clean_footprints()
    io.save_city(cleaned, tmp_path / "cleaned.json")
    loaded = io.load_city(tmp_path / "cleaned.json")
    assert all(SOURCE_ID_KEY not in p.attributes for p in loaded.buildings[0].building_parts)
    loaded.add_terrain(city_of().terrain)
    mesh = meshes.build_city_surface_mesh(loaded, report_mesh_quality=False)
    assert mesh.vertices[:, 2].max() == pytest.approx(20.0)
    again = loaded.clean_footprints().buildings[0]
    assert sorted(p.lod0.zmax for p in again.building_parts) == [10.0, 20.0]
    assert all(source_footprint(p) is not None for p in again.building_parts)


def test_new_roofs_keep_walls_unless_they_change_the_winner_of_an_overlap():
    cleaned = city_of(
        building("a", box(0, 0, 10, 10), 10.0),
        building("b", box(10, 0, 20, 10), 20.0),
        building("c", box(19.5, 0, 30, 10), 15.0),
    ).clean_footprints()
    block = cleaned.buildings[0]
    parts = parts_by_id(block)
    assert block.attributes[WALLS_KEY]["order"] == ["b", "c", "a"]
    # a overlaps nothing: its new roof changes the order but not the walls.
    set_part_roof(parts["a"], 25.0)
    assert block_walls(block).current
    # c above b changes the winner of their overlap.
    set_part_roof(parts["c"], 22.0)
    walls = block_walls(block)
    assert not walls.current
    pieces = {block.building_parts[i].id: cell for cell, i in walls.cells}
    assert pieces["c"].equals(box(19.5, 0, 30, 10))


def test_small_overlap_order_change_revalidates_and_persists_cells(tmp_path, monkeypatch):
    from dtcc_core.builder.cleaning import blocks
    from dtcc_core.builder.cleaning.contract import check_labelled_contract

    city = city_of(
        building("a", box(0, 0, 10, 10), 10, height=10, ground=0),
        building("b", box(9.98, 0, 20, 10), 20, height=20, ground=0),
    ).clean_footprints()
    block = city.buildings[0]
    before = cells_by_id(block)[1]
    set_part_roof(parts_by_id(block)["a"], 30)

    def unexpected_reconstruction(*args, **kwargs):
        pytest.fail("Faithful cells should be revalidated without reconstruction")

    monkeypatch.setattr(blocks, "wall_cells", unexpected_reconstruction)
    assert not update_walls(block)
    assert block.attributes[WALLS_KEY]["order"] == ["a", "b"]
    roofs, after = cells_by_id(block)
    assert roofs == {"a": 30, "b": 20}
    assert all(after[id].equals_exact(before[id], 0) for id in before)
    record = block.attributes[WALLS_KEY]
    walls = block_walls(block)
    checked = check_labelled_contract(
        blocks.part_sources(block), walls.order, blocks.block_region(block), walls.cells,
        delta=record["min_building_detail"], epsilon=record["fidelity_budget"],
    )
    assert checked["status"] != "fail"
    assert all(value == "pass" for value in checked["fidelity"].values())
    io.save_model(city, tmp_path / "revalidated.dtcc")
    loaded = io.load_model(tmp_path / "revalidated.dtcc").buildings[0]
    # Saving the new certificate avoids repeating validation on every consumer.
    monkeypatch.setattr(blocks, "condition_walls", unexpected_reconstruction)
    assert block_walls(loaded).current
    assert not update_walls(loaded)


def test_meshing_reports_the_outcomes_of_recomputed_walls(monkeypatch):
    from dtcc_core.builder.cleaning import blocks as blocks_module
    from dtcc_core.builder.cleaning.walls import LabelledRegion, label_region

    cleaned = raw_tower_city().clean_footprints()
    parts = parts_by_id(cleaned.buildings[0])
    set_part_roof(parts["low"], 6.0)
    set_part_roof(parts["tower"], 30.0)

    def warn(self):
        labelled = _original_xy(self)
        return LabelledRegion(labelled.region, labelled.cells, "warning", labelled.report)

    monkeypatch.setattr(walls_module._WallRepair, "_xy", warn)
    footprints = meshes._meshing_footprints(
        cleaned, lod=None, min_building_detail=0.5, min_building_area=0.0, max_mesh_size=5.0
    )
    walls = footprints.diagnostics["walls"]
    assert cleaned.attributes[RECORD_KEY]["walls"]["warning"] == 0
    assert (walls["warning"], walls["recomputed"]) == (1, 1)
    assert footprints.diagnostics["before_selection_contract"]["wall_warnings"] == 1


def total_area(cells):
    return sum(polygon.area for polygon, _ in cells)


def test_walls_follow_a_changed_region_and_changed_sources():
    def adjacent():
        return city_of(
            building("a", box(0, 0, 10, 10), 10.0), building("b", box(10, 0, 20, 10), 20.0)
        ).clean_footprints()

    # A larger region: the old cells no longer cover it.
    for region in (box(0, 0, 30, 10), box(0, 0, 15, 10)):
        cleaned = adjacent()
        block = cleaned.buildings[0]
        surface = Surface()
        surface.from_polygon(region, 20.0)
        block.add_geometry(surface, GeometryType.LOD0)
        walls = block_walls(block)
        assert not walls.current
        assert walls.outcome == "fallback" or total_area(walls.cells) == pytest.approx(region.area)
        again = cleaned.clean_footprints().buildings[0]
        cells, _ = footprint_cells(again)
        assert total_area(cells) == pytest.approx(region.area)
        footprints = meshes._meshing_footprints(
            cleaned, lod=None, min_building_detail=0.5, min_building_area=0.0, max_mesh_size=5.0
        )
        assert sum(s.to_polygon(simplify=0.0).area for s in footprints.surfaces) == pytest.approx(
            region.area
        )

    # A changed source footprint: the cells are made again from it.
    cleaned = city_of(
        building("a", box(0, 0, 10, 10), 10.0), building("b", box(8, 0, 20, 10), 20.0)
    ).clean_footprints()
    block = cleaned.buildings[0]
    parts = parts_by_id(block)
    assert polygon_of(parts["b"].lod0).equals(box(8, 0, 20, 10))
    set_source_footprint(parts["b"], box(9, 0, 20, 10), 20.0)
    walls = block_walls(block)
    assert not walls.current
    cells = {block.building_parts[i].id: cell for cell, i in walls.cells}
    assert cells["b"].equals(box(9, 0, 20, 10))
    assert update_walls(block)
    assert polygon_of(parts["b"].lod0).equals(box(9, 0, 20, 10))
    assert block_walls(block).current


def tower_block():
    """A 50 m tower inside a 6 m building, cleaned with its heights."""
    city = city_of(
        building("low", LOW, 6.0, height=6, ground=0),
        building("tower", TOWER, 50.0, height=50, ground=0),
    ).clean_footprints()
    return city, city.buildings[0]


@pytest.mark.parametrize("rebuild", [True, False])
def test_lod1_follows_a_new_partition(rebuild):
    city, block = tower_block()
    build_lod1_buildings(city.buildings)
    parts = parts_by_id(block)
    assert parts["tower"].lod1.bounds.zmax == pytest.approx(50.0)
    # New heights: the low building rises above the tower, which loses its cell.
    parts["low"].attributes["review_height"] = 40.0
    parts["tower"].attributes["review_height"] = 20.0
    set_building_heights_from_attribute(
        city.buildings,
        city.terrain.get_geometry(GeometryType.RASTER),
        height_attribute="review_height",
        always_use_default_ground=True,
    )
    assert parts["tower"].lod0 is None
    assert parts["tower"].lod1 is None
    build_lod1_buildings(city.buildings, rebuild=rebuild)
    assert parts["tower"].lod1 is None
    assert parts["low"].lod1.bounds.zmax == pytest.approx(40.0)
    lod1 = city.flatten_geometry(GeometryType.LOD1)
    assert len(lod1.surfaces) == len(parts["low"].lod1.surfaces)
    assert lod1.bounds.zmax == pytest.approx(40.0)


def test_lod1_follows_a_switch_to_and_from_one_footprint(monkeypatch):
    from dtcc_core.builder.cleaning import blocks as blocks_module

    city, block = tower_block()
    build_lod1_buildings(city.buildings)
    parts = parts_by_id(block)
    real = walls_module._WallRepair._xy

    def failure(*args, **kwargs):
        raise shapely.errors.GEOSException("controlled construction failure")

    # The walls cannot be made for the new roofs: the block keeps one
    # footprint and one prism, and the parts lose theirs.
    monkeypatch.setattr(walls_module._WallRepair, "_xy", failure)
    set_part_roof(parts["tower"], 4.0)
    build_lod1_buildings(city.buildings, rebuild=False)
    assert block.attributes[WALLS_KEY]["outcome"] == "fallback"
    assert all(part.lod1 is None for part in block.building_parts)
    assert block.lod1 is not None
    assert len(city.flatten_geometry(GeometryType.LOD1).surfaces) == len(block.lod1.surfaces)

    # Walls again: the parts have prisms and the block none.
    monkeypatch.setattr(walls_module._WallRepair, "_xy", real)
    set_part_roof(parts["tower"], 50.0)
    build_lod1_buildings(city.buildings, rebuild=False)
    assert block.lod1 is None
    assert parts["tower"].lod1.bounds.zmax == pytest.approx(50.0)
    assert len(city.flatten_geometry(GeometryType.LOD1).surfaces) == sum(
        len(part.lod1.surfaces) for part in block.building_parts
    )


def model_ids(obj):
    ids = [obj.id]
    for children in obj.children.values():
        for child in children:
            ids.extend(model_ids(child))
    return ids


def test_ids_are_unique_in_the_whole_model(tmp_path):
    # Descendants of a split source get an ID per occurrence and keep the ID
    # of what they copy.
    source, _ = split_source()
    child = BuildingPart(id="original-child")
    child.add_geometry(source.lod0.copy(), GeometryType.LOD0)
    source.add_child(child)
    cleaned = city_of(source).clean_footprints()
    ids = model_ids(cleaned)
    assert len(ids) == len(set(ids))
    copies = [c for b in cleaned.buildings for p in b.building_parts for c in p.building_parts]
    assert [c.id for c in copies] == ["original-child#1", "original-child#2"]
    assert all(c.attributes[SOURCE_ID_KEY] == "original-child" for c in copies)
    io.save_model(cleaned, tmp_path / "nested.dtcc")
    assert model_ids(io.load_model(tmp_path / "nested.dtcc")) == ids
    assert model_ids(cleaned.clean_footprints()) == ids

    # Split and merged with a neighbour.
    neighbour = building("neighbour", box(536780, 6474050, 536790, 6474056), 12.0, height=12, ground=0)
    merged = city_of(source, neighbour).clean_footprints()
    ids = model_ids(merged)
    assert len(ids) == len(set(ids))
    assert sorted(c.id for b in merged.buildings for p in b.building_parts for c in p.building_parts) == [
        "original-child#1",
        "original-child#2",
    ]
    io.save_model(merged, tmp_path / "merged.dtcc")
    assert model_ids(io.load_model(tmp_path / "merged.dtcc")) == ids

    # New IDs avoid the IDs of every object of the model, the city's too.
    source, _ = split_source()
    city = city_of(source)
    city.id = "source#1"
    cleaned = city.clean_footprints()
    ids = model_ids(cleaned)
    assert len(ids) == len(set(ids))
    assert [p.id for b in cleaned.buildings for p in b.building_parts] == ["source#2", "source#3"]
    io.save_model(cleaned, tmp_path / "city-id.dtcc")


def test_a_part_without_source_or_cell_keeps_the_labels_of_the_others():
    cleaned = city_of(
        building("low", box(0, 0, 20, 20), 30.0),
        building("hidden", box(5, 5, 10, 10), 10.0),
        building("side", box(20, 0, 30, 20), 20.0),
    ).clean_footprints()
    block = cleaned.buildings[0]
    assert [p.lod0 is not None for p in block.building_parts] == [True, False, True]
    expected = {p.id: polygon_of(p.lod0) for p in block.building_parts if p.lod0 is not None}
    for part in block.building_parts:
        part.attributes.pop(SOURCE_FOOTPRINT_KEY)
    again = cleaned.clean_footprints().buildings[0]
    got = {p.id: polygon_of(p.lod0) for p in again.building_parts if p.lod0 is not None}
    assert set(got) == set(expected)
    assert all(got[id].equals(expected[id]) for id in got)


def test_the_cause_of_a_fallback_is_kept_with_the_block(monkeypatch, tmp_path, caplog):
    from dtcc_core.builder.cleaning import blocks as blocks_module

    def failure(*args, **kwargs):
        raise shapely.errors.GEOSException("controlled construction failure")

    monkeypatch.setattr(walls_module._WallRepair, "_xy", failure)
    cleaned = city_of(
        building("a", box(0, 0, 10, 10), 10.0), building("b", box(10, 0, 20, 10), 20.0)
    ).clean_footprints()
    expected = {"reason": "construction_error", "error": "controlled construction failure"}
    assert cleaned.attributes[RECORD_KEY]["walls"]["fallback_reasons"] == [expected]
    assert cleaned.buildings[0].attributes[WALLS_KEY]["detail"] == expected
    monkeypatch.undo()
    # Saved, loaded and cleaned again, the walls are reused with their cause.
    io.save_model(cleaned, tmp_path / "fallback.dtcc")
    loaded = io.load_model(tmp_path / "fallback.dtcc")
    assert block_walls(loaded.buildings[0]).detail == expected
    again = loaded.clean_footprints()
    assert again.buildings[0].attributes[WALLS_KEY]["detail"] == expected
    assert again.attributes[RECORD_KEY]["walls"]["fallback_reasons"] == [expected]

    caplog.clear()
    stage = meshes._meshing_footprints(
        loaded, lod=GeometryType.LOD0, min_building_detail=0.5,
        min_building_area=15, max_mesh_size=10,
    )
    assert stage.contract["status"] == "warn"
    assert all(stage.contract["requirements"].values())
    assert any("construction_error" in message for message in stage.contract["warnings"])
    assert stage.diagnostics["kept_wall_fallbacks"][0]["ids"] == [loaded.buildings[0].id]
    assert loaded.buildings[0].id in caplog.text
    assert "wall fallback (construction_error)" in caplog.text
    caplog.clear()
    excluded = meshes._meshing_footprints(
        loaded, lod=GeometryType.LOD0, min_building_detail=0.5,
        min_building_area=1000, max_mesh_size=10,
    )
    assert excluded.diagnostics["kept_wall_fallbacks"] == []
    assert not any("wall fallback" in message for message in excluded.contract["warnings"])
    assert "wall fallback" not in caplog.text

    # A strict cleaning keeps the closest wall pair of the walls it refused.
    from dtcc_core.builder.cleaning.walls import LabelledRegion, label_region

    def warn(self):
        labelled = _original_xy(self)
        report = {**labelled.report, "wall_pairs": 1, "min_wall_distance": 0.01}
        return LabelledRegion(labelled.region, labelled.cells, "warning", report)

    monkeypatch.setattr(walls_module._WallRepair, "_xy", warn)
    cleaned = city_of(
        building("a", box(0, 0, 10, 10), 10.0), building("b", box(10, 0, 20, 10), 20.0)
    ).clean_footprints()
    record = cleaned.buildings[0].attributes[WALLS_KEY]
    assert record["outcome"] == "warning"
    assert record["detail"]["wall_pairs"] == 1
    assert record["detail"]["min_wall_distance"] == .01
    assert record["detail"]["contact_repair"]["state"] == "unneeded"
    strict = cleaned.clean_footprints(strict=True).buildings[0].attributes[WALLS_KEY]
    assert strict["outcome"] == "fallback"
    assert strict["detail"] == {"reason": "strict", **record["detail"]}


@pytest.mark.parametrize("owner", ["city", "terrain"])
def test_block_ids_avoid_surviving_model_objects(owner, tmp_path):
    def raw():
        return city_of(building("a", box(0, 0, 10, 10), 10),
                       building("b", box(10, 0, 20, 10), 20))

    expected = raw().clean_footprints().buildings[0].id
    city = raw()
    target = city if owner == "city" else city.terrain
    target.id = expected
    cleaned = city.clean_footprints()
    block_id = cleaned.buildings[0].id
    assert block_id == expected + "#1"
    io.save_model(cleaned, tmp_path / "unique.dtcc")
    loaded = io.load_model(tmp_path / "unique.dtcc")
    for _ in range(2):
        loaded = loaded.clean_footprints()
        assert loaded.buildings[0].id == block_id
        io.save_model(loaded, tmp_path / "unique.dtcc")


def test_height_only_changes_invalidate_generated_prisms_after_save_load(tmp_path):
    city = city_of(
        building("a", box(0, 0, 10, 10), 10, height=10, ground=0),
        building("b", box(10, 0, 20, 10), 20, height=20, ground=0),
    ).clean_footprints()
    build_lod1_buildings(city.buildings)
    io.save_model(city, tmp_path / "prisms.dtcc")
    city = io.load_model(tmp_path / "prisms.dtcc")
    for part in city.buildings[0].building_parts:
        part.attributes["new_height"] = part.estimated_height + 30
    set_building_heights_from_attribute(
        city.buildings, city.terrain.get_geometry(GeometryType.RASTER),
        height_attribute="new_height", always_use_default_ground=True,
    )
    block = city.buildings[0]
    assert block_walls(block).current
    assert all(part.lod1 is None for part in block.building_parts)
    build_lod1_buildings(city.buildings, rebuild=False)
    assert [part.lod1.bounds.zmax for part in block.building_parts] == [40, 50]


@pytest.mark.parametrize("edited_extrusion", [False, True])
def test_height_updates_preserve_caller_lod1(edited_extrusion):
    city = city_of(building("a", box(0, 0, 10, 10), 10, height=10, ground=0))
    obj = city.buildings[0]
    if edited_extrusion:
        build_lod1_buildings([obj])
    else:
        obj.add_geometry(extrude_building(obj), GeometryType.LOD1)
    # A caller edits the generated geometry before changing the modelling roof.
    obj.lod1.surfaces[0].vertices[0, 2] += 1
    caller_geometry = obj.lod1
    obj.attributes["new_height"] = 40
    set_building_heights_from_attribute(
        [obj], city.terrain.get_geometry(GeometryType.RASTER),
        height_attribute="new_height", always_use_default_ground=True,
    )
    assert obj.lod1 is caller_geometry
    build_lod1_buildings([obj], rebuild=False)
    assert obj.lod1 is caller_geometry


def test_ground_only_changes_invalidate_a_generated_prism():
    city = city_of(building("a", box(0, 0, 10, 10), 10, height=10, ground=0))
    obj = city.buildings[0]
    build_lod1_buildings([obj])
    obj.attributes["next_height"] = 5
    set_building_heights_from_attribute(
        [obj], city.terrain.get_geometry(GeometryType.RASTER),
        height_attribute="next_height", default_ground_height=5,
        always_use_default_ground=True,
    )
    assert obj.lod0.zmax == 10
    assert obj.lod1 is None
    build_lod1_buildings([obj], rebuild=False)
    assert obj.lod1.bounds.zmin == 5
    assert obj.lod1.bounds.zmax == 10


def test_compatible_cached_warning_cells_are_reconditioned_after_load(tmp_path):
    from dtcc_core.builder.cleaning import blocks
    from dtcc_core.model.mixins.mesh.quality import tri_element_quality

    payload = json.loads((DATA / "orebro-007-cached-warning.geojson").read_text())
    source_rows = [f for f in payload["features"] if f["properties"]["role"] == "source"]
    region = next(shape(f["geometry"]) for f in payload["features"] if f["properties"]["role"] == "region")
    city = city_of(*(building(str(f["properties"]["label"]), shape(f["geometry"]),
                             f["properties"]["roof"], ground=0) for f in source_rows)).clean_footprints()
    block = city.buildings[0]
    surface = Surface()
    surface.from_polygon(region, max(f["properties"]["roof"] for f in source_rows))
    block.add_geometry(surface, GeometryType.LOD0)
    cells = [(shape(f["geometry"]), f["properties"]["label"])
             for f in payload["features"] if f["properties"]["role"] == "cell"]
    blocks.assign_cells(block.building_parts, cells)
    block.attributes[WALLS_KEY] = blocks.walls_record(
        block, payload["metadata"]["order"], "roof", "warning", delta=0.5,
        epsilon=0.25, strict=False, detail={"wall_pairs": 3, "min_wall_distance": 0.00128793},
    )
    old_mesh = meshes.build_building_shell(block, max_mesh_size=10)
    assert np.count_nonzero(tri_element_quality(old_mesh.vertices, old_mesh.faces) < 0.02) > 0
    build_lod1_buildings([block])
    block.attributes[WALLS_KEY].pop("repair_version")
    io.save_model(city, tmp_path / "cached-warning.dtcc")
    loaded = io.load_model(tmp_path / "cached-warning.dtcc")
    block = loaded.buildings[0]
    assert update_walls(block)
    assert block.attributes[WALLS_KEY]["repair_version"] == blocks.WALLS_REPAIR_VERSION
    assert all(part.lod1 is None for part in block.building_parts)
    mesh = meshes.build_building_shell(block, max_mesh_size=10)
    assert np.count_nonzero(tri_element_quality(mesh.vertices, mesh.faces) < 0.02) == 0
    assert {part_roof for part_roof in cells_by_id(block)[0].values()} == {
        f["properties"]["roof"] for f in source_rows
    }


def test_failed_conditioning_upgrade_retains_valid_warning_cells_and_prisms(monkeypatch, tmp_path, caplog):
    from dtcc_core.builder.cleaning import blocks

    city = city_of(building("a", box(0, 0, 10, 10), 10, ground=0),
                   building("b", box(10, 0, 20, 10), 20, ground=0)).clean_footprints()
    block = city.buildings[0]
    build_lod1_buildings([block])
    prisms = [part.lod1 for part in block.building_parts]
    record = block.attributes[WALLS_KEY]
    record.update(outcome="warning", detail={"wall_pairs": 1, "min_wall_distance": 0.01})
    record.pop("repair_version")

    def failure(*args, **kwargs):
        raise shapely.errors.GEOSException("controlled conditioning failure")

    monkeypatch.setattr(walls_module._WallRepair, "_xy", failure)
    assert not update_walls(block)
    assert "retaining recorded warning cells" in caplog.text
    assert "it keeps one footprint" not in caplog.text
    assert all(part.lod1 is prism for part, prism in zip(block.building_parts, prisms))
    assert cells_by_id(block)[0] == {"a": 10, "b": 20}
    assert block.attributes[WALLS_KEY]["detail"]["repair_failure"] == {
        "reason": "construction_error", "error": "controlled conditioning failure",
    }
    io.save_model(city, tmp_path / "retained-warning.dtcc")
    loaded = io.load_model(tmp_path / "retained-warning.dtcc")
    assert loaded.buildings[0].attributes[WALLS_KEY]["detail"] == block.attributes[WALLS_KEY]["detail"]

    # Retention also reports contacts in an earlier warning partition accurately.
    pinched, *_ = _two_contact_block()
    blocks.walls_of(pinched).update(outcome="warning", detail={"wall_pairs": 1, "min_wall_distance": .01})
    result = block_walls(pinched)
    assert result.current and result.outcome == "warning"
    assert result.report["contact_repair"]["state"] == "unchanged"
    assert result.report["contact_repair"]["remaining_contacts"] == 2
    assert not update_walls(pinched)
    assert blocks.walls_of(pinched)["detail"]["contact_repair"]["remaining_contacts"] == 2
    assert not update_walls(pinched)


def test_changed_exposed_roof_groups_repair_contacts_even_with_a_hidden_source():
    city = city_of(
        building("high1", box(0, 0, 10, 10), 30, height=30, ground=0),
        building("high2", box(10, 10, 20, 20), 30, height=30, ground=0),
        building("low1", box(10, 0, 20, 10), 29.6, height=29.6, ground=0),
        building("low2", box(0, 10, 10, 20), 29.6, height=29.6, ground=0),
        building("hidden", box(2, 2, 4, 4), 29.2, height=29.2, ground=0),
    ).clean_footprints()
    block = city.buildings[0]
    assert parts_by_id(block)["hidden"].lod0 is None
    for part in block.building_parts:
        part.attributes["next_height"] = 29.4 if part.id.startswith("low") else part.estimated_height
    set_building_heights_from_attribute(
        [block], city.terrain.get_geometry(GeometryType.RASTER),
        height_attribute="next_height", always_use_default_ground=True,
    )
    assert block_walls(block).current
    stage = meshes._meshing_footprints(
        city, lod=None, min_building_detail=0.5, min_building_area=0,
        max_mesh_size=5, volume=True,
    )
    assert stage.diagnostics["levels"]["one_roof"] == []
    assert {surface.zmax for surface in stage.surfaces} == {29.4, 30}


def test_blocks_of_an_earlier_version_are_cleaned_again_from_their_buildings():
    cleaned = city_of(
        building("a", box(0, 0, 10, 10), 10.0), building("b", box(10, 0, 20, 10), 20.0)
    ).clean_footprints()
    block = cleaned.buildings[0]
    for record in ({k: v for k, v in block.attributes[WALLS_KEY].items() if k != "version"}, None):
        if record is None:
            del block.attributes[WALLS_KEY]
        else:
            block.attributes[WALLS_KEY] = record
        with pytest.raises(ValueError, match="earlier version of clean_footprints"):
            cleaned.clean_footprints()
        with pytest.raises(ValueError, match="Clean the buildings it was made from"):
            block_walls(block)


def test_compatible_warning_repairs_contacts_without_reconstructing_good_walls(tmp_path):
    from dtcc_core.builder.cleaning import blocks
    from dtcc_core.model.mixins.mesh.quality import tri_element_quality

    payload = json.loads((DATA / "uppsala-048-compatible-warning.geojson").read_text())
    sources = [f for f in payload["features"] if f["properties"]["role"] == "source"]
    region = next(shape(f["geometry"]) for f in payload["features"]
                  if f["properties"]["role"] == "region")
    block = Building(id="compatible", attributes={MERGED_FROM_KEY: list(map(str, range(len(sources))))})
    footprint = Surface()
    footprint.from_polygon(region, max(f["properties"]["roof"] for f in sources))
    block.add_geometry(footprint, GeometryType.LOD0)
    parts = []
    for feature in sources:
        part = BuildingPart(id=str(feature["properties"]["label"]))
        set_source_footprint(part, shape(feature["geometry"]), feature["properties"]["roof"])
        parts.append(part)
    block.add_children(parts)
    cells = [(shape(f["geometry"]), f["properties"]["label"])
             for f in payload["features"] if f["properties"]["role"] == "cell"]
    blocks.assign_cells(parts, cells)
    block.attributes[WALLS_KEY] = blocks.walls_record(
        block, payload["metadata"]["order"], "roof", "warning", delta=0.5,
        epsilon=0.25, strict=False,
        detail={"wall_pairs": 4, "min_wall_distance": 0.1562256861693886},
    )
    build_lod1_buildings([block])
    block.attributes[WALLS_KEY].pop("repair_version")
    city = city_of(block)
    io.save_model(city, tmp_path / "compatible.dtcc")
    loaded = io.load_model(tmp_path / "compatible.dtcc").buildings[0]
    update_walls(loaded)
    cells, _ = footprint_cells(loaded)
    from dtcc_core.builder.cleaning.contract import check_labelled_contract
    sources = blocks.part_sources(loaded)
    record = loaded.attributes[WALLS_KEY]
    result = check_labelled_contract(
        sources, payload["metadata"]["order"], region, cells,
        delta=0.5, epsilon=record["fidelity_budget"],
    )
    assert result["status"] in ("strict", "warning")
    mesh = meshes.build_building_shell(loaded, max_mesh_size=10)
    assert np.count_nonzero(tri_element_quality(mesh.vertices, mesh.faces) < 0.02) == 0


def _multiway_sources(spokes, transformed=False):
    from shapely.affinity import rotate, translate
    from shapely.geometry import Polygon

    vertices = np.column_stack([
        10 * np.cos(np.arange(spokes) * 2 * np.pi / spokes),
        10 * np.sin(np.arange(spokes) * 2 * np.pi / spokes),
    ])
    region = Polygon(vertices)
    sources = {
        i: Polygon([vertices[k] * .01, vertices[k], vertices[(k + 1) % spokes],
                    vertices[(k + 1) % spokes] * .01])
        for i, k in enumerate(range(0, spokes, 2))
    }
    roofs = {i: 20 + 5 * i for i in sources}
    sources[len(sources)] = region
    roofs[len(sources) - 1] = 6
    if transformed:
        sources = {i: translate(rotate(p, 23, origin=(0, 0)), 100000, 200000)
                   for i, p in sources.items()}
    return sources, roofs


@pytest.mark.parametrize("spokes", [6, 8])
@pytest.mark.parametrize("transformed", [False, True])
def test_cleaning_retains_nonpinched_multiway_roofs(spokes, transformed):
    from dtcc_core.builder.cleaning import blocks
    from dtcc_core.builder.cleaning.contract import check_labelled_contract, input_partition
    from dtcc_core.builder.cleaning.levels import level_regions, pinched_levels
    from dtcc_core.model.mixins.mesh.quality import tri_element_quality

    sources, roofs = _multiway_sources(spokes, transformed)
    city = city_of(*(building(str(i), p, roofs[i], height=roofs[i], ground=0)
                     for i, p in sources.items())).clean_footprints()
    block = city.buildings[0]
    walls = block_walls(block)
    current_sources = blocks.part_sources(block)
    current_roofs = {i: blocks.part_roof(p) for i, p in enumerate(block.building_parts)}
    record = blocks.walls_of(block)
    assert walls.outcome in {"strict", "warning"}
    assert {i for _, i in walls.cells} == set(current_sources)
    assert check_labelled_contract(current_sources, walls.order, blocks.block_region(block),
                                   walls.cells, delta=.5, epsilon=record["fidelity_budget"])["status"] != "fail"
    assert not pinched_levels([(p, z) for p, z, _ in level_regions(walls.cells, current_roofs, .5)])
    expected = input_partition(sources, sorted(sources, key=lambda i: -roofs[i]))
    original_solid = sum(p.area * roofs[i] for i, p in expected.items())
    assert abs(sum(p.area * current_roofs[i] for p, i in walls.cells) - original_solid) < 100
    shell = meshes.build_building_shell(block, max_mesh_size=5)
    assert np.all(tri_element_quality(shell.vertices, shell.faces) >= .02)
    stage = meshes._meshing_footprints(city, lod=GeometryType.LOD0, min_building_detail=.5,
                                      min_building_area=0, max_mesh_size=5, volume=True)
    assert not stage.diagnostics["levels"]["one_roof"]
    assert {s.zmax for s in stage.surfaces} == set(roofs.values())


def _two_contact_block(*, persisted_partial=False):
    from shapely.geometry import Polygon
    from dtcc_core.builder.cleaning import blocks, walls

    left = [Polygon([(0, 0), (20, 0), (20, 8), (10, 10), (0, 8)]),
            Polygon([(10, 10), (20, 8), (20, 12)]),
            Polygon([(0, 12), (10, 10), (20, 12), (20, 20), (0, 20)]),
            Polygon([(0, 8), (10, 10), (0, 12)])]
    region = box(0, 0, 40, 20)
    cells = list(zip(left + [box(20, 0, 30, 10), box(30, 0, 40, 10),
                             box(30, 10, 40, 20), box(20, 10, 30, 20)], range(8)))
    cells = walls._rebuild(region, walls._lines(shapely.union_all([p.boundary for p, _ in cells])), cells)
    sources = dict(enumerate(left + [box(20, 0, 30.17, 10.17), box(30.17, 0, 40, 10.17),
                                    box(30.17, 10.17, 40, 20), box(20, 10.17, 30.17, 20)]))
    roofs = {i: 20 if i % 2 == 0 else 6 for i in sources}
    order = [0, 2, 4, 6, 1, 3, 5, 7]
    block = building("two_contacts", region, 20)
    block.attributes["merged_from"] = list(map(str, sources))
    parts = []
    for i, p in sources.items():
        part = BuildingPart(id=str(i), attributes={"ground_height": 0})
        blocks.set_source_footprint(part, p, roofs[i])
        parts.append(part)
    block.add_children(parts)
    detail = {}
    if persisted_partial:
        result = walls.repair_roof_contacts(region, sources, order, cells, roofs,
                                            delta=.5, epsilon=.25)
        cells, report = result.cells, result.report
        assert result.remaining_contacts == 1
        detail = {"contact_repair": report}
    blocks.assign_cells(parts, cells)
    record = blocks.walls_record(block, order, "roof", "strict", delta=.5, epsilon=.25,
                                 strict=False, detail=detail)
    record["repair_version"] = 1 if persisted_partial else 0
    block.attributes[WALLS_KEY] = record
    return block, region, sources, order, cells, roofs


def test_partial_contact_repair_uses_permitted_retry():
    from dtcc_core.builder.cleaning import blocks
    from dtcc_core.builder.cleaning.contract import check_labelled_contract

    block, region, sources, order, cells, roofs = _two_contact_block()
    assert check_labelled_contract(sources, order, region, cells, delta=.5, epsilon=.25)["status"] != "fail"
    result = walls_module.condition_walls(region, sources, order, levels=roofs, delta=.5, epsilon=.25,
                                         previous=walls_module.LabelledRegion(region, cells, "strict"))
    assert result.contact_outcome == "complete" and result.remaining_contacts == 0
    assert result.accepted_budget == .375
    assert {i for _, i in result.cells} == set(sources)
    assert check_labelled_contract(sources, order, region, result.cells, delta=.5,
                                   epsilon=result.accepted_budget)["status"] != "fail"


@pytest.mark.parametrize("failure", ["no_improvement", "geometry_error"])
def test_failed_retry_retains_valid_partial_repair(monkeypatch, failure):
    from dtcc_core.builder.cleaning import blocks
    from dtcc_core.builder.cleaning.contract import check_labelled_contract

    block, region, sources, order, cells, roofs = _two_contact_block()
    original = walls_module.repair_roof_contacts

    def fail_later(*args, **kwargs):
        if kwargs["epsilon"] > .25:
            if failure == "geometry_error":
                raise shapely.errors.GEOSException("unusable later candidate")
            return walls_module.LabelledRegion(args[0], args[3], "warning", {}, "unchanged", 1, kwargs["epsilon"])
        return original(*args, **kwargs)

    monkeypatch.setattr(walls_module, "repair_roof_contacts", fail_later)
    result = walls_module.condition_walls(region, sources, order, levels=roofs, delta=.5, epsilon=.25,
                                         previous=walls_module.LabelledRegion(region, cells, "strict"))
    assert result.contact_outcome == "partial" and result.remaining_contacts == 1
    assert result.accepted_budget == .25
    assert result.report["contact_repair"]["changed_area"] > 0
    assert check_labelled_contract(sources, order, region, result.cells, delta=.5,
                                   epsilon=result.accepted_budget)["status"] != "fail"


@pytest.mark.parametrize("persisted_partial", [False, True])
def test_saved_partial_records_upgrade_once_and_preserve_roofs(tmp_path, monkeypatch, persisted_partial):
    from dtcc_core.builder.cleaning import blocks

    block, *_ = _two_contact_block(persisted_partial=persisted_partial)
    city = city_of(block)
    path = tmp_path / "contacts.dtcc"
    io.save_model(city, path)
    city = io.load_model(path)
    block = city.buildings[0]
    assert update_walls(block)
    record = blocks.walls_of(block)
    assert record["repair_version"] == blocks.WALLS_REPAIR_VERSION
    assert record["fidelity_budget"] == .375
    assert record["detail"]["contact_repair"]["remaining_contacts"] == 0
    stage = meshes._meshing_footprints(city, lod=GeometryType.LOD0, min_building_detail=.5,
                                      min_building_area=0, max_mesh_size=5, volume=True)
    assert not stage.diagnostics["levels"]["one_roof"]
    assert {s.zmax for s in stage.surfaces} == {6, 20}
    io.save_model(city, path)
    saved = io.load_model(path).buildings[0]
    monkeypatch.setattr(walls_module, "repair_roof_contacts", lambda *a, **k: pytest.fail("completed record repaired again"))
    assert block_walls(saved).current
    assert not update_walls(saved)
    assert blocks.walls_of(saved) == record


@pytest.mark.parametrize("persisted_partial", [False, True])
def test_recleaning_persists_the_budget_certifying_a_contact_retry(persisted_partial):
    from dtcc_core.builder.cleaning import blocks
    from dtcc_core.builder.cleaning.contract import check_labelled_contract

    block, *_ = _two_contact_block(persisted_partial=persisted_partial)
    cleaned = city_of(block).clean_footprints()
    block = cleaned.buildings[0]
    record = blocks.walls_of(block)
    walls = block_walls(block)
    assert record["fidelity_budget"] == .375
    assert record["detail"]["contact_repair"]["remaining_contacts"] == 0
    assert walls.current
    assert check_labelled_contract(blocks.part_sources(block), walls.order,
                                   blocks.block_region(block), walls.cells, delta=.5,
                                   epsilon=record["fidelity_budget"])["status"] != "fail"
    assert not update_walls(block)


@pytest.mark.parametrize("failure", ["construction_error", "contract_rejected", "unchanged"])
def test_failed_contact_diagnostics_persist_without_reassigning_cells(tmp_path, monkeypatch, failure):
    """Failed/no-progress attempts update the certificate without rebuilding geometry."""
    from dtcc_core.builder.cleaning import blocks, walls

    block, region, sources, order, cells, roofs = _two_contact_block()
    city = city_of(block)
    # Record generated prisms before the diagnostic-only upgrade.
    for part in block.building_parts:
        part.add_geometry(MultiSurface(surfaces=[Surface(vertices=np.array([
            [0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float))]), GeometryType.LOD1)
        blocks.record_lod1(part)
    footprints = [part.lod0 for part in block.building_parts]
    prisms = [part.lod1 for part in block.building_parts]
    budgets = []

    def failed(region, sources, order, cells, levels, *, delta, epsilon):
        budgets.append(epsilon)
        if failure == "construction_error":
            raise shapely.errors.GEOSException("controlled geometry failure")
        if failure == "contract_rejected":
            # A claimed improvement that removes all exposed labels must be rejected.
            return walls.LabelledRegion(region, [(region, order[0])], "strict", {},
                                        "complete", 0, epsilon, changed=True)
        return walls.LabelledRegion(region, cells, "strict", {}, "unchanged", 2, epsilon)

    monkeypatch.setattr(walls, "repair_roof_contacts", failed)
    view = block_walls(block)
    assert view.current and view.outcome == "strict"
    contact = view.report["contact_repair"]
    assert contact["state"] == "unchanged" and contact["remaining_contacts"] == 2
    assert contact["changed_area"] == 0
    assert [row["budget"] for row in contact["attempts"]] == [.25, .375]
    assert [row["state"] for row in contact["attempts"]] == [failure, failure]
    if failure == "construction_error":
        assert all(row["error"] == "controlled geometry failure" for row in contact["attempts"])
    assert not update_walls(block)
    assert all(part.lod0 is original for part, original in zip(block.building_parts, footprints))
    assert all(part.lod1 is original for part, original in zip(block.building_parts, prisms))
    record = json.loads(json.dumps(blocks.walls_of(block)))
    assert record["detail"]["contact_repair"] == contact
    assert record["fidelity_budget"] == .25
    path = tmp_path / "diagnostics.dtcc"
    io.save_model(city, path)
    loaded = io.load_model(path).buildings[0]
    calls = len(budgets)
    monkeypatch.setattr(walls, "repair_roof_contacts", lambda *a, **k: pytest.fail("unchanged record repaired again"))
    for _ in range(3):
        assert block_walls(loaded).current
        assert not update_walls(loaded)
        assert blocks.walls_of(loaded) == record
    assert len(budgets) == calls
    assert all(part.lod1 is not None for part in loaded.building_parts)


def test_later_success_replaces_failed_contact_diagnostics(monkeypatch):
    from dtcc_core.builder.cleaning import blocks, walls

    block, *_ = _two_contact_block()
    original = walls.repair_roof_contacts
    with monkeypatch.context() as failure:
        def fail(*args, **kwargs):
            raise shapely.errors.GEOSException("controlled geometry failure")
        failure.setattr(walls, "repair_roof_contacts", fail)
        assert not update_walls(block)
    assert blocks.walls_of(block)["detail"]["contact_repair"]["state"] == "unchanged"
    # A changed roof grouping requires a fresh contact check with current inputs.
    blocks.set_part_roof(block.building_parts[0], 22)
    assert update_walls(block)
    contact = blocks.walls_of(block)["detail"]["contact_repair"]
    assert contact["state"] == "complete" and contact["remaining_contacts"] == 0
    assert all(row["state"] in {"partial", "complete"} for row in contact["attempts"])
    assert "controlled geometry failure" not in json.dumps(contact)
    assert walls.repair_roof_contacts is original
    assert not update_walls(block)
