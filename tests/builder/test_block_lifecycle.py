"""Cleaned blocks through heights, LOD1, saving and cleaning again.

Regression tests for the review of issue #147: the order of cleaning and
height estimation, the lifecycle of parts, unique IDs with provenance, and
traversal of part geometry from a city.
"""

import json
from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import box, shape

from dtcc_core import io
from dtcc_core.builder import build_lod1_buildings
from dtcc_core.builder.cleaning.blocks import (
    MERGED_FROM_KEY,
    SOURCE_ID_KEY,
    WALLS_KEY,
    block_walls,
    set_part_roof,
    source_footprint,
    source_ids,
)
from dtcc_core.builder.cleaning.clean import RECORD_KEY, cleaning_record, footprint_cells
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.builder.geometry_builders.buildings import (
    compute_building_heights,
    extract_roof_points,
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
    Terrain,
)

DATA = Path(__file__).resolve().parents[1] / "data" / "cleaning"


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

    # The walls follow the heights: for meshing, LOD1 and cleaning again,
    # the same cells as when the heights came first.
    roofs, cells = cells_by_id(block)
    assert roofs == expected_roofs
    assert all(cells[id].equals(expected_cells[id]) for id in cells)
    assert not block_walls(block).current
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
