"""Walls between merged buildings: every part keeps a cleaned footprint."""

import json
import math
from pathlib import Path

import numpy as np
import pytest
import shapely
from shapely.geometry import LineString, Point, box, shape
from shapely.ops import unary_union

from dtcc_core.builder.cleaning import blocks as blocks_module
from dtcc_core.builder.cleaning import walls as walls_module
from dtcc_core.builder.cleaning.blocks import WALLS_KEY, set_part_roof
from dtcc_core.builder.cleaning.clean import MERGED_FROM_KEY, RECORD_KEY, clean_regions
from dtcc_core.builder.cleaning.contract import check_labelled_contract
from dtcc_core.builder.cleaning.walls import LabelledRegion, label_region
from dtcc_core.builder.cleaning.levels import level_regions, pinched_levels
from dtcc_core.builder.geometry_builders.meshes import build_city_flat_mesh_from_coverage
from dtcc_core.model.mixins.mesh.quality import tri_element_quality
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


_original_xy = walls_module._WallRepair._xy

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


def test_a_four_way_roof_contact_is_split_within_source_fidelity():
    sources = {
        0: box(0, 0, 10, 10), 1: box(10, 0, 20, 10),
        2: box(10, 10, 20, 20), 3: box(0, 10, 10, 20),
    }
    levels = {0: 20, 1: 6, 2: 12, 3: 6}
    order = [0, 2, 1, 3]
    region = box(0, 0, 20, 20)
    assert pinched_levels([(p, levels[i]) for i, p in sources.items()])
    result = blocks_module.wall_cells(
        region, sources, order, delta=0.5, epsilon=0.25,
        strict=False, levels=levels,
    )
    assert result.outcome in ("strict", "warning")
    assert result.report["contact_repair"]["remaining_contacts"] == 0
    assert result.report["contact_repair"]["changed_area"] < 10
    roofs = level_regions(result.cells, levels, 0.5)
    assert not pinched_levels([(p, z) for p, z, _ in roofs])
    assert {z for _, z, _ in roofs} == {6, 12, 20}
    assert {label for _, label in result.cells} == set(sources)
    contract = check_labelled_contract(
        sources, order, region, result.cells, delta=0.5,
        epsilon=result.report["fidelity_budget"],
    )
    assert contract["status"] in ("strict", "warning"), contract


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


def fallback(self):
    return LabelledRegion(self.region, [], "fallback", {"reason": "test"})


def test_failed_walls_fall_back_to_one_footprint(monkeypatch, caplog):
    monkeypatch.setattr(walls_module._WallRepair, "_xy", fallback)
    cleaned = city_of(
        building("a", box(0, 0, 10.2, 10), 10), building("b", box(10, 0, 20, 10), 11)
    ).clean_footprints()
    block = cleaned.buildings[0]
    assert all(p.get_geometry(GeometryType.LOD0) is None for p in block.building_parts)
    assert footprint(block).equals(box(0, 0, 20, 10))
    assert block.estimated_height == pytest.approx((10 * 102 + 11 * 100) / 202)
    walls = cleaned.attributes[RECORD_KEY]["walls"]
    assert walls["fallback"] == 1 and walls["fallback_ids"] == [["a", "b"]]
    assert block.attributes[WALLS_KEY]["outcome"] == "fallback"
    assert "keep one footprint" in caplog.text

    # Cleaning again changes nothing: the walls are current for the roofs.
    monkeypatch.undo()
    again = cleaned.clean_footprints()
    assert again.attributes[RECORD_KEY]["outcome"] == "unchanged"
    assert again.attributes[RECORD_KEY]["walls"]["fallback"] == 1
    assert again.buildings[0].attributes[MERGED_FROM_KEY] == ["a", "b"]

    # New roofs that do not change the winner of the overlap keep the walls;
    # new roofs that do recompute them from the sources.
    set_part_roof(parts_by_id(again.buildings[0])["b"], 13)
    assert again.clean_footprints().attributes[RECORD_KEY]["walls"]["fallback"] == 1
    set_part_roof(parts_by_id(again.buildings[0])["a"], 14)
    third = again.clean_footprints()
    assert third.attributes[RECORD_KEY]["walls"]["strict"] == 1
    parts = parts_by_id(third.buildings[0])
    assert footprint(parts["a"]).equals(box(0, 0, 10.2, 10))
    assert footprint(parts["b"]).equals(box(10.2, 0, 20, 10))
    assert parts["a"].get_geometry(GeometryType.LOD0).zmax == 14


def test_strict_cleaning_falls_back_instead_of_wall_warnings(monkeypatch):
    def warn(self):
        labelled = _original_xy(self)
        return LabelledRegion(labelled.region, labelled.cells, "warning", labelled.report)

    monkeypatch.setattr(walls_module._WallRepair, "_xy", warn)
    city = city_of(building("a", box(0, 0, 10, 10), 10), building("b", box(10, 0, 20, 10), 20))
    assert city.clean_footprints().attributes[RECORD_KEY]["walls"]["warning"] == 1
    walls = city.clean_footprints(strict=True).attributes[RECORD_KEY]["walls"]
    assert walls["warning"] == 0 and walls["fallback"] == 1


def test_a_geometry_failure_in_the_walls_leaves_only_that_block_one_footprint(
    monkeypatch, caplog
):
    def fail(self):
        sources = self.sources
        if 0 in sources:
            raise shapely.errors.GEOSException("TopologyException: test")
        return _original_xy(self)

    monkeypatch.setattr(walls_module._WallRepair, "_xy", fail)
    buildings = [
        building("a", box(0, 0, 10, 10), 10),
        building("b", box(10, 0, 20, 10), 20),
        building("c", box(40, 0, 50, 10), 10),
        building("d", box(50, 0, 60, 10), 20),
    ]
    regions = clean_regions(buildings, enable_logging=False)
    reports = {tuple(row["labels"]): row for row in regions.diagnostics["walls"]}
    assert reports[(0, 1)]["outcome"] == "fallback"
    assert reports[(0, 1)]["reason"] == "construction_error"
    assert reports[(2, 3)]["outcome"] == "strict"
    walls = city_of(*buildings).clean_footprints().attributes[RECORD_KEY]["walls"]
    assert walls["fallback_ids"] == [["a", "b"]] and walls["strict"] == 1
    assert "construction failed" in caplog.text


def test_walls_of_a_block_with_near_coincident_samples():
    # Helsingborg 007: Voronoi cells of samples a few nanometres apart did not
    # form an exact coverage, and their union raised a TopologyException.
    payload = json.loads((DATA / "helsingborg-wall-construction.geojson").read_text())
    meta = payload["metadata"]
    region = next(
        shape(f["geometry"]) for f in payload["features"] if f["properties"]["role"] == "region"
    )
    sources = {
        f["properties"]["label"]: shape(f["geometry"])
        for f in payload["features"]
        if f["properties"]["role"] == "source"
    }
    labelled = label_region(
        region,
        sources,
        meta["order"],
        delta=meta["delta"],
        epsilon=meta["epsilon"],
        retry_epsilon=meta["retry_epsilon"],
    )
    assert labelled.outcome == "strict"
    assert sorted(label for _, label in labelled.cells) == sorted(sources)
    union = unary_union([cell for cell, _ in labelled.cells])
    assert union.symmetric_difference(region).area < 1e-6


def _fixture(name):
    payload = json.loads((DATA / name).read_text())
    region = next(
        shape(f["geometry"]) for f in payload["features"] if f["properties"]["role"] == "region"
    )
    sources = {
        f["properties"]["label"]: shape(f["geometry"])
        for f in payload["features"]
        if f["properties"]["role"] == "source"
    }
    return payload["metadata"], region, sources


@pytest.mark.parametrize("name,from_sources", [
    ("orebro-007-boundary-junction.geojson", False),
    ("vasteras-057-boundary-junction.geojson", False),
    ("lund-056-boundary-contact.geojson", False),
    ("lund-003-boundary-opening.geojson", False),
    ("gothenburg-054-aligned-contact.geojson", False),
    ("uppsala-048-compatible-warning.geojson", True),
])
def test_local_wall_repairs_keep_stepped_roofs_and_good_closed_shells(name, from_sources):
    from dtcc_core.builder.cleaning.blocks import (
        assign_cells, set_source_footprint, walls_detail, walls_record,
    )
    from dtcc_core.builder.geometry_builders.meshes import build_building_shell
    from dtcc_core.model import BuildingPart

    meta, region, sources = _fixture(name)
    payload = json.loads((DATA / name).read_text())
    roofs = {f["properties"]["label"]: f["properties"]["roof"]
             for f in payload["features"] if f["properties"]["role"] == "source"}
    initial = [(shape(f["geometry"]), f["properties"]["label"])
               for f in payload["features"] if f["properties"]["role"] == "cell"]
    if from_sources:
        initial = []
    if initial:
        assert pinched_levels([(p, z) for p, z, _ in level_regions(initial, roofs, meta["delta"])])
        repair = walls_module.condition_walls(
            region, sources, meta["order"], levels=roofs,
            delta=meta["delta"], epsilon=meta["epsilon"],
            previous=LabelledRegion(region, initial, "warning"),
        )
        assert repair.contact_outcome == "complete" and repair.remaining_contacts == 0
        result = repair
        expected_area = sum(
            shapely.union_all([p for p, i in initial if i == label]).symmetric_difference(
                shapely.union_all([p for p, i in result.cells if i == label])
            ).area for label in sources
        ) / 2
        assert result.report["contact_repair"]["changed_area"] == expected_area
    else:
        result = blocks_module.wall_cells(
            region, sources, meta["order"], delta=meta["delta"],
            epsilon=meta["epsilon"], strict=False, levels=roofs,
        )
    assert result.outcome in ("strict", "warning")
    contract = check_labelled_contract(
        sources, meta["order"], region, result.cells, delta=meta["delta"],
        epsilon=result.report["fidelity_budget"],
    )
    assert contract["status"] in ("strict", "warning"), contract
    regions = level_regions(result.cells, roofs, meta["delta"])
    assert len({z for _, z, _ in regions}) > 1
    assert not pinched_levels([(p, z) for p, z, _ in regions])
    block = Building(id="block", attributes={MERGED_FROM_KEY: list(map(str, sources))})
    surface = Surface()
    surface.from_polygon(region, max(roofs.values()))
    block.add_geometry(surface, GeometryType.LOD0)
    parts = []
    for label, source in sources.items():
        part = BuildingPart(id=str(label))
        set_source_footprint(part, source, roofs[label])
        parts.append(part)
    block.add_children(parts)
    assign_cells(parts, result.cells)
    block.attributes[WALLS_KEY] = walls_record(
        block, meta["order"], "roof", result.outcome, delta=meta["delta"],
        epsilon=result.report["fidelity_budget"], strict=False,
        detail=walls_detail(result.outcome, result.report),
    )
    mesh = build_building_shell(block, max_mesh_size=10)
    quality = tri_element_quality(mesh.vertices, mesh.faces)
    assert np.count_nonzero(quality < 0.02) == 0
    edges = np.concatenate([mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]]])
    _, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    assert np.all(counts == 2)
    points = mesh.vertices[mesh.faces]
    volume = np.einsum("ij,ij->i", points[:, 0], np.cross(points[:, 1], points[:, 2])).sum() / 6
    assert volume == pytest.approx(sum(p.area * z for p, z, _ in regions), rel=1e-8)


@pytest.mark.parametrize(
    "name",
    [
        "stockholm-069-thin-piece.geojson",
        "stockholm-085-short-wall.geojson",
        "stockholm-085-junction-near-boundary.geojson",
        "orebro-085-thin-tail.geojson",
        "gothenburg-008-wall-near-corners.geojson",
        "linkoping-056-wall-near-vertex.geojson",
    ],
)
def test_walls_of_real_blocks_have_no_sub_delta_defects(name):
    # Real blocks of the surface survey whose walls left a piece of a cell of
    # about a square millimetre, two junctions 0.18 mm apart, a junction
    # 1.8 mm from a region vertex, a wall running millimetres from a region
    # edge, or a wall passing 2 to 10 cm from region vertices. The 2D mesher
    # failed on the first two (coverage seed conflicts), refined the fourth
    # into 27,000 triangles, and the last two gave tiny ground triangles,
    # which steep terrain lifts into slivers. Thin pieces go to a neighbour,
    # short walls contract, junctions move onto the boundary, and walls are
    # routed through region vertices next to them.
    meta, region, sources = _fixture(name)
    labelled = blocks_module.wall_cells(
        region, sources, meta["order"], delta=meta["delta"], epsilon=meta["epsilon"], strict=False
    )
    assert labelled.outcome == "strict"
    labels = [label for _, label in labelled.cells]
    assert len(labels) == len(set(labels))
    union = unary_union([cell for cell, _ in labelled.cells])
    assert union.symmetric_difference(region).area < 1e-6
    report = check_labelled_contract(
        sources,
        meta["order"],
        region,
        labelled.cells,
        delta=meta["delta"],
        epsilon=labelled.report["fidelity_budget"],
    )
    assert report["status"] == "strict", report
    x0, y0, x1, y1 = region.bounds
    polygons = [cell for cell, _ in labelled.cells]
    mesh = build_city_flat_mesh_from_coverage(
        region_polygons=polygons,
        region_markers=list(range(len(polygons))),
        bounds=(x0 - 10, y0 - 10, x1 + 10, y1 + 10),
        max_mesh_size=10,
        min_mesh_angle=25,
        backend="dtcc_mesher",
    )
    quality = tri_element_quality(mesh.vertices, mesh.faces)
    assert len(mesh.faces) < 1000
    assert quality.min() > 0.02


def test_site_offers_without_dispensable_labels_build_no_reach(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("site reach built without a dispensable label")

    region = box(0, 0, 20, 10)
    cells = [(box(0, 0, 10, 10), 0), (box(10, 0, 20, 10), 1)]
    monkeypatch.setattr(walls_module, "MultiPoint", forbidden)
    list(walls_module._site_offers(region, cells, [(10.0, 5.0)], 0.5, set()))
    monkeypatch.undo()
    calls = []
    multipoint = walls_module.MultiPoint
    monkeypatch.setattr(
        walls_module, "MultiPoint", lambda points: calls.append(points) or multipoint(points)
    )
    list(walls_module._site_offers(region, cells, [(10.0, 5.0)], 0.5, {1}))
    assert calls


def _greedy_sites(evaluation, radius):
    """The original quadratic grouping, as the reference for the grid index."""
    groups = [[tuple(vertex)] for vertex, _ in evaluation["sectors"]]
    for vertex, edge, _ in evaluation["pairs"]:
        line = LineString(edge)
        foot = line.interpolate(line.project(Point(vertex)))
        groups.append([tuple(vertex), (foot.x, foot.y)])
    sites = []
    for group in groups:
        for site in sites:
            if min(math.hypot(p[0] - q[0], p[1] - q[1]) for p in group for q in site) <= radius:
                site.extend(group)
                break
        else:
            sites.append(list(group))
    return [list(dict.fromkeys(site)) for site in sites]


@pytest.mark.parametrize("seed", range(5))
def test_sites_group_defects_as_the_greedy_rule(seed):
    rng = np.random.default_rng(seed)
    points = np.round(rng.uniform(0, 30, size=(120, 2)), 1)
    evaluation = {
        "sectors": [(tuple(p), 0.1) for p in points[:60]],
        "pairs": [
            (tuple(p), (tuple(q), tuple(q + rng.uniform(-2, 2, 2))), 0.2)
            for p, q in zip(points[60:90], points[90:])
        ],
    }
    assert walls_module._sites(evaluation, 1.5) == _greedy_sites(evaluation, 1.5)


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


def test_walls_closer_than_a_tenth_of_a_millimetre_fall_back(monkeypatch):
    def degenerate(self):
        labelled = _original_xy(self)
        return LabelledRegion(
            labelled.region, labelled.cells, "warning", {**labelled.report, "min_wall_distance": 5e-5}
        )

    monkeypatch.setattr(walls_module._WallRepair, "_xy", degenerate)
    regions = clean_regions(
        [building("a", box(0, 0, 10, 10), 10), building("b", box(10, 0, 20, 10), 20)],
        enable_logging=False,
    )
    (report,) = regions.diagnostics["walls"]
    assert (report["outcome"], report["reason"]) == ("fallback", "degenerate_walls")


def test_short_boundary_intervals_are_noded_before_unknown_roofs_are_set(monkeypatch):
    from dtcc_core.builder.cleaning.clean import footprint_cells
    from dtcc_core.builder.cleaning.blocks import (
        assign_cells, set_source_footprint, walls_record, update_walls,
    )
    from dtcc_core.builder.geometry_builders.meshes import build_building_shell
    from dtcc_core.model import BuildingPart

    name = "vasteras-053-world-noding.geojson"
    meta, region, sources = _fixture(name)
    payload = json.loads((DATA / name).read_text())
    roofs = {f["properties"]["label"]: f["properties"]["roof"]
             for f in payload["features"] if f["properties"]["role"] == "source"}
    recorded = [(shape(f["geometry"]), f["properties"]["label"])
                for f in payload["features"] if f["properties"]["role"] == "cell"]
    assert check_labelled_contract(
        sources, meta["order"], region, recorded, delta=meta["delta"],
        epsilon=meta["epsilon"],
    )["status"] == "strict"
    from dtcc_core.builder.cleaning import walls
    origin = region.bounds[:2]

    def recorded_construction(*args, **kwargs):
        # Retain the recorded arrangement to test world-coordinate noding
        # independently of which equivalent construction is selected.
        cells = [(walls._to_local(p, origin), label) for p, label in recorded]
        return cells, {"outcome": "strict", "small_sectors": 0,
                       "severe_wall_pairs": 0, "wall_pairs": 0, "wall_deficit": 0}

    monkeypatch.setattr(walls, "_repair", recorded_construction)
    result = blocks_module.wall_cells(
        region, sources, meta["order"], delta=meta["delta"],
        epsilon=meta["epsilon"], strict=False, levels=dict.fromkeys(sources, 0),
    )
    assert result.cells
    block = Building(id="block", attributes={MERGED_FROM_KEY: list(map(str, sources))})
    footprint = Surface()
    footprint.from_polygon(region, 0)
    block.add_geometry(footprint, GeometryType.LOD0)
    parts = []
    for label, source in sources.items():
        part = BuildingPart(id=str(label))
        set_source_footprint(part, source, 0)
        parts.append(part)
    block.add_children(parts)
    assign_cells(parts, result.cells)
    block.attributes[WALLS_KEY] = walls_record(
        block, meta["order"], "area", result.outcome, delta=meta["delta"],
        epsilon=result.report["fidelity_budget"], strict=False,
    )
    for label, part in enumerate(parts):
        set_part_roof(part, roofs[label])
    update_walls(block)
    cells, _ = footprint_cells(block)
    order = sorted(sources, key=lambda label: -roofs[label])
    contract = check_labelled_contract(
        sources, order, region, cells, delta=meta["delta"],
        epsilon=result.report["fidelity_budget"],
    )
    assert contract["status"] in ("strict", "warning"), contract
    mesh = build_building_shell(block, max_mesh_size=10)
    assert np.count_nonzero(tri_element_quality(mesh.vertices, mesh.faces) < 0.02) == 0


def test_failed_xy_retry_retains_the_valid_earlier_partition(monkeypatch):
    original = walls_module._repair
    sources = {0: box(0, 0, 10, 10), 1: box(10, 0, 20, 10)}
    region = box(0, 0, 20, 10)

    calls = []

    def failed_retry(region, exact, judge, delta, **kwargs):
        calls.append(judge)
        if len(calls) > 1:
            raise shapely.errors.GEOSException("controlled later XY failure")
        cells, report = original(region, exact, judge, delta, **kwargs)
        return cells, {**report, "outcome": "warning", "wall_pairs": 1}

    monkeypatch.setattr(walls_module, "_repair", failed_retry)
    result = label_region(region, sources, [1, 0], delta=.5, epsilon=.25, retry_epsilon=.375)
    assert result.outcome == "warning"
    assert result.accepted_budget == .25
    assert len(result.cells) == 2
    assert result.attempts[-1]["outcome"] == "construction_error"
    assert result.attempts[-1]["error"] == "controlled later XY failure"
    assert check_labelled_contract(sources, [1, 0], region, result.cells,
                                   delta=.5, epsilon=.25)["status"] == "strict"
    calls.append(None)
    failed = label_region(region, sources, [1, 0], delta=.5, epsilon=.25, retry_epsilon=.375)
    assert failed.outcome == "fallback" and failed.report["reason"] == "construction_error"
    assert failed.report["error"] == "controlled later XY failure"
    assert [a["fidelity_budget"] for a in failed.report["attempts"]] == [.25, .375]
