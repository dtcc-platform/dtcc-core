"""Geometric requirements, independent of repair branches and triangulations."""

import json
from pathlib import Path

import pytest
from shapely.affinity import translate
from shapely.geometry import Polygon, Point, box, shape
from shapely.ops import unary_union

from dtcc_core.builder.cleaning.contract import (
    admissibility,
    check_cleaning_contract,
    check_mesher_handoff_profile,
    fidelity,
    incident_sectors,
)
from dtcc_core.builder.geometry_builders.meshes import build_conditioned_footprints
from dtcc_core.model import Building, GeometryType, Surface


@pytest.mark.parametrize(
    "polygons,resolved",
    [
        ([box(0, 0, 2, 2), box(2, 0, 4, 2)], True),  # shared wall
        ([box(0, 0, 2, 2), box(2, 2, 4, 4)], False),  # point contact
        ([box(0, 0, 2, 2), box(1, 0, 3, 2)], False),  # overlapping interiors
        ([Polygon([(0, 0), (0.1, 0), (10, 0), (10, 6), (0, 6)])], True),
        ([Polygon([(0, 0), (0.1, 0.001), (10, 0), (10, 6), (0, 6)])], False),
        ([box(0, 0, 2, 2), box(2.2, 0, 4.2, 2)], False),  # narrow gap
        ([Polygon([(0, 0), (30, -2), (30, 2)])], True),  # genuine sharp corner
        ([], True),
    ],
)
def test_resolved_subdivision(polygons, resolved):
    assert admissibility(polygons, 0.5)["resolved"] is resolved


def test_separation_queries_count_pairs_across_batches_without_counting_incidence():
    # Only the last two rectangles are close; their facing vertices straddle
    # the query batch boundary. Two VV pairs and eight nonincident VE pairs.
    polygons = [box(3 * i, 0, 3 * i + 2, 2) for i in range(64)]
    polygons.append(box(191.25, 0, 193.25, 2))
    report = admissibility(polygons, 0.5)
    assert report["subscale_pairs"] == 10
    assert report["separation_capped_at_delta"] == 0.25
    assert report["closest_pair"] == [(191, 0), (191.25, 0)]


def test_fidelity_protects_both_occupied_and_open_cores():
    outer = box(0, 0, 10, 10)
    tiny_hole = Polygon(
        outer.exterior.coords, [box(4.9, 4.9, 5.1, 5.1).exterior.coords]
    )
    courtyard = Polygon(outer.exterior.coords, [box(3, 3, 7, 7).exterior.coords])
    assert fidelity([tiny_hole], [outer], 0.25)["status"] == "pass"
    assert fidelity([courtyard], [outer], 0.25)["status"] == "fail"
    assert fidelity([box(0, 0, 3, 3)], [], 0.25)["status"] == "fail"
    assert fidelity([box(0, 0, 0.2, 0.2)], [], 0.25)["status"] == "pass"
    assert fidelity([], [outer], 0.25)["status"] == "fail"


def test_incident_sector_profile_accounts_for_occupied_open_and_shared_sectors():
    acute = Polygon([(0, 0), (30, -0.2), (30, 0.2)])
    report = incident_sectors([acute], minimum_degrees=1.0)
    assert report["status"] == "fail"
    assert report["minimum_witness"]["kind"] == "occupied"
    assert report["minimum_angle_degrees"] < 1.0

    notched = Polygon(
        [(0, 0), (10, 0), (10, 10), (5.02, 10), (5, 0.1), (4.98, 10), (0, 10)]
    )
    report = incident_sectors([notched], minimum_degrees=1.0)
    assert report["status"] == "fail"
    assert report["minimum_witness"]["kind"] == "open"

    shared = incident_sectors([box(0, 0, 2, 2), box(2, 0, 4, 2)])
    assert shared["status"] == "pass"
    assert shared["occupied_sector_count"] > 0
    assert shared["open_sector_count"] > 0


def test_incident_sector_profile_is_stable_under_small_translation():
    polygon = Polygon([(0, 0), (30, -0.2), (30, 0.2)])
    origin = incident_sectors([polygon])
    shifted = incident_sectors([translate(polygon, 300000, 6000000)])
    assert shifted["status"] == origin["status"] == "fail"
    assert shifted["minimum_angle_degrees"] == pytest.approx(
        origin["minimum_angle_degrees"], abs=1e-7
    )


def test_full_survey_cusps_are_repaired_or_explicitly_removed_for_handoff():
    fixture = Path(__file__).parents[1] / "data/cleaning/full-survey-focused-groups.json"
    data = json.loads(fixture.read_text())
    from shapely import from_wkb
    from sandbox import cleaning_staged_probe as staged

    for case_id, group_id in (
        ("city_grid:gothenburg:006", 84),
        ("city_grid:lund:016", 8),
    ):
        case = next(row for row in data["cases"] if row["case_id"] == case_id)
        group = next(row for row in case["groups"] if row["group"] == group_id)
        raw = [from_wkb(bytes.fromhex(value)) for value in group["polygon_wkb_hex"]]
        output, construction = staged.construct(raw)
        assert construction["contract"]["status"] == "pass"
        assert check_mesher_handoff_profile(output)["status"] == "pass"
        assert construction["mesher_profile"]["operations"]
        if case_id == "city_grid:gothenburg:006":
            assert output == []
            assert construction["mesher_profile"][
                "empty_protected_core_removal"
            ] is True
        else:
            assert output


def test_fidelity_budget_is_independent_and_borderline_is_not_pass():
    raw = [box(0, 0, 10, 10)]
    cleaned = [box(-0.25, 0, 10, 10)]
    assert check_cleaning_contract(raw, cleaned, delta=0.5)["status"] == "borderline"
    assert (
        check_cleaning_contract(raw, cleaned, delta=0.5, epsilon=0.1)["status"]
        == "fail"
    )
    assert (
        check_cleaning_contract(raw, cleaned, delta=0.5, epsilon=0.5)["status"]
        == "pass"
    )


def test_fidelity_area_measurement_is_stable_at_large_map_coordinates():
    fixture = (
        Path(__file__).parents[1] / "data/cleaning/large-coordinate-overlay.geojson"
    )
    data = json.loads(fixture.read_text())
    protected, candidate = [shape(f["geometry"]) for f in data["features"]]
    world = fidelity([protected], [candidate], 0)
    x, y = data["metadata"]["local_origin"]
    local = fidelity([translate(protected, -x, -y)], [translate(candidate, -x, -y)], 0)
    # World-coordinate overlay previously reported zero occupied loss. The
    # reduced seven-vertex polygons retain the original narrow sliver.
    assert world["lost_protected_area"] > 1e-6
    for metric in ("lost_protected_area", "added_outside_budget_area"):
        assert world[metric] == pytest.approx(local[metric], abs=1e-11)
    assert world["status"] == local["status"] == "fail"


def test_invalid_inputs_are_explicit_and_invalid_output_cannot_pass():
    # GEOS interprets the retraced spur as a line, reported separately.
    raw = Polygon([(0, 0), (2, 0), (2, 2), (1, 2), (1, 3), (1, 2), (0, 2)])
    report = check_cleaning_contract([raw], [box(0, 0, 2, 2)], delta=0.5)
    assert report["status"] == "pass"
    assert report["fidelity"]["input_interpretation"]["nonarea_remnant_count"] == 1
    assert check_cleaning_contract([raw], [raw], delta=0.5)["status"] == "fail"
    with pytest.raises(ValueError, match="Polygon"):
        check_cleaning_contract([Point(0, 0)], [], delta=0.5)
    with pytest.raises(ValueError, match="delta"):
        check_cleaning_contract([], [], delta=0)
    with pytest.raises(ValueError, match="epsilon"):
        check_cleaning_contract([], [], delta=0.5, epsilon=float("nan"))


def test_helsingborg_courtyard_survives_final_cleaning():
    fixture = Path(__file__).parents[1] / "data/cleaning/helsingborg-courtyard.geojson"
    data = json.loads(fixture.read_text())
    raw = [shape(f["geometry"]) for f in data["features"]]
    occupied = unary_union(raw)
    courtyard = max(
        (Polygon(ring) for part in occupied.geoms for ring in part.interiors),
        key=lambda p: p.area,
    )
    buildings = []
    for polygon in raw:
        surface = Surface()
        surface.from_polygon(polygon)
        building = Building()
        building.add_geometry(surface, GeometryType.LOD0)
        buildings.append(building)
    result = build_conditioned_footprints(
        buildings,
        lod=GeometryType.LOD0,
        min_building_detail=0.5,
        min_building_area=15,
        merge_tolerance=0.5,
        merge_buildings=True,
        max_mesh_size=10,
        cleaning_diagnostics=False,
    )
    cleaned = [s.to_polygon(simplify=0) for s in result.surfaces]
    report = check_cleaning_contract(raw, cleaned, delta=0.5)
    assert report["admissibility"]["resolved"]
    # A protected interior point must remain open. This detects wholesale
    # courtyard deletion without blessing the remaining boundary deviations.
    open_core = courtyard.difference(occupied.buffer(report["epsilon"]))
    assert not unary_union(cleaned).covers(open_core.representative_point())
    assert sorted({i for group in result.source_map for i in group}) == list(
        range(len(raw))
    )


def test_public_cleaning_reports_fidelity_against_the_original_input():
    from dtcc_core.builder.cleaning import (
        ConditioningOptions,
        condition_polygon_coverage,
    )
    raw = box(0, 0, 10, 10)
    result = condition_polygon_coverage(
        [raw],
        options=ConditioningOptions(
            merge_distance=0,
            fidelity_tolerance=0.25,
            enable_logging=False,
        ),
    )
    assert unary_union(result.polygons).equals(raw)
    assert result.diagnostics["fidelity"]["status"] == "pass"
    assert result.diagnostics["before_selection_contract"]["status"] == "pass"


def test_cleaning_preserves_occupied_core_across_shared_walls():
    from dtcc_core.builder.cleaning import (
        ConditioningOptions,
        condition_polygon_coverage,
    )

    # Neither narrow parcel has its own protected core, but their union does.
    # Per-polygon opening must not delete the resolved building they form.
    raw = [box(0, 0, 0.3, 10), box(0.3, 0, 0.6, 10)]
    assert all(p.buffer(-0.25).is_empty for p in raw)
    result = condition_polygon_coverage(
        raw, options=ConditioningOptions(enable_logging=False)
    )
    assert check_cleaning_contract(raw, result.polygons, delta=0.5)["status"] == "pass"
    assert sorted({i for group in result.source_map for i in group}) == [0, 1]
