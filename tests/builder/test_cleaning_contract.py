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
    classify_candidate,
    fidelity,
    incident_sectors,
)
from dtcc_core.builder.cleaning.construction import construct
from dtcc_core.builder.geometry_builders.meshes import build_conditioned_footprints
from dtcc_core.model import Building, GeometryType, Surface
from tests.builder.cleaning_fixtures import survey_merge_group


@pytest.mark.parametrize(
    "polygons,admissible",
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
def test_admissible_subdivision(polygons, admissible):
    report = admissibility(polygons, 0.5)
    assert report["admissible"] is admissible
    assert report["topology_profile"] == "union"


def test_feature_size_counts_nonincident_vertex_edge_pairs_across_batches():
    # Only the last two rectangles are close; their facing vertices straddle
    # the query batch boundary: eight nonincident vertex-edge pairs define
    # mfs. The two vertex-vertex pairs are redundant (note, Observation 1) and
    # are counted only on request.
    polygons = [box(3 * i, 0, 3 * i + 2, 2) for i in range(64)]
    polygons.append(box(191.25, 0, 193.25, 2))
    report = admissibility(polygons, 0.5)
    assert report["subscale_pairs"] == 8
    assert "subscale_vertex_pairs" not in report
    assert report["mfs_capped_at_delta"] == 0.25
    assert report["closest_pair"] == [(191, 0), (191.25, 0)]
    finer = admissibility(polygons, 0.5, vertex_pairs=True)
    assert finer["subscale_pairs"] == 8 and finer["subscale_vertex_pairs"] == 2
    assert finer["mfs_capped_at_delta"] == report["mfs_capped_at_delta"]
    assert finer["admissible"] is report["admissible"] is False


def test_union_parts_shortcut_gives_the_same_report():
    # Parts of one union, including two that touch at a point (a pinch) and a
    # part with a hole, with and without the shortcut.
    union = unary_union(
        [
            box(0, 0, 2, 2),
            box(2, 2, 4, 4),
            box(4.2, 0, 6, 2),
            box(10, 0, 20, 10).difference(box(12, 2, 18, 8)),
        ]
    )
    parts = list(union.geoms)
    for vertex_pairs in (False, True):
        assert admissibility(
            parts, 0.5, vertex_pairs=vertex_pairs, union_parts=True
        ) == admissibility(parts, 0.5, vertex_pairs=vertex_pairs)


def test_achieved_budget_is_the_least_epsilon_the_output_satisfies():
    # Note, Figure 1: filling a 0.6 m wide inlet adds points at most 0.3 m
    # from the input, whatever the declared budget.
    inlet = box(0, 0, 4, 2.6).difference(box(1.7, 1.8, 2.3, 2.6))
    filled = box(0, 0, 4, 2.6)
    report = check_cleaning_contract(
        [inlet], [filled], delta=0.8, epsilon=0.4, achieved=True
    )["fidelity"]
    assert report["status"] == "pass"
    assert report["achieved_epsilon"] == pytest.approx(0.3, abs=0.005)
    assert 0 < report["achieved_epsilon_resolution"] < 0.005
    # Removing a 0.2 m band loses occupied points at most 0.2 m deep.
    report = fidelity([box(0, 0, 10, 10)], [box(0.2, 0, 10, 10)], 0.25, achieved=True)
    assert report["achieved_epsilon"] == pytest.approx(0.2, abs=0.004)
    unchanged = fidelity([box(0, 0, 10, 10)], [box(0, 0, 10, 10)], 0.25, achieved=True)
    assert unchanged["achieved_epsilon"] == 0.0
    failing = fidelity([box(0, 0, 10, 10)], [box(1, 0, 10, 10)], 0.25, achieved=True)
    assert failing["status"] == "fail" and failing["achieved_epsilon"] is None
    assert "achieved_epsilon" not in fidelity([box(0, 0, 1, 1)], [box(0, 0, 1, 1)], 0.25)


def test_candidates_are_classified_as_strict_warning_or_failed():
    raw = [box(0, 0, 5, 5), box(5.1, 0, 10, 5)]

    def classify(cleaned):
        return classify_candidate(
            check_cleaning_contract(raw, cleaned, delta=0.5),
            check_mesher_handoff_profile(cleaned),
        )

    assert classify([box(0, 0, 10, 5)]) == "strict"
    assert classify(raw) == "warning"
    assert classify([]) == "fail"
    touching = [box(0, 0, 5, 5), box(5, 5, 10, 10)]
    assert classify_candidate(
        check_cleaning_contract(touching, touching, delta=0.5),
        check_mesher_handoff_profile(touching),
    ) == "fail"


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
    for case_id, group_id in (
        ("city_grid:gothenburg:006", 84),
        ("city_grid:lund:016", 8),
    ):
        raw = survey_merge_group(case_id, group_id)
        output, construction = construct(raw)
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
    assert report["admissibility"]["admissible"]
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
            fidelity_budget=0.25,
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
