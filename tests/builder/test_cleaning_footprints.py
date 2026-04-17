import itertools
from types import SimpleNamespace

import pytest
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, box

import dtcc_core.builder.cleaning as cleaning
import dtcc_core.builder.cleaning.footprints as cleaning_footprints
from dtcc_core.model import Building, GeometryType, Surface


def make_building(polygon: Polygon, height: float = 10.0) -> Building:
    surface = Surface()
    surface.from_polygon(polygon, height)
    building = Building()
    building.add_geometry(surface, GeometryType.LOD0)
    building.attributes["height"] = height
    building.attributes["ground_height"] = 0.0
    return building


def exterior_vertex_set(polygon: Polygon) -> set[tuple[float, float]]:
    return {
        (round(float(x), 12), round(float(y), 12))
        for x, y, *_ in list(polygon.exterior.coords)[:-1]
    }


def assert_conditioning_invariants(
    result: cleaning.ConditioningResult,
    options: cleaning.ConditioningOptions,
) -> None:
    required_keys = {
        "input_count",
        "atomic_input_count",
        "output_count",
        "merged_group_count",
        "repaired_invalid_count",
        "collapsed_count",
        "dropped_small_count",
        "overlap_area_before",
        "overlap_area_after",
        "min_clearance_before",
        "min_clearance_after",
        "polygonize_cut_edges",
        "polygonize_dangles",
        "polygonize_invalid_rings",
        "geos_exception_count",
        "geos_exception_messages",
        "discarded_non_polygon_parts",
        "enable_logging",
        "precision_grid",
        "output_grid",
        "opening_radius",
        "closing_radius",
        "output_lattice_applied",
        "output_canonicalization_candidate_count",
        "output_canonicalization_output_grid_applied_count",
        "output_canonicalization_fallback_count",
        "output_canonicalization_reference_minus_candidate_area",
        "output_canonicalization_candidate_minus_reference_area",
        "output_canonicalization_signed_area_delta",
        "output_canonicalization_area_balance_budget",
        "final_min_area_filter_removed_count",
        "final_min_area_filter_removed_area",
        "global_reconstruction_applied",
        "global_reconstruction_reason",
        "global_reconstruction_overlap_threshold",
        "multi_component_group_count",
        "component_source_reassignment_count",
        "clearance_regularization_threshold",
        "clearance_regularization_applied",
        "clearance_regularization_candidate_count",
        "clearance_regularization_improved_count",
        "clearance_regularization_failed_count",
        "clearance_regularization_overlap_area",
        "source_coordinate_recovery_applied",
        "source_coordinate_recovery_candidate_count",
        "source_coordinate_recovery_applied_count",
        "source_coordinate_recovery_exact_count",
        "source_coordinate_recovery_vertex_count",
        "source_coordinate_recovery_contract_reject_count",
        "source_coordinate_recovery_short_edge_count_before",
        "source_coordinate_recovery_short_edge_count_after",
        "source_coordinate_recovery_pair_issue_count_before",
        "source_coordinate_recovery_pair_issue_count_after",
        "source_coordinate_recovery_ring_contact_count_before",
        "source_coordinate_recovery_ring_contact_count_after",
        "source_coordinate_recovery_min_clearance_before",
        "source_coordinate_recovery_min_clearance_after",
        "source_coordinate_recovery_reference_minus_candidate_area",
        "source_coordinate_recovery_candidate_minus_reference_area",
        "source_coordinate_recovery_signed_area_delta",
        "source_coordinate_recovery_support_reference_minus_candidate_area",
        "source_coordinate_recovery_support_candidate_minus_reference_area",
        "source_coordinate_recovery_support_signed_area_delta",
        "source_coordinate_recovery_rejected_overlap_count",
        "source_coordinate_recovery_rejected_non_improving_count",
        "source_coordinate_recovery_operator_attempts",
        "source_coordinate_recovery_operator_applied",
        "coverage_meshing_regularization_applied",
        "coverage_meshing_regularization_selected_branch",
        "coverage_meshing_regularization_short_edge_count_before",
        "coverage_meshing_regularization_short_edge_count_after",
        "coverage_meshing_regularization_pair_issue_count_before",
        "coverage_meshing_regularization_pair_issue_count_after",
        "coverage_meshing_regularization_ring_contact_count_before",
        "coverage_meshing_regularization_ring_contact_count_after",
        "coverage_meshing_regularization_ring_contact_polygon_count",
        "coverage_meshing_regularization_ring_contact_component_count",
        "coverage_meshing_regularization_ring_contact_failed_count",
        "coverage_meshing_regularization_ring_contact_area_delta",
        "coverage_meshing_regularization_reference_minus_candidate_area",
        "coverage_meshing_regularization_candidate_minus_reference_area",
        "coverage_meshing_regularization_signed_area_delta",
        "coverage_meshing_regularization_operator_attempts",
        "coverage_meshing_regularization_operator_applied",
        "small_component_absorb_applied",
        "small_component_absorb_candidate_count",
        "small_component_absorb_applied_count",
        "small_component_absorb_component_count",
        "small_component_absorb_edit_zone_area",
        "small_component_absorb_change_outside_edit_zone",
        "small_component_absorb_reference_minus_candidate_area",
        "small_component_absorb_candidate_minus_reference_area",
        "small_component_absorb_signed_area_delta",
        "small_component_absorb_rejected_nonlocal_count",
        "small_component_absorb_rejected_non_improving_count",
        "small_component_absorb_operator_attempts",
        "small_component_absorb_operator_applied",
        "coverage_simplify_short_edge_count_before",
        "coverage_simplify_short_edge_count_after",
        "coverage_simplify_edit_zone_area",
        "coverage_simplify_change_outside_edit_zone",
        "coverage_simplify_reference_minus_candidate_area",
        "coverage_simplify_candidate_minus_reference_area",
        "coverage_simplify_signed_area_delta",
        "coverage_simplify_area_balance_budget",
        "coverage_simplify_rejected_area_imbalance",
        "coverage_simplify_patch_count",
        "coverage_simplify_patch_applied_count",
        "coverage_simplify_operator_attempts",
        "coverage_simplify_operator_applied",
        "coverage_simplify_selected_branch",
        "local_defect_repair_candidate_count",
        "local_defect_repair_applied_count",
        "local_defect_repair_short_edge_count_before",
        "local_defect_repair_short_edge_count_after",
        "local_defect_repair_edit_zone_area",
        "local_defect_repair_change_outside_edit_zone",
        "local_defect_repair_reference_minus_candidate_area",
        "local_defect_repair_candidate_minus_reference_area",
        "local_defect_repair_signed_area_delta",
        "local_defect_repair_area_balance_budget",
        "local_defect_repair_rejected_area_imbalance_count",
        "source_reclaim_candidate_count",
        "source_reclaim_applied_count",
        "source_reclaim_component_count",
        "source_reclaim_short_edge_count_before",
        "source_reclaim_short_edge_count_after",
        "source_reclaim_edit_zone_area",
        "source_reclaim_change_outside_edit_zone",
        "source_reclaim_reference_minus_candidate_area",
        "source_reclaim_candidate_minus_reference_area",
        "source_reclaim_signed_area_delta",
        "polygon_simplify_candidate_count",
        "polygon_simplify_applied_count",
        "polygon_simplify_short_edge_count_before",
        "polygon_simplify_short_edge_count_after",
        "polygon_simplify_edit_zone_area",
        "polygon_simplify_change_outside_edit_zone",
        "polygon_simplify_reference_minus_candidate_area",
        "polygon_simplify_candidate_minus_reference_area",
        "polygon_simplify_signed_area_delta",
        "polygon_simplify_area_balance_budget",
        "polygon_simplify_rejected_area_imbalance_count",
        "stage_metrics",
    }
    assert required_keys.issubset(result.diagnostics)
    assert result.diagnostics["output_count"] == len(result.polygons)
    assert {
        "atomic_input",
        "opened",
        "regularized_groups",
        "reconstructed",
        "presimplify",
        "local_defect_repaired",
        "coverage_simplified",
        "source_reclaimed",
        "boundary_regularized",
        "clearance_regularized",
        "source_coordinate_recovered",
        "coverage_meshing_regularized",
        "final_output",
    }.issubset(result.diagnostics["stage_metrics"])
    assert result.diagnostics["output_grid"] >= result.diagnostics["precision_grid"]

    for polygon in result.polygons:
        assert isinstance(polygon, Polygon)
        assert polygon.is_valid
        assert polygon.area + 1e-9 >= options.min_area
        for hole in polygon.interiors:
            assert Polygon(hole).area + 1e-9 >= options.min_hole_area

    for left, right in itertools.combinations(result.polygons, 2):
        assert left.intersection(right).area == pytest.approx(0.0, abs=1e-9)

    for indices in result.source_map:
        assert indices == sorted(set(indices))


@pytest.mark.parametrize(
    ("geometries", "options", "expected_count"),
    [
        (
            [Polygon([(0, 0), (0, 4), (4, 4), (4, 0), (0, 0)][::-1])],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            1,
        ),
        (
            [Polygon(shell=[(0, 0), (6, 0), (6, 6), (0, 6)], holes=[[(4, 2), (2, 2), (2, 4), (4, 4), (4, 2)]])],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            1,
        ),
        (
            [Polygon([(0, 0), (3, 3), (0, 3), (3, 0), (0, 0)])],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            2,
        ),
        (
            [Polygon([(0, 0), (4, 0), (4, 0), (4, 4), (0, 4), (0, 0)])],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            1,
        ),
        (
            [Polygon([(0, 0), (4, 0), (4, 0), (4, 2), (4, 4), (0, 4), (0, 0)])],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            1,
        ),
        (
            [box(0, 0, 4, 4), box(3, 0, 7, 4)],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            1,
        ),
        (
            [box(0, 0, 4, 4), box(4.1, 0, 8.1, 4)],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.2,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            1,
        ),
        (
            [box(0, 0, 4, 4), box(4.4, 0, 8.4, 4)],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.2,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            2,
        ),
        (
            [Polygon([(0, 0), (4, 0), (4, 1.5), (6, 1.6), (4, 1.7), (4, 4), (0, 4)])],
            cleaning.ConditioningOptions(
                min_feature_size=0.5,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            1,
        ),
        (
            [Polygon([(0, 0), (2, 0), (2, 1.0), (4, 1.0), (4, 0), (6, 0), (6, 2), (4, 2), (4, 1.1), (2, 1.1), (2, 2), (0, 2)])],
            cleaning.ConditioningOptions(
                min_feature_size=0.4,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            2,
        ),
        (
            [Polygon([(0, 0), (4, 0), (4, 4), (0, 4), (0, 2.2), (-0.3, 2.05), (0, 1.9)])],
            cleaning.ConditioningOptions(
                min_feature_size=0.4,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            1,
        ),
        (
            [Polygon(shell=[(0, 0), (10, 0), (10, 10), (0, 10)], holes=[[(4.9, 4.9), (5.1, 4.9), (5.1, 5.1), (4.9, 5.1), (4.9, 4.9)]])],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.25,
            ),
            1,
        ),
        (
            [Polygon(shell=[(0, 0), (10, 0), (10, 10), (0, 10)], holes=[[(3, 3), (7, 3), (7, 7), (3, 7), (3, 3)]])],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.25,
            ),
            1,
        ),
        (
            [box(0, 0, 4, 4), box(4, 0, 8, 4)],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            2,
        ),
        (
            [MultiPolygon([box(0, 0, 2, 2), box(4, 0, 6, 2)])],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=0.0,
            ),
            2,
        ),
        (
            [box(0, 0, 1, 1)],
            cleaning.ConditioningOptions(
                min_feature_size=0.0,
                merge_distance=0.0,
                min_area=2.0,
                min_hole_area=0.0,
            ),
            0,
        ),
    ],
)
def test_condition_polygon_coverage_handles_pathologies(
    geometries,
    options,
    expected_count,
):
    result = cleaning.condition_polygon_coverage(geometries, options=options)
    assert len(result.polygons) == expected_count
    assert_conditioning_invariants(result, options)


def test_condition_polygon_coverage_accepts_geometry_collection_and_source_map():
    geometries = [
        GeometryCollection([box(0, 0, 2, 2), box(4, 0, 6, 2)]),
        box(8, 0, 10, 2),
    ]
    result = cleaning.condition_polygon_coverage(
        geometries,
        source_map=[[7, 3, 7], [2]],
        options=cleaning.ConditioningOptions(
            min_feature_size=0.0,
            merge_distance=0.0,
            min_area=0.0,
            min_hole_area=0.0,
        ),
    )
    assert len(result.polygons) == 3
    assert result.source_map[0] == [3, 7]
    assert_conditioning_invariants(
        result,
        cleaning.ConditioningOptions(
            min_feature_size=0.0,
            merge_distance=0.0,
            min_area=0.0,
            min_hole_area=0.0,
        ),
    )


def test_condition_polygon_coverage_regularizes_hole_touching_exterior():
    touching_hole = Polygon(
        [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)],
        [[(0, 5), (2, 4), (3, 5), (2, 6), (0, 5)]],
    )
    assert touching_hole.is_valid
    assert cleaning_footprints._polygon_has_ring_boundary_contacts(touching_hole)

    options = cleaning.ConditioningOptions(
        min_feature_size=0.0,
        merge_distance=0.0,
        min_area=0.0,
        min_hole_area=0.0,
        precision_grid=0.125,
    )
    result = cleaning.condition_polygon_coverage([touching_hole], options=options)

    assert len(result.polygons) == 1
    assert not cleaning_footprints._polygon_has_ring_boundary_contacts(
        result.polygons[0]
    )
    assert result.diagnostics["local_defect_repair_applied"] is True
    assert result.diagnostics["local_defect_repair_operator_applied"] == {
        "ring_contact_fill": 1
    }
    assert_conditioning_invariants(result, options)


def test_condition_polygon_coverage_repairs_hole_touching_exterior_upstream():
    touching_hole = Polygon(
        [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)],
        [[(0, 5), (2, 4), (3, 5), (2, 6), (0, 5)]],
    )
    assert touching_hole.is_valid
    assert cleaning_footprints._polygon_has_ring_boundary_contacts(touching_hole)

    options = cleaning.ConditioningOptions(
        min_feature_size=0.5,
        merge_distance=0.0,
        min_area=0.0,
        min_hole_area=0.0,
        precision_grid=0.03125,
    )
    result = cleaning.condition_polygon_coverage([touching_hole], options=options)

    assert len(result.polygons) == 1
    assert not cleaning_footprints._polygon_has_ring_boundary_contacts(
        result.polygons[0]
    )
    assert result.diagnostics["local_defect_repair_applied"] is True
    assert result.diagnostics["local_defect_repair_operator_applied"] == {
        "ring_contact_fill": 1
    }
    assert (
        result.diagnostics["coverage_meshing_regularization_ring_contact_count_before"]
        == 0
    )
    assert_conditioning_invariants(result, options)


def test_regularize_coverage_contacts_merges_point_touching_polygons():
    first = Polygon([(0, 0), (2, 0), (2, 2), (0, 2), (0, 0)])
    second = Polygon([(2, 2), (4, 2), (4, 4), (2, 4), (2, 2)])
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_coverage_contacts(
        [first, second],
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert len(polygons) == 1
    assert source_map == [[0, 1]]
    assert diagnostics["coverage_contact_regularization_pair_issue_count_before"] == 1
    assert diagnostics["coverage_contact_regularization_pair_issue_count_after"] == 0
    assert diagnostics["coverage_contact_regularization_applied"] is True
    assert any(
        key.startswith("coverage_pair_issue_point_")
        for key in diagnostics["coverage_contact_regularization_operator_applied"]
    )


def test_point_touch_bridge_operator_creates_bounded_bridge_neck():
    first = Polygon([(0, 0), (2, 0), (2, 2), (0, 2), (0, 0)])
    second = Polygon([(2, 2), (4, 2), (4, 4), (2, 4), (2, 2)])

    candidate = cleaning_footprints._apply_local_point_touch_bridge_operator(
        [first, second],
        [[0], [1]],
        radius=0.25,
        target_scale=0.5,
        grid=0.03125,
        diagnostics=cleaning_footprints._empty_diagnostics(2),
    )

    assert candidate is not None
    polygons, source_map = candidate
    assert len(polygons) == 1
    assert source_map == [[0, 1]]

    coords = list(polygons[0].exterior.coords)
    diagonal_segment_count = 0
    for start, end in zip(coords, coords[1:]):
        dx = abs(float(end[0] - start[0]))
        dy = abs(float(end[1] - start[1]))
        if dx > 1e-9 and dy > 1e-9:
            diagonal_segment_count += 1

    assert diagonal_segment_count >= 2
    assert polygons[0].area > first.area + second.area
    assert polygons[0].area < first.area + second.area + 0.25


def test_close_pair_bridge_operator_uses_same_bounded_neck_family():
    first = Polygon([(0, 0), (2, 0), (2, 2), (0, 2), (0, 0)])
    second = Polygon([(2.1, 2.1), (4.1, 2.1), (4.1, 4.1), (2.1, 4.1), (2.1, 2.1)])

    candidate = cleaning_footprints._apply_local_close_pair_bridge_operator(
        [first, second],
        [[0], [1]],
        radius=0.25,
        target_scale=0.5,
        grid=0.03125,
        diagnostics=cleaning_footprints._empty_diagnostics(2),
    )

    assert candidate is not None
    polygons, source_map = candidate
    assert len(polygons) == 1
    assert source_map == [[0, 1]]

    coords = list(polygons[0].exterior.coords)
    diagonal_segment_count = 0
    for start, end in zip(coords, coords[1:]):
        dx = abs(float(end[0] - start[0]))
        dy = abs(float(end[1] - start[1]))
        if dx > 1e-9 and dy > 1e-9:
            diagonal_segment_count += 1

    assert diagonal_segment_count >= 2
    assert polygons[0].area > first.area + second.area
    assert polygons[0].area < first.area + second.area + 0.25


def test_close_pair_bridge_operator_resolves_slanted_near_threshold_gap():
    left = Polygon(
        [
            (673632.179, 6581808.32),
            (673629.6246381734, 6581812.808878513),
            (673642.72, 6581820.183),
            (673639.603, 6581825.66),
            (673621.031, 6581815.061),
            (673624.157, 6581809.73),
            (673619.019, 6581806.814),
            (673621.586, 6581802.304),
            (673632.179, 6581808.32),
        ]
    )
    right = Polygon(
        [
            (673633.884, 6581814.606),
            (673634.669, 6581813.149),
            (673643.665, 6581818.21),
            (673642.821, 6581819.699),
            (673633.884, 6581814.606),
        ]
    )

    candidate = cleaning_footprints._apply_local_close_pair_bridge_operator(
        [left, right],
        [[0], [1]],
        radius=0.5,
        target_scale=0.5,
        grid=0.03125,
        diagnostics=cleaning_footprints._empty_diagnostics(2),
    )

    assert candidate is not None
    polygons, source_map = candidate
    assert source_map == [[0, 1]]
    signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
    assert signature.pair_issue_count == 0
    assert signature.short_edge_count == 0


def test_regularize_coverage_contacts_merges_courtyard_close_pair():
    outer = Polygon(
        shell=[(0, 0), (8, 0), (8, 8), (0, 8), (0, 0)],
        holes=[[(2, 2), (6, 2), (6, 6), (2, 6), (2, 2)]],
    )
    inner = Polygon(
        [(2.05, 3.0), (2.55, 3.0), (2.55, 3.5), (2.05, 3.5), (2.05, 3.0)]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_coverage_contacts(
        [outer, inner],
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert len(polygons) == 1
    assert source_map == [[0, 1]]
    signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
    assert signature.pair_issue_count == 0
    assert any(
        key.startswith("coverage_pair_issue_bridge_")
        for key in diagnostics["coverage_contact_regularization_operator_applied"]
    )


def test_regularize_coverage_contacts_records_post_contact_local_repair_metrics():
    first = Polygon([(0, 0), (2, 0), (2, 2), (0, 2), (0, 0)])
    second = Polygon([(2, 2), (4, 2), (4, 4), (2, 4), (2, 2)])
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    cleaning_footprints._regularize_coverage_contacts(
        [first, second],
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert "post_contact_local_defect_repair_short_edge_count_before" in diagnostics
    assert "post_contact_local_defect_repair_short_edge_count_after" in diagnostics
    assert diagnostics["post_contact_local_defect_repair_tolerance"] == 0.5


def test_regularize_coverage_contacts_keeps_low_drift_rescue_variant(
    monkeypatch,
):
    original = [box(0.0, 0.0, 1.0, 1.0), box(2.0, 0.0, 3.0, 1.0)]
    direct_candidate = [box(0.0, 0.0, 1.0, 1.0), box(1.75, 0.0, 3.0, 1.0)]
    rescue_candidate = [box(0.0, 0.0, 3.0, 1.0)]
    postprocessed_rescue_candidate = [box(0.0, 0.0, 2.5, 1.0)]
    original_sources = [[0], [1]]
    rescue_sources = [[0, 1]]

    original_key = cleaning_footprints._polygon_sequence_key(original)
    direct_key = cleaning_footprints._polygon_sequence_key(direct_candidate)
    rescue_key = cleaning_footprints._polygon_sequence_key(rescue_candidate)
    postprocessed_key = cleaning_footprints._polygon_sequence_key(
        postprocessed_rescue_candidate
    )

    signature_map = {
        original_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.1,
            pair_issue_count=2,
            point_touch_count=0,
            close_pair_count=2,
            min_pair_clearance=0.1,
            short_edge_count=0,
            min_edge_length=0.5,
            vertex_count=8,
            ring_contact_count=0,
        ),
        direct_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.15,
            pair_issue_count=1,
            point_touch_count=0,
            close_pair_count=1,
            min_pair_clearance=0.15,
            short_edge_count=0,
            min_edge_length=0.5,
            vertex_count=8,
            ring_contact_count=0,
        ),
        rescue_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.26,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=0.5,
            vertex_count=6,
            ring_contact_count=0,
        ),
        postprocessed_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.5,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=0.5,
            vertex_count=6,
            ring_contact_count=0,
        ),
    }
    difference_map = {
        (original_key, direct_key): {
            "reference_minus_candidate_area": 0.1,
            "candidate_minus_reference_area": 0.1,
            "symmetric_difference_area": 0.2,
            "union_area_delta": 0.0,
        },
        (original_key, rescue_key): {
            "reference_minus_candidate_area": 3.0,
            "candidate_minus_reference_area": 9.0,
            "symmetric_difference_area": 12.0,
            "union_area_delta": 6.0,
        },
        (original_key, postprocessed_key): {
            "reference_minus_candidate_area": 100.0,
            "candidate_minus_reference_area": 10.0,
            "symmetric_difference_area": 110.0,
            "union_area_delta": -90.0,
        },
    }

    monkeypatch.setattr(
        cleaning_footprints,
        "_direct_pair_issue_cluster_candidates",
        lambda *args, **kwargs: [
            ((0, 1), "coverage_pair_issue_bridge_0.500", direct_candidate, original_sources)
        ],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_repair_residual_pair_issues",
        lambda *args, **kwargs: (
            rescue_candidate,
            rescue_sources,
            {"coverage_pair_issue_bridge_residual_0.500": 1},
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_polygons_for_meshing",
        lambda polygons, source_map, **kwargs: (
            postprocessed_rescue_candidate
            if cleaning_footprints._polygon_sequence_key(polygons) == rescue_key
            else polygons,
            rescue_sources if cleaning_footprints._polygon_sequence_key(polygons) == rescue_key else source_map,
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_polygon_repairs",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_defect_signature",
        lambda cache, polygons, target_scale: signature_map[
            cleaning_footprints._polygon_sequence_key(polygons)
        ],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_difference_area_metrics",
        lambda cache, reference_polygons, candidate_polygons: difference_map.get(
            (
                cleaning_footprints._polygon_sequence_key(reference_polygons),
                cleaning_footprints._polygon_sequence_key(candidate_polygons),
            ),
            {
                "reference_minus_candidate_area": 0.0,
                "candidate_minus_reference_area": 0.0,
                "symmetric_difference_area": 0.0,
                "union_area_delta": 0.0,
            },
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_edit_zone",
        lambda *args, **kwargs: box(-1.0, -1.0, 4.0, 2.0),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_change_outside_edit_zone",
        lambda *args, **kwargs: 0.0,
    )

    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    polygons, source_map = cleaning_footprints._regularize_coverage_contacts(
        original,
        original_sources,
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert cleaning_footprints._polygon_sequence_key(polygons) == rescue_key
    assert source_map == rescue_sources
    assert diagnostics["coverage_contact_regularization_applied"] is True
    assert diagnostics["coverage_contact_regularization_pair_issue_count_after"] == 0


def test_regularize_coverage_contacts_prefers_small_shrink_candidate(
    monkeypatch,
):
    original = [box(0.0, 0.0, 1.0, 1.0), box(1.01, 0.0, 2.01, 1.0)]
    bridge_candidate = [box(0.0, 0.0, 2.01, 1.0)]
    shrink_candidate = [box(0.0, 0.0, 1.0, 1.0), box(1.26, 0.0, 1.76, 1.0)]
    original_sources = [[0], [1]]

    original_key = cleaning_footprints._polygon_sequence_key(original)
    bridge_key = cleaning_footprints._polygon_sequence_key(bridge_candidate)
    shrink_key = cleaning_footprints._polygon_sequence_key(shrink_candidate)

    signature_map = {
        original_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.01,
            pair_issue_count=1,
            point_touch_count=0,
            close_pair_count=1,
            min_pair_clearance=0.01,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=8,
            ring_contact_count=0,
        ),
        bridge_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.49,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=0.64,
            vertex_count=10,
            ring_contact_count=0,
        ),
        shrink_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=1.0,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=8,
            ring_contact_count=0,
        ),
    }
    difference_map = {
        (original_key, bridge_key): {
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 1.0,
            "symmetric_difference_area": 1.0,
            "union_area_delta": 1.0,
        },
        (original_key, shrink_key): {
            "reference_minus_candidate_area": 6.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 6.0,
            "union_area_delta": -6.0,
        },
    }

    monkeypatch.setattr(
        cleaning_footprints,
        "_direct_pair_issue_cluster_candidates",
        lambda *args, **kwargs: [
            ((0, 1), "coverage_pair_issue_bridge_0.500", bridge_candidate, [[0, 1]]),
            ((0, 1), "coverage_pair_issue_shrink_0.250", shrink_candidate, original_sources),
        ],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_polygons_for_meshing",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_polygon_repairs",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_defect_signature",
        lambda cache, polygons, target_scale: signature_map[
            cleaning_footprints._polygon_sequence_key(polygons)
        ],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_difference_area_metrics",
        lambda cache, reference_polygons, candidate_polygons: difference_map.get(
            (
                cleaning_footprints._polygon_sequence_key(reference_polygons),
                cleaning_footprints._polygon_sequence_key(candidate_polygons),
            ),
            {
                "reference_minus_candidate_area": 0.0,
                "candidate_minus_reference_area": 0.0,
                "symmetric_difference_area": 0.0,
                "union_area_delta": 0.0,
            },
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_edit_zone",
        lambda *args, **kwargs: box(-1.0, -1.0, 3.0, 2.0),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_change_outside_edit_zone",
        lambda *args, **kwargs: 0.0,
    )

    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    polygons, source_map = cleaning_footprints._regularize_coverage_contacts(
        original,
        original_sources,
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert cleaning_footprints._polygon_sequence_key(polygons) == shrink_key
    assert source_map == original_sources
    assert diagnostics["coverage_contact_regularization_operator_applied"] == {
        "coverage_pair_issue_shrink_0.250": 1
    }


def test_regularize_coverage_contacts_rejects_large_shrink_candidate(
    monkeypatch,
):
    original = [box(0.0, 0.0, 1.0, 1.0), box(1.01, 0.0, 2.01, 1.0)]
    bridge_candidate = [box(0.0, 0.0, 2.01, 1.0)]
    shrink_candidate = [box(0.0, 0.0, 1.0, 1.0), box(1.51, 0.0, 1.76, 1.0)]
    original_sources = [[0], [1]]

    original_key = cleaning_footprints._polygon_sequence_key(original)
    bridge_key = cleaning_footprints._polygon_sequence_key(bridge_candidate)
    shrink_key = cleaning_footprints._polygon_sequence_key(shrink_candidate)

    signature_map = {
        original_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.01,
            pair_issue_count=1,
            point_touch_count=0,
            close_pair_count=1,
            min_pair_clearance=0.01,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=8,
            ring_contact_count=0,
        ),
        bridge_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.49,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=0.64,
            vertex_count=10,
            ring_contact_count=0,
        ),
        shrink_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=1.0,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=8,
            ring_contact_count=0,
        ),
    }
    difference_map = {
        (original_key, bridge_key): {
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 1.0,
            "symmetric_difference_area": 1.0,
            "union_area_delta": 1.0,
        },
        (original_key, shrink_key): {
            "reference_minus_candidate_area": 12.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 12.0,
            "union_area_delta": -12.0,
        },
    }

    monkeypatch.setattr(
        cleaning_footprints,
        "_direct_pair_issue_cluster_candidates",
        lambda *args, **kwargs: [
            ((0, 1), "coverage_pair_issue_bridge_0.500", bridge_candidate, [[0, 1]]),
            ((0, 1), "coverage_pair_issue_shrink_0.250", shrink_candidate, original_sources),
        ],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_polygons_for_meshing",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_polygon_repairs",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_defect_signature",
        lambda cache, polygons, target_scale: signature_map[
            cleaning_footprints._polygon_sequence_key(polygons)
        ],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_difference_area_metrics",
        lambda cache, reference_polygons, candidate_polygons: difference_map.get(
            (
                cleaning_footprints._polygon_sequence_key(reference_polygons),
                cleaning_footprints._polygon_sequence_key(candidate_polygons),
            ),
            {
                "reference_minus_candidate_area": 0.0,
                "candidate_minus_reference_area": 0.0,
                "symmetric_difference_area": 0.0,
                "union_area_delta": 0.0,
            },
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_edit_zone",
        lambda *args, **kwargs: box(-1.0, -1.0, 3.0, 2.0),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_change_outside_edit_zone",
        lambda *args, **kwargs: 0.0,
    )

    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    polygons, source_map = cleaning_footprints._regularize_coverage_contacts(
        original,
        original_sources,
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert cleaning_footprints._polygon_sequence_key(polygons) == bridge_key
    assert source_map == [[0, 1]]
    assert diagnostics["coverage_contact_regularization_operator_applied"] == {
        "coverage_pair_issue_bridge_0.500": 1
    }


def test_regularize_coverage_contacts_uses_shrink_edit_zone_for_direct_shrink_candidates(
    monkeypatch,
):
    original = [box(0.0, 0.0, 1.0, 1.0), box(1.01, 0.0, 2.01, 1.0)]
    shrink_candidate = [box(0.0, 0.0, 1.0, 1.0), box(1.26, 0.0, 1.76, 1.0)]
    original_sources = [[0], [1]]

    original_key = cleaning_footprints._polygon_sequence_key(original)
    shrink_key = cleaning_footprints._polygon_sequence_key(shrink_candidate)

    signature_map = {
        original_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.01,
            pair_issue_count=1,
            point_touch_count=0,
            close_pair_count=1,
            min_pair_clearance=0.01,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=8,
            ring_contact_count=0,
        ),
        shrink_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.26,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=8,
            ring_contact_count=0,
        ),
    }
    difference_map = {
        (original_key, shrink_key): {
            "reference_minus_candidate_area": 5.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 5.0,
            "union_area_delta": -5.0,
        },
    }

    monkeypatch.setattr(
        cleaning_footprints,
        "_direct_pair_issue_cluster_candidates",
        lambda *args, **kwargs: [
            ((0, 1), "coverage_pair_issue_shrink_0.250", shrink_candidate, original_sources),
        ],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_polygons_for_meshing",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_polygon_repairs",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_defect_signature",
        lambda cache, polygons, target_scale: signature_map[
            cleaning_footprints._polygon_sequence_key(polygons)
        ],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_difference_area_metrics",
        lambda cache, reference_polygons, candidate_polygons: difference_map.get(
            (
                cleaning_footprints._polygon_sequence_key(reference_polygons),
                cleaning_footprints._polygon_sequence_key(candidate_polygons),
            ),
            {
                "reference_minus_candidate_area": 0.0,
                "candidate_minus_reference_area": 0.0,
                "symmetric_difference_area": 0.0,
                "union_area_delta": 0.0,
            },
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_edit_zone",
        lambda *args, **kwargs: box(-0.5, -0.5, 0.5, 0.5),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_pair_issue_shrink_budget_context",
        lambda *args, **kwargs: (box(-2.0, -2.0, 2.0, 2.0), 20.0),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_change_outside_edit_zone",
        lambda *args, **kwargs: 5.0,
    )

    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    polygons, source_map = cleaning_footprints._regularize_coverage_contacts(
        original,
        original_sources,
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert cleaning_footprints._polygon_sequence_key(polygons) == shrink_key
    assert source_map == original_sources
    assert diagnostics["coverage_contact_regularization_operator_applied"] == {
        "coverage_pair_issue_shrink_0.250": 1
    }


def test_post_contact_local_repair_accepts_subgrid_clearance_gain():
    reference_signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.014600733992761013,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=0,
        min_edge_length=0.5087131436281158,
        vertex_count=808,
        ring_contact_count=0,
    )
    candidate_signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.03125,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=0,
        min_edge_length=0.5087131436281158,
        vertex_count=812,
        ring_contact_count=0,
    )
    reference_difference_metrics = {
        "reference_minus_candidate_area": 0.0,
        "candidate_minus_reference_area": 0.0,
        "symmetric_difference_area": 0.0,
        "union_area_delta": 0.0,
    }
    candidate_difference_metrics = {
        "reference_minus_candidate_area": 0.16746342440630335,
        "candidate_minus_reference_area": 1.0185376431563034,
        "symmetric_difference_area": 1.1860010675626067,
        "union_area_delta": 0.85107421875,
    }

    assert (
        cleaning_footprints._coverage_signature_improves(
            reference_signature,
            candidate_signature,
            grid=0.03125,
            target_scale=0.5,
        )
        is False
    )
    assert (
        cleaning_footprints._should_accept_post_contact_local_repair(
            reference_signature,
            reference_difference_metrics,
            candidate_signature,
            candidate_difference_metrics,
            target_scale=0.5,
            grid=0.03125,
        )
        is True
    )


def test_regularize_coverage_contacts_accepts_local_point_bridge_on_case23_pair():
    first = Polygon(
        [
            (674248.707, 6579722.925),
            (674241.262, 6579748.708),
            (674235.883, 6579747.221),
            (674240.42, 6579729.477),
            (674233.568, 6579727.805),
            (674235.599, 6579719.836),
            (674248.707, 6579722.925),
        ]
    )
    second = Polygon(
        [
            (674241.261, 6579748.707),
            (674247.882, 6579750.29),
            (674245.251, 6579759.837),
            (674238.928, 6579758.325),
            (674241.261, 6579748.707),
        ]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_coverage_contacts(
        [first, second],
        [[20], [21]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert len(polygons) == 1
    assert source_map == [[20, 21]]
    assert diagnostics["coverage_contact_regularization_applied"] is True
    assert diagnostics["coverage_contact_regularization_pair_issue_count_before"] == 1
    assert diagnostics["coverage_contact_regularization_pair_issue_count_after"] == 0
    assert any(
        key.startswith("coverage_pair_issue_point_")
        for key in diagnostics["coverage_contact_regularization_operator_applied"]
    )


def test_regularize_coverage_contacts_uses_point_cluster_bridge_for_case62_pair():
    first = Polygon(
        [
            (673851.086, 6581961.473),
            (673844.989, 6581972.158),
            (673834.367, 6581966.078),
            (673867.295, 6581908.338),
            (673877.779, 6581914.315),
            (673865.408, 6581936.079),
            (673876.751, 6581942.56),
            (673886.957, 6581924.408),
            (673897.286, 6581930.534),
            (673867.97, 6581981.994),
            (673857.672, 6581975.823),
            (673862.203, 6581967.901),
            (673851.086, 6581961.473),
        ]
    )
    second = Polygon(
        [
            (673840.879, 6581979.362),
            (673844.989, 6581972.158),
            (673849.579, 6581974.784),
            (673845.277, 6581981.866),
            (673840.879, 6581979.362),
        ]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_coverage_contacts(
        [first, second],
        [[46], [47]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert len(polygons) == 1
    assert source_map == [[46, 47]]
    assert diagnostics["coverage_contact_regularization_pair_issue_count_before"] == 1
    assert diagnostics["coverage_contact_regularization_pair_issue_count_after"] == 0
    assert diagnostics["coverage_contact_regularization_short_edge_count_before"] == 0
    assert diagnostics["coverage_contact_regularization_short_edge_count_after"] == 0
    assert any(
        key.startswith("coverage_pair_issue_point_cluster_")
        for key in diagnostics["coverage_contact_regularization_operator_applied"]
    )

    signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
    assert signature.pair_issue_count == 0
    assert signature.short_edge_count == 0


def test_condition_polygon_coverage_exposes_contact_regularization_diagnostics(
    monkeypatch,
):
    def fake_regularize_coverage_contacts(*args, **kwargs):
        diagnostics = kwargs["diagnostics"]
        diagnostics["coverage_contact_regularization_applied"] = True
        diagnostics["coverage_contact_regularization_pair_issue_count_before"] = 2
        diagnostics["coverage_contact_regularization_pair_issue_count_after"] = 0
        diagnostics["coverage_contact_regularization_operator_applied"] = {
            "coverage_pair_issue_point_0.125": 1
        }
        return list(args[0]), [list(indices) for indices in args[1]]

    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_coverage_contacts",
        fake_regularize_coverage_contacts,
    )

    result = cleaning.condition_polygon_coverage(
        [box(0, 0, 2, 2)],
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.0,
            min_area=0.0,
            min_hole_area=0.0,
            precision_grid=0.03125,
            collect_stage_metrics=False,
            enable_logging=False,
        ),
    )

    assert result.diagnostics["coverage_contact_regularization_applied"] is True
    assert result.diagnostics["coverage_contact_regularization_pair_issue_count_before"] == 2
    assert result.diagnostics["coverage_contact_regularization_pair_issue_count_after"] == 0
    assert result.diagnostics["coverage_contact_regularization_operator_applied"] == {
        "coverage_pair_issue_point_0.125": 1
    }


def test_condition_polygon_coverage_is_deterministic():
    geometries = [box(10, 0, 14, 4), box(0, 0, 4, 4), box(4.1, 0, 8.1, 4)]
    options = cleaning.ConditioningOptions(
        min_feature_size=0.0,
        merge_distance=0.2,
        min_area=0.0,
        min_hole_area=0.0,
    )
    first = cleaning.condition_polygon_coverage(
        geometries,
        source_map=[[9], [2], [5, 2]],
        options=options,
    )
    second = cleaning.condition_polygon_coverage(
        geometries,
        source_map=[[9], [2], [5, 2]],
        options=options,
    )
    assert [polygon.wkt for polygon in first.polygons] == [
        polygon.wkt for polygon in second.polygons
    ]
    assert first.source_map == second.source_map
    assert_conditioning_invariants(first, options)


def test_condition_polygon_coverage_recovers_exact_source_coordinates_for_good_polygon():
    raw = Polygon(
        [
            (1.379, 6.173),
            (162.89, 54.754),
            (160.106, 63.836),
            (0.161, 15.261),
        ]
    )
    options = cleaning.ConditioningOptions(
        min_feature_size=0.5,
        merge_distance=0.5,
        min_area=0.0,
        min_hole_area=0.0,
    )

    result = cleaning.condition_polygon_coverage([raw], options=options)

    assert len(result.polygons) == 1
    assert exterior_vertex_set(result.polygons[0]) == exterior_vertex_set(raw)
    assert result.diagnostics["source_coordinate_recovery_applied"] is True
    assert result.diagnostics["source_coordinate_recovery_exact_count"] == 1


def test_recover_polygon_source_coordinates_skips_non_admissible_polygon():
    polygon = Polygon(
        [
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (3.75, 4.0),
            (3.75, 3.0),
            (0.0, 3.0),
        ]
    )
    support = Polygon(
        [
            (0.02, 0.0),
            (4.02, 0.0),
            (4.02, 4.0),
            (3.77, 4.0),
            (3.77, 3.0),
            (0.02, 3.0),
        ]
    )

    candidate, operator, metrics = (
        cleaning_footprints._recover_polygon_source_coordinates(
            polygon,
            support=support,
            min_segment_length=0.5,
            grid=0.125,
        )
    )

    assert candidate is None
    assert operator is None
    assert metrics is None


def test_recover_polygon_source_coordinates_normalizes_exact_support_polygon():
    polygon = Polygon(
        [
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
        ]
    )
    support = Polygon(
        [
            (0.0625, 0.0),
            (4.03125, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
            (0.0625, 0.0),
        ]
    )

    candidate, operator, metrics = (
        cleaning_footprints._recover_polygon_source_coordinates(
            polygon,
            support=support,
            min_segment_length=0.5,
            grid=0.03125,
        )
    )

    assert candidate is not None
    assert operator == "exact_source_polygon"
    assert metrics is not None
    lengths = [
        ((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5
        for a, b in zip(candidate.exterior.coords, list(candidate.exterior.coords)[1:])
    ]
    assert min(lengths) >= 0.125 - 1e-12


def test_source_coordinate_recovery_rejects_local_overlap_candidate(monkeypatch):
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    original = box(0.0, 0.0, 1.0, 1.0)
    neighbor = box(1.1, 0.0, 2.1, 1.0)
    overlap_candidate = box(0.8, 0.0, 1.8, 1.0)

    def fake_recover(polygon, *, support, min_segment_length, grid, **kwargs):
        if polygon.equals_exact(original, tolerance=0.0):
            return (
                overlap_candidate,
                "support_vertex_restore",
                {
                    "reference_minus_candidate_area": 0.0,
                    "candidate_minus_reference_area": 0.0,
                    "symmetric_difference_area": 0.0,
                    "union_area_delta": 0.0,
                },
            )
        return None, None, None

    monkeypatch.setattr(
        cleaning_footprints,
        "_recover_polygon_source_coordinates",
        fake_recover,
    )

    polygons, source_map = cleaning_footprints._recover_source_supported_coordinates(
        [original, neighbor],
        [[0], [1]],
        source_lookup={0: original, 1: neighbor},
        min_segment_length=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert polygons[0].equals_exact(original, tolerance=0.0)
    assert polygons[1].equals_exact(neighbor, tolerance=0.0)
    assert source_map == [[0], [1]]
    assert diagnostics["source_coordinate_recovery_applied_count"] == 0
    assert diagnostics["source_coordinate_recovery_rejected_overlap_count"] == 1


def test_source_coordinate_recovery_rejects_local_close_pair_candidate(monkeypatch):
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    original = box(0.0, 0.0, 1.0, 1.0)
    neighbor = box(2.0, 0.0, 3.0, 1.0)
    close_pair_candidate = box(1.4, 0.0, 2.4, 1.0)

    def fake_recover(polygon, *, support, min_segment_length, grid, **kwargs):
        if polygon.equals_exact(neighbor, tolerance=0.0):
            return (
                close_pair_candidate,
                "support_vertex_restore",
                {
                    "reference_minus_candidate_area": 0.0,
                    "candidate_minus_reference_area": 0.0,
                    "symmetric_difference_area": 0.0,
                    "union_area_delta": 0.0,
                },
            )
        return None, None, None

    monkeypatch.setattr(
        cleaning_footprints,
        "_recover_polygon_source_coordinates",
        fake_recover,
    )

    polygons, source_map = cleaning_footprints._recover_source_supported_coordinates(
        [original, neighbor],
        [[0], [1]],
        source_lookup={0: original, 1: neighbor},
        min_segment_length=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert polygons[0].equals_exact(original, tolerance=0.0)
    assert polygons[1].equals_exact(neighbor, tolerance=0.0)
    assert source_map == [[0], [1]]
    assert diagnostics["source_coordinate_recovery_applied_count"] == 0
    assert diagnostics["source_coordinate_recovery_rejected_non_improving_count"] >= 1


def test_source_coordinate_recovery_rejects_point_touch_to_close_pair_regression(
    monkeypatch,
):
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    original = box(0.0, 0.0, 1.0, 1.0)
    point_touch_neighbor = box(1.0, 1.0, 2.0, 2.0)
    close_pair_candidate = box(1.02, 1.02, 2.02, 2.02)

    def fake_recover(polygon, *, support, min_segment_length, grid, **kwargs):
        if polygon.equals_exact(point_touch_neighbor, tolerance=0.0):
            return (
                close_pair_candidate,
                "exact_source_polygon",
                {
                    "reference_minus_candidate_area": 0.0,
                    "candidate_minus_reference_area": 0.0,
                    "symmetric_difference_area": 0.0,
                    "union_area_delta": 0.0,
                },
            )
        return None, None, None

    monkeypatch.setattr(
        cleaning_footprints,
        "_recover_polygon_source_coordinates",
        fake_recover,
    )

    polygons, source_map = cleaning_footprints._recover_source_supported_coordinates(
        [original, point_touch_neighbor],
        [[0], [1]],
        source_lookup={0: original, 1: point_touch_neighbor},
        min_segment_length=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert polygons[0].equals_exact(original, tolerance=0.0)
    assert polygons[1].equals_exact(point_touch_neighbor, tolerance=0.0)
    assert source_map == [[0], [1]]
    assert diagnostics["source_coordinate_recovery_applied_count"] == 0
    assert diagnostics["source_coordinate_recovery_rejected_non_improving_count"] >= 1


def test_source_coordinate_recovery_reverts_full_contract_regression(monkeypatch):
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    original = box(0.0, 0.0, 4.0, 4.0)
    candidate = Polygon(
        [
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (2.0, 4.0),
            (2.0, 3.6),
            (1.6, 3.6),
            (1.6, 4.0),
            (0.0, 4.0),
        ]
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_recover_polygon_source_coordinates",
        lambda polygon, **kwargs: (
            (
                candidate,
                "support_vertex_restore",
                {
                    "reference_minus_candidate_area": 0.0,
                    "candidate_minus_reference_area": 0.0,
                    "symmetric_difference_area": 0.0,
                    "union_area_delta": 0.0,
                },
            )
            if polygon.equals_exact(original, tolerance=0.0)
            else (None, None, None)
        ),
    )

    original_key = cleaning_footprints._polygon_sequence_key([original])
    candidate_key = cleaning_footprints._polygon_sequence_key([candidate])
    signatures = {
        original_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.75,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=4.0,
            vertex_count=4,
            ring_contact_count=0,
        ),
        candidate_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.4,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=1,
            min_edge_length=0.4,
            vertex_count=8,
            ring_contact_count=0,
        ),
    }
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_defect_signature",
        lambda cache, polygons, **kwargs: signatures[
            cleaning_footprints._polygon_sequence_key(polygons)
        ],
    )

    polygons, source_map = cleaning_footprints._recover_source_supported_coordinates(
        [original],
        [[0]],
        source_lookup={0: original},
        min_segment_length=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert polygons[0].equals_exact(original, tolerance=0.0)
    assert source_map == [[0]]
    assert diagnostics["source_coordinate_recovery_applied"] is False
    assert diagnostics["source_coordinate_recovery_applied_count"] == 0
    assert diagnostics["source_coordinate_recovery_contract_reject_count"] == 1
    assert diagnostics["source_coordinate_recovery_short_edge_count_after"] == 0
    assert diagnostics["source_coordinate_recovery_pair_issue_count_after"] == 0
    assert diagnostics["source_coordinate_recovery_ring_contact_count_after"] == 0
    assert diagnostics["source_coordinate_recovery_min_clearance_after"] == pytest.approx(
        0.75
    )


def test_recover_polygon_source_coordinates_skips_vertex_restore_for_single_source_close_match(
    monkeypatch,
):
    polygon = Polygon(
        [
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (3.75, 4.0),
            (3.75, 3.0),
            (0.0, 3.0),
        ]
    )
    support = Polygon(
        [
            (0.02, 0.0),
            (4.02, 0.0),
            (4.0, 4.0),
            (3.77, 4.0),
            (3.77, 3.0),
            (0.02, 3.0),
        ]
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_recover_polygon_vertices_from_support",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("single-source close-match should not try vertex restore")
        ),
    )

    candidate, operator, metrics = cleaning_footprints._recover_polygon_source_coordinates(
        polygon,
        support=support,
        support_source_count=1,
        min_segment_length=0.5,
        grid=0.125,
    )

    assert candidate is None
    assert operator is None
    assert metrics is None


def test_condition_polygon_coverage_removes_meshing_hostile_short_edges():
    staircase = Polygon(
        [
            (0, 0),
            (10, 0),
            (10, 10),
            (9.875, 10),
            (9.875, 9.5),
            (9.75, 9.5),
            (9.75, 9.0),
            (9.625, 9.0),
            (9.625, 8.5),
            (9.5, 8.5),
            (9.5, 8.0),
            (9.375, 8.0),
            (9.375, 7.5),
            (9.25, 7.5),
            (9.25, 7.0),
            (9.125, 7.0),
            (9.125, 6.5),
            (9.0, 6.5),
            (9.0, 6.0),
            (0, 6.0),
        ]
    )
    options = cleaning.ConditioningOptions(
        min_feature_size=0.5,
        merge_distance=0.0,
        min_area=0.0,
        min_hole_area=0.0,
    )
    result = cleaning.condition_polygon_coverage([staircase], options=options)
    assert_conditioning_invariants(result, options)
    assert (
        result.diagnostics["stage_metrics"]["final_output"]["short_edge_count"] == 0
    )

    min_edge = min(
        ((bx - ax) ** 2 + (by - ay) ** 2) ** 0.5
        for polygon in result.polygons
        for ring in [polygon.exterior, *polygon.interiors]
        for (ax, ay), (bx, by) in zip(ring.coords, ring.coords[1:])
    )
    assert min_edge >= options.min_feature_size


def test_residual_scale_polygon_simplifier_enforces_contract():
    polygon = Polygon(
        [
            (0.0, 0.0),
            (20.0, 0.0),
            (20.0, 10.0),
            (10.2, 10.0),
            (10.2, 9.56),
            (9.7, 9.56),
            (9.7, 10.0),
            (0.0, 10.0),
        ]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    candidate = cleaning_footprints._simplify_coverage_residual_scale_polygons(
        [polygon],
        [[0]],
        tolerance=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
        cache=cleaning_footprints._CoverageEvalCache(),
    )

    assert candidate is not None
    assert candidate.operator_applied == {"coverage_residual_polygon_simplify": 1}
    assert cleaning_footprints._coverage_signature_satisfies_scale_contract(
        candidate.signature,
        target_scale=0.5,
        grid=0.03125,
    )
    assert diagnostics["coverage_residual_polygon_target_count"] == 1
    assert diagnostics["coverage_residual_polygon_attempt_count"] >= 1
    assert diagnostics["coverage_residual_polygon_viable_candidate_count"] >= 1
    assert diagnostics["coverage_residual_polygon_changed_count"] == 1


def test_condition_polygon_coverage_logs_stage_progress(monkeypatch):
    messages: list[str] = []
    monkeypatch.setattr(cleaning_footprints, "info", lambda message: messages.append(message))

    cleaning.condition_polygon_coverage(
        [box(0, 0, 4, 4), box(4.1, 0, 8.1, 4)],
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.2,
            min_area=0.0,
            min_hole_area=0.0,
        ),
    )

    joined = "\n".join(messages)
    assert "Footprint cleaning started" in joined
    assert "atomic_input" in joined
    assert "opened" in joined
    assert "regularized_groups" in joined
    assert "reconstructed" in joined
    assert "presimplify" in joined
    assert "local_defect_repaired" in joined
    assert "coverage_simplified" in joined
    assert "source_reclaimed" in joined
    assert "small_component_absorbed" in joined
    assert "boundary_regularized" in joined
    assert "clearance_regularized" in joined
    assert "source_coordinate_recovered" in joined
    assert "coverage_meshing_regularized" in joined
    assert "Footprint cleaning summary" in joined
    assert "local_defect_repair candidates=" in joined
    assert "coverage_simplify branch=" in joined
    assert "coverage_simplify operators attempted=" in joined
    assert "source_reclaim candidates=" in joined
    assert "small_component_absorb candidates=" in joined
    assert "source_coordinate_recovery candidates=" in joined
    assert "coverage_meshing_regularization branch=" in joined
    assert "boundary_regularization candidates=" in joined
    assert "final_min_area_filter removed=" in joined


def test_condition_polygon_coverage_skips_global_reconstruction_for_disjoint_coverage():
    result = cleaning.condition_polygon_coverage(
        [
            box(0, 0, 4, 4),
            box(6, 0, 10, 4),
            box(20, 0, 24, 4),
            box(26, 0, 30, 4),
        ],
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.5,
            min_area=0.0,
            min_hole_area=0.0,
        ),
    )

    assert result.diagnostics["global_reconstruction_applied"] is False
    assert (
        result.diagnostics["global_reconstruction_reason"]
        == "coverage_already_disjoint"
    )
    assert (
        result.diagnostics["stage_metrics"]["regularized_groups"]["union_area"]
        == pytest.approx(
            result.diagnostics["stage_metrics"]["reconstructed"]["union_area"],
            abs=1e-9,
        )
    )


def test_condition_polygon_coverage_can_disable_stage_metrics(monkeypatch):
    messages: list[str] = []
    monkeypatch.setattr(cleaning_footprints, "info", lambda message: messages.append(message))

    result = cleaning.condition_polygon_coverage(
        [box(0, 0, 4, 4), box(4.1, 0, 8.1, 4)],
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.2,
            min_area=0.0,
            min_hole_area=0.0,
            collect_stage_metrics=False,
        ),
    )

    assert result.diagnostics["collect_stage_metrics"] is False
    assert result.diagnostics["stage_metrics"] == {}

    joined = "\n".join(messages)
    assert "Footprint cleaning started" in joined
    assert "Footprint cleaning summary" in joined
    assert "stage_metrics=disabled" in joined
    assert "atomic_input:" not in joined


def test_condition_polygon_coverage_can_disable_logging(monkeypatch):
    messages: list[str] = []
    monkeypatch.setattr(cleaning_footprints, "info", lambda message: messages.append(message))

    result = cleaning.condition_polygon_coverage(
        [box(0, 0, 4, 4), box(4.1, 0, 8.1, 4)],
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.2,
            min_area=0.0,
            min_hole_area=0.0,
            enable_logging=False,
        ),
    )

    assert result.diagnostics["enable_logging"] is False
    assert messages == []


def test_condition_polygon_coverage_skips_zero_scale_source_recovery():
    result = cleaning.condition_polygon_coverage(
        [box(0, 0, 2, 2), box(2, 2, 4, 4)],
        options=cleaning.ConditioningOptions(
            min_feature_size=0.0,
            merge_distance=0.0,
            min_area=0.0,
            min_hole_area=0.0,
        ),
    )

    assert result.diagnostics["source_coordinate_recovery_applied"] is False
    assert result.diagnostics["source_coordinate_recovery_applied_count"] == 0


def test_condition_polygon_coverage_applies_min_area_at_final_output():
    result = cleaning.condition_polygon_coverage(
        [MultiPolygon([box(0, 0, 4, 4), box(10, 0, 11, 1)])],
        options=cleaning.ConditioningOptions(
            min_feature_size=0.0,
            merge_distance=0.0,
            min_area=2.0,
            min_hole_area=0.0,
        ),
    )

    assert len(result.polygons) == 1
    assert result.polygons[0].equals(box(0, 0, 4, 4))
    assert result.diagnostics["final_min_area_filter_removed_count"] == 1
    assert result.diagnostics["final_min_area_filter_removed_area"] == pytest.approx(
        1.0,
        abs=1e-12,
    )


def test_small_component_absorb_merges_supported_remnant_into_large_neighbor():
    large = box(0, 0, 4, 4)
    remnant = box(4.8, 0, 5.6, 1.0)
    support = large.union(remnant).union(box(4.0, 0.0, 4.8, 1.0))
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._absorb_small_supported_components(
        [large, remnant],
        [[1], [1]],
        source_lookup={1: support},
        raw_support_union=support,
        min_area=2.0,
        min_segment_length=0.75,
        grid=0.1,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert len(polygons) == 1
    assert source_map == [[1]]
    assert polygons[0].area > large.area
    assert diagnostics["small_component_absorb_applied"] is True
    assert diagnostics["small_component_absorb_applied_count"] == 1


def test_small_component_absorb_allows_scale_limited_raw_supported_bridge():
    large = box(0, 0, 4, 4)
    remnant = box(7.2, 0, 8.2, 1.0)
    support = large.union(remnant).union(box(4.0, 0.0, 7.2, 1.0))
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._absorb_small_supported_components(
        [large, remnant],
        [[1], [2]],
        source_lookup={1: large, 2: remnant},
        raw_support_union=support,
        min_area=16.0,
        min_segment_length=0.5,
        grid=0.1,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert len(polygons) == 1
    assert source_map == [[1, 2]]
    assert diagnostics["small_component_absorb_applied"] is True
    assert diagnostics["small_component_absorb_applied_count"] == 1


def test_assign_component_sources_splits_disconnected_components_by_overlap():
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    components = [box(0, 0, 1, 1), box(3, 0, 4, 1)]

    assigned = cleaning_footprints._assign_component_sources(
        components,
        components,
        [[10], [20]],
        grid=0.01,
        diagnostics=diagnostics,
    )

    assert assigned == [[10], [20]]
    assert diagnostics["multi_component_group_count"] == 1
    assert diagnostics["component_source_reassignment_count"] == 2


def test_low_clearance_regularization_improves_polygon_clearance():
    polygon = Polygon(
        [
            (0, 0),
            (0, 8),
            (1, 8),
            (1, 5),
            (4, 5),
            (4, 4.8),
            (1, 4.8),
            (1, 0),
        ]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_low_clearance_polygons(
        [polygon],
        [[0]],
        min_clearance=0.5,
        grid=0.125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert source_map == [[0]]
    assert len(polygons) == 1
    assert polygons[0].minimum_clearance > polygon.minimum_clearance
    assert diagnostics["clearance_regularization_candidate_count"] == 1
    assert diagnostics["clearance_regularization_improved_count"] == 1
    assert diagnostics["clearance_regularization_failed_count"] == 0


def test_polygon_simplify_only_touches_polygons_with_short_edges():
    staircase = Polygon(
        [
            (0, 0),
            (10, 0),
            (10, 10),
            (9.875, 10),
            (9.875, 9.5),
            (9.75, 9.5),
            (9.75, 9.0),
            (9.625, 9.0),
            (9.625, 8.5),
            (9.5, 8.5),
            (9.5, 8.0),
            (9.375, 8.0),
            (9.375, 7.5),
            (9.25, 7.5),
            (9.25, 7.0),
            (9.125, 7.0),
            (9.125, 6.5),
            (9.0, 6.5),
            (9.0, 6.0),
            (0, 6.0),
        ]
    )
    clean = box(20, 0, 28, 8)
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._simplify_polygons_for_meshing(
        [staircase, clean],
        [[1], [2]],
        min_segment_length=0.5,
        grid=0.125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    keyed = {
        tuple(indices): polygon
        for polygon, indices in zip(polygons, source_map)
    }
    assert keyed[(2,)].equals_exact(clean, tolerance=0.0)
    assert diagnostics["polygon_simplify_candidate_count"] == 1
    assert diagnostics["polygon_simplify_applied_count"] == 1
    assert diagnostics["polygon_simplify_short_edge_count_after"] < diagnostics[
        "polygon_simplify_short_edge_count_before"
    ]


def test_coverage_simplify_applies_local_shared_boundary_patches():
    shared_boundary = [
        (5.0, 0.0),
        (5.0, 1.0),
        (4.875, 1.0),
        (4.875, 1.5),
        (5.0, 1.5),
        (5.0, 2.0),
        (4.875, 2.0),
        (4.875, 2.5),
        (5.0, 2.5),
        (5.0, 6.0),
    ]
    left = Polygon([(0.0, 0.0), *shared_boundary, (0.0, 6.0)])
    right = Polygon([shared_boundary[0], (10.0, 0.0), (10.0, 6.0), shared_boundary[-1], *reversed(shared_boundary[1:-1])])
    distant = box(20.0, 0.0, 24.0, 4.0)

    diagnostics = cleaning_footprints._empty_diagnostics(3)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._simplify_coverage(
        [left, right, distant],
        [[0], [1], [2]],
        tolerance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert len(polygons) == 3
    assert source_map == [[0], [1], [2]]
    assert diagnostics["coverage_simplify_patch_count"] >= 1
    assert diagnostics["coverage_simplify_patch_applied_count"] >= 1
    assert diagnostics["coverage_simplify_short_edge_count_after"] < diagnostics[
        "coverage_simplify_short_edge_count_before"
    ]
    assert any(bounds[0] > 19.9 for bounds in (polygon.bounds for polygon in polygons))


def test_coverage_simplify_local_operator_family_includes_patch_morphology():
    shared_boundary = [
        (5.0, 0.0),
        (5.0, 1.0),
        (4.875, 1.0),
        (4.875, 1.5),
        (5.0, 1.5),
        (5.0, 2.0),
        (4.875, 2.0),
        (4.875, 2.5),
        (5.0, 2.5),
        (5.0, 6.0),
    ]
    left = Polygon([(0.0, 0.0), *shared_boundary, (0.0, 6.0)])
    right = Polygon(
        [
            shared_boundary[0],
            (10.0, 0.0),
            (10.0, 6.0),
            shared_boundary[-1],
            *reversed(shared_boundary[1:-1]),
        ]
    )
    distant = box(20.0, 0.0, 24.0, 4.0)

    diagnostics = cleaning_footprints._empty_diagnostics(3)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._simplify_coverage_locally(
        [left, right, distant],
        [[0], [1], [2]],
        tolerance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.patch_count >= 1
    assert candidate.patch_applied_count >= 1
    assert any(
        operator.startswith("coverage_patch_")
        or operator == "coverage_graph_short_edges"
        for operator in candidate.operator_attempts
    )


def test_nonpair_patch_fallback_configs_skip_light_short_edge_clusters():
    signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.45,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=1,
        min_edge_length=0.4,
        vertex_count=12,
    )

    configs = cleaning_footprints._iter_nonpair_patch_fallback_configs(
        signature,
        subset_size=3,
        patch_radii=(0.125, 0.1875, 0.25),
    )

    assert configs == ()


def test_nonpair_patch_fallback_configs_keep_small_deterministic_rescue():
    signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.3,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=2,
        min_edge_length=0.25,
        vertex_count=10,
    )

    configs = cleaning_footprints._iter_nonpair_patch_fallback_configs(
        signature,
        subset_size=2,
        patch_radii=(0.125, 0.1875, 0.25),
    )

    assert configs == (("open_close", 0.125),)


def test_local_coverage_simplify_tolerances_skip_fine_trial_for_light_nonpair_cluster():
    signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.45,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=1,
        min_edge_length=0.4,
        vertex_count=12,
    )

    tolerances = cleaning_footprints._iter_local_coverage_simplify_tolerances(
        signature,
        patch_tolerances=(0.5, 0.625),
    )

    assert tolerances == (0.625,)


def test_local_coverage_simplify_tolerances_keep_pair_cluster_trials():
    signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.2,
        pair_issue_count=1,
        point_touch_count=0,
        close_pair_count=1,
        min_pair_clearance=0.2,
        short_edge_count=1,
        min_edge_length=0.2,
        vertex_count=14,
    )

    tolerances = cleaning_footprints._iter_local_coverage_simplify_tolerances(
        signature,
        patch_tolerances=(0.5, 0.625),
    )

    assert tolerances == (0.5, 0.625)


def test_graph_short_edge_simplifier_is_noop_for_clean_rectangle():
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate, removed = cleaning_footprints._simplify_polygon_short_edge_graphically(
        box(0.0, 0.0, 10.0, 5.0),
        target_scale=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is None
    assert removed == 0


def test_coverage_simplify_local_uses_graph_short_edge_fast_path():
    shared_boundary = [
        (5.0, 0.0),
        (5.0, 1.0),
        (4.875, 1.0),
        (4.875, 1.5),
        (5.0, 1.5),
        (5.0, 2.0),
        (4.875, 2.0),
        (4.875, 2.5),
        (5.0, 2.5),
        (5.0, 6.0),
    ]
    left = Polygon([(0.0, 0.0), *shared_boundary, (0.0, 6.0)])
    right = Polygon(
        [
            shared_boundary[0],
            (10.0, 0.0),
            (10.0, 6.0),
            shared_boundary[-1],
            *reversed(shared_boundary[1:-1]),
        ]
    )

    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._simplify_coverage_short_edge_graphically(
        [left, right],
        [[0], [1]],
        tolerance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.signature.short_edge_count < 6
    assert candidate.operator_applied == {"coverage_graph_short_edges": 2}


def test_coverage_defect_clusters_group_local_shared_boundary_defects():
    shared_boundary = [
        (5.0, 0.0),
        (5.0, 1.0),
        (4.875, 1.0),
        (4.875, 1.5),
        (5.0, 1.5),
        (5.0, 2.0),
        (4.875, 2.0),
        (4.875, 2.5),
        (5.0, 2.5),
        (5.0, 6.0),
    ]
    left = Polygon([(0.0, 0.0), *shared_boundary, (0.0, 6.0)])
    right = Polygon(
        [
            shared_boundary[0],
            (10.0, 0.0),
            (10.0, 6.0),
            shared_boundary[-1],
            *reversed(shared_boundary[1:-1]),
        ]
    )
    distant = box(20.0, 0.0, 24.0, 4.0)

    clusters = cleaning_footprints._coverage_defect_clusters(
        [left, right, distant],
        target_scale=0.5,
        cluster_radius=1.5,
    )

    assert clusters == [[0, 1]]


def test_clean_ring_coords_removes_tiny_closing_segment():
    cleaned = cleaning_footprints._clean_ring_coords(
        [
            (0.0625, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
            (0.0625, 0.0),
        ],
        grid=0.03125,
    )

    assert cleaned is not None
    assert cleaned[0] == cleaned[-1]
    lengths = [
        ((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5
        for a, b in zip(cleaned, cleaned[1:])
    ]
    assert min(lengths) >= 0.125 - 1e-12


def test_coverage_defect_cluster_descriptors_classify_short_edge_cluster():
    shared_boundary = [
        (5.0, 0.0),
        (5.0, 1.0),
        (4.875, 1.0),
        (4.875, 1.5),
        (5.0, 1.5),
        (5.0, 2.0),
        (4.875, 2.0),
        (4.875, 2.5),
        (5.0, 2.5),
        (5.0, 6.0),
    ]
    left = Polygon([(0.0, 0.0), *shared_boundary, (0.0, 6.0)])
    right = Polygon(
        [
            shared_boundary[0],
            (10.0, 0.0),
            (10.0, 6.0),
            shared_boundary[-1],
            *reversed(shared_boundary[1:-1]),
        ]
    )

    descriptors = cleaning_footprints._cached_coverage_defect_cluster_descriptors(
        None,
        [left, right],
        target_scale=0.5,
        cluster_radius=1.5,
    )

    assert len(descriptors) == 1
    assert descriptors[0]["indices"] == [0, 1]
    assert descriptors[0]["kind"] in {"short_edge_only", "mixed_pair_short_edge"}


def test_cpp_cluster_rewrite_removes_short_edge_chain():
    polygon = Polygon(
        [
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (2.0, 4.0),
            (2.0, 4.2),
            (1.8, 4.2),
            (1.8, 4.0),
            (0.0, 4.0),
        ]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    rewrite = cleaning_footprints._rewrite_coverage_cluster_graphically(
        [polygon],
        [[0]],
        cluster={"indices": [0], "kind": "short_edge_only"},
        tolerance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert rewrite is not None
    operator, candidate_polygons, candidate_sources = rewrite
    assert operator == "coverage_cpp_cluster_rewrite"
    assert candidate_sources == [[0]]
    signature = cleaning_footprints._coverage_defect_signature(
        candidate_polygons,
        target_scale=0.5,
    )
    assert signature.short_edge_count == 0


def test_coverage_simplify_local_operator_resolves_point_touch_pair():
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._simplify_coverage_locally(
        [box(0.0, 0.0, 1.0, 1.0), box(1.0, 1.0, 2.0, 2.0)],
        [[0], [1]],
        tolerance=0.01,
        grid=0.01,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.signature.pair_issue_count == 0
    assert candidate.operator_applied == {"coverage_pair_issue_point_0.010": 1}


def test_coverage_simplify_local_skips_generic_patch_fallback_for_simple_single_polygon(
    monkeypatch,
):
    polygon = Polygon(
        [
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (2.0, 4.0),
            (2.0, 4.2),
            (1.8, 4.2),
            (1.8, 4.0),
            (0.0, 4.0),
        ]
    )
    signature = cleaning_footprints._coverage_defect_signature(
        [polygon],
        target_scale=0.5,
    )
    assert signature.pair_issue_count == 0
    assert signature.short_edge_count <= 4

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_pair_issues_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints.shapely,
        "coverage_simplify",
        lambda geometries, tolerance, simplify_boundary=True: geometries,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_patch_union_operator",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("generic patch fallback should have been skipped")
        ),
    )

    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._simplify_coverage_locally(
        [polygon],
        [[0]],
        tolerance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    if candidate is not None:
        assert candidate.operator_applied == {"coverage_cpp_cluster_rewrite": 1}


def test_coverage_simplify_local_skips_pair_patch_fallback_for_simple_pair_cluster(
    monkeypatch,
):
    polygons = [box(0.0, 0.0, 1.0, 1.0), box(1.0, 1.0, 2.0, 2.0)]
    signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
    assert signature.pair_issue_count > 0
    assert signature.short_edge_count <= 4

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_pair_issues_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints.shapely,
        "coverage_simplify",
        lambda geometries, tolerance, simplify_boundary=True: geometries,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_point_touch_bridge_operator",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_pair_merge_operator",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_patch_union_operator",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("pair patch fallback should have been skipped")
        ),
    )

    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._simplify_coverage_locally(
        polygons,
        [[0], [1]],
        tolerance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is None


def test_coverage_simplify_local_can_disable_patch_union_fallback(monkeypatch):
    shared_boundary = [
        (5.0, 0.0),
        (5.0, 1.0),
        (4.875, 1.0),
        (4.875, 1.5),
        (5.0, 1.5),
        (5.0, 2.0),
        (4.875, 2.0),
        (4.875, 2.5),
        (5.0, 2.5),
        (5.0, 6.0),
    ]
    left = Polygon([(0.0, 0.0), *shared_boundary, (0.0, 6.0)])
    right = Polygon(
        [
            shared_boundary[0],
            (10.0, 0.0),
            (10.0, 6.0),
            shared_boundary[-1],
            *reversed(shared_boundary[1:-1]),
        ]
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_pair_issues_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints.shapely,
        "coverage_simplify",
        lambda geometries, tolerance, simplify_boundary=True: geometries,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_patch_union_operator",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("patch union fallback should have been disabled")
        ),
    )

    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._simplify_coverage_locally(
        [left, right],
        [[0], [1]],
        tolerance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
        enable_patch_union_fallback=False,
    )

    if candidate is not None:
        assert candidate.operator_applied == {"coverage_cpp_cluster_rewrite": 1}


def test_coverage_simplify_local_can_disable_pair_cluster_rescue(monkeypatch):
    polygons = [box(0.0, 0.0, 1.0, 1.0), box(1.0, 1.0, 2.0, 2.0)]

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_pair_issues_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints.shapely,
        "coverage_simplify",
        lambda geometries, tolerance, simplify_boundary=True: geometries,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_point_touch_bridge_operator",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("pair cluster rescue should have been disabled")
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_pair_merge_operator",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("pair cluster rescue should have been disabled")
        ),
    )

    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._simplify_coverage_locally(
        polygons,
        [[0], [1]],
        tolerance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
        enable_patch_union_fallback=False,
        enable_pair_cluster_rescue=False,
    )

    assert candidate is None


def test_coverage_simplify_local_attempts_direct_pair_rewrite_before_simplify(
    monkeypatch,
):
    class _StopCoverageSimplify(RuntimeError):
        pass

    polygons = [box(0.0, 0.0, 1.0, 1.0), box(1.0, 1.0, 2.0, 2.0)]
    events: list[str] = []

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_pair_issues_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_coverage_signature_improves",
        lambda *args, **kwargs: True,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_coverage_signature_satisfies_scale_contract",
        lambda *args, **kwargs: True,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_difference_metrics_have_small_fidelity_drift",
        lambda *args, **kwargs: True,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_direct_pair_issue_cluster_candidates",
        lambda *args, **kwargs: (events.append("direct_pair"), [])[1],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_defect_cluster_descriptors",
        lambda *args, **kwargs: [{"indices": [0, 1], "kind": "pair_issue"}],
    )

    def record_coverage_simplify(*args, **kwargs):
        events.append("coverage_simplify")
        raise _StopCoverageSimplify

    monkeypatch.setattr(
        cleaning_footprints.shapely,
        "coverage_simplify",
        record_coverage_simplify,
    )
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    with pytest.raises(_StopCoverageSimplify):
        cleaning_footprints._simplify_coverage_locally(
            polygons,
            [[0], [1]],
            tolerance=0.5,
            grid=0.03125,
            diagnostics=diagnostics,
            enable_patch_union_fallback=False,
            enable_pair_cluster_rescue=True,
        )

    assert events[0] == "direct_pair"
    assert "coverage_simplify" in events


def test_coverage_simplify_local_skips_pair_cluster_rescue_for_light_close_pair(
    monkeypatch,
):
    polygons = [box(0.0, 0.0, 1.0, 1.0), box(1.1, 0.0, 2.1, 1.0)]

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_pair_issues_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints.shapely,
        "coverage_simplify",
        lambda geometries, tolerance, simplify_boundary=True: geometries,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_point_touch_bridge_operator",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("pair cluster rescue should have been skipped")
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_pair_merge_operator",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("pair cluster rescue should have been skipped")
        ),
    )

    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._simplify_coverage_locally(
        polygons,
        [[0], [1]],
        tolerance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
        enable_patch_union_fallback=False,
        enable_pair_cluster_rescue=True,
    )

    assert candidate is None


def test_coverage_simplify_pair_issue_fast_path_resolves_close_gap_pair():
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._simplify_coverage_pair_issues_graphically(
        [box(0.0, 0.0, 1.0, 1.0), box(1.01, 0.0, 2.01, 1.0)],
        [[0], [1]],
        tolerance=0.5,
        grid=0.01,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.signature.pair_issue_count == 0
    assert any(
        operator.startswith("coverage_pair_issue_bridge_")
        for operator in candidate.operator_applied
    )


def test_coverage_simplify_pair_issue_fast_path_prefers_best_full_candidate(
    monkeypatch,
):
    original = [box(0.0, 0.0, 1.0, 1.0), box(1.0, 1.0, 2.0, 2.0)]
    bridge_candidate = [box(0.0, 0.0, 2.0, 2.0)]
    shrink_candidate = [box(0.0, 0.0, 1.0, 1.0), box(1.3, 1.3, 2.3, 2.3)]
    source_map = [[0], [1]]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    original_key = cleaning_footprints._polygon_sequence_key(original)
    bridge_key = cleaning_footprints._polygon_sequence_key(bridge_candidate)
    shrink_key = cleaning_footprints._polygon_sequence_key(shrink_candidate)

    signature_map = {
        original_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.0,
            pair_issue_count=1,
            point_touch_count=1,
            close_pair_count=0,
            min_pair_clearance=0.0,
            short_edge_count=4,
            min_edge_length=0.1,
            vertex_count=8,
            ring_contact_count=0,
        ),
        bridge_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.5,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=6,
            min_edge_length=0.08,
            vertex_count=12,
            ring_contact_count=0,
        ),
        shrink_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.5,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=4,
            min_edge_length=0.2,
            vertex_count=8,
            ring_contact_count=0,
        ),
    }

    difference_map = {
        (original_key, bridge_key): {
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.2,
            "symmetric_difference_area": 0.2,
            "union_area_delta": 0.2,
        },
        (original_key, shrink_key): {
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.05,
            "symmetric_difference_area": 0.05,
            "union_area_delta": 0.05,
        },
    }

    def fake_pair_issue_candidates(cache, polygons, **kwargs):
        key = cleaning_footprints._polygon_sequence_key(polygons)
        if key == original_key:
            return [(0, 1, 0.0, "point")]
        return []

    def fake_signature(cache, polygons, **kwargs):
        return signature_map[cleaning_footprints._polygon_sequence_key(polygons)]

    def fake_difference(cache, reference, candidate):
        reference_key = cleaning_footprints._polygon_sequence_key(reference)
        candidate_key = cleaning_footprints._polygon_sequence_key(candidate)
        if (reference_key, candidate_key) in difference_map:
            return difference_map[(reference_key, candidate_key)]
        return {
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        }

    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_pair_issue_candidates",
        fake_pair_issue_candidates,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_defect_signature",
        fake_signature,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_difference_area_metrics",
        fake_difference,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_edit_zone",
        lambda *args, **kwargs: box(-10.0, -10.0, 10.0, 10.0),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_change_outside_edit_zone",
        lambda *args, **kwargs: 0.0,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_point_touch_bridge_operator",
        lambda *args, **kwargs: (bridge_candidate, [[0, 1]]),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iter_pair_issue_shrink_radii",
        lambda **kwargs: (0.25,),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_smaller_polygon_shrink_operator",
        lambda *args, **kwargs: (shrink_candidate, source_map),
    )

    candidate = cleaning_footprints._simplify_coverage_pair_issues_graphically(
        original,
        source_map,
        tolerance=0.5,
        grid=0.01,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert cleaning_footprints._polygon_sequence_key(
        candidate.polygons
    ) == cleaning_footprints._polygon_sequence_key(shrink_candidate)
    assert candidate.operator_applied == {"coverage_pair_issue_shrink_0.250": 1}


def test_direct_pair_issue_cluster_candidates_rewrite_close_gap_pair():
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    candidates = cleaning_footprints._direct_pair_issue_cluster_candidates(
        [box(0.0, 0.0, 1.0, 1.0), box(1.01, 0.0, 2.01, 1.0)],
        [[0], [1]],
        tolerance=0.5,
        grid=0.01,
        diagnostics=diagnostics,
    )

    bridge_candidates = [
        (affected_indices, operator, polygons, source_map)
        for affected_indices, operator, polygons, source_map in candidates
        if operator.startswith("coverage_pair_issue_bridge_")
    ]
    assert bridge_candidates
    affected_indices, operator, polygons, source_map = bridge_candidates[0]
    signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
    assert tuple(affected_indices) == (0, 1)
    assert signature.pair_issue_count == 0
    assert source_map == [[0, 1]]


def test_apply_local_smaller_polygon_shrink_operator_separates_point_touch_pair():
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    polygons = [box(0.0, 0.0, 1.0, 1.0), box(1.0, 1.0, 2.0, 2.0)]
    candidate = cleaning_footprints._apply_local_smaller_polygon_shrink_operator(
        polygons,
        [[0], [1]],
        radius=0.375,
        grid=0.01,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    polygons, source_map = candidate
    signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
    assert signature.pair_issue_count == 0
    assert source_map == [[0], [1]]


def test_direct_pair_issue_cluster_candidates_include_shrink_candidate_for_point_touch_pair():
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    candidates = cleaning_footprints._direct_pair_issue_cluster_candidates(
        [box(0.0, 0.0, 1.0, 1.0), box(1.0, 1.0, 2.0, 2.0)],
        [[0], [1]],
        tolerance=0.5,
        grid=0.01,
        diagnostics=diagnostics,
    )

    shrink_candidates = [
        (affected_indices, operator, polygons, source_map)
        for affected_indices, operator, polygons, source_map in candidates
        if operator.startswith("coverage_pair_issue_shrink_")
    ]

    assert shrink_candidates
    assert any(
        tuple(affected_indices) == (0, 1)
        and
        cleaning_footprints._coverage_defect_signature(
            polygons,
            target_scale=0.5,
        ).pair_issue_count
        == 0
        and source_map == [[0], [1]]
        for affected_indices, _, polygons, source_map in shrink_candidates
    )


def test_should_attempt_local_coverage_candidate_skips_when_global_is_already_good():
    reference_signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.2,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=4,
        min_edge_length=0.2,
        vertex_count=8,
    )
    global_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=[box(0.0, 0.0, 10.0, 10.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.7,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=0.7,
            vertex_count=4,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.4,
            "candidate_minus_reference_area": 0.3,
            "symmetric_difference_area": 0.6,
            "union_area_delta": -0.1,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=50.0,
        area_balance_budget=10.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_simplify_global": 1},
        operator_applied={"coverage_simplify_global": 1},
    )

    assert (
        cleaning_footprints._should_attempt_local_coverage_candidate(
            reference_signature,
            global_candidate,
            target_scale=0.5,
            grid=0.03125,
            purpose="coverage",
        )
        is False
    )


def test_should_attempt_meshing_local_candidate_keeps_mixed_residual_cleanup_when_contract_fails():
    reference_signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.0,
        pair_issue_count=3,
        point_touch_count=2,
        close_pair_count=1,
        min_pair_clearance=0.0,
        short_edge_count=5,
        min_edge_length=0.1,
        vertex_count=24,
    )
    global_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=[box(0.0, 0.0, 10.0, 10.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.0,
            pair_issue_count=3,
            point_touch_count=2,
            close_pair_count=1,
            min_pair_clearance=0.0,
            short_edge_count=4,
            min_edge_length=0.15,
            vertex_count=18,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.2,
            "candidate_minus_reference_area": 0.1,
            "symmetric_difference_area": 0.3,
            "union_area_delta": -0.1,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=10.0,
        area_balance_budget=1.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_simplify_global": 1},
        operator_applied={"coverage_simplify_global": 1},
    )

    assert (
        cleaning_footprints._should_attempt_meshing_local_candidate(
            reference_signature,
            global_candidate,
            target_scale=0.5,
            grid=0.03125,
        )
        is True
    )


def test_should_attempt_meshing_local_candidate_keeps_close_gap_cleanup():
    reference_signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.11,
        pair_issue_count=1,
        point_touch_count=0,
        close_pair_count=1,
        min_pair_clearance=0.11,
        short_edge_count=3,
        min_edge_length=0.14,
        vertex_count=18,
    )
    global_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=[box(0.0, 0.0, 10.0, 10.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.12,
            pair_issue_count=1,
            point_touch_count=0,
            close_pair_count=1,
            min_pair_clearance=0.12,
            short_edge_count=3,
            min_edge_length=0.14,
            vertex_count=18,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.2,
            "candidate_minus_reference_area": 0.1,
            "symmetric_difference_area": 0.3,
            "union_area_delta": -0.1,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=10.0,
        area_balance_budget=1.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_simplify_global": 1},
        operator_applied={"coverage_simplify_global": 1},
    )

    assert (
        cleaning_footprints._should_attempt_meshing_local_candidate(
            reference_signature,
            global_candidate,
            target_scale=0.5,
            grid=0.03125,
        )
        is True
    )
    assert (
        cleaning_footprints._should_attempt_local_coverage_candidate(
            reference_signature,
            global_candidate,
            target_scale=0.5,
            grid=0.03125,
            purpose="meshing",
        )
        is True
    )


def test_should_attempt_meshing_local_candidate_for_residual_short_edges():
    reference_signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.43,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=5,
        min_edge_length=0.43,
        vertex_count=24,
    )
    global_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=[box(0.0, 0.0, 10.0, 10.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.43,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=5,
            min_edge_length=0.43,
            vertex_count=18,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.2,
            "candidate_minus_reference_area": 0.1,
            "symmetric_difference_area": 0.3,
            "union_area_delta": -0.1,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=10.0,
        area_balance_budget=1.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_simplify_global": 1},
        operator_applied={"coverage_simplify_global": 1},
    )

    assert (
        cleaning_footprints._should_attempt_meshing_local_candidate(
            reference_signature,
            global_candidate,
            target_scale=0.5,
            grid=0.03125,
        )
        is True
    )
    assert (
        cleaning_footprints._should_attempt_local_coverage_candidate(
            reference_signature,
            global_candidate,
            target_scale=0.5,
            grid=0.03125,
            purpose="meshing",
        )
        is True
    )


def test_select_post_coverage_candidates_for_evaluation_skips_dominated_runner_up():
    identity_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="identity",
        polygons=[box(0.0, 0.0, 10.0, 10.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.2,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=4,
            min_edge_length=0.2,
            vertex_count=8,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=0.0,
        area_balance_budget=0.0,
        patch_count=0,
        patch_applied_count=0,
        operator_attempts={},
        operator_applied={},
    )
    global_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=[box(0.0, 0.0, 10.0, 10.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.7,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=0.7,
            vertex_count=4,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.4,
            "candidate_minus_reference_area": 0.3,
            "symmetric_difference_area": 0.6,
            "union_area_delta": -0.1,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=50.0,
        area_balance_budget=10.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_simplify_global": 1},
        operator_applied={"coverage_simplify_global": 1},
    )
    local_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="local",
        polygons=[box(0.0, 0.0, 10.0, 10.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.68,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=0.68,
            vertex_count=5,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.8,
            "candidate_minus_reference_area": 0.4,
            "symmetric_difference_area": 1.2,
            "union_area_delta": -0.4,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=50.0,
        area_balance_budget=10.0,
        patch_count=3,
        patch_applied_count=2,
        operator_attempts={"coverage_patch_open_close_0.250": 3},
        operator_applied={"coverage_patch_open_close_0.250": 2},
    )

    selected = cleaning_footprints._select_post_coverage_candidates_for_evaluation(
        identity_candidate,
        [global_candidate, local_candidate],
        target_scale=0.5,
        grid=0.03125,
        output_min_area=15.0,
    )

    assert [candidate.label for candidate in selected] == ["identity", "global"]


def test_select_post_coverage_candidates_for_evaluation_skips_identity_when_global_strongly_dominates():
    identity_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="identity",
        polygons=[box(0.0, 0.0, 10.0, 10.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.2,
            pair_issue_count=1,
            point_touch_count=1,
            close_pair_count=0,
            min_pair_clearance=0.2,
            short_edge_count=24,
            min_edge_length=0.2,
            vertex_count=24,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=50.0,
        area_balance_budget=10.0,
        patch_count=0,
        patch_applied_count=0,
        operator_attempts={},
        operator_applied={},
    )
    global_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=[box(0.0, 0.0, 10.0, 10.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.8,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=0.8,
            vertex_count=4,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.2,
            "candidate_minus_reference_area": 0.1,
            "symmetric_difference_area": 0.3,
            "union_area_delta": -0.1,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=50.0,
        area_balance_budget=10.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_simplify_global": 1},
        operator_applied={"coverage_simplify_global": 1},
    )

    selected = cleaning_footprints._select_post_coverage_candidates_for_evaluation(
        identity_candidate,
        [global_candidate],
        target_scale=0.5,
        grid=0.03125,
        output_min_area=15.0,
    )

    assert [candidate.label for candidate in selected] == ["global"]


def test_evaluate_post_coverage_branch_skips_noop_stages(monkeypatch):
    polygon = box(0.0, 0.0, 4.0, 4.0)
    candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=[polygon],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=4.0,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=4.0,
            vertex_count=4,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=0.0,
        area_balance_budget=0.0,
        patch_count=0,
        patch_applied_count=0,
        operator_attempts={},
        operator_applied={},
    )

    for name in (
        "_reclaim_source_supported_area",
        "_absorb_small_supported_components",
        "_simplify_polygons_for_meshing",
        "_regularize_low_clearance_polygons",
    ):
        monkeypatch.setattr(
            cleaning_footprints,
            name,
            lambda *args, **kwargs: (_ for _ in ()).throw(
                AssertionError(f"{name} should have been skipped")
            ),
        )

    monkeypatch.setattr(
        cleaning_footprints,
        "_recover_source_supported_coordinates",
        lambda polygons, source_map, **kwargs: (polygons, source_map),
    )

    branch = cleaning_footprints._evaluate_post_coverage_branch(
        candidate,
        reference_union=polygon,
        source_lookup={0: polygon},
        raw_support_union=polygon,
        min_feature_size=0.5,
        source_recovery_scale=0.5,
        grid=0.03125,
        min_area=0.0,
        output_min_area=15.0,
        min_hole_area=0.0,
        cache=cleaning_footprints._CoverageEvalCache(),
    )

    assert branch.final_output_polygons[0].equals_exact(polygon, tolerance=0.0)


def test_regularize_coverage_for_meshing_skips_local_search_when_global_is_sufficient(
    monkeypatch,
):
    polygon = Polygon(
        [
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (2.0, 4.0),
            (2.0, 4.2),
            (1.8, 4.2),
            (1.8, 4.0),
            (0.0, 4.0),
        ]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    global_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=[box(0.0, 0.0, 4.0, 4.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=4.0,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=4.0,
            vertex_count=4,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.04,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.04,
            "union_area_delta": -0.04,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=20.0,
        area_balance_budget=1.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_simplify_global": 1},
        operator_applied={"coverage_simplify_global": 1},
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: global_candidate,
    )

    local_called = {"value": False}

    def fail_if_called(*args, **kwargs):
        local_called["value"] = True
        raise AssertionError("local coverage search should have been skipped")

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        fail_if_called,
    )

    polygons, source_map = cleaning_footprints._regularize_coverage_for_meshing(
        [polygon],
        [[0]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert local_called["value"] is False
    assert diagnostics["coverage_meshing_regularization_local_candidate_attempted"] is False
    assert len(polygons) == 1
    assert source_map == [[0]]


def test_regularize_coverage_for_meshing_uses_graphical_fast_path_when_sufficient(
    monkeypatch,
):
    polygon = Polygon(
        [
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (2.0, 4.0),
            (2.0, 4.2),
            (1.8, 4.2),
            (1.8, 4.0),
            (0.0, 4.0),
        ]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    graphical_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="local",
        polygons=[box(0.0, 0.0, 4.0, 4.0)],
        source_map=[[0]],
        signature=cleaning_footprints._CoverageDefectSignature(
            min_clearance=4.0,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=4.0,
            vertex_count=4,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.04,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.04,
            "union_area_delta": -0.04,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=20.0,
        area_balance_budget=1.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_graph_short_edges": 1},
        operator_applied={"coverage_graph_short_edges": 1},
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: graphical_candidate,
    )

    def fail_if_called(*args, **kwargs):
        raise AssertionError("heavier meshing regularization should have been skipped")

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        fail_if_called,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        fail_if_called,
    )

    polygons, source_map = cleaning_footprints._regularize_coverage_for_meshing(
        [polygon],
        [[0]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert diagnostics["coverage_meshing_regularization_selected_branch"] == "local"
    assert diagnostics["coverage_meshing_regularization_local_candidate_attempted"] is False
    assert diagnostics["coverage_meshing_regularization_operator_applied"] == {
        "coverage_graph_short_edges": 1
    }
    assert len(polygons) == 1
    assert source_map == [[0]]


def test_regularize_coverage_for_meshing_returns_early_when_no_candidate_exists(
    monkeypatch,
):
    polygons = [box(0.0, 0.0, 1.0, 1.0), box(3.0, 3.0, 4.0, 4.0)]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_polygons_for_meshing",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "unary_union",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("unary_union should not run when no candidate exists")
        ),
    )

    result_polygons, result_sources = cleaning_footprints._regularize_coverage_for_meshing(
        polygons,
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert result_polygons == polygons
    assert result_sources == [[0], [1]]
    assert diagnostics["coverage_meshing_regularization_selected_branch"] == "identity"
    assert diagnostics["coverage_meshing_regularization_applied"] is False


def test_regularize_coverage_for_meshing_applies_residual_pair_issue_rescue(
    monkeypatch,
):
    polygons = [box(0.0, 0.0, 10.0, 10.0), box(10.0, 10.0, 20.0, 20.0)]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    identity_signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
    identity_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="local",
        polygons=polygons,
        source_map=[[0], [1]],
        signature=identity_signature,
        difference_metrics={
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=0.0,
        area_balance_budget=0.0,
        patch_count=0,
        patch_applied_count=0,
        operator_attempts={},
        operator_applied={},
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: identity_candidate,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_residual_scale_polygons",
        lambda *args, **kwargs: None,
    )

    result_polygons, result_sources = cleaning_footprints._regularize_coverage_for_meshing(
        polygons,
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    signature = cleaning_footprints._coverage_defect_signature(
        result_polygons,
        target_scale=0.5,
    )
    assert signature.pair_issue_count == 0
    assert sorted(sorted(indices) for indices in result_sources) in (
        [[0], [1]],
        [[0, 1]],
    )
    assert diagnostics["coverage_meshing_regularization_selected_branch"] == (
        "residual_pair_issue_rescue"
    )
    operator_applied = diagnostics["coverage_meshing_regularization_operator_applied"]
    assert any(
        key.startswith("coverage_pair_issue_") for key in operator_applied
    )


def test_regularize_coverage_for_meshing_applies_residual_point_bridge_rescue(
    monkeypatch,
):
    polygons = [box(0.0, 0.0, 10.0, 10.0), box(10.0, 10.0, 11.0, 11.0)]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    identity_signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=1.0,
    )
    identity_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="local",
        polygons=polygons,
        source_map=[[0], [1]],
        signature=identity_signature,
        difference_metrics={
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=0.0,
        area_balance_budget=0.0,
        patch_count=0,
        patch_applied_count=0,
        operator_attempts={},
        operator_applied={},
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: identity_candidate,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_residual_scale_polygons",
        lambda *args, **kwargs: None,
    )

    result_polygons, result_sources = cleaning_footprints._regularize_coverage_for_meshing(
        polygons,
        [[0], [1]],
        min_segment_length=1.0,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    signature = cleaning_footprints._coverage_defect_signature(
        result_polygons,
        target_scale=1.0,
    )
    operator_applied = diagnostics["coverage_meshing_regularization_operator_applied"]
    assert signature.pair_issue_count == 0
    assert len(result_polygons) == 1
    assert result_sources == [[0, 1]]
    assert diagnostics["coverage_meshing_regularization_selected_branch"] == (
        "residual_pair_issue_rescue"
    )
    assert any(
        key.startswith("coverage_pair_issue_point_")
        for key in operator_applied
    )


def test_regularize_coverage_for_meshing_prefers_pair_rescue_with_better_scale_contract(
    monkeypatch,
):
    left = box(0.0, 0.0, 10.0, 10.0)
    right = box(10.0, 10.0, 20.0, 20.0)
    polygons = [left, right]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    identity_signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
    identity_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="local",
        polygons=polygons,
        source_map=[[0], [1]],
        signature=identity_signature,
        difference_metrics={
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=0.0,
        area_balance_budget=0.0,
        patch_count=0,
        patch_applied_count=0,
        operator_attempts={},
        operator_applied={},
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: identity_candidate,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_residual_scale_polygons",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iter_pair_issue_bridge_radii",
        lambda **kwargs: (0.03125, 0.75),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iter_pair_issue_shrink_radii",
        lambda **kwargs: (),
    )

    def fake_point_bridge_operator(
        pair_polygons,
        pair_sources,
        *,
        radius,
        target_scale=None,
        grid,
        diagnostics,
    ):
        if radius <= 0.05:
            adjusted_right = Polygon(
                [
                    (10.03125, 10.0),
                    (20.0, 10.0),
                    (20.0, 20.0),
                    (10.0, 20.0),
                    (10.0, 10.03125),
                ]
            )
        else:
            adjusted_right = Polygon(
                [
                    (10.75, 10.0),
                    (20.0, 10.0),
                    (20.0, 20.0),
                    (10.0, 20.0),
                    (10.0, 10.75),
                ]
            )
        return [pair_polygons[0], adjusted_right], [list(pair_sources[0]), list(pair_sources[1])]

    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_point_touch_bridge_operator",
        fake_point_bridge_operator,
    )

    result_polygons, result_sources = cleaning_footprints._regularize_coverage_for_meshing(
        polygons,
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    signature = cleaning_footprints._coverage_defect_signature(
        result_polygons,
        target_scale=0.5,
    )
    assert signature.pair_issue_count == 0
    assert signature.short_edge_count == 0
    assert signature.min_edge_length >= 0.75
    assert result_sources == [[0], [1]]
    operator_applied = diagnostics["coverage_meshing_regularization_operator_applied"]
    assert operator_applied == {"coverage_pair_issue_point_residual_0.750": 1}


def test_regularize_coverage_for_meshing_runs_direct_pair_rescue_on_current_best(
    monkeypatch,
):
    original_polygons = [box(0.0, 0.0, 10.0, 10.0), box(10.0, 10.0, 20.0, 20.0)]
    best_polygons = [
        box(0.0, 0.0, 10.0, 10.0),
        Polygon(
            [
                (10.5, 10.0),
                (20.0, 10.0),
                (20.0, 20.0),
                (10.0, 20.0),
                (10.0, 10.5),
            ]
        ),
    ]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    global_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=best_polygons,
        source_map=[[0], [1]],
        signature=cleaning_footprints._coverage_defect_signature(
            best_polygons,
            target_scale=0.5,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=25.0,
        area_balance_budget=1.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_simplify_global": 1},
        operator_applied={"coverage_simplify_global": 1},
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: global_candidate,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_should_attempt_meshing_local_candidate",
        lambda *args, **kwargs: False,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_residual_scale_polygons",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_polygons_for_meshing",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iter_pair_issue_bridge_radii",
        lambda **kwargs: (),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iter_pair_issue_shrink_radii",
        lambda **kwargs: (),
    )

    captured_calls: list[tuple[tuple[int, ...], list[list[int]]]] = []

    def fake_direct_pair_issue_cluster_candidates(
        subset_polygons,
        subset_sources,
        **kwargs,
    ):
        captured_calls.append(
            (
                cleaning_footprints._polygon_sequence_key(subset_polygons),
                subset_sources,
            )
        )
        return []

    monkeypatch.setattr(
        cleaning_footprints,
        "_direct_pair_issue_cluster_candidates",
        fake_direct_pair_issue_cluster_candidates,
    )

    cleaning_footprints._regularize_coverage_for_meshing(
        original_polygons,
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert captured_calls[0][0] == cleaning_footprints._polygon_sequence_key(best_polygons)
    assert captured_calls[0][1] == [[0], [1]]


def test_regularize_coverage_for_meshing_can_use_original_direct_pair_rescue_candidates(
    monkeypatch,
):
    original_polygons = [box(0.0, 0.0, 10.0, 10.0), box(10.0, 10.0, 20.0, 20.0)]
    best_polygons = [
        box(0.0, 0.0, 10.0, 10.0),
        Polygon(
            [
                (10.5, 10.0),
                (20.0, 10.0),
                (20.0, 20.0),
                (10.0, 20.0),
                (10.0, 10.5),
            ]
        ),
    ]
    original_rescue_polygons = [box(0.0, 0.0, 10.0, 10.0), box(10.75, 10.75, 20.0, 20.0)]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    global_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=best_polygons,
        source_map=[[0], [1]],
        signature=cleaning_footprints._coverage_defect_signature(
            best_polygons,
            target_scale=0.5,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=25.0,
        area_balance_budget=1.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_simplify_global": 1},
        operator_applied={"coverage_simplify_global": 1},
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: global_candidate,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_should_attempt_meshing_local_candidate",
        lambda *args, **kwargs: False,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_residual_scale_polygons",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_polygons_for_meshing",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iter_pair_issue_bridge_radii",
        lambda **kwargs: (),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iter_pair_issue_shrink_radii",
        lambda **kwargs: (),
    )

    def fake_direct_pair_issue_cluster_candidates(
        subset_polygons,
        subset_sources,
        **kwargs,
    ):
        if (
            cleaning_footprints._polygon_sequence_key(subset_polygons)
            == cleaning_footprints._polygon_sequence_key(best_polygons)
        ):
            return []
        if (
            cleaning_footprints._polygon_sequence_key(subset_polygons)
            == cleaning_footprints._polygon_sequence_key(original_polygons)
        ):
            return [
                (
                    (0, 1),
                    "coverage_pair_issue_bridge_residual_0.250",
                    original_rescue_polygons,
                    [[0], [1]],
                )
            ]
        return []

    monkeypatch.setattr(
        cleaning_footprints,
        "_direct_pair_issue_cluster_candidates",
        fake_direct_pair_issue_cluster_candidates,
    )

    result_polygons, result_sources = cleaning_footprints._regularize_coverage_for_meshing(
        original_polygons,
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert cleaning_footprints._polygon_sequence_key(
        result_polygons
    ) == cleaning_footprints._polygon_sequence_key(original_rescue_polygons)
    assert result_sources == [[0], [1]]
    operator_applied = diagnostics["coverage_meshing_regularization_operator_applied"]
    assert operator_applied["coverage_simplify_global"] == 1
    assert operator_applied["coverage_pair_issue_bridge_residual_0.250"] == 1


def test_regularize_coverage_for_meshing_rejects_original_direct_pair_rescue_with_residual_pairs(
    monkeypatch,
):
    original_polygons = [box(0.0, 0.0, 10.0, 10.0), box(10.0, 10.0, 20.0, 20.0)]
    best_polygons = [
        box(0.0, 0.0, 10.0, 10.0),
        Polygon(
            [
                (10.5, 10.0),
                (20.0, 10.0),
                (20.0, 20.0),
                (10.0, 20.0),
                (10.0, 10.5),
            ]
        ),
    ]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    global_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="global",
        polygons=best_polygons,
        source_map=[[0], [1]],
        signature=cleaning_footprints._coverage_defect_signature(
            best_polygons,
            target_scale=0.5,
        ),
        difference_metrics={
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=25.0,
        area_balance_budget=1.0,
        patch_count=1,
        patch_applied_count=1,
        operator_attempts={"coverage_simplify_global": 1},
        operator_applied={"coverage_simplify_global": 1},
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: global_candidate,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_should_attempt_meshing_local_candidate",
        lambda *args, **kwargs: False,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_residual_scale_polygons",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_polygons_for_meshing",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iter_pair_issue_bridge_radii",
        lambda **kwargs: (),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iter_pair_issue_shrink_radii",
        lambda **kwargs: (),
    )

    def fake_direct_pair_issue_cluster_candidates(
        subset_polygons,
        subset_sources,
        **kwargs,
    ):
        if (
            cleaning_footprints._polygon_sequence_key(subset_polygons)
            == cleaning_footprints._polygon_sequence_key(best_polygons)
        ):
            return []
        if (
            cleaning_footprints._polygon_sequence_key(subset_polygons)
            == cleaning_footprints._polygon_sequence_key(original_polygons)
        ):
            return [
                (
                    (0, 1),
                    "coverage_pair_issue_bridge_residual_0.250",
                    original_polygons,
                    [[0], [1]],
                )
            ]
        return []

    monkeypatch.setattr(
        cleaning_footprints,
        "_direct_pair_issue_cluster_candidates",
        fake_direct_pair_issue_cluster_candidates,
    )

    result_polygons, result_sources = cleaning_footprints._regularize_coverage_for_meshing(
        original_polygons,
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert cleaning_footprints._polygon_sequence_key(
        result_polygons
    ) == cleaning_footprints._polygon_sequence_key(best_polygons)
    assert result_sources == [[0], [1]]
    assert diagnostics["coverage_meshing_regularization_operator_applied"] == {
        "coverage_simplify_global": 1
    }


def test_regularize_coverage_for_meshing_keeps_raw_rescue_when_postprocess_regresses(
    monkeypatch,
):
    polygons = [box(0.0, 0.0, 10.0, 10.0), box(10.0, 10.0, 20.0, 20.0)]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    identity_signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
    identity_candidate = cleaning_footprints._CoverageSimplifyCandidate(
        label="local",
        polygons=polygons,
        source_map=[[0], [1]],
        signature=identity_signature,
        difference_metrics={
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
        change_outside_edit_zone=0.0,
        edit_zone_area=0.0,
        area_balance_budget=0.0,
        patch_count=0,
        patch_applied_count=0,
        operator_attempts={},
        operator_applied={},
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: identity_candidate,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_residual_scale_polygons",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_polygons_for_meshing",
        lambda polygons, sources, **kwargs: ([box(0.0, 0.0, 10.0, 10.0), box(10.0, 10.0, 20.0, 20.0)], [[0], [1]]),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )

    result_polygons, result_sources = cleaning_footprints._regularize_coverage_for_meshing(
        polygons,
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    signature = cleaning_footprints._coverage_defect_signature(
        result_polygons,
        target_scale=0.5,
    )
    assert signature.pair_issue_count == 0
    assert sorted(sorted(indices) for indices in result_sources) in (
        [[0], [1]],
        [[0, 1]],
    )
    assert diagnostics["coverage_meshing_regularization_selected_branch"] == (
        "residual_pair_issue_rescue"
    )


def test_regularize_coverage_for_meshing_runs_pair_rescue_even_without_primary_candidates(
    monkeypatch,
):
    polygons = [box(0.0, 0.0, 10.0, 10.0), box(10.0, 10.0, 20.0, 20.0)]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_residual_scale_polygons",
        lambda *args, **kwargs: None,
    )

    result_polygons, result_sources = cleaning_footprints._regularize_coverage_for_meshing(
        polygons,
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    signature = cleaning_footprints._coverage_defect_signature(
        result_polygons,
        target_scale=0.5,
    )
    assert signature.pair_issue_count == 0
    assert diagnostics["coverage_meshing_regularization_selected_branch"] == (
        "residual_pair_issue_rescue"
    )


def test_regularize_coverage_for_meshing_retries_ring_contacts_after_pair_rescue(
    monkeypatch,
):
    original_polygons = [box(0.0, 0.0, 10.0, 10.0), box(10.0, 10.0, 20.0, 20.0)]
    rescue_polygons = [
        box(0.0, 0.0, 10.0, 10.0),
        Polygon(
            [
                (10.0, 10.0),
                (15.0, 10.0),
                (20.0, 10.0),
                (20.0, 20.0),
                (10.0, 20.0),
                (10.0, 10.0),
            ]
        ),
    ]
    fixed_polygons = [box(0.0, 0.0, 10.0, 10.0), box(10.5, 10.5, 20.0, 20.0)]
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False
    diagnostics["enable_logging"] = False

    original_signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.25,
        pair_issue_count=1,
        point_touch_count=0,
        close_pair_count=1,
        min_pair_clearance=0.25,
        short_edge_count=1,
        min_edge_length=0.25,
        vertex_count=8,
        ring_contact_count=1,
    )
    rescue_signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.25,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=1,
        min_edge_length=0.25,
        vertex_count=11,
        ring_contact_count=1,
    )
    fixed_signature = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.5,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=0,
        min_edge_length=0.5,
        vertex_count=8,
        ring_contact_count=0,
    )

    def _signature_for(polygons):
        if len(polygons) != 2:
            raise KeyError("unexpected polygon count")
        if polygons[0].equals_exact(original_polygons[0], 0.0) and polygons[
            1
        ].equals_exact(original_polygons[1], 0.0):
            return original_signature
        if polygons[0].equals_exact(rescue_polygons[0], 0.0) and polygons[
            1
        ].equals_exact(rescue_polygons[1], 0.0):
            return rescue_signature
        if polygons[0].equals_exact(fixed_polygons[0], 0.0) and polygons[
            1
        ].equals_exact(fixed_polygons[1], 0.0):
            return fixed_signature
        raise KeyError("unexpected polygon sequence")

    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_defect_signature",
        lambda cache, polygons, **kwargs: _signature_for(polygons),
    )

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_residual_scale_polygons",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_polygons_for_meshing",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, sources, **kwargs: (polygons, sources),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_pair_issue_candidates",
        lambda *args, **kwargs: [],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_direct_pair_issue_cluster_candidates",
        lambda *args, **kwargs: [
            (
                (0, 1),
                "coverage_pair_issue_bridge_residual_0.250",
                rescue_polygons,
                [[0], [1]],
            )
        ],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_difference_area_metrics",
        lambda *args, **kwargs: {
            "reference_minus_candidate_area": 0.0,
            "candidate_minus_reference_area": 0.0,
            "symmetric_difference_area": 0.0,
            "union_area_delta": 0.0,
        },
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_polygon_ring_boundary_contact_count",
        lambda polygon: 1 if polygon is rescue_polygons[1] else 0,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_ring_contact_polygon",
        lambda polygon, **kwargs: [fixed_polygons[1]]
        if polygon is rescue_polygons[1]
        else [polygon],
    )

    result_polygons, result_sources = cleaning_footprints._regularize_coverage_for_meshing(
        original_polygons,
        [[0], [1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert len(result_polygons) == len(fixed_polygons)
    assert result_polygons[0].equals(fixed_polygons[0])
    assert result_polygons[1].equals(fixed_polygons[1])
    assert result_sources == [[0], [1]]
    assert diagnostics["coverage_meshing_regularization_selected_branch"] == (
        "residual_pair_issue_rescue"
    )
    assert diagnostics["coverage_meshing_regularization_pair_issue_count_after"] == 0
    assert diagnostics["coverage_meshing_regularization_ring_contact_count_after"] == 0


def test_builder_boundary_defect_clusters_sorts_native_clusters(monkeypatch):
    polygons = [box(float(index), 0.0, float(index) + 0.5, 0.5) for index in range(64)]

    monkeypatch.setattr(
        cleaning_footprints,
        "_dtcc_builder",
        SimpleNamespace(
            boundary_defect_clusters=lambda *args, **kwargs: [
                {
                    "indices": [9, 4, 7],
                    "kind": "close_pair",
                    "short_edge_count": 2,
                    "pair_issue_count": 1,
                },
                {
                    "indices": [3, 1],
                    "kind": "short_edge_only",
                    "short_edge_count": 1,
                    "pair_issue_count": 0,
                },
                {
                    "indices": [8, 6],
                    "kind": "mixed_pair_short_edge",
                    "short_edge_count": 3,
                    "pair_issue_count": 1,
                },
            ]
        ),
    )

    descriptors = cleaning_footprints._builder_boundary_defect_clusters(
        polygons,
        target_scale=0.5,
        pair_tolerance=0.75,
    )

    assert descriptors == [
        {
            "indices": [1, 3],
            "kind": "short_edge_only",
            "short_edge_count": 1,
            "pair_issue_count": 0,
        },
        {
            "indices": [4, 7, 9],
            "kind": "close_pair",
            "short_edge_count": 2,
            "pair_issue_count": 1,
        },
        {
            "indices": [6, 8],
            "kind": "mixed_pair_short_edge",
            "short_edge_count": 3,
            "pair_issue_count": 1,
        },
    ]


def test_polygon_simplify_rejects_nonlocal_changes(monkeypatch):
    polygon = Polygon(
        [
            (0, 0),
            (8, 0),
            (8, 4),
            (8.2, 4),
            (8.2, 3.8),
            (8.4, 3.8),
            (8.4, 8),
            (0, 8),
        ]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False
    monkeypatch.setattr(
        cleaning_footprints,
        "_try_close_courtyard_passage",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_try_polygon_clearance_opening",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iteratively_open_polygon_short_edges",
        lambda *args, **kwargs: None,
    )

    monkeypatch.setattr(
        cleaning_footprints.shapely,
        "simplify",
        lambda geometry, tolerance, preserve_topology=True: box(-2, -2, 10, 10),
    )

    polygons, source_map = cleaning_footprints._simplify_polygons_for_meshing(
        [polygon],
        [[7]],
        min_segment_length=0.5,
        grid=0.125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert source_map == [[7]]
    assert len(polygons) == 1
    assert polygons[0].equals_exact(polygon, tolerance=0.0)
    assert diagnostics["polygon_simplify_applied"] is False
    assert diagnostics["polygon_simplify_rejected_nonlocal_count"] >= 1
    assert diagnostics["polygon_simplify_change_outside_edit_zone"] > 0.0


def test_polygon_simplify_rejects_area_imbalanced_changes(monkeypatch):
    polygon = Polygon(
        [
            (0, 0),
            (8, 0),
            (8, 4),
            (8.4, 4),
            (8.4, 5),
            (8, 5),
            (8, 8),
            (0, 8),
        ]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False
    monkeypatch.setattr(
        cleaning_footprints,
        "_try_close_courtyard_passage",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_try_polygon_clearance_opening",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_iteratively_open_polygon_short_edges",
        lambda *args, **kwargs: None,
    )

    monkeypatch.setattr(
        cleaning_footprints.shapely,
        "simplify",
        lambda geometry, tolerance, preserve_topology=True: box(0, 0, 8, 8),
    )

    polygons, source_map = cleaning_footprints._simplify_polygons_for_meshing(
        [polygon],
        [[7]],
        min_segment_length=0.5,
        grid=0.125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert source_map == [[7]]
    assert len(polygons) == 1
    assert polygons[0].equals_exact(polygon, tolerance=0.0)
    assert diagnostics["polygon_simplify_applied"] is False
    assert diagnostics["polygon_simplify_rejected_area_imbalance_count"] >= 1
    assert diagnostics["polygon_simplify_reference_minus_candidate_area"] == pytest.approx(
        0.0,
        abs=1e-12,
    )


def test_courtyard_passage_operator_returns_hole_candidate():
    outer = box(0, 0, 10, 10)
    courtyard = box(3, 3, 7, 7)
    passage = box(0, 4.9, 3.2, 5.1)
    polygon = outer.difference(courtyard.union(passage))
    assert isinstance(polygon, Polygon)
    assert len(polygon.interiors) == 0

    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._try_close_courtyard_passage(
        polygon,
        min_clearance=0.5,
        grid=0.125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.operator == "courtyard_passage"
    assert len(candidate.polygon.interiors) == 1
    assert candidate.edit_zone.area > 0.0


def test_ring_contact_fill_operator_returns_local_single_polygon():
    touching_hole = Polygon(
        [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)],
        [[(0, 5), (2, 4), (3, 5), (2, 6), (0, 5)]],
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._try_fill_ring_contact_vertices(
        touching_hole,
        min_clearance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.operator == "ring_contact_fill"
    assert candidate.polygon.area > touching_hole.area
    assert not cleaning_footprints._polygon_has_ring_boundary_contacts(
        candidate.polygon
    )
    assert candidate.edit_zone.area > 0.0
    assert candidate.area_balance_budget_override is not None


def test_ring_contact_connector_operator_returns_local_single_polygon():
    touching_hole = Polygon(
        [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)],
        [[(0, 5), (2, 4), (3, 5), (2, 6), (0, 5)]],
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._try_ring_contact_connector_fill(
        touching_hole,
        min_clearance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.operator == "ring_contact_connector"
    signature = cleaning_footprints._polygon_defect_signature(
        candidate.polygon,
        target_scale=0.5,
    )
    assert signature.ring_contact_count == 0
    assert candidate.polygon.area > touching_hole.area
    assert candidate.edit_zone.area > 0.0


def test_self_clearance_connector_resolves_case23_style_vertex_edge_junction():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674233.5625 6579727.8125, 674235.59375 6579719.84375, "
        "674248.71875 6579722.9375, 674241.25 6579748.75, "
        "674247.875 6579750.28125, 674245.25 6579759.84375, "
        "674238.9375 6579758.3125, 674241.25 6579748.6875, "
        "674235.875 6579747.21875, 674240.40625 6579729.46875, "
        "674233.5625 6579727.8125))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._try_polygon_self_clearance_connector_fill(
        polygon,
        min_clearance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.operator == "self_clearance_connector"
    signature = cleaning_footprints._polygon_defect_signature(
        candidate.polygon,
        target_scale=0.5,
    )
    assert cleaning_footprints._signature_satisfies_scale_contract(
        signature,
        target_scale=0.5,
        grid=0.03125,
    )
    assert candidate.edit_zone.area > 0.0
    assert candidate.area_balance_budget_override is not None


def test_self_clearance_connector_cuts_case23_style_hole_hole_wall():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674348.09375 6579639.875, 674351.9375 6579626.1875, "
        "674340.28125 6579622.59375, 674342.3125 6579614.8125, "
        "674354.21875 6579618.09375, 674355.75 6579612.75, "
        "674323.1875 6579602.96875, 674321.375 6579609.1875, "
        "674335.28125 6579613.15625, 674333.3125 6579619.90625, "
        "674321.46875 6579616.8125, 674309.5 6579614.59375, "
        "674316.71875 6579590.6875, 674357.90625 6579602.78125, "
        "674425.0625 6579622.5625, 674410.75 6579670.5625, "
        "674381.1875 6579660.8125, 674362.71875 6579655.5625, "
        "674345.3125 6579650.96875, 674325.375 6579644.4375, "
        "674302.4375 6579638.03125, 674306.28125 6579625.25, "
        "674323.625 6579630.34375, 674324.71875 6579626.125, "
        "674307.59375 6579620.96875, 674308.25 6579618.75, "
        "674322.6875 6579621.34375, 674345.4375 6579628.15625, "
        "674343.9375 6579631.9375, 674334.6875 6579629.15625, "
        "674333.78125 6579632.1875, 674329.03125 6579630.9375, "
        "674328.09375 6579634.4375, 674348.09375 6579639.875), "
        "(674407.21875 6579635.71875, 674405.84375 6579635.28125, "
        "674407.8125 6579629.15625, 674403.28125 6579627.75, "
        "674401.03125 6579635.21875, 674413.09375 6579638.625, "
        "674415.28125 6579631.375, 674410.09375 6579629.875, "
        "674408.96875 6579634.09375, 674407.875 6579633.71875, "
        "674407.21875 6579635.71875), "
        "(674412.40625 6579643.6875, 674408.59375 6579642.59375, "
        "674407.4375 6579646.21875, 674400.6875 6579644.1875, "
        "674398.3125 6579651.75, 674400.6875 6579652.375, "
        "674398.84375 6579659.5625, 674396.1875 6579658.875, "
        "674396.9375 6579656.78125, 674388.46875 6579654.46875, "
        "674387.53125 6579657.46875, 674402.375 6579662.28125, "
        "674405.0625 6579653.6875, 674409.125 6579655.09375, "
        "674412.40625 6579643.6875), "
        "(674386.53125 6579643.6875, 674390.34375 6579644.75, "
        "674392.34375 6579638.34375, 674388.53125 6579637.15625, "
        "674386.53125 6579643.6875), "
        "(674376.125 6579618.21875, 674373.0625 6579628.0625, "
        "674367.71875 6579626.53125, 674366.28125 6579631.3125, "
        "674368.96875 6579632.03125, 674366.1875 6579642.40625, "
        "674385.1875 6579647.9375, 674386.5 6579643.6875, "
        "674379.71875 6579641.6875, 674383.1875 6579630.75, "
        "674376.3125 6579628.875, 674379.46875 6579619.1875, "
        "674376.125 6579618.21875), "
        "(674361.5 6579614.125, 674353.90625 6579641.6875, "
        "674359.09375 6579643.28125, 674367 6579615.75, "
        "674361.5 6579614.125))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._try_polygon_self_clearance_connector_fill(
        polygon,
        min_clearance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.operator in {
        "self_clearance_connector_cut",
        "hole_pair_clearance_merge",
    }
    signature = cleaning_footprints._polygon_defect_signature(
        candidate.polygon,
        target_scale=0.5,
    )
    assert signature.clearance is not None
    assert signature.clearance > 0.03125
    assert len(candidate.polygon.interiors) == len(polygon.interiors) - 1
    assert candidate.edit_zone.area > 0.0


def test_self_clearance_connector_fills_case54_style_shell_hole_wall():
    polygon = Polygon(
        [(0, 0), (12, 0), (12, 10), (0, 10), (0, 0)],
        [[(0.25, 3), (5, 3), (5, 7), (0.25, 7), (0.25, 3)]],
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._try_polygon_self_clearance_connector_fill(
        polygon,
        min_clearance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.operator == "self_clearance_connector_fill"
    signature = cleaning_footprints._polygon_defect_signature(
        candidate.polygon,
        target_scale=0.5,
    )
    assert signature.clearance is not None
    assert signature.clearance >= 0.5 - 1.0e-9
    assert len(candidate.polygon.interiors) == len(polygon.interiors)
    assert candidate.polygon.area > polygon.area
    assert candidate.edit_zone.area > 0.0


def test_self_clearance_connector_fills_case54_style_same_ring_exterior_slit():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674533.59375 6581130.875, 674534.75 6581130.59375, "
        "674534.03125 6581130.5, 674538.03125 6581098.75, "
        "674562.3125 6581101.8125, 674562.40625 6581101.3125, "
        "674537.96875 6581098.03125, 674543.75 6581054.9375, 674544.28125 6581055, "
        "674545.5625 6581049.65625, 674544.46875 6581048.5625, "
        "674546.25 6581046.75, 674546.59375 6581045.21875, "
        "674559.59375 6581047, 674583.21875 6581050.8125, "
        "674601.03125 6581055.03125, 674610.9375 6581057.09375, "
        "674604.6875 6581087.28125, 674602.09375 6581106.625, "
        "674562.9375 6581101.375, 674562.875 6581101.875, "
        "674575.71875 6581103.5, 674571.71875 6581135.25, "
        "674555.84375 6581133.25, 674546.65625 6581132.40625, "
        "674570.78125 6581135.71875, 674628 6581143.875, "
        "674708.375 6581155.1875, 674701.34375 6581205.84375, "
        "674696.625 6581205.15625, 674695.5625 6581212.3125, "
        "674672.40625 6581209.0625, 674669.8125 6581226.96875, "
        "674697.78125 6581231.03125, 674694.34375 6581261.90625, "
        "674690.90625 6581286.34375, 674591.375 6581272.09375, "
        "674597.6875 6581224.40625, 674598.21875 6581224.46875, "
        "674599.28125 6581216.90625, 674640.9375 6581222.875, "
        "674643.53125 6581205.03125, 674638.53125 6581204.34375, "
        "674637.59375 6581211.0625, 674615.6875 6581207.96875, "
        "674615.28125 6581210.78125, 674595.96875 6581208.28125, "
        "674593.84375 6581224.0625, 674582.875 6581222.6875, "
        "674529.40625 6581215.78125, 674529.90625 6581212.15625, "
        "674523.1875 6581206.84375, 674531 6581203.46875, "
        "674531.46875 6581199.875, 674524.4375 6581198.9375, "
        "674533.59375 6581130.875), "
        "(674670.53125 6581260.34375, 674669.125 6581269.59375, "
        "674678.03125 6581270.875, 674682 6581240.875, "
        "674671.1875 6581239.28125, 674670.03125 6581240.53125, "
        "674667.28125 6581259.84375, 674670.53125 6581260.34375), "
        "(674623.375 6581176.9375, 674627.8125 6581144.84375, "
        "674622.0625 6581176.84375, 674623.375 6581176.9375))"
    )
    before_signature = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False
    candidate = cleaning_footprints._try_polygon_self_clearance_connector_fill(
        polygon,
        min_clearance=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    assert candidate.operator == "self_clearance_connector_segment"
    signature = cleaning_footprints._polygon_defect_signature(
        candidate.polygon,
        target_scale=0.5,
    )
    assert signature.clearance is not None
    assert signature.clearance > before_signature.clearance
    assert candidate.polygon.area > polygon.area
    assert candidate.edit_zone.area > 0.0


def test_remove_meshing_hostile_holes_drops_case54_style_sliver_hole_only():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674533.59375 6581130.875, 674534.75 6581130.59375, "
        "674534.03125 6581130.5, 674538.03125 6581098.75, "
        "674562.3125 6581101.8125, 674562.40625 6581101.3125, "
        "674537.96875 6581098.03125, 674543.75 6581054.9375, 674544.28125 6581055, "
        "674545.5625 6581049.65625, 674544.46875 6581048.5625, "
        "674546.25 6581046.75, 674546.59375 6581045.21875, "
        "674559.59375 6581047, 674583.21875 6581050.8125, "
        "674601.03125 6581055.03125, 674610.9375 6581057.09375, "
        "674604.6875 6581087.28125, 674602.09375 6581106.625, "
        "674562.9375 6581101.375, 674562.875 6581101.875, "
        "674575.71875 6581103.5, 674571.71875 6581135.25, "
        "674555.84375 6581133.25, 674546.65625 6581132.40625, "
        "674570.78125 6581135.71875, 674628 6581143.875, "
        "674708.375 6581155.1875, 674701.34375 6581205.84375, "
        "674696.625 6581205.15625, 674695.5625 6581212.3125, "
        "674672.40625 6581209.0625, 674669.8125 6581226.96875, "
        "674697.78125 6581231.03125, 674694.34375 6581261.90625, "
        "674690.90625 6581286.34375, 674591.375 6581272.09375, "
        "674597.6875 6581224.40625, 674598.21875 6581224.46875, "
        "674599.28125 6581216.90625, 674640.9375 6581222.875, "
        "674643.53125 6581205.03125, 674638.53125 6581204.34375, "
        "674637.59375 6581211.0625, 674615.6875 6581207.96875, "
        "674615.28125 6581210.78125, 674595.96875 6581208.28125, "
        "674593.84375 6581224.0625, 674582.875 6581222.6875, "
        "674529.40625 6581215.78125, 674529.90625 6581212.15625, "
        "674523.1875 6581206.84375, 674531 6581203.46875, "
        "674531.46875 6581199.875, 674524.4375 6581198.9375, "
        "674533.59375 6581130.875), "
        "(674670.53125 6581260.34375, 674669.125 6581269.59375, "
        "674678.03125 6581270.875, 674682 6581240.875, "
        "674671.1875 6581239.28125, 674670.03125 6581240.53125, "
        "674667.28125 6581259.84375, 674670.53125 6581260.34375), "
        "(674623.375 6581176.9375, 674627.8125 6581144.84375, "
        "674622.0625 6581176.84375, 674623.375 6581176.9375))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    cleaned = cleaning_footprints._remove_meshing_hostile_holes(
        polygon,
        min_hole_area=0.25,
        min_clearance=0.5,
        diagnostics=diagnostics,
    )

    assert len(cleaned.interiors) == 1
    assert diagnostics.get("dropped_meshing_hole_count", 0) == 1
    kept_hole = Polygon(cleaned.interiors[0])
    assert kept_hole.area > 300.0


def test_regularize_meshing_hostile_voids_keeps_case54_style_wedge_pair_separate():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674593.84375 6581224.0625, 674582.875 6581222.6875, "
        "674529.40625 6581215.78125, 674529.90625 6581212.15625, "
        "674523.1875 6581206.84375, 674531 6581203.46875, 674531.46875 6581199.875, "
        "674524.4375 6581198.9375, 674533.59375 6581130.875, 674534.25 6581130.71875, "
        "674570.78125 6581135.71875, 674628 6581143.875, 674708.375 6581155.1875, "
        "674701.34375 6581205.84375, 674696.625 6581205.15625, 674695.5625 6581212.3125, "
        "674672.40625 6581209.0625, 674669.8125 6581226.96875, 674697.78125 6581231.03125, "
        "674694.34375 6581261.90625, 674690.90625 6581286.34375, 674591.375 6581272.09375, "
        "674597.6875 6581224.40625, 674598.21875 6581224.46875, 674599.28125 6581216.90625, "
        "674640.9375 6581222.875, 674643.53125 6581205.03125, 674638.53125 6581204.34375, "
        "674637.59375 6581211.0625, 674615.6875 6581207.96875, 674615.28125 6581210.78125, "
        "674595.96875 6581208.28125, 674593.84375 6581224.0625), "
        "(674623.375 6581176.9375, 674627.8125 6581144.84375, 674622.0625 6581176.84375, "
        "674623.375 6581176.9375), "
        "(674682 6581240.875, 674671.1875 6581239.28125, 674670.03125 6581240.53125, "
        "674667.28125 6581259.84375, 674670.53125 6581260.34375, 674669.125 6581269.59375, "
        "674678.03125 6581270.875, 674682 6581240.875))"
    )
    neighbor = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674546.84375 6581045.5625, 674583.15625 6581051.125, "
        "674610.5625 6581057.34375, 674604.375 6581087.21875, 674601.8125 6581106.28125, "
        "674575.96875 6581102.8125, 674575.28125 6581104.40625, 674571.4375 6581134.90625, "
        "674534.375 6581130.21875, 674538.3125 6581099.09375, 674562.5625 6581102.15625, "
        "674562.78125 6581101.03125, 674538.3125 6581097.75, 674544.03125 6581055.28125, "
        "674544.53125 6581055.34375, 674545.90625 6581049.5625, 674544.90625 6581048.5625, "
        "674546.53125 6581046.90625, 674546.84375 6581045.5625))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_meshing_hostile_voids(
        [polygon, neighbor],
        [[2], [4]],
        min_segment_length=0.5,
        grid=0.03125,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    assert source_map == [[2], [4]]
    assert len(polygons) == 2
    assert diagnostics["coverage_void_regularization_applied"] is True
    assert diagnostics["coverage_void_regularization_hole_cleanup_count"] >= 1
    assert diagnostics["coverage_void_regularization_gap_patch_applied_count"] == 0
    assert diagnostics["coverage_void_regularization_operator_applied"] == {}


def test_regularize_meshing_hostile_voids_rejects_case54_style_long_bridge_wedge():
    left = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674595.96875 6581208.28125, 674593.84375 6581224.0625, "
        "674582.875 6581222.6875, 674529.40625 6581215.78125, "
        "674529.90625 6581212.15625, 674523.1875 6581206.84375, "
        "674531 6581203.46875, 674531.46875 6581199.875, 674524.4375 6581198.9375, "
        "674533.59375 6581130.875, 674534.25 6581130.71875, 674570.78125 6581135.71875, "
        "674628 6581143.875, 674708.375 6581155.1875, 674701.34375 6581205.84375, "
        "674696.625 6581205.15625, 674695.5625 6581212.3125, 674672.40625 6581209.0625, "
        "674669.8125 6581226.96875, 674697.78125 6581231.03125, 674694.34375 6581261.90625, "
        "674690.90625 6581286.34375, 674591.375 6581272.09375, 674597.6875 6581224.40625, "
        "674598.21875 6581224.46875, 674599.28125 6581216.90625, 674640.9375 6581222.875, "
        "674643.53125 6581205.03125, 674638.53125 6581204.34375, 674637.59375 6581211.0625, "
        "674615.6875 6581207.96875, 674615.28125 6581210.78125, 674595.96875 6581208.28125), "
        "(674678.03125 6581270.875, 674682 6581240.875, 674671.1875 6581239.28125, "
        "674670.03125 6581240.53125, 674667.28125 6581259.84375, 674670.53125 6581260.34375, "
        "674669.125 6581269.59375, 674678.03125 6581270.875))"
    )
    right = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674546.84375 6581045.5625, 674583.15625 6581051.125, "
        "674610.5625 6581057.34375, 674604.375 6581087.21875, 674601.8125 6581106.28125, "
        "674575.875 6581102.8125, 674571.4375 6581134.90625, 674534.375 6581130.21875, "
        "674538.3125 6581099.09375, 674562.5625 6581102.15625, 674562.78125 6581101.03125, "
        "674538.3125 6581097.75, 674544.03125 6581055.28125, 674544.53125 6581055.34375, "
        "674545.90625 6581049.59375, 674544.96875 6581047.75, 674546.53125 6581047, "
        "674546.84375 6581045.5625))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_meshing_hostile_voids(
        [left, right],
        [[3, 4, 5, 6, 7, 8, 9], [0, 1]],
        min_segment_length=0.5,
        grid=0.03125,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    assert source_map == [[3, 4, 5, 6, 7, 8, 9], [0, 1]]
    assert len(polygons) == 2
    assert diagnostics["coverage_void_regularization_gap_patch_applied_count"] == 0
    assert diagnostics["coverage_void_regularization_operator_applied"] == {}
    assert polygons[0].bounds[1] == pytest.approx(6581130.71875)


def test_regularize_meshing_hostile_voids_rejects_case15_style_notch_wedge():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((675045.15625 6579425.21875, 675047.65625 6579416.90625, "
        "675035.21875 6579413.03125, 675029 6579434.3125, 675056.59375 6579442.34375, "
        "675055.09375 6579448.3125, 675055.71875 6579449.375, 675063 6579451.34375, "
        "675064.3125 6579446.71875, 675068.34375 6579448.5, 675092.21875 6579452.46875, "
        "675085.84375 6579483.90625, 675061.09375 6579479.1875, 675043.65625 6579474.59375, "
        "675045.0625 6579469.25, 675049.5 6579470.4375, 675050.5 6579466.5, "
        "675063.34375 6579469.96875, 675064.1875 6579466, 675077.15625 6579468.28125, "
        "675079.25 6579457.1875, 675066.875 6579454.9375, 675065.78125 6579459.34375, "
        "675053.125 6579456.03125, 675054.75 6579449.625, 675054.875 6579449.09375, "
        "675033.09375 6579443.15625, 675030.375 6579452.40625, 675020.8125 6579449.96875, "
        "675023.78125 6579440.5, 675014.40625 6579437.84375, 675013.71875 6579440.375, "
        "675015.625 6579440.90625, 675014.8125 6579443.78125, 675012.1875 6579445.53125, "
        "675011.4375 6579448.03125, 675005.53125 6579446.46875, 675008.75 6579435.21875, "
        "675010.71875 6579429.0625, 675014.5625 6579430.15625, 675025.15625 6579395.46875, "
        "675074.34375 6579410.1875, 675138.90625 6579429.5625, 675134.46875 6579443.9375, "
        "675159.1875 6579451.625, 675155.59375 6579463.75, 675124.90625 6579454.21875, "
        "675128.6875 6579442.125, 675098.84375 6579432.875, 675060.1875 6579420.8125, "
        "675057.65625 6579428.96875, 675045.15625 6579425.21875))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_meshing_hostile_voids(
        [polygon],
        [[114, 116, 117, 118, 119, 121, 122, 125, 128, 131, 186]],
        min_segment_length=0.5,
        grid=0.03125,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    assert source_map == [[114, 116, 117, 118, 119, 121, 122, 125, 128, 131, 186]]
    assert len(polygons) == 1
    assert polygons[0].equals_exact(polygon, tolerance=0.0)
    assert diagnostics["coverage_void_regularization_notch_simplify_count"] == 0
    assert diagnostics["coverage_void_regularization_operator_applied"] == {}


def test_micro_detour_chain_simplify_removes_case54_style_boundary_wedge():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674524.4375 6581198.9375, 674533.59375 6581130.875, "
        "674534.15625 6581130.75, 674534.40625 6581130.0625, 674538.3125 6581099.09375, "
        "674562.5625 6581102.15625, 674562.78125 6581101.03125, 674538.3125 6581097.75, "
        "674544.03125 6581055.28125, 674545.90625 6581049.5625, 674544.90625 6581048.5625, "
        "674546.84375 6581045.5625, 674583.15625 6581051.125, 674610.5625 6581057.34375, "
        "674604.375 6581087.21875, 674601.8125 6581106.28125, 674575.96875 6581102.8125, "
        "674571.3125 6581135.78125, 674628 6581143.875, 674708.375 6581155.1875, "
        "674701.34375 6581205.84375, 674696.625 6581205.15625, 674695.5625 6581212.3125, "
        "674672.40625 6581209.0625, 674669.8125 6581226.96875, 674697.78125 6581231.03125, "
        "674694.34375 6581261.90625, 674690.90625 6581286.34375, 674591.375 6581272.09375, "
        "674597.6875 6581224.40625, 674598.21875 6581224.46875, 674599.28125 6581216.90625, "
        "674640.9375 6581222.875, 674643.53125 6581205.03125, 674638.53125 6581204.34375, "
        "674637.59375 6581211.0625, 674615.6875 6581207.96875, 674615.28125 6581210.78125, "
        "674595.96875 6581208.28125, 674593.84375 6581224.0625, 674582.875 6581222.6875, "
        "674529.40625 6581215.78125, 674529.90625 6581212.15625, 674523.1875 6581206.84375, "
        "674531 6581203.46875, 674531.46875 6581199.875, 674524.4375 6581198.9375), "
        "(674671.1875 6581239.28125, 674670.03125 6581240.53125, 674667.28125 6581259.84375, "
        "674670.53125 6581260.34375, 674669.125 6581269.59375, 674678.03125 6581270.875, "
        "674682 6581240.875, 674671.1875 6581239.28125))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._try_polygon_micro_detour_chain_simplify(
        polygon,
        target_scale=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        candidate.polygon,
        target_scale=0.5,
    )
    assert len(candidate.polygon.interiors) == 1
    assert signature_after.clearance is not None
    assert signature_before.clearance is not None
    assert signature_after.clearance >= signature_before.clearance
    assert signature_after.vertex_count < signature_before.vertex_count


def test_regularize_final_polygon_shapes_removes_case54_visible_wedge():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674595.96875 6581208.28125, 674593.84375 6581224.0625, "
        "674582.875 6581222.6875, 674529.40625 6581215.78125, "
        "674529.90625 6581212.15625, 674523.1875 6581206.84375, "
        "674531 6581203.46875, 674531.46875 6581199.875, 674524.4375 6581198.9375, "
        "674533.59375 6581130.875, 674534.25 6581130.71875, 674570.78125 6581135.71875, "
        "674628 6581143.875, 674708.375 6581155.1875, 674701.34375 6581205.84375, "
        "674696.625 6581205.15625, 674695.5625 6581212.3125, 674672.40625 6581209.0625, "
        "674669.8125 6581226.96875, 674697.78125 6581231.03125, 674694.34375 6581261.90625, "
        "674690.90625 6581286.34375, 674591.375 6581272.09375, 674597.6875 6581224.40625, "
        "674598.21875 6581224.46875, 674599.28125 6581216.90625, 674640.9375 6581222.875, "
        "674643.53125 6581205.03125, 674638.53125 6581204.34375, 674637.59375 6581211.0625, "
        "674615.6875 6581207.96875, 674615.28125 6581210.78125, 674595.96875 6581208.28125), "
        "(674678.03125 6581270.875, 674682 6581240.875, 674671.1875 6581239.28125, "
        "674670.03125 6581240.53125, 674667.28125 6581259.84375, 674670.53125 6581260.34375, "
        "674669.125 6581269.59375, 674678.03125 6581270.875))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_final_polygon_shapes(
        [polygon],
        [[0]],
        min_segment_length=0.5,
        grid=0.03125,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    assert source_map == [[0]]
    assert len(polygons) == 1
    assert not polygons[0].equals_exact(polygon, tolerance=0.0)
    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    difference_metrics = cleaning_footprints._difference_area_metrics(
        polygon,
        polygons[0],
    )
    assert difference_metrics["reference_minus_candidate_area"] == pytest.approx(
        0.0,
        abs=1.0e-9,
    )
    assert difference_metrics["candidate_minus_reference_area"] > 0.0
    assert len(list(polygons[0].exterior.coords)) < len(list(polygon.exterior.coords))
    assert signature_after.clearance is not None
    assert signature_before.clearance is not None
    assert signature_after.clearance + 1.0e-9 >= signature_before.clearance
    assert diagnostics["final_shape_regularization_applied"] is True
    assert any(
        operator.startswith(
            (
                "final_shape_local_simplify_",
                "final_shape_micro_detour_chain",
                "final_shape_same_turn_short_walk",
                "final_shape_bevel_corner_collapse",
            )
        )
        for operator in diagnostics["final_shape_regularization_operator_applied"]
    )


def test_regularize_final_polygon_shapes_removes_case15_visible_wedge():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((675045.15625 6579425.21875, 675047.65625 6579416.90625, "
        "675035.21875 6579413.03125, 675029 6579434.3125, 675056.59375 6579442.34375, "
        "675055.09375 6579448.3125, 675055.71875 6579449.375, 675063 6579451.34375, "
        "675064.3125 6579446.71875, 675068.34375 6579448.5, 675092.21875 6579452.46875, "
        "675085.84375 6579483.90625, 675061.09375 6579479.1875, 675043.65625 6579474.59375, "
        "675045.0625 6579469.25, 675049.5 6579470.4375, 675050.5 6579466.5, "
        "675063.34375 6579469.96875, 675064.1875 6579466, 675077.15625 6579468.28125, "
        "675079.25 6579457.1875, 675066.875 6579454.9375, 675065.78125 6579459.34375, "
        "675053.125 6579456.03125, 675054.75 6579449.625, 675054.875 6579449.09375, "
        "675033.09375 6579443.15625, 675030.375 6579452.40625, 675020.8125 6579449.96875, "
        "675023.78125 6579440.5, 675014.40625 6579437.84375, 675013.71875 6579440.375, "
        "675015.625 6579440.90625, 675014.8125 6579443.78125, 675012.1875 6579445.53125, "
        "675011.4375 6579448.03125, 675005.53125 6579446.46875, 675008.75 6579435.21875, "
        "675010.71875 6579429.0625, 675014.5625 6579430.15625, 675025.15625 6579395.46875, "
        "675074.34375 6579410.1875, 675138.90625 6579429.5625, 675134.46875 6579443.9375, "
        "675159.1875 6579451.625, 675155.59375 6579463.75, 675124.90625 6579454.21875, "
        "675128.6875 6579442.125, 675098.84375 6579432.875, 675060.1875 6579420.8125, "
        "675057.65625 6579428.96875, 675045.15625 6579425.21875))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_final_polygon_shapes(
        [polygon],
        [[0]],
        min_segment_length=0.5,
        grid=0.03125,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    assert source_map == [[0]]
    assert len(polygons) == 1
    assert not polygons[0].equals_exact(polygon, tolerance=0.0)
    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    difference_metrics = cleaning_footprints._difference_area_metrics(
        polygon,
        polygons[0],
    )
    assert difference_metrics["reference_minus_candidate_area"] == pytest.approx(
        0.0,
        abs=1.0e-9,
    )
    assert difference_metrics["candidate_minus_reference_area"] > 0.0
    assert len(list(polygons[0].exterior.coords)) < len(list(polygon.exterior.coords))
    assert signature_after.clearance is not None
    assert signature_before.clearance is not None
    assert signature_after.clearance + 1.0e-9 >= signature_before.clearance
    assert diagnostics["final_shape_regularization_applied"] is True
    assert any(
        operator.startswith(
            (
                "final_shape_local_simplify_",
                "final_shape_micro_detour_chain",
                "final_shape_same_turn_short_walk",
                "final_shape_bevel_corner_collapse",
            )
        )
        for operator in diagnostics["final_shape_regularization_operator_applied"]
    )


def test_regularize_final_polygon_shapes_fills_case15_full_building_wedge():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((675338.03125 6579481.15625, 675330.6875 6579478.96875, "
        "675332.3125 6579473.1875, 675341.0625 6579475.5625, 675341.375 6579473.59375, "
        "675343.34375 6579473.96875, 675343.03125 6579475.78125, 675345.25 6579476.3125, "
        "675346.3125 6579471.96875, 675344.5625 6579471.4375, 675346.46875 6579462.46875, "
        "675354.5625 6579464.8125, 675352 6579476.90625, 675352.4375 6579477.53125, "
        "675355.78125 6579478.25, 675357.90625 6579465.8125, 675373.9375 6579470.53125, "
        "675373 6579474.40625, 675368.3125 6579473.34375, 675366.875 6579482.375, "
        "675359.28125 6579481.46875, 675358.4375 6579486.59375, 675350.53125 6579484.5, "
        "675351.4375 6579477.34375, 675347.96875 6579476.59375, 675346.71875 6579483.34375, "
        "675338.5625 6579481.0625, 675338.71875 6579477.25, 675338.03125 6579481.15625))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_final_polygon_shapes(
        [polygon],
        [[0]],
        min_segment_length=0.5,
        grid=0.03125,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    assert source_map == [[0]]
    assert len(polygons) == 1
    assert not polygons[0].equals_exact(polygon, tolerance=0.0)
    assert polygons[0].area > polygon.area
    assert len(list(polygons[0].exterior.coords)) < len(list(polygon.exterior.coords))
    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    assert signature_after.min_edge_length > signature_before.min_edge_length
    assert diagnostics["final_shape_regularization_applied"] is True
    assert any(
        operator.startswith(
            (
                "final_shape_local_simplify_",
                "final_shape_micro_detour_chain",
                "final_shape_same_turn_short_walk",
            )
        )
        for operator in diagnostics["final_shape_regularization_operator_applied"]
    )


def test_regularize_final_polygon_shapes_fills_case54_full_building_wedge():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674546.84375 6581045.5625, 674583.15625 6581051.125, "
        "674610.5625 6581057.34375, 674604.375 6581087.21875, 674601.8125 6581106.28125, "
        "674575.875 6581102.8125, 674571.4375 6581134.90625, 674534.375 6581130.21875, "
        "674538.3125 6581099.09375, 674562.5625 6581102.15625, 674562.78125 6581101.03125, "
        "674538.3125 6581097.75, 674544.03125 6581055.28125, 674544.53125 6581055.34375, "
        "674545.90625 6581049.5625, 674544.90625 6581048.5625, 674546.53125 6581046.90625, "
        "674546.84375 6581045.5625))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_final_polygon_shapes(
        [polygon],
        [[0]],
        min_segment_length=0.5,
        grid=0.03125,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    assert source_map == [[0]]
    assert len(polygons) == 1
    assert not polygons[0].equals_exact(polygon, tolerance=0.0)
    assert polygons[0].area > polygon.area
    assert len(list(polygons[0].exterior.coords)) < len(list(polygon.exterior.coords))
    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    assert signature_after.min_edge_length > signature_before.min_edge_length
    assert diagnostics["final_shape_regularization_applied"] is True
    assert any(
        operator.startswith(
            (
                "final_shape_local_simplify_",
                "final_shape_micro_detour_chain",
                "final_shape_same_turn_short_walk",
            )
        )
        for operator in diagnostics["final_shape_regularization_operator_applied"]
    )


def test_same_turn_short_walk_collapse_removes_case54_residual_boundary_wedge():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674602.09375 6581106.625, 674562.9375 6581101.375, "
        "674562.875 6581101.875, 674575.71875 6581103.5, "
        "674571.71875 6581135.25, 674570.78125 6581135.71875, "
        "674628 6581143.875, 674708.375 6581155.1875, "
        "674701.34375 6581205.84375, 674696.625 6581205.15625, "
        "674695.5625 6581212.3125, 674672.40625 6581209.0625, "
        "674669.8125 6581226.96875, 674697.78125 6581231.03125, "
        "674694.34375 6581261.90625, 674690.90625 6581286.34375, "
        "674591.375 6581272.09375, 674597.6875 6581224.40625, "
        "674598.21875 6581224.46875, 674599.28125 6581216.90625, "
        "674640.9375 6581222.875, 674643.53125 6581205.03125, "
        "674638.53125 6581204.34375, 674637.59375 6581211.0625, "
        "674615.6875 6581207.96875, 674615.28125 6581210.78125, "
        "674595.96875 6581208.28125, 674593.84375 6581224.0625, "
        "674582.875 6581222.6875, 674529.40625 6581215.78125, "
        "674529.90625 6581212.15625, 674523.1875 6581206.84375, "
        "674531 6581203.46875, 674531.46875 6581199.875, "
        "674524.4375 6581198.9375, 674533.59375 6581130.875, "
        "674534.25 6581130.71875, 674536.6875 6581131.0625, "
        "674539.34375 6581131.15625, 674534.03125 6581130.5, "
        "674538.03125 6581098.75, 674537.96875 6581098.03125, "
        "674543.75 6581054.9375, 674544.28125 6581055, "
        "674545.5625 6581049.65625, 674544.46875 6581048.5625, "
        "674546.25 6581046.75, 674546.59375 6581045.21875, "
        "674559.59375 6581047, 674583.21875 6581050.8125, "
        "674601.03125 6581055.03125, 674610.9375 6581057.09375, "
        "674604.6875 6581087.28125, 674602.09375 6581106.625), "
        "(674670.03125 6581240.53125, 674667.28125 6581259.84375, "
        "674670.53125 6581260.34375, 674669.125 6581269.59375, "
        "674678.03125 6581270.875, 674682 6581240.875, "
        "674671.1875 6581239.28125, 674670.03125 6581240.53125))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._try_polygon_same_turn_short_walk_collapse(
        polygon,
        target_scale=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        candidate.polygon,
        target_scale=0.5,
    )
    assert signature_after.clearance + 1.0e-9 >= signature_before.clearance
    assert signature_after.min_edge_length + 1.0e-9 >= signature_before.min_edge_length
    assert signature_after.vertex_count < signature_before.vertex_count
    diff = cleaning_footprints._difference_area_metrics(polygon, candidate.polygon)
    assert diff["reference_minus_candidate_area"] == pytest.approx(0.0, abs=0.01)
    assert diff["candidate_minus_reference_area"] < 16.0


def test_fill_chain_collapse_removes_case54_residual_upper_notch():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674583.21875 6581050.8125, 674601.03125 6581055.03125, "
        "674610.9375 6581057.09375, 674604.6875 6581087.28125, "
        "674602.09375 6581106.625, 674575.71875 6581103.5, "
        "674571.71875 6581135.25, 674570.78125 6581135.71875, "
        "674628 6581143.875, 674708.375 6581155.1875, "
        "674701.34375 6581205.84375, 674696.625 6581205.15625, "
        "674695.5625 6581212.3125, 674672.40625 6581209.0625, "
        "674669.8125 6581226.96875, 674697.78125 6581231.03125, "
        "674694.34375 6581261.90625, 674690.90625 6581286.34375, "
        "674591.375 6581272.09375, 674597.6875 6581224.40625, "
        "674598.21875 6581224.46875, 674599.28125 6581216.90625, "
        "674640.9375 6581222.875, 674637.59375 6581211.0625, "
        "674615.6875 6581207.96875, 674615.28125 6581210.78125, "
        "674595.96875 6581208.28125, 674593.84375 6581224.0625, "
        "674582.875 6581222.6875, 674529.40625 6581215.78125, "
        "674529.90625 6581212.15625, 674523.1875 6581206.84375, "
        "674524.4375 6581198.9375, 674533.59375 6581130.875, "
        "674534.25 6581130.71875, 674536.6875 6581131.0625, "
        "674539.34375 6581131.15625, 674534.03125 6581130.5, "
        "674538.03125 6581098.75, 674537.96875 6581098.03125, "
        "674543.75 6581054.9375, 674544.46875 6581048.5625, "
        "674546.25 6581046.75, 674546.59375 6581045.21875, "
        "674559.59375 6581047, 674583.21875 6581050.8125), "
        "(674669.125 6581269.59375, 674678.03125 6581270.875, "
        "674682 6581240.875, 674671.1875 6581239.28125, "
        "674670.03125 6581240.53125, 674667.28125 6581259.84375, "
        "674670.53125 6581260.34375, 674669.125 6581269.59375))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._try_polygon_fill_chain_collapse(
        polygon,
        target_scale=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        candidate.polygon,
        target_scale=0.5,
    )
    diff = cleaning_footprints._difference_area_metrics(polygon, candidate.polygon)
    assert signature_after.clearance > signature_before.clearance
    assert signature_after.clearance_deficit == pytest.approx(0.0)
    assert signature_after.min_edge_length + 1.0e-9 >= signature_before.min_edge_length
    assert signature_after.vertex_count < signature_before.vertex_count
    assert diff["reference_minus_candidate_area"] == pytest.approx(0.0, abs=0.01)
    assert diff["candidate_minus_reference_area"] < 16.0


def test_bevel_corner_collapse_removes_case54_tiny_remnant_corner():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674559.59375 6581047, 674583.21875 6581050.8125, "
        "674601.03125 6581055.03125, 674610.9375 6581057.09375, "
        "674604.6875 6581087.28125, 674602.09375 6581106.625, "
        "674575.71875 6581103.5, 674571.71875 6581135.25, "
        "674570.78125 6581135.71875, 674628 6581143.875, "
        "674708.375 6581155.1875, 674701.34375 6581205.84375, "
        "674696.625 6581205.15625, 674695.5625 6581212.3125, "
        "674672.40625 6581209.0625, 674669.8125 6581226.96875, "
        "674697.78125 6581231.03125, 674694.34375 6581261.90625, "
        "674690.90625 6581286.34375, 674591.375 6581272.09375, "
        "674597.6875 6581224.40625, 674598.21875 6581224.46875, "
        "674599.28125 6581216.90625, 674640.9375 6581222.875, "
        "674637.59375 6581211.0625, 674615.6875 6581207.96875, "
        "674615.28125 6581210.78125, 674595.96875 6581208.28125, "
        "674593.84375 6581224.0625, 674582.875 6581222.6875, "
        "674529.40625 6581215.78125, 674529.90625 6581212.15625, "
        "674523.1875 6581206.84375, 674524.4375 6581198.9375, "
        "674533.59375 6581130.875, 674538.03125 6581098.75, "
        "674537.96875 6581098.03125, 674543.75 6581054.9375, "
        "674544.46875 6581048.5625, 674546.25 6581046.75, "
        "674546.59375 6581045.21875, 674559.59375 6581047), "
        "(674669.125 6581269.59375, 674678.03125 6581270.875, "
        "674682 6581240.875, 674671.1875 6581239.28125, "
        "674670.03125 6581240.53125, 674667.28125 6581259.84375, "
        "674670.53125 6581260.34375, 674669.125 6581269.59375))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._try_polygon_bevel_corner_collapse(
        polygon,
        target_scale=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
    )

    assert candidate is not None
    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        candidate.polygon,
        target_scale=0.5,
    )
    diff = cleaning_footprints._difference_area_metrics(polygon, candidate.polygon)
    assert signature_after.clearance + 1.0e-9 >= signature_before.clearance
    assert signature_after.min_edge_length + 1.0e-9 >= signature_before.min_edge_length
    assert signature_after.vertex_count < signature_before.vertex_count
    assert diff["reference_minus_candidate_area"] == pytest.approx(0.0, abs=1e-6)
    assert diff["candidate_minus_reference_area"] < 1.0


def test_regularize_final_polygon_shapes_fills_case54_contact_cluster_wedges():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674562.875 6581101.875, 674575.71875 6581103.5, "
        "674571.71875 6581135.25, 674543.5625 6581131.6875, "
        "674540.375 6581131.5625, 674570.78125 6581135.71875, "
        "674628 6581143.875, 674708.375 6581155.1875, "
        "674701.34375 6581205.84375, 674696.625 6581205.15625, "
        "674695.5625 6581212.3125, 674672.40625 6581209.0625, "
        "674669.8125 6581226.96875, 674697.78125 6581231.03125, "
        "674694.34375 6581261.90625, 674690.90625 6581286.34375, "
        "674591.375 6581272.09375, 674597.6875 6581224.40625, "
        "674598.21875 6581224.46875, 674599.28125 6581216.90625, "
        "674640.9375 6581222.875, 674643.53125 6581205.03125, "
        "674638.53125 6581204.34375, 674637.59375 6581211.0625, "
        "674615.6875 6581207.96875, 674615.28125 6581210.78125, "
        "674595.96875 6581208.28125, 674593.84375 6581224.0625, "
        "674582.875 6581222.6875, 674529.40625 6581215.78125, "
        "674529.90625 6581212.15625, 674523.1875 6581206.84375, "
        "674531 6581203.46875, 674531.46875 6581199.875, "
        "674524.4375 6581198.9375, 674533.59375 6581130.875, "
        "674534.25 6581130.71875, 674536.6875 6581131.0625, "
        "674539.34375 6581131.15625, 674534.03125 6581130.5, "
        "674538.03125 6581098.75, 674562.3125 6581101.8125, "
        "674562.40625 6581101.3125, 674537.96875 6581098.03125, "
        "674543.75 6581054.9375, 674544.28125 6581055, "
        "674545.5625 6581049.65625, 674544.46875 6581048.5625, "
        "674546.25 6581046.75, 674546.59375 6581045.21875, "
        "674559.59375 6581047, 674583.21875 6581050.8125, "
        "674601.03125 6581055.03125, 674610.9375 6581057.09375, "
        "674604.6875 6581087.28125, 674602.09375 6581106.625, "
        "674562.9375 6581101.375, 674562.875 6581101.875), "
        "(674670.53125 6581260.34375, 674669.125 6581269.59375, "
        "674678.03125 6581270.875, 674682 6581240.875, "
        "674671.1875 6581239.28125, 674670.03125 6581240.53125, "
        "674667.28125 6581259.84375, 674670.53125 6581260.34375), "
        "(674623.375 6581176.9375, 674627.8125 6581144.84375, "
        "674622.0625 6581176.84375, 674623.375 6581176.9375))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_final_polygon_shapes(
        [polygon],
        [[0]],
        min_segment_length=0.5,
        grid=0.03125,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )

    assert source_map == [[0]]
    assert len(polygons) == 1
    assert not polygons[0].equals_exact(polygon, tolerance=0.0)
    assert polygons[0].area > polygon.area
    assert len(list(polygons[0].exterior.coords)) < len(list(polygon.exterior.coords))
    assert signature_after.clearance + 1.0e-9 >= signature_before.clearance
    assert signature_after.min_edge_length > signature_before.min_edge_length
    assert diagnostics["final_shape_regularization_applied"] is True
    assert diagnostics["final_shape_regularization_operator_applied"].get(
        "final_shape_micro_detour_chain", 0
    ) >= 1
    assert diagnostics["final_shape_regularization_operator_applied"].get(
        "final_shape_same_turn_short_walk", 0
    ) >= 1


def test_regularize_final_polygon_shapes_closes_case54_residual_same_turn_wedges():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674610.9375 6581057.09375, 674604.6875 6581087.28125, "
        "674602.09375 6581106.625, 674575.71875 6581103.5, "
        "674571.71875 6581135.25, 674570.78125 6581135.71875, "
        "674628 6581143.875, 674708.375 6581155.1875, "
        "674701.34375 6581205.84375, 674696.625 6581205.15625, "
        "674695.5625 6581212.3125, 674672.40625 6581209.0625, "
        "674669.8125 6581226.96875, 674697.78125 6581231.03125, "
        "674694.34375 6581261.90625, 674690.90625 6581286.34375, "
        "674591.375 6581272.09375, 674597.6875 6581224.40625, "
        "674598.21875 6581224.46875, 674599.28125 6581216.90625, "
        "674640.9375 6581222.875, 674643.53125 6581205.03125, "
        "674638.53125 6581204.34375, 674637.59375 6581211.0625, "
        "674615.6875 6581207.96875, 674615.28125 6581210.78125, "
        "674595.96875 6581208.28125, 674593.84375 6581224.0625, "
        "674582.875 6581222.6875, 674529.40625 6581215.78125, "
        "674529.90625 6581212.15625, 674523.1875 6581206.84375, "
        "674531 6581203.46875, 674531.46875 6581199.875, "
        "674524.4375 6581198.9375, 674533.59375 6581130.875, "
        "674534.25 6581130.71875, 674536.6875 6581131.0625, "
        "674539.34375 6581131.15625, 674534.03125 6581130.5, "
        "674538.03125 6581098.75, 674537.96875 6581098.03125, "
        "674543.75 6581054.9375, 674544.46875 6581048.5625, "
        "674546.25 6581046.75, 674546.59375 6581045.21875, "
        "674559.59375 6581047, 674583.21875 6581050.8125, "
        "674601.03125 6581055.03125, 674610.9375 6581057.09375), "
        "(674682 6581240.875, 674671.1875 6581239.28125, "
        "674670.03125 6581240.53125, 674667.28125 6581259.84375, "
        "674670.53125 6581260.34375, 674669.125 6581269.59375, "
        "674678.03125 6581270.875, 674682 6581240.875))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_final_polygon_shapes(
        [polygon],
        [[0]],
        min_segment_length=0.5,
        grid=0.03125,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    diff = cleaning_footprints._difference_area_metrics(polygon, polygons[0])

    assert source_map == [[0]]
    assert len(polygons) == 1
    assert not polygons[0].equals_exact(polygon, tolerance=0.0)
    assert diff["reference_minus_candidate_area"] == pytest.approx(0.0, abs=0.01)
    assert diff["candidate_minus_reference_area"] > 100.0
    assert signature_after.clearance > signature_before.clearance
    assert signature_after.clearance_deficit == pytest.approx(0.0)
    assert signature_after.min_edge_length + 1.0e-9 >= signature_before.min_edge_length
    assert signature_after.vertex_count < signature_before.vertex_count
    assert diagnostics["final_shape_regularization_operator_applied"].get(
        "final_shape_same_turn_short_walk", 0
    ) >= 2
    assert diagnostics["final_shape_regularization_operator_applied"].get(
        "final_shape_fill_chain_collapse", 0
    ) >= 1
    assert diagnostics["final_shape_regularization_operator_applied"].get(
        "final_shape_bevel_corner_collapse", 0
    ) >= 1


def test_regularize_coverage_for_meshing_repairs_case54_residual_self_clearance():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674533.59375 6581130.875, 674534.75 6581130.59375, "
        "674534.03125 6581130.5, 674538.03125 6581098.75, "
        "674562.3125 6581101.8125, 674562.40625 6581101.3125, "
        "674537.96875 6581098.03125, 674543.75 6581054.9375, 674544.28125 6581055, "
        "674545.5625 6581049.65625, 674544.46875 6581048.5625, "
        "674546.25 6581046.75, 674546.59375 6581045.21875, "
        "674559.59375 6581047, 674583.21875 6581050.8125, "
        "674601.03125 6581055.03125, 674610.9375 6581057.09375, "
        "674604.6875 6581087.28125, 674602.09375 6581106.625, "
        "674562.9375 6581101.375, 674562.875 6581101.875, "
        "674575.71875 6581103.5, 674571.71875 6581135.25, "
        "674555.84375 6581133.25, 674546.65625 6581132.40625, "
        "674570.78125 6581135.71875, 674628 6581143.875, "
        "674708.375 6581155.1875, 674701.34375 6581205.84375, "
        "674696.625 6581205.15625, 674695.5625 6581212.3125, "
        "674672.40625 6581209.0625, 674669.8125 6581226.96875, "
        "674697.78125 6581231.03125, 674694.34375 6581261.90625, "
        "674690.90625 6581286.34375, 674591.375 6581272.09375, "
        "674597.6875 6581224.40625, 674598.21875 6581224.46875, "
        "674599.28125 6581216.90625, 674640.9375 6581222.875, "
        "674643.53125 6581205.03125, 674638.53125 6581204.34375, "
        "674637.59375 6581211.0625, 674615.6875 6581207.96875, "
        "674615.28125 6581210.78125, 674595.96875 6581208.28125, "
        "674593.84375 6581224.0625, 674582.875 6581222.6875, "
        "674529.40625 6581215.78125, 674529.90625 6581212.15625, "
        "674523.1875 6581206.84375, 674531 6581203.46875, "
        "674531.46875 6581199.875, 674524.4375 6581198.9375, "
        "674533.59375 6581130.875), "
        "(674670.53125 6581260.34375, 674669.125 6581269.59375, "
        "674678.03125 6581270.875, 674682 6581240.875, "
        "674671.1875 6581239.28125, 674670.03125 6581240.53125, "
        "674667.28125 6581259.84375, 674670.53125 6581260.34375), "
        "(674623.375 6581176.9375, 674627.8125 6581144.84375, "
        "674622.0625 6581176.84375, 674623.375 6581176.9375))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_coverage_for_meshing(
        [polygon],
        [[0]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=15.0,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
    assert source_map == [[0]]
    assert diagnostics["coverage_meshing_regularization_selected_branch"] == (
        "residual_self_clearance_repair"
    )
    assert signature.min_clearance is not None
    assert signature.min_clearance > 0.2609715353086986


def test_regularize_coverage_for_meshing_prefers_lower_drift_direct_clearance_repair(
    monkeypatch,
):
    original = [box(0.0, 0.0, 4.0, 4.0)]
    local_repair = [box(0.0, 0.0, 4.1, 4.0)]
    direct_clearance = [box(0.0, 0.0, 4.2, 4.0)]
    local_then_clearance = [box(0.0, 0.0, 4.6, 4.0)]
    source_map = [[0]]

    original_key = cleaning_footprints._polygon_sequence_key(original)
    local_key = cleaning_footprints._polygon_sequence_key(local_repair)
    direct_key = cleaning_footprints._polygon_sequence_key(direct_clearance)
    local_then_clearance_key = cleaning_footprints._polygon_sequence_key(
        local_then_clearance
    )

    signature_map = {
        original_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.26,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=4,
            ring_contact_count=0,
        ),
        local_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.41,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=4,
            ring_contact_count=0,
        ),
        direct_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.51,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=4,
            ring_contact_count=0,
        ),
        local_then_clearance_key: cleaning_footprints._CoverageDefectSignature(
            min_clearance=0.51,
            pair_issue_count=0,
            point_touch_count=0,
            close_pair_count=0,
            min_pair_clearance=None,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=4,
            ring_contact_count=0,
        ),
    }
    difference_map = {
        (original_key, original_key): cleaning_footprints._zero_difference_metrics(),
        (original_key, local_key): {
            "reference_minus_candidate_area": 0.1,
            "candidate_minus_reference_area": 0.4,
            "symmetric_difference_area": 0.5,
            "union_area_delta": 0.3,
        },
        (original_key, direct_key): {
            "reference_minus_candidate_area": 0.2,
            "candidate_minus_reference_area": 0.6,
            "symmetric_difference_area": 0.8,
            "union_area_delta": 0.4,
        },
        (original_key, local_then_clearance_key): {
            "reference_minus_candidate_area": 0.6,
            "candidate_minus_reference_area": 1.8,
            "symmetric_difference_area": 2.4,
            "union_area_delta": 1.2,
        },
    }

    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_short_edge_graphically",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_globally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_locally",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_simplify_coverage_residual_scale_polygons",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_polygon_repairs",
        lambda polygons, source_map, **kwargs: (
            (local_repair, source_map)
            if kwargs.get("stage_prefix") == "meshing_contract_repair"
            else (polygons, source_map)
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_regularize_low_clearance_polygons",
        lambda polygons, source_map, **kwargs: (
            (direct_clearance, source_map)
            if cleaning_footprints._polygon_sequence_key(polygons) == original_key
            else (
                (local_then_clearance, source_map)
                if cleaning_footprints._polygon_sequence_key(polygons) == local_key
                else (polygons, source_map)
            )
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_coverage_defect_signature",
        lambda cache, polygons, target_scale: signature_map[
            cleaning_footprints._polygon_sequence_key(polygons)
        ],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_cached_difference_area_metrics",
        lambda cache, reference_polygons, candidate_polygons: difference_map[
            (
                cleaning_footprints._polygon_sequence_key(reference_polygons),
                cleaning_footprints._polygon_sequence_key(candidate_polygons),
            )
        ],
    )

    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, repaired_sources = cleaning_footprints._regularize_coverage_for_meshing(
        original,
        source_map,
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    assert cleaning_footprints._polygon_sequence_key(polygons) == direct_key
    assert repaired_sources == source_map
    assert diagnostics["coverage_meshing_regularization_selected_branch"] == (
        "residual_self_clearance_repair"
    )
    assert diagnostics["coverage_meshing_regularization_operator_applied"] == {
        "meshing_contract_direct_clearance_regularization": 1
    }


def test_contact_resolution_candidate_score_prefers_lower_extra_until_clearance_contract_is_met():
    lower_extra = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.1902816349840361,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=0,
        min_edge_length=0.5038911092686593,
        vertex_count=64,
        ring_contact_count=0,
    )
    higher_extra = cleaning_footprints._CoverageDefectSignature(
        min_clearance=0.2609715353086986,
        pair_issue_count=0,
        point_touch_count=0,
        close_pair_count=0,
        min_pair_clearance=None,
        short_edge_count=0,
        min_edge_length=0.5038911092686593,
        vertex_count=62,
        ring_contact_count=0,
    )

    lower_extra_diff = {
        "reference_minus_candidate_area": 3.1576304008481744,
        "candidate_minus_reference_area": 5.087948325380748,
        "symmetric_difference_area": 8.245578726228922,
        "union_area_delta": 1.9303179245325737,
    }
    higher_extra_diff = {
        "reference_minus_candidate_area": 2.9101616926682032,
        "candidate_minus_reference_area": 8.921046023450778,
        "symmetric_difference_area": 11.83120771611898,
        "union_area_delta": 6.010884330782574,
    }

    assert cleaning_footprints._contact_resolution_candidate_score(
        lower_extra,
        lower_extra_diff,
        target_scale=0.5,
    ) < cleaning_footprints._contact_resolution_candidate_score(
        higher_extra,
        higher_extra_diff,
        target_scale=0.5,
    )


def test_regularize_low_clearance_polygons_prefers_lower_spike_growth_until_contract_met(
    monkeypatch,
):
    original = box(0.0, 0.0, 10.0, 10.0)
    smooth = box(0.0, 0.0, 10.2, 10.0)
    spiky = Polygon(
        [
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (6.0, 10.0),
            (5.0, 25.0),
            (4.0, 10.0),
            (0.0, 10.0),
        ]
    )

    signature_map = {
        original.wkt: cleaning_footprints._PolygonDefectSignature(
            clearance=0.26,
            clearance_deficit=0.24,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=4,
            ring_contact_count=0,
        ),
        smooth.wkt: cleaning_footprints._PolygonDefectSignature(
            clearance=0.39,
            clearance_deficit=0.11,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=4,
            ring_contact_count=0,
        ),
        spiky.wkt: cleaning_footprints._PolygonDefectSignature(
            clearance=0.41,
            clearance_deficit=0.09,
            short_edge_count=0,
            min_edge_length=1.0,
            vertex_count=7,
            ring_contact_count=0,
        ),
    }
    clearance_map = {
        original.wkt: 0.26,
        smooth.wkt: 0.39,
        spiky.wkt: 0.41,
    }

    monkeypatch.setattr(
        cleaning_footprints.shapely,
        "minimum_clearance",
        lambda polygon: clearance_map[polygon.wkt],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_polygon_defect_signature",
        lambda polygon, target_scale: signature_map[polygon.wkt],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_remove_meshing_hostile_holes",
        lambda polygon, **kwargs: polygon,
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_canonicalize",
        lambda geometry, grid, diagnostics: [geometry]
        if isinstance(geometry, Polygon)
        else [],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_opening",
        lambda polygon, radius, grid, diagnostics: [],
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_try_polygon_self_clearance_connector_fill",
        lambda polygon, **kwargs: cleaning_footprints._RepairCandidate(
            polygon=spiky,
            edit_zone=spiky.buffer(0.0),
            operator="spike_candidate",
        ),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_apply_local_polygon_repairs",
        lambda polygons, source_map, **kwargs: ([smooth], source_map),
    )
    monkeypatch.setattr(
        cleaning_footprints,
        "_try_polygon_local_simplify",
        lambda *args, **kwargs: None,
    )

    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._regularize_low_clearance_polygons(
        [original],
        [[0]],
        min_clearance=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
    )

    assert len(polygons) == 1
    assert source_map == [[0]]
    assert polygons[0].equals_exact(smooth, tolerance=0.0)


def test_regularize_low_clearance_polygons_removes_case54_style_sliver_hole():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674533.59375 6581130.875, 674534.75 6581130.59375, "
        "674534.03125 6581130.5, 674538.03125 6581098.75, "
        "674562.3125 6581101.8125, 674562.40625 6581101.3125, "
        "674537.96875 6581098.03125, 674543.75 6581054.9375, 674544.28125 6581055, "
        "674545.5625 6581049.65625, 674544.46875 6581048.5625, "
        "674546.25 6581046.75, 674546.59375 6581045.21875, "
        "674559.59375 6581047, 674583.21875 6581050.8125, "
        "674601.03125 6581055.03125, 674610.9375 6581057.09375, "
        "674604.6875 6581087.28125, 674602.09375 6581106.625, "
        "674562.9375 6581101.375, 674562.875 6581101.875, "
        "674575.71875 6581103.5, 674571.71875 6581135.25, "
        "674555.84375 6581133.25, 674546.65625 6581132.40625, "
        "674570.78125 6581135.71875, 674628 6581143.875, "
        "674708.375 6581155.1875, 674701.34375 6581205.84375, "
        "674696.625 6581205.15625, 674695.5625 6581212.3125, "
        "674672.40625 6581209.0625, 674669.8125 6581226.96875, "
        "674697.78125 6581231.03125, 674694.34375 6581261.90625, "
        "674690.90625 6581286.34375, 674591.375 6581272.09375, "
        "674597.6875 6581224.40625, 674598.21875 6581224.46875, "
        "674599.28125 6581216.90625, 674640.9375 6581222.875, "
        "674643.53125 6581205.03125, 674638.53125 6581204.34375, "
        "674637.59375 6581211.0625, 674615.6875 6581207.96875, "
        "674615.28125 6581210.78125, 674595.96875 6581208.28125, "
        "674593.84375 6581224.0625, 674582.875 6581222.6875, "
        "674529.40625 6581215.78125, 674529.90625 6581212.15625, "
        "674523.1875 6581206.84375, 674531 6581203.46875, "
        "674531.46875 6581199.875, 674524.4375 6581198.9375, "
        "674533.59375 6581130.875), "
        "(674670.53125 6581260.34375, 674669.125 6581269.59375, "
        "674678.03125 6581270.875, 674682 6581240.875, "
        "674671.1875 6581239.28125, 674670.03125 6581240.53125, "
        "674667.28125 6581259.84375, 674670.53125 6581260.34375), "
        "(674623.375 6581176.9375, 674627.8125 6581144.84375, "
        "674622.0625 6581176.84375, 674623.375 6581176.9375))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, _ = cleaning_footprints._regularize_low_clearance_polygons(
        [polygon],
        [[0]],
        min_clearance=0.5,
        grid=0.03125,
        min_area=15.0,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    assert len(polygons) == 1
    assert len(polygons[0].interiors) == 1
    signature_before = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    signature_after = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    assert signature_after.clearance is not None
    assert signature_before.clearance is not None
    assert signature_after.clearance + 1.0e-9 >= signature_before.clearance
    assert diagnostics["clearance_regularization_applied_count"] == 1
    assert diagnostics.get("dropped_meshing_hole_count", 0) >= 1


def test_regularize_low_clearance_polygons_closes_case54_style_same_ring_slit():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674533.59375 6581130.875, 674534.75 6581130.59375, "
        "674534.03125 6581130.5, 674538.03125 6581098.75, "
        "674562.3125 6581101.8125, 674562.40625 6581101.3125, "
        "674537.96875 6581098.03125, 674543.75 6581054.9375, 674544.28125 6581055, "
        "674545.5625 6581049.65625, 674544.46875 6581048.5625, "
        "674546.25 6581046.75, 674546.59375 6581045.21875, "
        "674559.59375 6581047, 674583.21875 6581050.8125, "
        "674601.03125 6581055.03125, 674610.9375 6581057.09375, "
        "674604.6875 6581087.28125, 674602.09375 6581106.625, "
        "674562.9375 6581101.375, 674562.875 6581101.875, "
        "674575.71875 6581103.5, 674571.71875 6581135.25, "
        "674555.84375 6581133.25, 674546.65625 6581132.40625, "
        "674570.78125 6581135.71875, 674628 6581143.875, "
        "674708.375 6581155.1875, 674701.34375 6581205.84375, "
        "674696.625 6581205.15625, 674695.5625 6581212.3125, "
        "674672.40625 6581209.0625, 674669.8125 6581226.96875, "
        "674697.78125 6581231.03125, 674694.34375 6581261.90625, "
        "674690.90625 6581286.34375, 674591.375 6581272.09375, "
        "674597.6875 6581224.40625, 674598.21875 6581224.46875, "
        "674599.28125 6581216.90625, 674640.9375 6581222.875, "
        "674643.53125 6581205.03125, 674638.53125 6581204.34375, "
        "674637.59375 6581211.0625, 674615.6875 6581207.96875, "
        "674615.28125 6581210.78125, 674595.96875 6581208.28125, "
        "674593.84375 6581224.0625, 674582.875 6581222.6875, "
        "674529.40625 6581215.78125, 674529.90625 6581212.15625, "
        "674523.1875 6581206.84375, 674531 6581203.46875, "
        "674531.46875 6581199.875, 674524.4375 6581198.9375, "
        "674533.59375 6581130.875), "
        "(674670.53125 6581260.34375, 674669.125 6581269.59375, "
        "674678.03125 6581270.875, 674682 6581240.875, "
        "674671.1875 6581239.28125, 674670.03125 6581240.53125, "
        "674667.28125 6581259.84375, 674670.53125 6581260.34375), "
        "(674623.375 6581176.9375, 674627.8125 6581144.84375, "
        "674622.0625 6581176.84375, 674623.375 6581176.9375))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False
    polygon = cleaning_footprints._remove_meshing_hostile_holes(
        polygon,
        min_hole_area=0.25,
        min_clearance=0.5,
        diagnostics=diagnostics,
    )

    polygons, _ = cleaning_footprints._regularize_low_clearance_polygons(
        [polygon],
        [[0]],
        min_clearance=0.5,
        grid=0.03125,
        min_area=15.0,
        min_hole_area=0.25,
        diagnostics=diagnostics,
    )

    signature = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    assert len(polygons) == 1
    assert len(polygons[0].interiors) == 1
    assert diagnostics["clearance_regularization_applied_count"] == 1
    assert signature.clearance is not None
    assert signature.short_edge_count == 0
    assert signature.min_edge_length is not None
    assert signature.min_edge_length >= 0.5 - 1.0e-9


def test_regularize_coverage_ring_contacts_repairs_ring_contacts_without_clearance_deficit():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674539.8125 6579688.65625, 674541.4375 6579683.03125, "
        "674520.59375 6579677.25, 674520.0625 6579679.125, "
        "674514.65625 6579676.90625, 674512.6875 6579674.5625, "
        "674511.5 6579670.46875, 674512.65625 6579666.09375, "
        "674517.9375 6579667.875, 674520.53125 6579666.8125, "
        "674521.9375 6579662.125, 674514.3125 6579659.78125, "
        "674517.03125 6579650.03125, 674542.15625 6579657.34375, "
        "674570.5625 6579665.46875, 674557.84375 6579712.90625, "
        "674534.6875 6579706.59375, 674536.4375 6579700.4375, "
        "674540.375 6579701.65625, 674541.65625 6579697.09375, "
        "674546.4375 6579698.46875, 674545.34375 6579702.3125, "
        "674553.15625 6579704.46875, 674555.25 6579696.8125, "
        "674547.5 6579694.6875, 674548.375 6579691.1875, "
        "674539.8125 6579688.65625), "
        "(674545.40625 6579669, 674541.4375 6579683.03125, "
        "674546.40625 6579684.59375, 674547.4375 6579680.34375, "
        "674550.40625 6579678.5, 674558.71875 6579680.75, "
        "674560.75 6579673.96875, 674556.9375 6579672.78125, "
        "674556.40625 6579674.46875, 674553.9375 6579675.46875, "
        "674548.40625 6579673.9375, 674549.5 6579670.09375, "
        "674545.40625 6579669))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    before_signature = cleaning_footprints._polygon_defect_signature(
        polygon,
        target_scale=0.5,
    )
    assert before_signature.ring_contact_count == 1
    assert before_signature.clearance is not None
    assert before_signature.clearance > 0.5

    candidate = cleaning_footprints._regularize_coverage_ring_contacts(
        [polygon],
        [[0]],
        min_segment_length=0.5,
        grid=0.03125,
    )
    assert candidate is not None
    (polygons, _), stats = candidate

    after_signature = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    assert len(polygons) == 1
    assert after_signature.ring_contact_count == 0
    assert stats["changed_count"] == 1
    assert stats["failed_count"] == 0


def test_iteratively_open_polygon_short_edges_resolves_case62_style_exterior_stub():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((673576.1875 6581592.21875, 673527.65625 6581564.28125, "
        "673527.5625 6581564.0625, 673534.34375 6581552, 673540.625 6581555.59375, "
        "673541.8125 6581553.625, 673539.4375 6581552.21875, 673546.53125 6581539.9375, "
        "673542.6875 6581537.6875, 673550.03125 6581524.96875, 673565.90625 6581534.21875, "
        "673550.4375 6581561.1875, 673568.65625 6581571.65625, 673570.1875 6581575.90625, "
        "673574.125 6581574.46875, 673573 6581570.375, 673597.03125 6581526.8125, "
        "673585.34375 6581520.125, 673592.40625 6581508.09375, 673616.03125 6581521.71875, "
        "673576.1875 6581592.21875))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    candidate = cleaning_footprints._iteratively_open_polygon_short_edges(
        polygon,
        target_scale=0.5,
        grid=0.03125,
        diagnostics=diagnostics,
        operator_prefix="short_edge_angle_open",
    )

    assert candidate is not None
    signature = cleaning_footprints._polygon_defect_signature(
        candidate.polygon,
        target_scale=0.5,
    )
    assert signature.short_edge_count == 0
    assert signature.min_edge_length is not None
    assert signature.min_edge_length >= 0.5 - 1e-12


def test_apply_local_polygon_repairs_resolves_case62_style_self_clearance_junction():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((673578.34375 6581552.40625, 673569.71875 6581568, "
        "673563.4375 6581564.21875, 673570.09375 6581549.90625, "
        "673565.09375 6581547.46875, 673570.875 6581537.34375, "
        "673574.875 6581539.625, 673577.5 6581533.46875, "
        "673585.34375 6581537.09375, 673582.5 6581543.03125, "
        "673574.96875 6581539.5, 673570.90625 6581547.90625, "
        "673578.34375 6581552.40625))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._apply_local_polygon_repairs(
        [polygon],
        [[10]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
        stage_prefix="post_contact_local_defect_repair",
        enable_defect_operators=True,
        enable_simplify_operators=True,
    )

    assert source_map == [[10]]
    signature = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    assert cleaning_footprints._signature_satisfies_scale_contract(
        signature,
        target_scale=0.5,
        grid=0.03125,
    )
    assert diagnostics["post_contact_local_defect_repair_applied"] is True
    assert any(
        operator.startswith("self_clearance_connector")
        or operator.startswith("hole_pair_clearance_merge")
        for operator in diagnostics["post_contact_local_defect_repair_operator_applied"]
    )


def test_apply_local_polygon_repairs_resolves_case23_style_hole_hole_wall():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((674348.09375 6579639.875, 674351.9375 6579626.1875, "
        "674340.28125 6579622.59375, 674342.3125 6579614.8125, "
        "674354.21875 6579618.09375, 674355.75 6579612.75, "
        "674323.1875 6579602.96875, 674321.375 6579609.1875, "
        "674335.28125 6579613.15625, 674333.3125 6579619.90625, "
        "674321.46875 6579616.8125, 674309.5 6579614.59375, "
        "674316.71875 6579590.6875, 674357.90625 6579602.78125, "
        "674425.0625 6579622.5625, 674410.75 6579670.5625, "
        "674381.1875 6579660.8125, 674362.71875 6579655.5625, "
        "674345.3125 6579650.96875, 674325.375 6579644.4375, "
        "674302.4375 6579638.03125, 674306.28125 6579625.25, "
        "674323.625 6579630.34375, 674324.71875 6579626.125, "
        "674307.59375 6579620.96875, 674308.25 6579618.75, "
        "674322.6875 6579621.34375, 674345.4375 6579628.15625, "
        "674343.9375 6579631.9375, 674334.6875 6579629.15625, "
        "674333.78125 6579632.1875, 674329.03125 6579630.9375, "
        "674328.09375 6579634.4375, 674348.09375 6579639.875), "
        "(674407.21875 6579635.71875, 674405.84375 6579635.28125, "
        "674407.8125 6579629.15625, 674403.28125 6579627.75, "
        "674401.03125 6579635.21875, 674413.09375 6579638.625, "
        "674415.28125 6579631.375, 674410.09375 6579629.875, "
        "674408.96875 6579634.09375, 674407.875 6579633.71875, "
        "674407.21875 6579635.71875), "
        "(674412.40625 6579643.6875, 674408.59375 6579642.59375, "
        "674407.4375 6579646.21875, 674400.6875 6579644.1875, "
        "674398.3125 6579651.75, 674400.6875 6579652.375, "
        "674398.84375 6579659.5625, 674396.1875 6579658.875, "
        "674396.9375 6579656.78125, 674388.46875 6579654.46875, "
        "674387.53125 6579657.46875, 674402.375 6579662.28125, "
        "674405.0625 6579653.6875, 674409.125 6579655.09375, "
        "674412.40625 6579643.6875), "
        "(674386.53125 6579643.6875, 674390.34375 6579644.75, "
        "674392.34375 6579638.34375, 674388.53125 6579637.15625, "
        "674386.53125 6579643.6875), "
        "(674376.125 6579618.21875, 674373.0625 6579628.0625, "
        "674367.71875 6579626.53125, 674366.28125 6579631.3125, "
        "674368.96875 6579632.03125, 674366.1875 6579642.40625, "
        "674385.1875 6579647.9375, 674386.5 6579643.6875, "
        "674379.71875 6579641.6875, 674383.1875 6579630.75, "
        "674376.3125 6579628.875, 674379.46875 6579619.1875, "
        "674376.125 6579618.21875), "
        "(674361.5 6579614.125, 674353.90625 6579641.6875, "
        "674359.09375 6579643.28125, 674367 6579615.75, "
        "674361.5 6579614.125))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._apply_local_polygon_repairs(
        [polygon],
        [[29]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
        stage_prefix="post_contact_local_defect_repair",
        enable_defect_operators=True,
        enable_simplify_operators=True,
    )

    assert source_map == [[29]]
    signature = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    assert cleaning_footprints._signature_satisfies_scale_contract(
        signature,
        target_scale=0.5,
        grid=0.03125,
    )
    assert diagnostics["post_contact_local_defect_repair_applied"] is True
    assert any(
        operator.startswith("self_clearance_connector")
        or operator.startswith("hole_pair_clearance_merge")
        for operator in diagnostics["post_contact_local_defect_repair_operator_applied"]
    )


def test_apply_local_polygon_repairs_resolves_case62_hole_hole_ring_contact():
    polygon = cleaning_footprints.shapely.from_wkt(
        "POLYGON ((673804.96875 6581615.5625, 673765.21875 6581592.6875, "
        "673745.46875 6581581.65625, 673719.90625 6581567.15625, 673752.53125 6581509.96875, "
        "673792.53125 6581532, 673858 6581567.40625, 673923.4375 6581603.71875, "
        "673889.375 6581664.21875, 673864.0625 6581649.09375, 673838.15625 6581634.3125, "
        "673804.96875 6581615.5625), (673850.8125 6581624, 673872.125 6581635.4375, "
        "673879.8125 6581618.875, 673885 6581609.78125, 673899.46875 6581617.96875, "
        "673904.46875 6581609.09375, 673877.0625 6581593.90625, 673867.1875 6581611.59375, "
        "673855.0625 6581604.03125, 673854.84375 6581604.4375, 673842.59375 6581596.71875, "
        "673844 6581594, 673829.78125 6581586.28125, 673827.78125 6581589.78125, "
        "673821.90625 6581586.75, 673812.15625 6581603.3125, 673818.84375 6581606.875, "
        "673823.875 6581597.96875, 673835.9375 6581604.46875, 673830.8125 6581613.3125, "
        "673840.125 6581618.28125, 673845 6581609.34375, 673855.6875 6581615.125, "
        "673850.8125 6581624), (673885.21875 6581643.125, 673872.125 6581635.4375, "
        "673869.09375 6581640.5625, 673882.375 6581648.15625, 673885.21875 6581643.125), "
        "(673855.6875 6581593.25, 673863.9375 6581597.8125, 673862.0625 6581601.40625, "
        "673866.6875 6581604.03125, 673873.28125 6581591.8125, 673855.59375 6581582.03125, "
        "673849.125 6581594.125, 673853.84375 6581596.78125, 673855.6875 6581593.25), "
        "(673750.09375 6581541.03125, 673753.8125 6581543.15625, 673756.125 6581539.03125, "
        "673770.8125 6581547.375, 673766.3125 6581555.65625, 673775.5625 6581560.90625, "
        "673772.5625 6581565.71875, 673763.65625 6581560.59375, 673756.625 6581573.5625, "
        "673764.90625 6581578, 673760.5 6581585.09375, 673763.8125 6581586.90625, "
        "673773.21875 6581570.71875, 673777.125 6581572.96875, 673771.84375 6581581.71875, "
        "673775.8125 6581583.84375, 673781 6581575.21875, 673793.375 6581582.34375, "
        "673788.375 6581590.5625, 673808.3125 6581601.25, 673811.71875 6581595.65625, "
        "673795.875 6581585.6875, 673801.8125 6581575.5, 673814.125 6581582.3125, "
        "673818.6875 6581574.28125, 673840.0625 6581586.15625, 673845.59375 6581576.5, "
        "673833.90625 6581570, 673831.6875 6581573.625, 673825.65625 6581570.40625, "
        "673827.75 6581566.59375, 673785.84375 6581543.40625, 673783.75 6581547, "
        "673773.9375 6581541.5625, 673775.90625 6581537.90625, 673757.6875 6581527.78125, "
        "673750.09375 6581541.03125), (673762.875 6581553.71875, 673747.90625 6581544.90625, "
        "673737.40625 6581563.25, 673752.0625 6581571.09375, 673762.875 6581553.71875))"
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    polygons, source_map = cleaning_footprints._apply_local_polygon_repairs(
        [polygon],
        [[26]],
        min_segment_length=0.5,
        grid=0.03125,
        min_area=0.0,
        min_hole_area=0.0,
        diagnostics=diagnostics,
        stage_prefix="post_contact_local_defect_repair",
        enable_defect_operators=True,
        enable_simplify_operators=True,
    )

    assert source_map == [[26]]
    signature = cleaning_footprints._polygon_defect_signature(
        polygons[0],
        target_scale=0.5,
    )
    assert signature.ring_contact_count == 0
    assert signature.short_edge_count == 0
    assert diagnostics["post_contact_local_defect_repair_applied"] is True
    assert any(
        operator.startswith("ring_contact_connector")
        for operator in diagnostics["post_contact_local_defect_repair_operator_applied"]
    )
    assert any(
        "short_edge_angle_open" in operator
        for operator in diagnostics["post_contact_local_defect_repair_operator_applied"]
    )


def test_output_canonicalization_falls_back_on_area_imbalanced_snap(monkeypatch):
    polygon = Polygon(
        [
            (0, 0),
            (8, 0),
            (8, 4),
            (8.4, 4),
            (8.4, 5),
            (8, 5),
            (8, 8),
            (0, 8),
        ]
    )
    diagnostics = cleaning_footprints._empty_diagnostics(1)
    diagnostics["collect_stage_metrics"] = False

    def fake_canonicalize(geom, grid, diagnostics):
        if grid == 0.125:
            return [box(0, 0, 8, 8)]
        return [geom]

    monkeypatch.setattr(cleaning_footprints, "_canonicalize", fake_canonicalize)

    polygons = cleaning_footprints._canonicalize_for_output(
        polygon,
        output_grid=0.125,
        min_area=0.0,
        diagnostics=diagnostics,
    )

    assert len(polygons) == 1
    assert polygons[0].equals_exact(polygon, tolerance=0.0)
    assert diagnostics["output_canonicalization_output_grid_applied_count"] == 0
    assert diagnostics["output_canonicalization_fallback_count"] == 1


def test_condition_polygon_coverage_rejects_bad_arguments():
    with pytest.raises(ValueError):
        cleaning.condition_polygon_coverage(
            [box(0, 0, 1, 1)],
            options=cleaning.ConditioningOptions(min_feature_size=-1.0),
        )
    with pytest.raises(ValueError):
        cleaning.condition_polygon_coverage(
            [box(0, 0, 1, 1)],
            source_map=[[0], [1]],
            options=cleaning.ConditioningOptions(),
        )


def test_cleaning_public_api_exports_are_callable():
    assert set(cleaning.__all__) == {
        "ConditioningOptions",
        "ConditioningResult",
        "condition_polygon_coverage",
        "condition_building_footprints",
    }

    options = cleaning.ConditioningOptions(
        min_feature_size=0.0,
        merge_distance=0.0,
        min_area=0.0,
        min_hole_area=0.0,
    )
    polygon_result = cleaning.condition_polygon_coverage([box(0, 0, 2, 2)], options=options)
    assert isinstance(polygon_result, cleaning.ConditioningResult)

    buildings = [
        make_building(box(0, 0, 2, 2)),
        make_building(box(2.1, 0, 4.1, 2)),
    ]
    building_result = cleaning.condition_building_footprints(
        buildings,
        lod=GeometryType.LOD0,
        options=cleaning.ConditioningOptions(
            min_feature_size=0.0,
            merge_distance=0.2,
            min_area=0.0,
            min_hole_area=0.0,
        ),
    )
    assert isinstance(building_result, cleaning.ConditioningResult)
    assert_conditioning_invariants(
        building_result,
        cleaning.ConditioningOptions(
            min_feature_size=0.0,
            merge_distance=0.2,
            min_area=0.0,
            min_hole_area=0.0,
        ),
    )
