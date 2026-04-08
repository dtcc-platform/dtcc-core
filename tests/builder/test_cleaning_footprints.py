import itertools

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
        "source_coordinate_recovery_short_edge_count_before",
        "source_coordinate_recovery_short_edge_count_after",
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
    assert (
        result.diagnostics["coverage_meshing_regularization_ring_contact_count_before"]
        > 0
    )
    assert (
        result.diagnostics["coverage_meshing_regularization_ring_contact_count_after"]
        == 0
    )
    assert (
        result.diagnostics["coverage_meshing_regularization_ring_contact_polygon_count"]
        == 1
    )
    assert_conditioning_invariants(result, options)


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
    original_coverage_simplify = cleaning_footprints.shapely.coverage_simplify

    def record_coverage_simplify(*args, **kwargs):
        events.append("coverage_simplify")
        return original_coverage_simplify(*args, **kwargs)

    monkeypatch.setattr(
        cleaning_footprints.shapely,
        "coverage_simplify",
        record_coverage_simplify,
    )
    diagnostics = cleaning_footprints._empty_diagnostics(2)
    diagnostics["collect_stage_metrics"] = False

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
        (operator, polygons, source_map)
        for operator, polygons, source_map in candidates
        if operator.startswith("coverage_pair_issue_bridge_")
    ]
    assert bridge_candidates
    operator, polygons, source_map = bridge_candidates[0]
    signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=0.5,
    )
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
        (operator, polygons, source_map)
        for operator, polygons, source_map in candidates
        if operator.startswith("coverage_pair_issue_shrink_")
    ]

    assert shrink_candidates
    assert any(
        cleaning_footprints._coverage_defect_signature(
            polygons,
            target_scale=0.5,
        ).pair_issue_count
        == 0
        and source_map == [[0], [1]]
        for _, polygons, source_map in shrink_candidates
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
