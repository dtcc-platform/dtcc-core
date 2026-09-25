"""Public-boundary tests for the authoritative footprint constructor."""

import pytest
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, box
from shapely.ops import unary_union

import dtcc_core.builder.cleaning as cleaning
from dtcc_core.builder.cleaning import construction
from dtcc_core.builder.cleaning.contract import (
    check_cleaning_contract,
    check_mesher_handoff_profile,
)
from dtcc_core.model import Building, GeometryType, Surface
from dtcc_core.plotting.style import DTCC_THEMES


def make_building(polygon, height=10.0):
    surface = Surface()
    surface.from_polygon(polygon, height)
    building = Building()
    building.add_geometry(surface, GeometryType.LOD0)
    return building


@pytest.mark.parametrize("show_changes", [False, True])
def test_plot_footprint_cleaning_comparison_returns_figure(show_changes):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg", force=True)
    to_rgba = pytest.importorskip("matplotlib.colors").to_rgba
    plt = pytest.importorskip("matplotlib.pyplot")

    figure, axes = cleaning.plot_footprint_cleaning_comparison(
        [box(0.0, 0.0, 4.0, 4.0)],
        [box(0.0, 0.0, 3.5, 3.5)],
        show=False,
        show_changes=show_changes,
    )

    assert len(axes) == (3 if show_changes else 2)
    assert axes[0].get_title() == "Input footprints"
    assert axes[1].get_title() == "Conditioned footprints"
    assert axes[0].get_facecolor() == to_rgba(DTCC_THEMES["dark"]["axes"])
    if show_changes:
        assert [text.get_text() for text in axes[2].get_legend().get_texts()] == [
            "Added: 0.00 m²",
            "Removed: 3.75 m²",
        ]
    plt.close(figure)


def test_unscaled_identity_preserves_parts_sources_and_is_deterministic():
    raw = [MultiPolygon([box(4, 0, 6, 2), box(0, 0, 2, 2)]), box(8, 0, 10, 2)]
    options = cleaning.ConditioningOptions(
        min_feature_size=0.0,
        merge_distance=0.0,
        min_hole_area=0.0,
        enable_logging=False,
    )

    first = cleaning.condition_polygon_coverage(
        raw, source_map=[[7, 3, 7], [2]], options=options
    )
    second = cleaning.condition_polygon_coverage(
        raw, source_map=[[7, 3, 7], [2]], options=options
    )

    assert [polygon.wkb for polygon in first.polygons] == [
        polygon.wkb for polygon in second.polygons
    ]
    assert first.source_map == second.source_map
    assert sorted(first.source_map) == [[2], [3, 7], [3, 7]]
    assert first.diagnostics["outcome"] == "unscaled_identity"
    assert first.diagnostics["mesher_profile"]["status"] == "pass"


def test_raw_boundary_rejects_non_polygonal_container():
    with pytest.raises(ValueError, match="Polygon or MultiPolygon"):
        cleaning.condition_polygon_coverage(
            [GeometryCollection([box(0, 0, 2, 2)])],
            options=cleaning.ConditioningOptions(
                min_feature_size=0.0,
                enable_logging=False,
            ),
        )


def test_positive_scale_uses_authoritative_constructor_and_reports_contract(
    monkeypatch,
):
    polygon = box(0, 0, 10, 10)
    calls = []

    def fake_construct(
        raw,
        source_map,
        *,
        delta,
        epsilon,
        merge_distance,
        allow_source_merging,
        allow_residual_separation,
    ):
        assert allow_residual_separation is True
        calls.append(
            (
                raw,
                source_map,
                delta,
                epsilon,
                merge_distance,
                allow_source_merging,
            )
        )
        return [polygon], [[4]], {
            "outcome": "unchanged",
            "before_selection_contract": {
                "status": "pass",
                "fidelity": {"status": "pass"},
            },
            "mesher_profile": {"status": "pass"},
        }

    monkeypatch.setattr(construction, "construct_coverage", fake_construct)
    result = cleaning.condition_polygon_coverage(
        [polygon],
        source_map=[[4]],
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.0,
            fidelity_tolerance=0.2,
            enable_logging=False,
        ),
    )

    assert calls == [([polygon], [[4]], 0.5, 0.2, 0.0, False)]
    assert result.polygons == [polygon]
    assert result.source_map == [[4]]
    assert result.diagnostics["fidelity"]["status"] == "pass"


def test_unresolved_result_is_a_typed_failure(monkeypatch):
    diagnostics = {
        "outcome": "unresolved",
        "unresolved_groups": [{"group": 12}],
        "before_selection_contract": {"status": "fail"},
    }
    monkeypatch.setattr(
        construction,
        "construct_coverage",
        lambda *args, **kwargs: (None, None, diagnostics),
    )

    with pytest.raises(
        cleaning.UnresolvedFootprintCleaningError, match="unresolved groups: 12"
    ) as error:
        cleaning.condition_polygon_coverage(
            [box(0, 0, 1, 1)],
            options=cleaning.ConditioningOptions(enable_logging=False),
        )
    assert error.value.diagnostics["outcome"] == "unresolved"


def test_residual_separation_warns_by_default_without_relaxing_fidelity(monkeypatch):
    from dtcc_core.builder.cleaning import footprints

    messages = []
    monkeypatch.setattr(footprints, "warning", messages.append)
    raw = [box(0, 0, 5, 5), box(5.1, 0, 10, 5)]
    options = cleaning.ConditioningOptions(fidelity_tolerance=0, enable_logging=False)
    result = cleaning.condition_polygon_coverage(raw, options=options)
    contract = result.diagnostics["before_selection_contract"]
    assert result.diagnostics["outcome"] == "warning"
    assert contract["status"] == "fail"
    assert contract["fidelity"]["status"] == "pass"
    assert contract["admissibility"]["topology_ok"] is True
    assert result.diagnostics["mesher_profile"]["status"] == "pass"
    assert unary_union(result.polygons).equals(unary_union(raw))
    assert len(messages) == 1 and "CONTRACT NOT SATISFIED" in messages[0]

    options.allow_residual_separation = False
    with pytest.raises(cleaning.UnresolvedFootprintCleaningError):
        cleaning.condition_polygon_coverage(raw, options=options)


@pytest.mark.parametrize("failure", ["topology", "sector", "fidelity", "borderline"])
def test_warning_policy_does_not_admit_other_contract_failures(monkeypatch, failure):
    raw = [box(0, 0, 5, 5), box(5.1, 0, 10, 5)]
    if failure == "topology":
        raw = [box(0, 0, 5, 5), box(5, 5, 10, 10)]
    elif failure == "sector":
        raw = [Polygon([(0, 0), (10, 0), (0, 0.01)])]
    else:
        original = construction.check_cleaning_contract

        def failed_fidelity(*args, **kwargs):
            report = original(*args, **kwargs)
            report["status"] = "fail"
            report["fidelity"]["status"] = "fail" if failure == "fidelity" else "borderline"
            return report

        monkeypatch.setattr(construction, "check_cleaning_contract", failed_fidelity)
    with pytest.raises(cleaning.UnresolvedFootprintCleaningError):
        cleaning.condition_polygon_coverage(
            raw, options=cleaning.ConditioningOptions(
                fidelity_tolerance=0, enable_logging=False,
            ),
        )


def test_simple_positive_scale_result_satisfies_shared_contracts():
    raw = [box(0, 0, 10, 10)]
    result = cleaning.condition_polygon_coverage(
        raw,
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.0,
            enable_logging=False,
        ),
    )

    assert check_cleaning_contract(raw, result.polygons, delta=0.5)["status"] == "pass"
    assert check_mesher_handoff_profile(result.polygons)["status"] == "pass"
    assert result.diagnostics["before_selection_contract"]["status"] == "pass"


def test_area_selection_is_explicit_and_does_not_change_cleaned_geometry():
    result = cleaning.condition_polygon_coverage(
        [box(0, 0, 1, 1), box(3, 0, 6, 3)],
        source_map=[[7], [9]],
        options=cleaning.ConditioningOptions(
            min_feature_size=0.0,
            merge_distance=0.0,
            enable_logging=False,
        ),
    )
    original_wkb = [polygon.wkb for polygon in result.polygons]

    selected = cleaning.select_footprints(result, min_area=2.0)

    assert [polygon.wkb for polygon in result.polygons] == original_wkb
    assert selected.source_map == [[9]]
    assert selected.diagnostics["selection"]["unrepresented_source_indices"] == [7]
    assert selected.diagnostics["selection"]["removed_area"] == pytest.approx(1.0)


def test_building_adapter_keeps_original_building_indices():
    buildings = [
        make_building(box(0, 0, 2, 2)),
        Building(),
        make_building(box(4, 0, 6, 2)),
    ]
    result = cleaning.condition_building_footprints(
        buildings,
        lod=GeometryType.LOD0,
        options=cleaning.ConditioningOptions(
            min_feature_size=0.0,
            merge_distance=0.0,
            enable_logging=False,
        ),
    )
    assert sorted(result.source_map) == [[0], [2]]


@pytest.mark.parametrize(
    "options",
    [
        cleaning.ConditioningOptions(min_feature_size=-1.0),
        cleaning.ConditioningOptions(merge_distance=float("nan")),
        cleaning.ConditioningOptions(precision_grid=0.0),
        cleaning.ConditioningOptions(fidelity_tolerance=float("inf")),
    ],
)
def test_condition_polygon_coverage_rejects_bad_options(options):
    with pytest.raises(ValueError):
        cleaning.condition_polygon_coverage([box(0, 0, 1, 1)], options=options)


def test_condition_polygon_coverage_rejects_bad_source_maps():
    options = cleaning.ConditioningOptions(enable_logging=False)
    with pytest.raises(ValueError, match="length"):
        cleaning.condition_polygon_coverage(
            [box(0, 0, 1, 1)], source_map=[[0], [1]], options=options
        )
    with pytest.raises(ValueError, match="must not be empty"):
        cleaning.condition_polygon_coverage(
            [box(0, 0, 1, 1)], source_map=[[]], options=options
        )
    with pytest.raises(TypeError, match="integers"):
        cleaning.condition_polygon_coverage(
            [box(0, 0, 1, 1)], source_map=[["0"]], options=options
        )
    with pytest.raises(TypeError, match="integers"):
        cleaning.condition_polygon_coverage(
            [box(0, 0, 1, 1)], source_map=[[True]], options=options
        )
    with pytest.raises(ValueError, match="non-negative"):
        cleaning.condition_polygon_coverage(
            [box(0, 0, 1, 1)], source_map=[[-1]], options=options
        )


@pytest.mark.parametrize("name", ["allow_source_merging", "allow_residual_separation"])
def test_condition_polygon_coverage_rejects_nonboolean_permissions(name):
    with pytest.raises(TypeError, match="boolean"):
        cleaning.condition_polygon_coverage(
            [box(0, 0, 1, 1)],
            options=cleaning.ConditioningOptions(**{name: "false"}),
        )


def test_cleaning_public_api_contains_only_supported_entry_points():
    assert set(cleaning.__all__) == {
        "ConditioningOptions",
        "ConditioningResult",
        "UnresolvedFootprintCleaningError",
        "condition_polygon_coverage",
        "condition_building_footprints",
        "plot_footprint_cleaning_comparison",
        "select_footprints",
    }


def test_merge_distance_controls_original_source_merge_eligibility():
    raw = [box(0, 0, 4, 4), box(4.2, 0, 8.2, 4)]

    separate = cleaning.condition_polygon_coverage(
        raw,
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.01,
            collect_stage_metrics=False,
            enable_logging=False,
        ),
    )
    merged = cleaning.condition_polygon_coverage(
        raw,
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.5,
            collect_stage_metrics=False,
            enable_logging=False,
        ),
    )

    assert len(separate.polygons) == 2
    assert sorted(separate.source_map) == [[0], [1]]
    assert len(merged.polygons) == 1
    assert merged.source_map == [[0, 1]]


def test_merge_eligibility_does_not_assign_disconnected_sources():
    raw = [box(0, 0, 4, 4), box(5, 0, 9, 4)]
    result = cleaning.condition_polygon_coverage(
        raw, source_map=[[3], [7]],
        options=cleaning.ConditioningOptions(merge_distance=2, enable_logging=False),
    )
    assert result.diagnostics["outcome"] == "unchanged"
    assert result.diagnostics["groups"] == 1
    assert {polygon.bounds: owners for polygon, owners in zip(
        result.polygons, result.source_map
    )} == {raw[0].bounds: [3], raw[1].bounds: [7]}


def test_source_attribution_retains_resolved_small_scale_input():
    result = cleaning.condition_polygon_coverage(
        [box(0, 0, 1e-6, 1e-6)], source_map=[[7]],
        options=cleaning.ConditioningOptions(
            min_feature_size=1e-7, enable_logging=False,
        ),
    )
    assert result.source_map == [[7]]


def test_merge_policy_handles_zero_touching_and_explicit_disable():
    raw = [box(0, 0, 4, 4), box(4, 0, 8, 4)]

    touching = cleaning.condition_polygon_coverage(
        raw,
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.0,
            allow_source_merging=True,
            collect_stage_metrics=False,
            enable_logging=False,
        ),
    )
    disabled = cleaning.condition_polygon_coverage(
        raw,
        options=cleaning.ConditioningOptions(
            min_feature_size=0.5,
            merge_distance=0.5,
            allow_source_merging=False,
            collect_stage_metrics=False,
            enable_logging=False,
        ),
    )

    assert len(touching.polygons) == 1
    assert touching.source_map == [[0, 1]]
    assert touching.diagnostics["before_selection_contract"] == check_cleaning_contract(
        raw, touching.polygons, delta=0.5
    )
    assert len(disabled.polygons) == 2
    assert sorted(disabled.source_map) == [[0], [1]]


def test_merge_distance_is_inclusive_and_transitive():
    raw = [
        box(0.0, 0.0, 4.0, 4.0),
        box(4.5, 0.0, 8.5, 4.0),
        box(9.0, 0.0, 13.0, 4.0),
    ]

    result = cleaning.condition_polygon_coverage(
        raw,
        options=cleaning.ConditioningOptions(
            min_feature_size=0.6,
            merge_distance=0.5,
            collect_stage_metrics=False,
            enable_logging=False,
        ),
    )

    assert len(result.polygons) == 1
    assert result.source_map == [[0, 1, 2]]
    assert result.diagnostics["merge_eligibility_groups"] == [[0, 1, 2]]
