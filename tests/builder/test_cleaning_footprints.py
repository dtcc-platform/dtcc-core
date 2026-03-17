import itertools

import pytest
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, box

import dtcc_core.builder.cleaning as cleaning
from dtcc_core.model import Building, GeometryType, Surface


def make_building(polygon: Polygon, height: float = 10.0) -> Building:
    surface = Surface()
    surface.from_polygon(polygon, height)
    building = Building()
    building.add_geometry(surface, GeometryType.LOD0)
    building.attributes["height"] = height
    building.attributes["ground_height"] = 0.0
    return building


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
    }
    assert required_keys.issubset(result.diagnostics)
    assert result.diagnostics["output_count"] == len(result.polygons)

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
