import dtcc_core
from dtcc_core.builder.cleaning.footprints import (
    condition_polygon_coverage,
    ConditioningOptions,
)
import shapely
from shapely.geometry import Polygon

from polyforge import fix_clearance

test_polygons = []

with open("regression_polys.txt", "r") as f:
    for line in f:
        wkt_string = line.strip()
        if len(wkt_string) > 0:
            try:
                test_polygons.append(shapely.from_wkt(wkt_string))
            except:
                print(f"failed to parse: {wkt_string}")
                continue


results = []
opts = ConditioningOptions()
opts.min_feature_size = 1.0

for idx, test_poly in enumerate(test_polygons):
    initial_clearance = test_poly.minimum_clearance
    initial_area = test_poly.area

    pf_fix = fix_clearance(test_poly, min_clearance=1.0)
    pf_fix_area = pf_fix.area
    pf_fix_clearance = pf_fix.minimum_clearance

    dtcc_fix = condition_polygon_coverage([test_poly], options=opts)
    dtcc_fix_polygon = dtcc_fix.polygons[0]
    dtcc_fix_area = dtcc_fix_polygon.area
    dtcc_fix_clearance = dtcc_fix_polygon.minimum_clearance

    result = {
        "initial_area": initial_area,
        "initial_clearance": initial_clearance,
        "pf_fix_area": pf_fix_area,
        "pf_fix_clearance": pf_fix_clearance,
        "dtcc_fix_area": dtcc_fix_area,
        "dtcc_fix_clearance": dtcc_fix_clearance,
    }
    results.append(result)


for idx, result in enumerate(results):
    print("-------")
    print(f"polygon {idx}:")
    print(
        f"  initial area: {result['initial_area']}, initial clearance: {result['initial_clearance']}"
    )
    print(
        f"  polyforge fix area: {result['pf_fix_area']}, polyforge fix clearance: {result['pf_fix_clearance']}"
    )
    print(
        f"  dtcc fix area: {result['dtcc_fix_area']}, dtcc fix clearance: {result['dtcc_fix_clearance']}"
    )
