import pytest
from shapely.geometry import Polygon

from dtcc_core.io.data.deso import filter_deso_geodataframe
from dtcc_core.model import Bounds

try:
    import geopandas as gpd
except ImportError:
    gpd = None


pytestmark = pytest.mark.skipif(gpd is None, reason="GeoPandas is required")


def test_filter_deso_geodataframe_by_bounds():
    inside = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
    outside = Polygon([(10, 10), (11, 10), (11, 11), (10, 11), (10, 10)])
    gdf = gpd.GeoDataFrame(
        {"desokod": ["inside", "outside"]},
        geometry=[inside, outside],
        crs="EPSG:3006",
    )

    filtered = filter_deso_geodataframe(gdf, Bounds(-1.0, -1.0, 2.0, 2.0))

    assert filtered["desokod"].tolist() == ["inside"]
