import matplotlib
import numpy as np
import pytest
from shapely.geometry import Polygon

from dtcc_core.plotting.style import DTCC_THEMES
from dtcc_core.io.data.deso import deso_from_geodataframe
from dtcc_core.model import Bounds, DeSO, Field, GeometryType, Object

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import geopandas as gpd
except ImportError:
    gpd = None


pytestmark = pytest.mark.skipif(gpd is None, reason="GeoPandas is required")


def _deso_gdf():
    polygon = Polygon(
        [
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.0, 0.0),
        ]
    )
    return gpd.GeoDataFrame(
        {
            "desokod": ["1480C1970"],
            "kommunkod": ["1480"],
            "version": ["2025_v2"],
        },
        geometry=[polygon],
        crs="EPSG:3006",
    )


def test_deso_from_geodataframe():
    deso = deso_from_geodataframe(
        _deso_gdf(),
        bounds=Bounds(0.0, 0.0, 2.0, 2.0),
        year=2025,
    )

    assert isinstance(deso, DeSO)
    assert len(deso) == 1
    assert deso.codes == ["1480C1970"]
    assert deso.transform.srs == "EPSG:3006"
    assert deso.attributes["source"] == "SCB"
    assert deso.attributes["year"] == 2025
    assert deso.bounds.xmax == 2.0

    area = deso.areas[0]
    assert isinstance(area, Object)
    assert area.id == "1480C1970"
    assert area.attributes["kommunkod"] == "1480"
    assert GeometryType.LOD0 in area.geometry
    assert len(area.geometry[GeometryType.LOD0].surfaces) == 1


def test_deso_info_arrays_dataframe_and_plot():
    deso = deso_from_geodataframe(_deso_gdf(), year=2025)
    deso.attach_field(
        Field(
            name="population_total",
            unit="persons",
            description="Total population",
            values=np.array([1729.0]),
            dim=1,
        )
    )

    info = deso.info(print=False)
    arrays = deso.to_arrays()
    dataframe = deso.to_dataframe()
    ax = deso.plot(show=False)

    assert "DTCC DeSO" in info
    assert "Areas: 1" in info
    assert "population_total" in info
    assert arrays["codes"].tolist() == ["1480C1970"]
    assert arrays["centroids"].shape == (1, 3)
    assert np.allclose(arrays["centroids"][0, :2], [0.5, 0.5])
    assert arrays["fields"]["population_total"].tolist() == [[1729.0]]
    assert dataframe["desokod"].tolist() == ["1480C1970"]
    assert dataframe["population_total"].tolist() == [1729.0]
    assert ax is not None
    assert ax.get_title() == "DTCC DeSO 2025"
    assert "Areas: 1" in ax.texts[0].get_text()
    assert "Bounds:" in ax.texts[0].get_text()
    plt.close(ax.figure)

    legend_ax = deso.plot(column="population_total", show=False)
    assert len(legend_ax.figure.axes) == 2
    assert legend_ax.figure.axes[1].yaxis.label.get_color() == DTCC_THEMES["dark"]["text"]
    plt.close(legend_ax.figure)


def test_deso_attach_field_validates_area_count():
    deso = deso_from_geodataframe(_deso_gdf(), year=2025)

    with pytest.raises(ValueError):
        deso.attach_field(
            Field(
                name="population_total",
                unit="persons",
                values=np.array([1.0, 2.0]),
                dim=1,
            )
        )


def test_deso_protobuf_roundtrip():
    deso = deso_from_geodataframe(_deso_gdf(), year=2025)
    deso.attach_field(
        Field(
            name="population_total",
            unit="persons",
            description="Total population",
            values=np.array([1729.0]),
            dim=1,
        )
    )

    payload = deso.to_proto().SerializeToString()
    restored = DeSO()
    restored.from_proto(payload)

    assert len(restored) == 1
    assert restored.codes == ["1480C1970"]
    assert restored.attributes["year"] == 2025
    assert restored.fields["population_total"].values.tolist() == [[1729.0]]
