import geopandas as gpd
import pytest
from shapely.geometry import box

from dtcc_core.io.data import urban_atlas
from dtcc_core.model import Landuse, LanduseClasses
from dtcc_core.model.geometry import MultiSurface, Surface
from dtcc_core.model.object.landuse import (
    LANDUSE_CLASS_COLORS,
    URBAN_ATLAS_CLASS_COLORS,
    URBAN_ATLAS_SHORT_NAMES,
)


@pytest.fixture
def plt():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    yield plt
    plt.close("all")


def _urban_atlas():
    rows = [
        ("14100", "Green urban areas", box(0, 0, 60, 100)),
        ("11100", "Continuous urban fabric (S.L. : > 80%)", box(60, 0, 100, 90)),
        ("12220", "Other roads and associated land", box(60, 90, 100, 100)),
    ]
    gdf = gpd.GeoDataFrame(
        [
            {"ua_code": c, "ua_class": n, "fua_name": "Göteborg", "geometry": g}
            for c, n, g in rows
        ],
        crs="EPSG:3006",
    )
    landuse = urban_atlas.landuse_from_urban_atlas(gdf)
    landuse.attributes["requested_bounds"] = [0.0, 0.0, 100.0, 100.0]
    return landuse


def _texts(ax):
    return [text.get_text() for text in ax.texts]


def test_every_urban_atlas_class_has_a_distinct_color_and_short_name():
    assert set(URBAN_ATLAS_CLASS_COLORS) == set(urban_atlas.URBAN_ATLAS_LANDUSE_MAP)
    assert set(URBAN_ATLAS_SHORT_NAMES) == set(urban_atlas.URBAN_ATLAS_LANDUSE_MAP)
    assert len(set(URBAN_ATLAS_CLASS_COLORS.values())) == 27
    assert set(LANDUSE_CLASS_COLORS) == {value.name for value in LanduseClasses}


def test_plot_colors_urban_atlas_classes_and_ranks_legend_by_area(plt):
    from matplotlib.colors import to_hex

    ax = _urban_atlas().plot(show=False)

    faces = [to_hex(patch.get_facecolor()) for patch in ax.patches]
    assert faces == [
        URBAN_ATLAS_CLASS_COLORS[c].lower() for c in ("14100", "11100", "12220")
    ]
    legend = [text.get_text() for text in ax.get_legend().get_texts()]
    assert legend == [
        "Green urban areas (60%)",
        "Continuous urban fabric (>80% sealed) (36%)",
        "Other roads (4%)",
    ]
    assert ax.get_title() == "Copernicus Urban Atlas 2018: Göteborg"
    assert any("Copernicus Land Monitoring Service" in t for t in _texts(ax))
    assert ax.get_aspect() == 1.0


def test_plot_groups_classes_beyond_the_legend_limit(plt, monkeypatch):
    import dtcc_core.model.object.landuse as landuse_module

    monkeypatch.setattr(landuse_module, "MAX_LEGEND_CLASSES", 2)

    ax = _urban_atlas().plot(show=False)

    legend = [text.get_text() for text in ax.get_legend().get_texts()]
    assert legend[-1] == "1 other classes (4%)"


def test_plot_outlines_requested_bounds_and_draws_scale_bar(plt):
    ax = _urban_atlas().plot(show=False)

    dashed = [line for line in ax.lines if line.get_linestyle() == "--"]
    assert list(dashed[0].get_xdata()) == [0.0, 100.0, 100.0, 0.0, 0.0]
    assert any(t.endswith(" m") for t in _texts(ax))


def test_plot_without_urban_atlas_codes_uses_landuse_classes(plt):
    from matplotlib.colors import to_hex

    landuse = Landuse(landuses=[LanduseClasses.WATER, LanduseClasses.FOREST])
    surfaces = MultiSurface()
    surfaces.surfaces = [
        Surface().from_polygon(box(0, 0, 1, 1), 0),
        Surface().from_polygon(box(1, 0, 3, 1), 0),
    ]
    landuse.add_geometry(surfaces, "MULTISURFACE")

    ax = landuse.plot(show=False)

    assert [to_hex(p.get_facecolor()) for p in ax.patches] == ["#3b82c4", "#2e7d32"]
    assert [t.get_text() for t in ax.get_legend().get_texts()] == [
        "Forest (67%)",
        "Water (33%)",
    ]


def test_plot_empty_urban_atlas_result_explains_coverage(plt):
    landuse = urban_atlas.landuse_from_urban_atlas(
        urban_atlas.empty_urban_atlas_geodataframe()
    )
    landuse.attributes["requested_bounds"] = [0.0, 0.0, 100.0, 100.0]

    ax = landuse.plot(show=False)

    assert any("Functional Urban Areas only" in t for t in _texts(ax))
    assert ax.get_legend() is None


def test_plot_failed_download_says_so(plt):
    landuse = urban_atlas.landuse_from_urban_atlas(
        urban_atlas.empty_urban_atlas_geodataframe()
    )
    landuse.attributes["partial_result"] = True

    ax = landuse.plot(show=False)

    assert any("download failed" in t for t in _texts(ax))
