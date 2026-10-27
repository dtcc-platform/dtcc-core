import numpy as np
import pytest
from affine import Affine
from shapely.geometry import LineString

from dtcc_core.model import Bounds, PopulationGrid, Raster
from dtcc_core.model.values.population_grid import POPULATION_GRID_CREDIT, ZERO_COLOR


def _grid(data, **kwargs):
    data = np.array(data, dtype=np.float32)
    return PopulationGrid(
        data=data,
        georef=Affine(1000.0, 0, 319000.0, 0, -1000.0, 6399000.0),
        nodata=np.nan,
        crs="EPSG:3006",
        year=kwargs.pop("year", 2025),
        variable=kwargs.pop("variable", "total"),
        requested_bounds=kwargs.pop(
            "requested_bounds", Bounds(319500, 6397500, 320500, 6398500)
        ),
        **kwargs,
    )


@pytest.fixture
def plt():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    yield plt
    plt.close("all")


def _texts(ax):
    return [text.get_text() for text in ax.texts]


def test_plot_draws_cells_with_unit_title_credit_and_note(plt):
    ax = _grid([[100, 0], [np.nan, 300]]).plot(show=False)

    assert ax.get_title() == "Total population, 31 December 2025"
    colorbar_ax = ax.child_axes[0]
    assert colorbar_ax.get_ylabel() == "Persons per 1 km cell"
    note = next(t for t in _texts(ax) if POPULATION_GRID_CREDIT in t)
    assert "whole 1 km cells" in note and "beyond the requested bounds" in note
    assert ax.get_aspect() == 1.0
    # Cells are drawn without smoothing on the 1 km grid.
    image = ax.images[0]
    assert image.get_interpolation() == "nearest"
    assert image.get_extent() == [319000.0, 321000.0, 6397000.0, 6399000.0]


def test_plot_distinguishes_zero_and_nodata(plt):
    from matplotlib.colors import to_hex

    ax = _grid([[100, 0], [np.nan, 300]]).plot(show=False)

    inhabited, zero = ax.images[0].get_array(), ax.images[1].get_array()
    assert inhabited.mask.tolist() == [[False, True], [True, False]]
    assert zero.mask.tolist() == [[True, False], [True, True]]
    hatched = [p for p in ax.patches if p.get_hatch() == "///"]
    assert len(hatched) == 1 and hatched[0].get_xy() == (319000.0, 6397000.0)
    legend = [t.get_text() for t in ax.get_legend().get_texts()]
    assert legend == [
        "No registered residents (0)",
        "Outside Sweden (no data)",
        "Requested bounds",
    ]
    assert (
        to_hex(ax.get_legend().legend_handles[0].get_facecolor()) == ZERO_COLOR.lower()
    )
    # Counts are written in cells, except outside Sweden.
    assert {"100", "0", "300"} <= set(_texts(ax)) and "nan" not in _texts(ax)


def test_plot_outlines_requested_bounds(plt):
    ax = _grid([[1, 2], [3, 4]]).plot(show=False)

    dashed = [line for line in ax.lines if line.get_linestyle() == "--"]
    assert len(dashed) == 1
    assert list(dashed[0].get_xdata()) == [319500, 320500, 320500, 319500, 319500]


def test_plot_all_zero_grid_explains_itself(plt):
    ax = _grid([[0, 0], [0, 0]]).plot(show=False)

    assert "No registered residents in these cells" in _texts(ax)
    assert ax.images[0].get_array().mask.all()
    assert ax.child_axes == []  # no colour scale without residents


def test_plot_all_nodata_grid_explains_itself(plt):
    ax = _grid([[np.nan, np.nan]]).plot(show=False)

    assert "All cells are outside Sweden (no data)" in _texts(ax)
    assert len([p for p in ax.patches if p.get_hatch() == "///"]) == 2


def test_plot_titles_age_bands_and_skips_labels_for_large_grids(plt):
    ax = _grid(np.ones((9, 9)), variable="age_20_24", year=2020).plot(show=False)

    assert ax.get_title() == "Aged 20–24, 31 December 2020"
    assert "1" not in _texts(ax)


def test_plot_draws_area_boundaries(plt):
    boundaries = LineString([(319000, 6398000), (321000, 6398000)])

    ax = _grid([[1, 2], [3, 4]], area_boundaries=boundaries).plot(show=False)

    assert any(list(line.get_xdata()) == [319000.0, 321000.0] for line in ax.lines)


def test_info_counts_cells_and_persons():
    rows = dict(_grid([[100, 0], [np.nan, 300]])._info_sections()[0][2])

    assert rows["Cells with residents"] == 2
    assert rows["Cells outside Sweden (nodata)"] == 1
    assert rows["Persons in cells"] == 400
    assert (rows["Variable"], rows["Year"]) == ("total", 2025)


def test_to_raster_returns_plain_raster_that_serializes():
    grid = _grid([[100, 0], [np.nan, 300]])

    raster = grid.to_raster()

    assert type(raster) is Raster
    assert raster.georef == grid.georef
    assert np.array_equal(raster.data, grid.data, equal_nan=True)
    assert raster.to_proto() is not None


def test_to_raster_can_replace_nodata():
    raster = _grid([[100, 0], [np.nan, 300]]).to_raster(nodata=-1.0)

    assert raster.nodata == -1.0
    assert raster.data.tolist() == [[100.0, 0.0], [-1.0, 300.0]]
