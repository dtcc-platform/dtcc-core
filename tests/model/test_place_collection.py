import numpy as np
import pytest

from dtcc_core.model import Object, PlaceCollection
from dtcc_core.model.geometry import Point
from dtcc_core.model.object.place_collection import PLACE_CATEGORY_COLORS


def _place(name, category, x, y):
    place = Object()
    place.attributes = {"name": name, "category": category}
    place.add_geometry(Point(x=x, y=y, z=0.0), "location")
    return place


def _collection():
    collection = PlaceCollection()
    collection.add_place(_place("Panos", "food_and_drink", 1.0, 2.0))
    collection.add_place(_place("Salong Marie", "shopping", 3.0, 4.0))
    collection.add_place(_place("Café Linsen", "food_and_drink", 5.0, 6.0))
    return collection


def test_places_returns_all_children():
    assert [p.attributes["name"] for p in _collection().places()] == [
        "Panos",
        "Salong Marie",
        "Café Linsen",
    ]


def test_add_place_rejects_non_objects():
    with pytest.raises(TypeError, match="Place must be an Object"):
        PlaceCollection().add_place("Panos")


def test_to_geodataframe():
    gdf = _collection().to_geodataframe()

    assert gdf["name"].tolist() == ["Panos", "Salong Marie", "Café Linsen"]
    assert gdf.crs.to_epsg() == 3006
    assert [(p.x, p.y) for p in gdf.geometry] == [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]


def test_summary_and_info_sections():
    collection = _collection()

    assert collection._summary_items() == [("num_places", 3)]
    sections = dict((title, rows) for title, _, rows in collection._info_sections())
    assert sections["Places by category"] == [("food_and_drink", 2), ("shopping", 1)]
    assert sections["Sample places (first 3)"][0][0] == "Panos"


@pytest.fixture
def plt():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    yield plt
    plt.close("all")


def _texts(ax):
    return [text.get_text() for text in ax.texts]


def test_plot_uses_stable_category_colors_and_counted_legend(plt):
    from matplotlib.colors import to_hex

    ax = _collection().plot(show=False)

    scatter = ax.collections[-1]
    assert [to_hex(c) for c in scatter.get_facecolors()] == [
        PLACE_CATEGORY_COLORS["food_and_drink"].lower(),
        PLACE_CATEGORY_COLORS["shopping"].lower(),
        PLACE_CATEGORY_COLORS["food_and_drink"].lower(),
    ]
    legend = [text.get_text() for text in ax.get_legend().get_texts()]
    assert legend == ["food and drink (2)", "shopping (1)"]
    assert ax.get_aspect() == 1.0


@pytest.mark.parametrize("theme", ["dark", "light"])
def test_marker_outlines_differ_from_the_background(plt, theme):
    from matplotlib.colors import to_hex

    from dtcc_core.plotting.style import get_theme

    ax = _collection().plot(show=False, theme=theme)

    background = get_theme(theme)["axes"].lower()
    outlines = {to_hex(c) for c in ax.collections[-1].get_edgecolors()}
    assert background not in outlines


def test_other_category_differs_from_the_bounds_frame():
    from dtcc_core.plotting.style import get_theme

    for theme in ("dark", "light"):
        assert (
            PLACE_CATEGORY_COLORS["other"].lower() != get_theme(theme)["muted"].lower()
        )


def _lab(hex_color):
    """CIELAB coordinates of an sRGB hex colour (D65)."""
    rgb = [int(hex_color.lstrip("#")[i : i + 2], 16) / 255 for i in (0, 2, 4)]
    r, g, b = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
    x = (0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.95047
    y = 0.2126 * r + 0.7152 * g + 0.0722 * b
    z = (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.08883

    def f(t):
        return t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116

    return 116 * f(y) - 16, 500 * (f(x) - f(y)), 200 * (f(y) - f(z))


def test_category_colors_are_easy_to_tell_apart():
    from itertools import combinations

    def distance(a, b):
        return sum((p - q) ** 2 for p, q in zip(_lab(a), _lab(b), strict=True)) ** 0.5

    for (name_a, a), (name_b, b) in combinations(PLACE_CATEGORY_COLORS.items(), 2):
        assert distance(a, b) > 20, (name_a, name_b)


def test_every_dataset_category_has_a_distinct_color():
    from dtcc_core.io.data.osm_places import PLACE_CATEGORY_NAMES

    assert set(PLACE_CATEGORY_NAMES) == set(PLACE_CATEGORY_COLORS)
    assert len(set(PLACE_CATEGORY_COLORS.values())) == len(PLACE_CATEGORY_COLORS)


def test_plot_shows_attribution_labels_and_requested_bounds(plt):
    collection = _collection()
    collection.attributes = {
        "attribution": "© OpenStreetMap contributors, ODbL 1.0",
        "requested_bounds": [0.0, 0.0, 100.0, 100.0],
    }

    ax = collection.plot(show=False)

    assert "© OpenStreetMap contributors, ODbL 1.0" in _texts(ax)
    # Places lie within a few metres of each other, so only one name fits.
    assert (
        sum(name in _texts(ax) for name in ("Panos", "Salong Marie", "Café Linsen"))
        == 1
    )
    xmin, xmax = ax.get_xlim()
    assert xmin < 0.0 and xmax > 100.0
    assert any(text.endswith(" m") for text in _texts(ax))


def test_plot_labels_can_be_disabled(plt):
    ax = _collection().plot(labels=0, show=False)

    assert not any(
        name in _texts(ax) for name in ("Panos", "Salong Marie", "Café Linsen")
    )


def test_plot_hover_shows_place_details(plt):
    from matplotlib.backend_bases import MouseEvent

    collection = PlaceCollection()
    place = _place("Panos", "food_and_drink", 50.0, 50.0)
    place.attributes.update(osm_key="amenity", osm_value="restaurant")
    collection.add_place(place)
    collection.attributes = {"requested_bounds": [0.0, 0.0, 100.0, 100.0]}
    ax = collection.plot(labels=0, show=False)
    ax.figure.canvas.draw()

    x, y = ax.transData.transform((50.0, 50.0))
    ax.figure.canvas.callbacks.process(
        "motion_notify_event", MouseEvent("motion_notify_event", ax.figure.canvas, x, y)
    )

    tooltip = next(t for t in ax.texts if t.get_visible() and "Panos" in t.get_text())
    assert "amenity=restaurant" in tooltip.get_text()


def _roads():
    from dtcc_core.model import RoadNetwork

    roads = RoadNetwork()
    roads.vertices = np.array([[0.0, 0.0], [10.0, 10.0]])
    roads.edges = np.array([[0, 1]])
    return roads


def test_plot_draws_roads_underneath(plt):
    ax = _collection().plot(roads=_roads(), show=False)

    assert len(ax.collections) == 2


def test_plot_draws_context_roads_by_default(plt):
    collection = _collection()
    collection.context_roads = _roads()

    assert len(collection.plot(show=False).collections) == 2
    assert len(collection.plot(roads=False, show=False).collections) == 1
    assert len(collection.places()) == 3


def test_plot_empty_collection_explains_itself(plt):
    collection = PlaceCollection()
    collection.attributes = {"requested_bounds": [0.0, 0.0, 100.0, 100.0]}

    ax = collection.plot(show=False)

    assert "No places in the requested bounds" in _texts(ax)
    assert ax.get_legend() is None
    assert ax.get_xlim()[0] < 0.0


def test_plot_title_shows_with_presentation_panel(plt, monkeypatch):
    collection = _collection()
    context = object()
    monkeypatch.setattr(
        PlaceCollection, "dataset_context", property(lambda self: context)
    )
    monkeypatch.setattr(
        "dtcc_core.datasets.presentation.draw_presentation_panel",
        lambda *args, **kwargs: None,
    )

    ax = collection.plot(title="Food and drink", show=False)

    assert ax.get_title() == "Food and drink"
