"""Real loader boundaries and end-to-end prism conversion without network IO."""

import meshio
import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core import io
from dtcc_core.builder import (
    build_lod1_buildings,
    clean_footprints,
    footprints_from_lod1,
)
from tests.builder.test_footprints_from_lod1 import _source, _volume
from tests.builder.test_prism_decomposition import (
    BASE,
    COURTYARD,
    TOWER,
    _input,
    _prism,
    _stack,
)


@pytest.mark.parametrize("merge", [False, True])
@pytest.mark.parametrize("case", ["box", "step", "courtyard", "two_boxes"])
def test_binary_stl_loader_conversion_cleaning_and_rebuild(tmp_path, case, merge):
    fixtures = {
        "box": (_prism(box(0, 0, 10, 10), 2, 12), 1, 1, 1000),
        "step": (_stack(BASE, TOWER), 1, 2, 3264),
        "courtyard": (_prism(COURTYARD), 1, 1, 3640),
        "two_boxes": (
            _prism(box(0, 0, 10, 10), 2, 12) + _prism(box(30, 0, 35, 5), 3, 11),
            2,
            2,
            1200,
        ),
    }
    faces, loaded_count, piece_count, volume = fixtures[case]
    soup = _input(faces, "soup")
    path = tmp_path / "source.stl"
    meshio.write(
        path,
        meshio.Mesh(soup.vertices, [("triangle", soup.faces)]),
        file_format="stl",
        binary=True,
    )
    loaded = io.load_mesh_as_city(path, merge_coplanar_surfaces=merge)
    assert len(loaded.buildings) == loaded_count
    converted = footprints_from_lod1(loaded)
    assert len(converted.buildings) == piece_count
    assert _volume(converted.buildings) == pytest.approx(volume, abs=1e-7, rel=0)
    cleaned = clean_footprints(converted, min_building_detail=0.1)
    build_lod1_buildings(cleaned.buildings, rebuild=False)
    assert _volume(cleaned.buildings) == pytest.approx(volume, abs=1e-7, rel=0)
    target = tmp_path / "rebuilt.dtcc"
    io.save_model(cleaned, target)
    restored = io.load_model(target)
    assert _volume(restored.buildings) == pytest.approx(volume, abs=1e-7, rel=0)


@pytest.mark.parametrize("merge", [False, True])
def test_noisy_stl_loader_fragmentation_is_explicit_and_never_recombined(
    tmp_path, merge
):
    soup = _input(_stack(BASE, TOWER), "soup")
    soup.vertices += np.random.default_rng(53).uniform(
        -0.0002, 0.0002, soup.vertices.shape
    )
    path = tmp_path / "noise.stl"
    meshio.write(
        path,
        meshio.Mesh(soup.vertices, [("triangle", soup.faces)]),
        file_format="stl",
        binary=True,
    )
    direct = _source("direct", [])
    direct.geometry.clear()
    direct.add_geometry(io.load_mesh(path), "lod1")
    assert _volume(footprints_from_lod1([direct])) == pytest.approx(3264, abs=1, rel=0)
    loaded = io.load_mesh_as_city(path, merge_coplanar_surfaces=merge)
    assert len(loaded.buildings) == len(soup.faces)
    with pytest.raises(ValueError):
        footprints_from_lod1(loaded)
    preserved = footprints_from_lod1(loaded, on_error="warn")
    report = preserved.attributes["footprints_from_lod1"]
    assert report["converted"] == 0 and report["preserved"] == len(loaded.buildings)
    assert len(report["failures"]) == len(loaded.buildings)
    assert [b.id for b in preserved.buildings] == [b.id for b in loaded.buildings]
