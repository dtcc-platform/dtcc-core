"""
Guard the model methods that builder modules attach via @register_model_method.

Those methods only exist as a side effect of importing the module that defines
them, so moving a builder submodule behind a lazy import silently removes
public API from the model classes. This happened once already: deferring
`dtcc_core.builder.raster` made Raster.slope_aspect() and eleven sibling
methods appear or vanish depending on unrelated import order (issue #87).

Each check runs in a subprocess doing nothing but `import dtcc_core`, so the
result cannot be masked by modules another test happened to import first.
"""

import subprocess
import sys

import pytest

# Methods that must be present on a model class after a bare `import dtcc_core`.
# Add to these lists when registering new model methods; do not trim a list to
# make a failure go away, since a missing entry is the regression itself.
EXPECTED_MODEL_METHODS = {
    "Raster": [
        "TPI",
        "TRI",
        "VRM",
        "burn_polygons",
        "erode_small_lines",
        "fill_holes",
        "fill_small_holes",
        "remove_small_masks",
        "resample",
        "slope_aspect",
        "stats",
        "to_pointcloud",
    ],
    "City": [
        "add_flat_terrain",
        "fix_building_clearance",
        "merge_buildings",
        "remove_small_buildings",
        "simplify_buildings",
    ],
    "Building": ["get_footprint"],
    "RoadNetwork": ["to_matrix", "to_surfaces"],
    "Surface": ["mesh", "ray_intersection"],
    "MultiSurface": ["mesh", "to_polygon"],
}


def _missing_methods(class_name: str, methods: list[str]) -> list[str]:
    """Return the expected methods absent from a freshly imported dtcc_core."""
    script = (
        "import dtcc_core\n"
        f"from dtcc_core.model import {class_name} as cls\n"
        f"missing = [m for m in {methods!r} if not hasattr(cls, m)]\n"
        "print(','.join(missing))\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"import dtcc_core failed:\n{result.stderr}"
    return [name for name in result.stdout.strip().split(",") if name]


@pytest.mark.parametrize(
    "class_name,methods",
    sorted(EXPECTED_MODEL_METHODS.items()),
)
def test_model_methods_registered_on_bare_import(class_name, methods):
    missing = _missing_methods(class_name, methods)
    assert not missing, (
        f"{class_name} is missing {missing} after a bare `import dtcc_core`. "
        "A builder module that registers these methods is probably behind a "
        "lazy import; a module using @register_model_method must be imported "
        "eagerly. See issue #87."
    )


def test_register_model_method_rejects_string_annotations():
    """
    register_model_method resolves the first parameter's annotation with
    issubclass(), so a module using `from __future__ import annotations` turns
    that annotation into a string and fails at import. This test pins the
    behaviour so the failure stays loud rather than becoming a silent skip.
    """
    from dtcc_core.builder.register import register_model_method

    def to_matrix(roadnetwork: "RoadNetwork"):  # noqa: F821 - deliberately a string
        return roadnetwork

    with pytest.raises(TypeError):
        register_model_method(to_matrix)
