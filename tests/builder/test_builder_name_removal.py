"""
Verify that legacy meshing builder names have been removed.

These tests ensure the old ``build_terrain_mesh`` and ``build_city_mesh``
symbols are gone from the public API.  Any re-introduction of those names
will cause these tests to fail.
"""

import importlib

import pytest

import dtcc_core.builder


class TestLegacyNamesRemoved:
    """Old builder names must not exist anywhere in the public API."""

    def test_build_terrain_mesh_not_in_builder(self):
        assert not hasattr(dtcc_core.builder, "build_terrain_mesh")

    def test_build_city_mesh_not_in_builder(self):
        assert not hasattr(dtcc_core.builder, "build_city_mesh")

    def test_import_build_terrain_mesh_fails(self):
        with pytest.raises(ImportError):
            from dtcc_core.builder import build_terrain_mesh  # noqa: F401

    def test_import_build_city_mesh_fails(self):
        with pytest.raises(ImportError):
            from dtcc_core.builder import build_city_mesh  # noqa: F401


class TestCanonicalNamesExist:
    """Canonical builder names must be importable."""

    def test_build_terrain_surface_mesh_exists(self):
        assert hasattr(dtcc_core.builder, "build_terrain_surface_mesh")

    def test_build_city_surface_mesh_exists(self):
        assert hasattr(dtcc_core.builder, "build_city_surface_mesh")

    def test_build_city_volume_mesh_exists(self):
        assert hasattr(dtcc_core.builder, "build_city_volume_mesh")

    def test_build_city_flat_mesh_exists(self):
        assert hasattr(dtcc_core.builder, "build_city_flat_mesh")

    def test_import_build_terrain_surface_mesh(self):
        from dtcc_core.builder import build_terrain_surface_mesh  # noqa: F401

    def test_import_build_city_surface_mesh(self):
        from dtcc_core.builder import build_city_surface_mesh  # noqa: F401

    def test_import_build_city_flat_mesh(self):
        from dtcc_core.builder import build_city_flat_mesh  # noqa: F401
