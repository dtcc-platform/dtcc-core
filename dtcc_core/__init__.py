from __future__ import annotations

import importlib

from dtcc_core.common.dtcc_logging import (
    init_logging,
    set_log_level,
    get_logger,
    get_python_logger,
)
from dtcc_core import common
from dtcc_core import model
from dtcc_core import builder

register_model_method = builder.register_model_method
Bounds = model.Bounds

_LAZY_SUBMODULES = {
    "io": "dtcc_core.io",
    "datasets": "dtcc_core.datasets",
    "plotting": "dtcc_core.plotting",
    "reproject": "dtcc_core.reproject",
}

_LAZY_LOADERS = (
    "load_3dbag", "load_model", "load_city", "load_cityjson",
    "load_mesh", "load_volume_mesh", "load_mesh_as_city",
    "load_pointcloud", "load_pointcloud_directory", "load_raster",
    "load_footprints", "load_landuse", "load_roadnetwork",
)


def __getattr__(name: str):
    if name in _LAZY_SUBMODULES:
        module = importlib.import_module(_LAZY_SUBMODULES[name])
        globals()[name] = module
        return module
    if name in _LAZY_LOADERS:
        loader = getattr(importlib.import_module("dtcc_core.io"), name)
        globals()[name] = loader
        return loader
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


__all__ = [
    "init_logging",
    "set_log_level",
    "get_logger",
    "get_python_logger",
    "common",
    "model",
    "builder",
    "Bounds",
    *_LAZY_LOADERS,
    "register_model_method",
    "io",
    "datasets",
    "plotting",
    "reproject",
]
